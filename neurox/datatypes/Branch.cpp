#include "neurox/Neurox.h"
#include <cstring>
#include <algorithm>
#include <numeric>

using namespace NeuroX;

Branch::~Branch()
{
    delete [] b;
    delete [] d;
    delete [] a;
    delete [] v;
    delete [] rhs;
    delete [] area;
    delete [] branches;
    for (int m=0; m<mechanismsCount; m++)
    {
        delete [] mechsInstances[m].data;
        delete [] mechsInstances[m].pdata;
        delete [] mechsInstances[m].nodesIndices;
    }
    delete [] mechsInstances;
}

hpx_action_t Branch::initMechanismsInstances = 0;
int Branch::initMechanismsInstances_handler(const int nargs, const void *args[], const size_t sizes[])
{
    /**
     * nargs = 4 where
     * args[0] = arrays of number of instances per mechanism
     * args[1] = array of data for all mechanisms
     * args[2] = array of pdata for all mechanisms
     * args[3] = array of compartment/node index where the mechanisms are applied to
     */

    neurox_hpx_pin(Branch);
    assert(local==DEBUG_BRANCH_DELETE);
    assert(nargs==4);
    const int recvMechanismsCount = sizes[0]/sizeof(int);
    assert (recvMechanismsCount == mechanismsCount);

    const int * instancesCount = (const int*) args[0];
    const double * data =  (const double*) args[1];
    const int * pdata = (const int*) args[2];
    const int * nodesIndices = (const int*) args[3];

    int dataOffset=0;
    int pdataOffset=0;
    int instancesOffset=0;

    local->mechsInstances = new MechanismInstance[mechanismsCount];
    for (int m=0; m<mechanismsCount; m++)
    {
        MechanismInstance & instance = local->mechsInstances[m];
        instance.instancesCount = instancesCount[m];

        instance.nodesIndices = instance.instancesCount>0 ? new int[instance.instancesCount] : nullptr;
        memcpy(instance.nodesIndices, &nodesIndices[instancesOffset], sizeof(int)*instance.instancesCount);
        instancesOffset += instance.instancesCount;

        int totalDataSize = mechanisms[m]->dataSize * instance.instancesCount;
        instance.data = totalDataSize>0 ?  new double[totalDataSize]: nullptr;
        memcpy(instance.data, &data[dataOffset], sizeof(double)*totalDataSize);
        dataOffset += totalDataSize;

        int totalPdataSize = mechanisms[m]->pdataSize * instance.instancesCount;
        instance.pdata = totalPdataSize>0 ? new int[totalPdataSize] : nullptr;
        memcpy(instance.pdata, &pdata[pdataOffset], sizeof(int)*totalPdataSize);
        pdataOffset += totalPdataSize;
    }
    neurox_hpx_unpin;
}

hpx_action_t Branch::initNetCons = 0;
int Branch::initNetCons_handler(const int nargs, const void *args[], const size_t sizes[])
{
    /**
     * nargs = 2 where
     * args[0] = array of netcons
     * args[1] = array of args per netcon
     */

    neurox_hpx_pin(Branch);
    assert(local==DEBUG_BRANCH_DELETE);
    assert(nargs==7);

    //inform pre-synaptic neurons that we connect (my hpx address is stored in variable "target")
    hpx_addr_t lco =  local->netcons.size() ?  local->netcons.size() : HPX_NULL;
    for (auto nc = local->netcons.begin(); nc != local->netcons.end(); nc++)
    {
        int preNeuronId = getNeuronAddr(nc->first); //false, he may not be in the network
        int e = hpx_call(preNeuronId, Neuron::addSynapseTarget, lco, &target, sizeof(target)) ;
        assert(e==HPX_SUCCESS);
    }

    neurox_hpx_unpin;
}

Branch::MechanismInstance & Branch::getMechanismInstanceFromType(int type)
{
    return mechsInstances[mechanismsMap[type]];
}

hpx_action_t Branch::init = 0;
int Branch::init_handler(const int nargs, const void *args[], const size_t sizes[])
{
    /**
     * nargs = 9 where
     * args[0] = isSoma (1 or 0)
     * args[1] = a, vector A per compartment
     * args[2] = b, vector B per compartment
     * args[3] = d, vector D per compartment
     * args[4] = v, vector V per compartment
     * args[5] = rhs, vector RHS per compartment
     * args[6] = area, vector 'area' per compartment
     * args[7] = branches, children branches (if any)
     * args[8] = parent compartment indices (if any)
     * NOTE: args[7] or args[8] or both must have size>0
     */

    neurox_hpx_pin(Branch);
    assert(nargs==9);

    DEBUG_BRANCH_DELETE = local;

    local->isSoma = *(const char*) args[0];
    local->n = sizes[1]/sizeof(double);
    local->a = new double[local->n];
    local->b = new double[local->n];
    local->d = new double[local->n];
    local->v = new double[local->n];
    local->rhs = new double[local->n];
    local->area = new double[local->n];

    memcpy(local->a, args[1], sizes[1]);
    memcpy(local->b, args[2], sizes[2]);
    memcpy(local->d, args[3], sizes[3]);
    memcpy(local->v, args[4], sizes[4]);
    memcpy(local->rhs, args[5], sizes[5]);
    memcpy(local->area, args[6], sizes[6]);

    if (sizes[7]>0)
    {
        local->branchesCount = sizes[7]/sizeof(hpx_t);
        local->branches = new hpx_t[local->branchesCount];
        memcpy(local->branches, args[7], sizes[7]);
    }
    else
    {
        local->branchesCount = 0;
        local->branches = nullptr;
    }

    if (sizes[8]>0)
    {
        local->p = new int[local->n];
        memcpy(local->p, args[8], sizes[8]);
    }
    else
    {
        local->p = NULL;
    }

    neurox_hpx_unpin;;
}

void Branch::setupTreeMatrixMinimal()
{
    //////// setup_tree_minimal() --> nrn_rhs() //////

    for (int i=0; i<n; i++)
    {
        rhs[i]=0;
        d[i]=0;
    }

    Mechanism::ModFunction func;
    func = Mechanism::ModFunction::before_breakpoint;
    callModFunction_handler(&func, sizeof(func));
    func = Mechanism::ModFunction::current;
    callModFunction_handler(&func, sizeof(func));

    /* now the internal axial currents.
    The extracellular mechanism contribution is already done.
        rhs += ai_j*(vi_j - vi)
    */
    if (p!=NULL)
    {
        double dv=0;
        for (int i=1; i<n; i++)
        {
            dv = v[p[i]] - v[i];
            rhs[i] -= b[i]*dv;
            rhs[p[i]] += a[i]*dv;
        }
    }
    else
    {
        assert(0); //BRUNO TODO: insert Hines here
    }

    //////// setup_tree_minimal() --> nrn_lhs() //////
    // calculate left hand side of
    //cm*dvm/dt = -i(vm) + is(vi) + ai_j*(vi_j - vi)
    //cx*dvx/dt - cm*dvm/dt = -gx*(vx - ex) + i(vm) + ax_j*(vx_j - vx)
    //with a matrix so that the solution is of the form [dvm+dvx,dvx] on the right
    //hand side after solving.
    //This is a common operation for fixed step, cvode, and daspk methods

    func = Mechanism::ModFunction::jacob;
    callModFunction_handler(&func, sizeof(func));

    /* now the cap current can be computed because any change to cm
     * by another model has taken effect. */
    func = Mechanism::ModFunction::jacobCapacitance;
    callModFunction_handler(&func, sizeof(func));

    /* now add the axial currents */
    if (p!=NULL)
    {
      for (int i=1; i<n; i++)
      {
          d[i] -= b[i];
          d[p[i]] -= a[i];
      }
    }
    else
    {
        assert(0);
    }
}

hpx_action_t Branch::finitialize = 0;
int Branch::finitialize_handler(const double * v, const size_t v_size)
{
    neurox_hpx_pin(Branch);
    assert(local== DEBUG_BRANCH_DELETE);
    neurox_hpx_recursive_branch_async_call(Branch::finitialize, v, v_size);
    //set up by finitialize.c:nrn_finitialize(): if (setv)
    for (int n=0; n<local->n; n++)
        local->v[n]=*v;

    // the INITIAL blocks are ordered so that mechanisms that write
    // concentrations are after ions and before mechanisms that read
    // concentrations.
    Mechanism::ModFunction func;
    func = Mechanism::ModFunction::before_initialize;
    callModFunction_handler(&func, sizeof(func));
    func = Mechanism::ModFunction::initialize;
    callModFunction_handler(&func, sizeof(func));
    func = Mechanism::ModFunction::after_initialize;
    callModFunction_handler(&func, sizeof(func));

    local->setupTreeMatrixMinimal();

    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t Branch::updateV = 0;
int Branch::updateV_handler(const int * secondOrder, size_t secondOrder_size)
{
    neurox_hpx_pin(Branch);
    assert(local==DEBUG_BRANCH_DELETE);
    neurox_hpx_recursive_branch_async_call(Branch::updateV, &secondOrder, secondOrder_size);
    for (int i=0; i<local->n; i++)
        local->v[i] += (*secondOrder ? 2 : 1) * local->rhs[i];
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t Branch::getSomaVoltage=0;
int Branch::getSomaVoltage_handler()
{
    neurox_hpx_pin(Branch);
    assert(local==DEBUG_BRANCH_DELETE);
    neurox_hpx_unpin_continue(local->v[0]);
}

hpx_action_t Branch::callNetReceiveFunction = 0;
int Branch::callNetReceiveFunction_handler(const int nargs, const void *args[], const size_t sizes[])
{
    neurox_hpx_pin(Branch);

    assert(nargs==3);
    neurox_hpx_recursive_branch_sync(Branch::callNetReceiveFunction,
        args[0], sizes[0], args[1], sizes[1], args[2], sizes[2]);

    /** nargs=3 where:
     * args[0] = function flag: NetReceiveInit (1) or NetReceive(0)
     * args[1] = actual time
     * args[2] = timestep size
     */

    const char isInitFunction = *(const char*) args[0];
    const double t = *(const double *) args[1];
    const double dt = *(const double *) args[2];

    //*sequential* delivery of received spikes
    while (!local->spikesQueue.empty() &&
           local->spikesQueue.top().deliveryTime < t+dt)
    {
        Spike spike = local->spikesQueue.top();
        if (spike.netcon->active)
        {
             int mechType = spike.netcon->mechType;
             getMechanismFromType(mechType)->callNetReceiveFunction
                     (local, &spike, isInitFunction, t, dt);
        }
        local->spikesQueue.pop();
    }
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t Branch::callModFunction = 0;
int Branch::callModFunction_handler(const Mechanism::ModFunction * functionId_ptr, const size_t functionId_size)
{
    neurox_hpx_pin(Branch);
    assert(local==DEBUG_BRANCH_DELETE);
    neurox_hpx_recursive_branch_async_call(Branch::callModFunction, functionId_ptr, functionId_size);

    //only for capacitance mechanism
    if (*functionId_ptr == Mechanism::ModFunction::currentCapacitance
     || *functionId_ptr == Mechanism::ModFunction::jacobCapacitance)
    {
        int mechType = CAP;
        hpx_call_sync(HPX_HERE, Mechanism::callModFunction, NULL, 0,
                      &local, sizeof(local),
                      functionId_ptr, functionId_size,
                      &mechType, sizeof(mechType));
    }
    else
    {
      //*parallel* execution of independent mechanisms
      int topDependenciesCount = 0;
      for (int m=0; m<mechanismsCount; m++)
        if (mechanisms[m]->isTopMechanism)
            topDependenciesCount++;
      assert(topDependenciesCount>0);

      hpx_t lco_mechs = hpx_lco_and_new(topDependenciesCount);
      for (int m=0; m<mechanismsCount; m++)
      {
        Mechanism * mech = mechanisms[m];
        if (mechanisms[m]->isTopMechanism)
            hpx_call(HPX_HERE, Mechanism::callModFunction, lco_mechs,
                     &local, sizeof(local),
                     functionId_ptr, functionId_size,
                     &(mech->type), sizeof(mech->type));
      }
      hpx_lco_wait(lco_mechs);
      hpx_lco_delete(lco_mechs, HPX_NULL);
    }

    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

//from eion.c
#define nparm 5
#define cur	3
#define dcurdv 4

hpx_action_t Branch::secondOrderCurrent = 0;
int Branch::secondOrderCurrent_handler()
{
    neurox_hpx_pin(Branch);
    assert(local==DEBUG_BRANCH_DELETE);
    neurox_hpx_recursive_branch_async_call(Branch::secondOrderCurrent);
    MechanismInstance * mechInstances= local->mechsInstances;
    for (int m=0; m<mechanismsCount; m++)
    {
        Mechanism * mech = mechanisms[m];
        if (!mech->isIon) continue;
        assert(nparm==mech->dataSize); //see '#define nparm 5' in eion.c

        for (int i=0; i<mechInstances[m].instancesCount; i++)
        {
            double * data = &mechInstances[m].data[i*mech->dataSize];
            int & nodeIndex = mechInstances[m].nodesIndices[i];
            data[cur] += data[dcurdv] * local->rhs[nodeIndex];
        }
    }
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t Branch::queueSpikes = 0;
int Branch::queueSpikes_handler(const int nargs, const void *args[], const size_t sizes[])
{
    neurox_hpx_pin(Branch);
    assert(nargs==2);
    const int * preNeuronId = (const int *) args[0];
    const double * deliveryTime = (const double*) args[1];

    //netcvode::PreSyn::send()
    for (auto nc = local->netcons.at(*preNeuronId).begin();
         nc!=local->netcons.at(*preNeuronId).end(); nc++)
    {
      hpx_lco_sema_p(local->spikesQueueMutex);
      local->spikesQueue.push( Spike(*deliveryTime, &(*nc)) );
      hpx_lco_sema_v_sync(local->spikesQueueMutex);
    }
    neurox_hpx_unpin;
}

void Branch::registerHpxActions()
{
    neurox_hpx_register_action(1, Branch::updateV);
    neurox_hpx_register_action(1, Branch::finitialize);
    neurox_hpx_register_action(1, Branch::callModFunction);
    neurox_hpx_register_action(2, Branch::callNetReceiveFunction);
    neurox_hpx_register_action(2, Branch::init);
    neurox_hpx_register_action(2, Branch::initMechanismsInstances);
    neurox_hpx_register_action(2, Branch::initNetCons);
    neurox_hpx_register_action(2, Branch::queueSpikes);
    neurox_hpx_register_action(0, Branch::getSomaVoltage);
}
