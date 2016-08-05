#include "neurox/Neurox.h"
#include <cstring>
#include <algorithm>
#include <numeric>

using namespace Neurox;

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
    assert(nargs==4);

    const int m = sizes[0]/sizeof(int);
    assert (m==mechanismsCount);

    const int * instancesCount = (const int*) args[0];
    const double * data =  (const double*) args[1];
    const int * pdata = (const int*) args[2];
    const int * nodesIndices = (const int*) args[3];

    int dataOffset=0;
    int pdataOffset=0;
    int nodesIndicesOffset=0;

    local->mechsInstances = new MechanismInstance[m];
    for (int m=0; m<mechanismsCount; m++)
    {
        MechanismInstance & instance = local->mechsInstances[m];
        instance.instancesCount = instancesCount[m];
        int totalDataSize = mechanisms[m]->dataSize * instance.instancesCount;
        int totalPdataSize = mechanisms[m]->pdataSize * instance.instancesCount;
        instance.data = new double[totalDataSize];
        instance.pdata = new int[totalPdataSize];
        instance.nodesIndices = new int[instance.instancesCount];
        memcpy(instance.data, &data[dataOffset], sizeof(double)*totalDataSize);
        memcpy(instance.pdata, &pdata[pdataOffset], sizeof(int)*totalPdataSize);
        memcpy(instance.nodesIndices, &nodesIndices[nodesIndicesOffset], sizeof(int)*instance.instancesCount);
        dataOffset += totalDataSize;
        pdataOffset += totalPdataSize;
        nodesIndicesOffset += instance.instancesCount;
    }

    neurox_hpx_unpin;;
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
     * nargs = 8 where
     * args[0] = isSoma (1 or 0)
     * args[1] = a, vector A per compartment
     * args[2] = b, vector B per compartment
     * args[3] = d, vector D per compartment
     * args[4] = v, vector V per compartment
     * args[5] = rhs, vector RHS per compartment
     * args[6] = area, vector 'area' per compartment
     * args[7] = branches, children branches
     */

    neurox_hpx_pin(Branch);
    assert(nargs==8);

    local->isSoma = *(const char*) args[0];
    local->n = sizes[1]/sizeof(double);
    local->branchesCount = sizes[7]/sizeof(hpx_t);

    const double * a = (const double*) args[1];
    const double * b = (const double*) args[2];
    const double * d = (const double*) args[3];
    const double * v = (const double*) args[4];
    const double * rhs = (const double*) args[5];
    const double * area = (const double*) args[6];
    const hpx_t * branches = (const hpx_t*) args[7];

    local->a = new double[local->n];
    local->b = new double[local->n];
    local->d = new double[local->n];
    local->v = new double[local->n];
    local->rhs = new double[local->n];
    local->area = new double[local->n];
    local->branches = new hpx_t[local->branchesCount];

    memcpy(local->a, a, sizes[1]);
    memcpy(local->b, b, sizes[2]);
    memcpy(local->d, d, sizes[3]);
    memcpy(local->v, v, sizes[4]);
    memcpy(local->rhs, rhs, sizes[5]);
    memcpy(local->area, area, sizes[6]);
    memcpy(local->branches, branches, sizes[7]);

    neurox_hpx_unpin;;
}

hpx_action_t Branch::setV = 0;
int Branch::setV_handler(const double * v, const size_t v_size)
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::setV, v, v_size);
    for (int n=0; n<local->n; n++)
        local->v[n]=*v;
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t Branch::updateV = 0;
int Branch::updateV_handler(const int * secondOrder, size_t secondOrder_size)
{
    neurox_hpx_pin(Branch);
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
    neurox_hpx_unpin_continue(local->v[0]);
}

hpx_action_t Branch::setupMatrixInitValues = 0;
int Branch::setupMatrixInitValues_handler()
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::setupMatrixInitValues);
    for (int n=0; n<local->n; n++)
    {
        local->rhs[n]=0;
        local->d[n]=0;
    }
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t Branch::setupMatrixRHS = 0;
int Branch::setupMatrixRHS_handler(const double * parentV_ptr, const size_t)
{
    neurox_hpx_pin(Branch);
    const int n = local->n;
    double *a   = local->a;
    double *b   = local->b;
    double *v   = local->v;
    double *rhs = local->rhs;
    int branchesCount = local->branchesCount;

    double returnValue = -1; //contribution to upper branch
    double dv=-1;

    //for (i = i2; i < i3; ++i))
    if (!local->isSoma)
    {
        dv = *parentV_ptr-v[0];
        rhs[0] -= b[0]*dv;
        returnValue = a[0]*dv;
    }
    for (int i=1; i<local->n; i++)
    {
        dv = v[i-1]-v[i];
        rhs[i] -= b[i]*dv;
        rhs[i-1] += a[i]*dv;
    }

    //send/receive contribution to/from branches
    hpx_t * futures = new hpx_t[branchesCount];
    void  ** addrs  = new void*[branchesCount];
    size_t * sizes  = new size_t[branchesCount];
    double * values = new double[branchesCount];
    for (int c = 0; c < local->branchesCount; c++)
    {
        futures[c] = hpx_lco_future_new(sizeof(double));
        addrs[c]   = &values[c];
        sizes[c]   = sizeof(double);
        hpx_call(local->branches[c], Branch::setupMatrixRHS, futures[c],
                 &v[n-1], sizeof(v[n-1]));
    }

    if (branchesCount > 0) //required or fails
        hpx_lco_get_all(branchesCount, futures, sizes, addrs, NULL);

    //received contributions, can now update value
    for (int c = 0; c < branchesCount; c++)
    {
        rhs[n-1] += values[c];
        hpx_lco_delete(futures[c], HPX_NULL);
    }

    delete [] futures;
    delete [] addrs;
    delete [] sizes;
    delete [] values;

    if (!local->isSoma)
        neurox_hpx_unpin_continue(returnValue);
    neurox_hpx_unpin;
}

struct BackTriangFuture
{
    double rhs;
    double b;
}; ///> future value of the back-triangulation method

hpx_action_t Branch::gaussianBackTriangulation = 0;
int Branch::gaussianBackTriangulation_handler()
{
    neurox_hpx_pin(Branch);
    int n = local->n;
    double *a   = local->a;
    double *b   = local->b;
    double *d   = local->d;
    double *rhs = local->rhs;
    int branchesCount = local->branchesCount;

    hpx_t * futures = branchesCount ? new hpx_t[branchesCount]  : nullptr;
    void  ** addrs  = branchesCount ? new void*[branchesCount]  : nullptr;
    size_t * sizes  = branchesCount ? new size_t[branchesCount] : nullptr;

    BackTriangFuture * values = branchesCount ? new BackTriangFuture[branchesCount] : nullptr;

    for (int c = 0; c < branchesCount; c++)
    {
        futures[c] = hpx_lco_future_new(sizeof (BackTriangFuture));
        addrs[c]   = &values[c];
        sizes[c]   = sizeof(BackTriangFuture);
        hpx_call(local->branches[c], Branch::gaussianBackTriangulation, futures[c]);
    }

    if (branchesCount > 0) //required or fails
        hpx_lco_get_all(branchesCount, futures, sizes, addrs, NULL);

    //bottom compartment can now be updated with children's contribution
    for (int c = 0; c < branchesCount; c++)
    {
        d[n-1]   -= values[c].b;
        rhs[n-1] -= values[c].rhs;
        hpx_lco_delete(futures[c], HPX_NULL);
    }

    double q;
    for (int i = n - 1; i >= 1; --i)
    {
        q = a[i]/d[i];
        d[i-1]   -= q * b[i];
        rhs[i-1] -= q * rhs[i];
    }

    delete [] futures;
    delete [] addrs;
    delete [] sizes;
    delete [] values;

    //value to be decremented will be sent to parent branch (except soma)
    if (!local->isSoma)
    {
        BackTriangFuture futureData;
        q = a[0] / d[0];
        futureData.b   = q * b[0];
        futureData.rhs = q * rhs[0];
        neurox_hpx_unpin_continue(futureData);
    }
    neurox_hpx_unpin;
}


hpx_action_t Branch::gaussianFwdSubstitution = 0;
int Branch::gaussianFwdSubstitution_handler(const double * parentRHS_ptr, const size_t)
{
    neurox_hpx_pin(Branch);
    double *b   = local->b;
    double *d   = local->d;
    double *rhs = local->rhs;
    int n = local->n;

    if(local->isSoma)
    {
        rhs[0] /= d[0];
    }
    else for (int i=0; i<n; i++)
    {
        rhs[i] -= b[i] * (i==0 ? *parentRHS_ptr : rhs[i-1]);
        rhs[i] /= d[i];
    }

    double childrenRHS=rhs[n-1];
    neurox_hpx_recursive_branch_sync(Branch::gaussianFwdSubstitution, &childrenRHS, sizeof(childrenRHS));
    neurox_hpx_unpin;
}

hpx_action_t Branch::setupMatrixLHS = 0;
int Branch::setupMatrixLHS_handler()
{
    neurox_hpx_pin(Branch);
    int n = local->n;
    double *a   = local->a;
    double *b   = local->b;
    double *d   = local->d;
    int branchesCount = local->branchesCount;

    double returnValue = -1; //contribution to upper branch

    //for (i = i2; i < i3; ++i))
    if (!local->isSoma)
    {
        d[0] -= b[0];
        returnValue = a[0];
    }
    for (int i=1; i<local->n; i++)
    {
        d[i]   -= b[i];
        d[i-1] -= a[i];
    }

    //send/receive contribution to/from branches
    hpx_t * futures = new hpx_t[branchesCount];
    void  ** addrs  = new void*[branchesCount];
    size_t * sizes  = new size_t[branchesCount];
    double * values = new double[branchesCount];
    for (int c = 0; c<local->branchesCount; c++)
    {
        futures[c] = hpx_lco_future_new(sizeof(double));
        addrs[c]   = &values[c];
        sizes[c]   = sizeof(double);
        hpx_call(local->branches[c], Branch::setupMatrixLHS, futures[c]);
    }

    if (branchesCount > 0) //required or fails
        hpx_lco_get_all(branchesCount, futures, sizes, addrs, NULL);

    //received contributions, can now update value
    for (int c = 0; c < branchesCount; c++)
    {
        d[n-1] -= values[c];
        hpx_lco_delete(futures[c], HPX_NULL);
    }

    delete [] futures;
    delete [] addrs;
    delete [] sizes;
    delete [] values;

    if (!local->isSoma)
        neurox_hpx_unpin_continue(returnValue);
    neurox_hpx_unpin;
}

hpx_action_t Branch::callNetReceiveFunction = 0;
int Branch::callNetReceiveFunction_handler(const int nargs, const void *args[], const size_t sizes[])
{
    neurox_hpx_pin(Branch);

    assert(nargs==3);
    /** nargs=3 where:
     * args[0] = function flag: NetReceiveInit (1) or NetReceive(0)
     * args[1] = actual time
     * args[2] = timestep size
     */
    neurox_hpx_recursive_branch_sync(Branch::callNetReceiveFunction,
        args[0], sizes[0], args[1], sizes[1], args[2], sizes[2]);

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
    neurox_hpx_recursive_branch_async_call(Branch::callModFunction, functionId_ptr, functionId_size);

    printf ("FLAG1 : function id is %d\n", *functionId_ptr);

    //only for capacitance mechanism
    if (*functionId_ptr == Mechanism::ModFunction::capacitanceCurrent
     || *functionId_ptr == Mechanism::ModFunction::capacitanceJacob)
    {
        int mechType = CAP;
        hpx_call_sync(HPX_HERE, Mechanism::callModFunction, NULL, 0,
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
                              functionId_ptr, functionId_size,
                              &(mech->type), sizeof(mech->type));
      }
      hpx_lco_wait(lco_mechs);
      hpx_lco_delete(lco_mechs, HPX_NULL);
    }

    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

//eion.c:
#define	nparm 5
#define cur 3
#define dcurdv 4
static void ion_alloc() { assert(0); } //used in secondOrderCurrent (below)

hpx_action_t Branch::secondOrderCurrent = 0;
int Branch::secondOrderCurrent_handler()
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::secondOrderCurrent);
    MechanismInstance * mechInstances= local->mechsInstances;
    for (int m=0; m<mechanismsCount; m++)
    {
        Mechanism * mech = mechanisms[m];
        if (mech->BAfunctions[Mechanism::ModFunction::alloc] ) //TODO used to be if == ion_alloc()
        {
            for (int i=0; i<mechInstances[m].instancesCount; i++)
            {
                double * data = &mechInstances[m].data[i*nparm];
                int & nodeIndex = mechInstances[m].nodesIndices[i];
                data[cur] += data[dcurdv] * local->rhs[nodeIndex]; // cur += dcurdv * rhs(ni[i])
            }
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
    neurox_hpx_register_action(1, Branch::setupMatrixRHS);
    neurox_hpx_register_action(0, Branch::setupMatrixLHS);
    neurox_hpx_register_action(1, Branch::updateV);
    neurox_hpx_register_action(0, Branch::gaussianBackTriangulation);
    neurox_hpx_register_action(1, Branch::gaussianFwdSubstitution);
    neurox_hpx_register_action(1, Branch::setV);
    neurox_hpx_register_action(1, Branch::callModFunction);
    neurox_hpx_register_action(2, Branch::callNetReceiveFunction);
    neurox_hpx_register_action(0, Branch::setupMatrixInitValues);
    neurox_hpx_register_action(2, Branch::init);
    neurox_hpx_register_action(2, Branch::initMechanismsInstances);
    neurox_hpx_register_action(2, Branch::initNetCons);
    neurox_hpx_register_action(2, Branch::queueSpikes);
    neurox_hpx_register_action(0, Branch::getSomaVoltage);
}
