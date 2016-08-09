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
    delete [] p;
    delete [] rhs;
    delete [] area;
    delete [] branches;
    delete [] data;

    for (int m=0; m<mechanismsCount; m++)
    {
        delete [] mechsInstances[m].nodesIndices;
        delete [] mechsInstances[m].pdata;
    }
    delete [] mechsInstances;
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
     * nargs = 8 where
     * args[0] = number of compartments
     * args[1] = isSoma (char, 1 or 0)
     * args[2] = data vector (RHS, D, A, V, B, area, and mechs...)
     * args[3] = pdata vector
     * args[4] = arrays of number of instances per mechanism
     * args[5] = array of compartment/node index where the mechanisms are applied to
     * args[6] = branches, children branches (if any)
     * args[7] = parent compartment indices (if any)
     * NOTE: either args[6] or args[7] (or both) must have size>0
     */

    neurox_hpx_pin(Branch);
    assert(nargs==8);

    DEBUG_BRANCH_DELETE = local;

    //reconstruct and and topology data;
    local->n = *(const int*) args[0];
    local->isSoma = *(const char*) args[1];
    local->data = new double[sizes[2]/sizeof(double)];
    memcpy(local->data, args[2], sizes[2]);
    local->rhs = local->data + local->n*0;
    local->d = local->data + local->n*1;
    local->a = local->data + local->n*2;
    local->b = local->data + local->n*3;
    local->v = local->data + local->n*4;
    local->area = local->data + local->n*5;

    // reconstruct mechanisms
    const int * pdata = (const int*) args[3];
    const int * instancesCount = (const int*) args[4];
    const int * nodesIndices = (const int*) args[5];
    const int recvMechanismsCount = sizes[4]/sizeof(int);
    assert (recvMechanismsCount == mechanismsCount);

    int dataOffset=6*local->n;
    int pdataOffset=0;
    int instancesOffset=0;
    local->mechsInstances = new MechanismInstance[mechanismsCount];

    //allocate space for vdata
    int vdataTotalSize=0;
    for (int m=0; m<mechanismsCount; m++)
        if (mechanisms[m]->pntMap>0)
            vdataTotalSize += mechanisms[m]->vdataSize;
    local->vdata = new void*[vdataTotalSize];

    int vdataOffset=0;
    for (int m=0; m<mechanismsCount; m++)
    {
        MechanismInstance & instance = local->mechsInstances[m];
        Mechanism * mech = mechanisms[m];
        instance.instancesCount = instancesCount[m];

        instance.nodesIndices = instance.instancesCount>0 ? new int[instance.instancesCount] : nullptr;
        memcpy(instance.nodesIndices, &nodesIndices[instancesOffset], sizeof(int)*instance.instancesCount);
        instancesOffset += instance.instancesCount;

        int totalDataSize = mech->dataSize * instance.instancesCount;
        instance.data = local->data+dataOffset;

        int totalPdataSize = mech->pdataSize * instance.instancesCount;
        instance.pdata = totalPdataSize>0 ? new int[totalPdataSize] : nullptr;
        memcpy(instance.pdata, &pdata[pdataOffset], sizeof(int)*totalPdataSize);

        //vdata: if is point process we need to allocate the vdata (by calling bbcore_reg in mod file)
        //and assign the correct offset in pdata (offset of vdata is in pdata[1])
        if (mech->pntMap>0)
        {
            assert((mech->vdataSize == 1 && mech->pdataSize == 2)   //IClamp
                || (mech->vdataSize == 2 && mech->pdataSize == 3)); //ProbAMPA_EMS, ProbGABA_EMD

            //initialize vdata
            //for (int i=0; i<mech->vdataSize; i++)
            //    local->vdata[vdataOffset+i] = NULL;

            //point pdata to the correct offset, and allocate vdata
            int * pointProcPdata = (int*) &pdata[pdataOffset];
            for (int i=0; i<instance.instancesCount; i++)
            {
                int * dataProcOffsetInPdata = (int*) &pointProcPdata[i*mech->pdataSize+0];
                printf("BEFORE %d.%d :: pdata[0]==%d\n", mech->type, i, *dataProcOffsetInPdata);
                (void) *dataProcOffsetInPdata ;//TODO we will not change this for now
                printf("AFTER  %d.%d :: pdata[0]==%d\n", mech->type, i, *dataProcOffsetInPdata);

                int * pointProcOffsetInPdata = (int*) &pointProcPdata[i*mech->pdataSize+1];
                printf("BEFORE %d.%d :: pdata[1]==%d\n", mech->type, i, *pointProcOffsetInPdata);
                *pointProcOffsetInPdata = vdataOffset;
                printf("AFTER  %d.%d :: pdata[1]==%d\n", mech->type, i, *pointProcOffsetInPdata);
                //local->vdata[*pointProcOffsetInPdata] = NULL; //never used? new Point_process();

                if (mech->vdataSize==2)
                {
                    int * vdataOffsetInPdata = (int*) &pointProcPdata[i*mech->pdataSize+2]; //if any
                    printf("BEFORE %d.%d :: pdata[2]==%d\n", mech->type, i, *vdataOffsetInPdata);
                    *vdataOffsetInPdata = vdataOffset+1; //i.e. the one in vdata after the point proc
                    printf("AFTER  %d.%d :: pdata[2]==%d\n", mech->type, i, *vdataOffsetInPdata);
                    //local->vdata[*pointProcOffsetInPdata] = (pointer to Random Gen from bbcore_reg);
                    //TODO cant be called because darray is read from file!
                    //(*nrn_bbcore_read_[type])(dArray, iArray, &dk, &ik, 0, aln_cntml, d, pd, ml->_thread, &nt, 0.0);

                }
                vdataOffset+=mech->vdataSize;

                //We will call bbcore_red on the correct vdata offset
                //local->vdata[vdataOffset];
            }
        }

        dataOffset += totalDataSize;
        pdataOffset += totalPdataSize;
    }

    //reconstructs parents or branching
    if (sizes[6]>0)
    {
        local->branchesCount = sizes[6]/sizeof(hpx_t);
        local->branches = new hpx_t[local->branchesCount];
        memcpy(local->branches, args[6], sizes[6]);
    }
    else
    {
        local->branchesCount = 0;
        local->branches = nullptr;
    }

    if (sizes[7]>0)
    {
        local->p = new int[local->n];
        memcpy(local->p, args[7], sizes[7]);
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
    //for all others except capacitance
    else
    {
        for (int m=0; m<mechanismsCount; m++)
        {
            if (mechanisms[m]->type==CAP) continue;
            int mechType = mechanisms[m]->type;
            hpx_call_sync(HPX_HERE, Mechanism::callModFunction, NULL, 0,
                      &local, sizeof(local),
                      functionId_ptr, functionId_size,
                      &mechType, sizeof(mechType));
        }

        /*
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
            hpx_call(target, Mechanism::callModFunction, lco_mechs,
                     &local, sizeof(local),
                     functionId_ptr, functionId_size,
                     &(mech->type), sizeof(mech->type));
      }
      hpx_lco_wait(lco_mechs);
      hpx_lco_delete(lco_mechs, HPX_NULL);
    */
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
    neurox_hpx_register_action(2, Branch::initNetCons);
    neurox_hpx_register_action(2, Branch::queueSpikes);
    neurox_hpx_register_action(0, Branch::getSomaVoltage);
}
