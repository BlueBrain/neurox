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
     * args[8] = vdata //TODO to be removed, should be instanteated here
     * NOTE: either args[6] or args[7] (or both) must have size>0
     */

    neurox_hpx_pin(Branch);
    assert(nargs==8 || nargs==9);

    //reconstruct and and topology data;
    local->n = *(const int*) args[0];
    local->isSoma = *(const char*) args[1];
    local->data = sizes[2]==0 ? nullptr : new double[sizes[2]/sizeof(double)];
    memcpy(local->data, args[2], sizes[2]);
    local->vdata = sizes[8]==0 ? nullptr :  new void*[sizes[8]/sizeof(void*)];
    memcpy(local->vdata, args[8], sizes[8]);
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
    int vdataOffset=0;
    int instancesOffset=0;
    local->mechsInstances = new MechanismInstance[mechanismsCount];

    for (int m=0; m<mechanismsCount; m++)
    {
        MechanismInstance & instance = local->mechsInstances[m];
        Mechanism * mech = mechanisms[m];
        instance.count = instancesCount[m];

        //data, pdata, and nodesIndices arrays
        instance.data = mech->dataSize ==0 ? nullptr : local->data+dataOffset;
        instance.pdata = mech->pdataSize==0 ? nullptr : new int[mech->pdataSize * instance.count];
        memcpy(instance.pdata, &pdata[pdataOffset], sizeof(int)*(mech->pdataSize * instance.count));
        instance.nodesIndices = instance.count>0 ? new int[instance.count] : nullptr;
        memcpy(instance.nodesIndices, &nodesIndices[instancesOffset], sizeof(int)*instance.count);

        //vdata: if is point process we need to allocate the vdata (by calling bbcore_reg in mod file)
        //and assign the correct offset in pdata (offset of vdata is in pdata[1])
        for (int i=0; i<instance.count; i++)
        {
            //point pdata to the correct offset, and allocate vdata
            assert(dataOffset  <= sizes[2]/sizeof(double));
            assert(vdataOffset <= sizes[8]/sizeof(void*) );
            double * instanceData  = (double*) &local->data[dataOffset ];
            double * instanceData2 = (double*) &instance.data [i*mech->dataSize];
            int *    instancePdata = (int   *) &instance.pdata[i*mech->pdataSize];
            assert (instanceData = instanceData2); //Make sure data offsets are good so far
            if (mech->pntMap>0)
            {
                assert( (mech->type == IClamp && mech->vdataSize == 1 && mech->pdataSize == 2)
                    || ((mech->type == ProbAMPANMDA_EMS || mech->type == ProbGABAAB_EMS)
                        && mech->vdataSize == 2 && mech->pdataSize == 3));

                void** instanceVdata = (void**) &local->vdata[vdataOffset];
                Point_process * pp = (Point_process *) instanceVdata[0];
                (void) pp;
                (void) instanceVdata;
                if (mech->vdataSize==2)
                {
                    int offsetPP  = instancePdata[1];
                    assert(offsetPP  < sizes[8]/sizeof(void*));
                    int offsetRNG = instancePdata[2];
                    assert(offsetRNG < sizes[8]/sizeof(void*));
                    void * rng = instanceVdata[1];
                }
                //We will call bbcore_red on the correct vdata offset
                //local->vdata[vdataOffset];
            }
            else
            { assert (mech->vdataSize==0);}
            dataOffset  += mech->dataSize;
            pdataOffset += mech->pdataSize;
            vdataOffset += mech->vdataSize;
            instancesOffset++;
        }
    }
    assert( dataOffset == sizes[2]/sizeof(double));
    assert(pdataOffset == sizes[3]/sizeof(int));
    assert(vdataOffset == sizes[8]/sizeof(void*));
    assert(instancesOffset == sizes[5]/sizeof(int));

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

    //initializes mechanisms graphs (capacitance is excluded from graph)
    local->mechsGraphAllNodesLCO = hpx_lco_and_new(mechanismsCount-1); //excludes 'capacitance'
    local->mechsGraphLCOs = new hpx_t[mechanismsCount];
    local->mechsGraphLCOs[mechanismsMap[capacitance]] = HPX_NULL;
    int terminalMechanismsCount=0;
    for (int m=0; m<mechanismsCount; m++)
    {
        if (mechanisms[m]->type == capacitance) continue; //exclude capacitance
        local->mechsGraphLCOs[m] = hpx_lco_reduce_new(max((short) 1,mechanisms[m]->dependenciesCount),
                      sizeof(Mechanism::ModFunction), Mechanism::initModFunction, Mechanism::reduceModFunction);
        if (mechanisms[m]->successorsCount==0) //bottom of mechs graph
            terminalMechanismsCount++;
        hpx_call(target,  Branch::mechsGraphNodeFunction, local->mechsGraphAllNodesLCO, &mechanisms[m]->type, sizeof(int));
    }
    local->mechsGraphEndLCO = hpx_lco_and_new(terminalMechanismsCount);
    neurox_hpx_unpin;
}

hpx_action_t Branch::clear=0;
int Branch::clear_handler()
{
    neurox_hpx_pin(Branch);
    //hpx_lco_wait(-local->mechsGraphLco);
    hpx_lco_delete_sync(local->mechsGraphAllNodesLCO);
    local->mechsGraphAllNodesLCO=HPX_NULL;
    neurox_hpx_unpin;
}

hpx_action_t Branch::finitialize = 0;
int Branch::finitialize_handler(const double * v, const size_t v_size)
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::finitialize, v, v_size);
    //set up by finitialize.c:nrn_finitialize(): if (setv)
    for (int n=0; n<local->n; n++)
        local->v[n]=*v;

    // the INITIAL blocks are ordered so that mechanisms that write
    // concentrations are after ions and before mechanisms that read
    // concentrations.
    local->callModFunction2(Mechanism::ModFunction::before_initialize);
    local->callModFunction2(Mechanism::ModFunction::initialize);
    local->callModFunction2(Mechanism::ModFunction::after_initialize);
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

void Branch::callModFunction2(const Mechanism::ModFunction functionId)
{
    callModFunction_handler(&functionId, sizeof(functionId));
}

hpx_action_t Branch::callModFunction = 0;
int Branch::callModFunction_handler(const Mechanism::ModFunction * functionId_ptr, const size_t functionId_size)
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::callModFunction, functionId_ptr, functionId_size);

    //only for capacitance mechanism
    if (*functionId_ptr == Mechanism::ModFunction::currentCapacitance
     || *functionId_ptr == Mechanism::ModFunction::jacobCapacitance)
    {
        mechanisms[mechanismsMap[capacitance]]->callModFunction(local, *functionId_ptr);
    }
    //for all others except capacitance (mechanisms graph)
    else
    {
        //launch execution on top nodes of the branch
        for (int m=0; m<mechanismsCount;m++)
        {
            if (mechanisms[m]->type == capacitance) continue; //not capacitance
            if (mechanisms[m]->dependenciesCount > 0) continue; //not top branch
            hpx_lco_set(local->mechsGraphLCOs[m], functionId_size, functionId_ptr, HPX_NULL, HPX_NULL);
        }
        //wait for the completion of the graph by waiting at the 'end node' lco
        hpx_lco_wait_reset(local->mechsGraphEndLCO);
    }
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t Branch::mechsGraphNodeFunction = 0;
int Branch::mechsGraphNodeFunction_handler(const int * mechType_ptr, const size_t)
{
    neurox_hpx_pin(Branch);
    int type = *mechType_ptr;
    assert(type!=capacitance); //capacitance should be outside mechanisms graph
    assert(local->mechsGraphLCOs[mechanismsMap[type]] != HPX_NULL);
    Mechanism * mech = getMechanismFromType(type);

    Mechanism::ModFunction functionId;
    while (local->mechsGraphAllNodesLCO != HPX_NULL)
    {
        //wait until all dependencies have completed, and get the argument (function id) from the hpx_lco_set
        hpx_lco_get_reset(local->mechsGraphLCOs[mechanismsMap[type]], sizeof(Mechanism::ModFunction), &functionId);

        mech->callModFunction(local, functionId);

        if (mech->successorsCount==0) //bottom mechanism
            hpx_lco_set(local->mechsGraphEndLCO, 0, NULL, HPX_NULL, HPX_NULL);
        else
            for (int c=0; c<mech->successorsCount; c++)
              hpx_lco_set(local->mechsGraphLCOs[mechanismsMap[mech->successors[c]]],
                sizeof(functionId), &functionId, HPX_NULL, HPX_NULL);
    }
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
    neurox_hpx_recursive_branch_async_call(Branch::secondOrderCurrent);
    MechanismInstance * mechInstances= local->mechsInstances;
    for (int m=0; m<mechanismsCount; m++)
    {
        Mechanism * mech = mechanisms[m];
        if (!mech->isIon) continue;
        assert(nparm==mech->dataSize); //see '#define nparm 5' in eion.c

        for (int i=0; i<mechInstances[m].count; i++)
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
    neurox_hpx_register_action(1, Branch::mechsGraphNodeFunction);
    neurox_hpx_register_action(2, Branch::callNetReceiveFunction);
    neurox_hpx_register_action(2, Branch::init);
    neurox_hpx_register_action(2, Branch::initNetCons);
    neurox_hpx_register_action(2, Branch::queueSpikes);
    neurox_hpx_register_action(0, Branch::getSomaVoltage);
    neurox_hpx_register_action(0, Branch::clear);
}
