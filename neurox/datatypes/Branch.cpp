#include "neurox/Neurox.h"
#include <cstring>
#include <algorithm>
#include <numeric>
#include <set>

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

hpx_action_t Branch::broadcastNetCons = 0;
int Branch::broadcastNetCons_handler(const int nargs, const void *args[], const size_t sizes[])
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::broadcastNetCons, args[0], sizes[0], args[1], sizes[1]);
    assert(nargs==2);
    const int   * neuronsIds  = (const int*)   args[0]; //list of all existing neurons ids
    const hpx_t * neuronsAddr = (const hpx_t*) args[1]; //list of all existing neurons addrs (same order as ids)

    //create index of unique pre-neurons ids
    std::set<int> netconsPreNeuronIds;
    for (auto nc : local->netcons)
        netconsPreNeuronIds.insert(nc.first);

    //inform pre-synaptic neurons (once) that we connect
    for (int i=0; i<sizes[0]/sizeof(int); i++) //for all existing neurons
        //if I'm connected to it (ie is not artificial)
        if (netconsPreNeuronIds.find(neuronsIds[i]) != netconsPreNeuronIds.end())
            //tell the neuron to add the synapse to this branch
            hpx_call_sync(neuronsAddr[i], Neuron::addSynapseTarget, &target, sizeof(target)) ;
    neurox_hpx_recursive_branch_async_wait;
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
     * args[8] = vecplay T-data
     * args[9] = vecplay Y-data
     * args[10] = vecplay Info
     * args[11] = netcons information
     * args[12] = netcons preneuron ids
     * args[13] = netcons args (size of each netcon args array in netcons info)
     * args[14] = vdata //TODO to be removed, should be instanteated here
     * NOTE: either args[6] or args[7] (or both) must have size>0
     */

    neurox_hpx_pin(Branch);
    assert(nargs==14 || nargs==15);

    local->t = inputParams->t;
    local->dt = inputParams->dt;
    local->eventsQueueMutex = hpx_lco_sema_new(1);

    //reconstruct and and topology data;
    local->n = *(const int*) args[0];
    local->isSoma = *(const char*) args[1];
    local->data = sizes[2]==0 ? nullptr : new double[sizes[2]/sizeof(double)];
    memcpy(local->data, args[2], sizes[2]);
    local->vdata = sizes[14]==0 ? nullptr :  new void*[sizes[14]/sizeof(void*)];
    memcpy(local->vdata, args[14], sizes[14]);
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
            assert(vdataOffset <= sizes[14]/sizeof(void*) );
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
                    int offsetRNG = instancePdata[2];
                    assert(offsetPP  < sizes[14]/sizeof(void*));
                    assert(offsetRNG < sizes[14]/sizeof(void*));
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
    assert(vdataOffset == sizes[14]/sizeof(void*));
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
#if PARALLEL_MECHS_DEPENDENCY==true
    local->mechsGraph = new MechanismsExecutionGraph;
    local->mechsGraph->graphLCO = hpx_lco_and_new(mechanismsCount-1); //excludes 'capacitance'
    local->mechsGraph->mechsLCOs = new hpx_t[mechanismsCount];
    local->mechsGraph->mechsLCOs[mechanismsMap[capacitance]] = HPX_NULL;
    int terminalMechanismsCount=0;
    for (int m=0; m<mechanismsCount; m++)
    {
        if (mechanisms[m]->type == capacitance) continue; //exclude capacitance
        local->mechsGraph->mechsLCOs[m] = hpx_lco_reduce_new(max((short) 1,mechanisms[m]->dependenciesCount),
                      sizeof(Mechanism::ModFunction), Mechanism::initModFunction, Mechanism::reduceModFunction);
        if (mechanisms[m]->successorsCount==0) //bottom of mechs graph
            terminalMechanismsCount++;
        hpx_call(target,  Branch::MechanismsExecutionGraph::nodeFunction,
                 local->mechsGraph->graphLCO, &mechanisms[m]->type, sizeof(int));
    }
    local->mechsGraph->endLCO = hpx_lco_and_new(terminalMechanismsCount);
#else
    local->mechsGraph = NULL;
#endif

    //reconstructs vecplay
    //* args[8] = vecplay T-data
    //* args[9] = vecplay Y-data
    //* args[10] = vecplay Info
    local->vecplayCount = sizes[10]/sizeof(Input::Coreneuron::PointProcInfo);
    local->vecplay = new VecPlayContinuouX*[local->vecplayCount];
    const double * vecplayTdata = (const double *) args[8];
    const double * vecplayYdata = (const double *) args[9];
    const Input::Coreneuron::PointProcInfo * ppis = (const Input::Coreneuron::PointProcInfo*) args[10];

    int vOffset=0;
    for (int v=0; v<local->vecplayCount; v++)
    {
        int m = mechanismsMap[ppis[v].mechType];
        double *instancesData = local->mechsInstances[m].data;
        double *pd = &(instancesData[ppis->mechInstance*mechanisms[m]->dataSize+ppis->instanceDataOffset]);
        local->vecplay[v] = new VecPlayContinuouX(pd, &vecplayTdata[vOffset], &vecplayYdata[vOffset], ppis[v].size);
        vOffset += ppis[v].size;
    }

    //reconstructs netcons
    //* args[11] = netcons information
    //* args[12] = netcons preneuron ids
    //* args[13] = netcons args (size of each netcon args array in netcons info)
    int netconsCount = sizes[11]/sizeof(NetConX);
    NetConX *netcons    = (NetConX*) args[11];
    int *netConsPreId   = (int *)    args[12];
    double *netConsArgs = (double*)  args[13];

    int argsOffset=0;
    for (int nc=0; nc<netconsCount; nc++)
    {
        local->netcons = map<int, vector<NetConX*> > (); //initialize map
        local->netcons[netConsPreId[nc]] = vector<NetConX*> (); //initialize vector
        vector<NetConX*> & vecNetCons = local->netcons[netConsPreId[nc]];
        vecNetCons.push_back(new NetConX(netcons[nc].mechType, netcons[nc].mechInstance, netcons[nc].delay,
                                         &netConsArgs[argsOffset], netcons[nc].argsCount, netcons[nc].active));
        //(int mechType, int mechInstance, double delay, double * args, short int argsCount, bool active);
        argsOffset += netcons[nc].argsCount;
    }

    neurox_hpx_unpin;
}

hpx_action_t Branch::clear=0;
int Branch::clear_handler()
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::callModFunction);
    //hpx_lco_wait(-local->mechsGraphLco);
    hpx_lco_delete_sync(local->mechsGraph->graphLCO);
    local->mechsGraph->graphLCO=HPX_NULL;
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

void Branch::initEventsQueue()
{
    hpx_lco_sema_p(this->eventsQueueMutex);
    for (int v=0; v<this->vecplayCount; v++)
        eventsQueue.push(make_pair(this->vecplay[v]->getFirstInstant(),
                                   (Event*) this->vecplay[v]));
    hpx_lco_sema_v_sync(this->eventsQueueMutex);
}

hpx_action_t Branch::fixedPlayContinuous = 0;
int Branch::fixedPlayContinuous_handler()
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::callModFunction);
    for (int v=0; v<local->vecplayCount; v++)
        local->vecplay[v]->continuous(local->t);
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t Branch::deliverEvents = 0;
int Branch::deliverEvents_handler(const double *t_ptr, const size_t size)
{
    neurox_hpx_pin(Branch);
    double t = t_ptr ? *t_ptr : local->t;
    neurox_hpx_recursive_branch_async_call(Branch::deliverEvents, t_ptr, size);
    hpx_lco_sema_p(local->eventsQueueMutex);
    while (!local->eventsQueue.empty() &&
           local->eventsQueue.top().first <= t)
    {
        Event * e = local->eventsQueue.top().second;
        e->deliver(t, local);
        local->eventsQueue.pop();
    }
    hpx_lco_sema_v_sync(local->eventsQueueMutex);
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
        if (local->mechsGraph!=NULL) //parallel
        {
          //launch execution on top nodes of the branch
          for (int m=0; m<mechanismsCount;m++)
          {
            if (mechanisms[m]->type == capacitance) continue; //not capacitance
            if (mechanisms[m]->dependenciesCount > 0) continue; //not top branch
            hpx_lco_set(local->mechsGraph->mechsLCOs[m], functionId_size, functionId_ptr, HPX_NULL, HPX_NULL);
          }
          //wait for the completion of the graph by waiting at the 'end node' lco
          hpx_lco_wait_reset(local->mechsGraph->endLCO);
        }
        else //serial
        {
            for (int m=0; m<mechanismsCount; m++)
                if ( mechanisms[m]->type == capacitance
                   && (*functionId_ptr == Mechanism::ModFunction::current
                      || *functionId_ptr == Mechanism::ModFunction::jacob))
                    continue;
                else
                    mechanisms[m]->callModFunction(local, *functionId_ptr);
        }
    }
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t Branch::MechanismsExecutionGraph::nodeFunction = 0;
int Branch::MechanismsExecutionGraph::nodeFunction_handler(const int * mechType_ptr, const size_t)
{
    neurox_hpx_pin(Branch);
    int type = *mechType_ptr;
    assert(type!=capacitance); //capacitance should be outside mechanisms graph
    assert(local->mechsGraph->mechsLCOs[mechanismsMap[type]] != HPX_NULL);
    Mechanism * mech = getMechanismFromType(type);

    Mechanism::ModFunction functionId;
    while (local->mechsGraph->graphLCO != HPX_NULL)
    {
        //wait until all dependencies have completed, and get the argument (function id) from the hpx_lco_set
        hpx_lco_get_reset(local->mechsGraph->mechsLCOs[mechanismsMap[type]], sizeof(Mechanism::ModFunction), &functionId);
        assert(functionId!=Mechanism::ModFunction::jacobCapacitance);
        assert(functionId!=Mechanism::ModFunction::currentCapacitance);
        mech->callModFunction(local, functionId);

        if (mech->successorsCount==0) //bottom mechanism
            hpx_lco_set(local->mechsGraph->endLCO, 0, NULL, HPX_NULL, HPX_NULL);
        else
            for (int c=0; c<mech->successorsCount; c++)
              hpx_lco_set(local->mechsGraph->mechsLCOs[mechanismsMap[mech->successors[c]]],
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
    hpx_lco_sema_p(local->eventsQueueMutex);
    for (auto nc : local->netcons.at(*preNeuronId))
        local->eventsQueue.push( make_pair(*deliveryTime, (Event*) &nc) );
    hpx_lco_sema_v_sync(local->eventsQueueMutex);
    neurox_hpx_unpin;
}

void Branch::registerHpxActions()
{
    neurox_hpx_register_action(1, Branch::deliverEvents);
    neurox_hpx_register_action(1, Branch::updateV);
    neurox_hpx_register_action(1, Branch::callModFunction);
    neurox_hpx_register_action(2, Branch::init);
    neurox_hpx_register_action(2, Branch::broadcastNetCons);
    neurox_hpx_register_action(2, Branch::queueSpikes);
    neurox_hpx_register_action(0, Branch::getSomaVoltage);
    neurox_hpx_register_action(0, Branch::clear);
    neurox_hpx_register_action(0, Branch::fixedPlayContinuous);
    neurox_hpx_register_action(1, Branch::MechanismsExecutionGraph::nodeFunction);
}
