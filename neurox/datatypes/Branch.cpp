#include "neurox/neurox.h"
#include <cstring>
#include <algorithm>
#include <numeric>
#include <set>

using namespace neurox;
using namespace neurox::Solver;

void* Branch::operator new(size_t bytes, void* addr) {
  return addr;
}

void Branch::operator delete(void* worker) {
  free(worker);
}

Branch::Branch(offset_t n,
               hpx_t branchHpxAddr,
               floble_t * data, size_t dataCount,
               offset_t *pdata, size_t pdataCount,
               offset_t * instancesCount, size_t recvMechanismsCount,
               offset_t * nodesIndices, size_t nodesIndicesCount,
               hpx_t * branches, size_t branchesCount,
               offset_t * p, size_t pCount,
               floble_t * vecplayT, size_t vecplayTCount,
               floble_t * vecplayY, size_t vecplayYCount,
               PointProcInfo * ppis, size_t vecplayCount,
               NetConX * netcons, size_t netconsCount,
               neuron_id_t * netConsPreId, size_t netConsPreIdsCount,
               floble_t *netConsWeights, size_t netConsWeightsCount,
               void** vdata, size_t vdataCount):
    soma(nullptr)
{
    NrnThread & nt = this->nt;

    //all non usable values
    nt.tml = NULL;
    nt._ml_list = NULL;
    nt.pntprocs = NULL;
    nt.presyns = NULL;
    nt.presyns_helper = NULL;
    nt.pnt2presyn_ix = NULL;
    nt.netcons = NULL;
    nt.n_pntproc=-1;
    nt.n_presyn=-1;
    nt.n_input_presyn = -1;
    nt.n_netcon = -1;
    nt.ncell=-1;
    nt.id = this->soma ? this->soma->id : -1;
    nt._stop_stepping = -1;
    nt._permute = NULL;
    nt._sp13mat = NULL;
    nt._ecell_memb_list = NULL;
    nt._ctime = -1;
    for (int i=0; i<BEFORE_AFTER_SIZE; i++)
      nt.tbl[i] = NULL;
    nt.shadow_rhs_cnt=-1;
    nt.compute_gpu=0;
    nt._net_send_buffer_size=-1;
    nt._net_send_buffer_cnt=-1;
    nt._net_send_buffer = NULL;
    nt.mapping = NULL;
    nt._idata = NULL;
    nt._nidata = -1;

    //assignemnts start here
    nt._t = inputParams->t;
    nt._dt = inputParams->dt;
    nt.cj = (inputParams->secondorder ? 2.0 : 1.0 ) / inputParams->dt;
    nt.end = n;

    nt._data = dataCount==0 ? nullptr : new floble_t[dataCount];
    memcpy(nt._data, data, dataCount*sizeof(floble_t));
    nt._ndata = dataCount;

    nt._vdata = vdataCount==0 ? nullptr :  new void*[vdataCount];
    memcpy(nt._vdata, vdata, vdataCount*sizeof(void*));
    nt._nvdata = vdataCount;

    nt.weights = netConsWeightsCount==0 ? nullptr : new floble_t[netConsWeightsCount];
    memcpy(nt.weights, netConsWeights, sizeof(floble_t)*netConsWeightsCount);
    nt.n_weight = netConsWeightsCount;

    nt._actual_rhs  = nt._data + n*0;
    nt._actual_d    = nt._data + n*1;
    nt._actual_a    = nt._data + n*2;
    nt._actual_b    = nt._data + n*3;
    nt._actual_v    = nt._data + n*4;
    nt._actual_area = nt._data + n*5;

    this->eventsQueueMutex = hpx_lco_sema_new(1);
    this->initMechanismsGraph(branchHpxAddr);

    //parent index
    if (pCount>0)
    {
        assert(pCount==n);
        nt._v_parent_index = new offset_t[pCount];
        memcpy(nt._v_parent_index, p, n*sizeof(offset_t));
    }
    else
    {
        this->nt._v_parent_index = nullptr;
    }

    //vecplay
    nt.n_vecplay = vecplayCount;
    nt._vecplay = vecplayCount == 0 ? nullptr : new void*[vecplayCount];

    offset_t vOffset=0;
    for (size_t v=0; v < nt.n_vecplay; v++)
    {
        int m = mechanismsMap[ppis[v].mechType];
        floble_t *instancesData = this->mechsInstances[m].data;
        floble_t *pd = &(instancesData[ppis->mechInstance*mechanisms[m]->dataSize+ppis->instanceDataOffset]);
        nt._vecplay[v] = new VecPlayContinuouX(pd, &vecplayT[vOffset], &vecplayY[vOffset], ppis[v].size);
        vOffset += ppis[v].size;
    }

    // reconstruct mechanisms
    assert (recvMechanismsCount == mechanismsCount);
    offset_t dataOffset=6*n;
    offset_t pdataOffset=0;
    offset_t vdataOffset=0;
    offset_t instancesOffset=0;
    this->mechsInstances = new Memb_list[mechanismsCount];

    for (offset_t m=0; m<mechanismsCount; m++)
    {
        Memb_list & instance = this->mechsInstances[m];
        Mechanism * mech = mechanisms[m];
        instance.nodecount = instancesCount[m];

        //data, pdata, and nodesIndices arrays
        instance.data = mech->dataSize ==0 ? nullptr : this->nt._data+dataOffset;
        instance.pdata = mech->pdataSize==0 ? nullptr : new offset_t[mech->pdataSize * instance.nodecount];
        memcpy(instance.pdata, &pdata[pdataOffset], sizeof(offset_t)*(mech->pdataSize * instance.nodecount));
        instance.nodeindices = instance.nodecount>0 ? new offset_t[instance.nodecount] : nullptr;
        memcpy(instance.nodeindices, &nodesIndices[instancesOffset], sizeof(offset_t)*instance.nodecount);

        //vdata: if is point process we need to allocate the vdata (by calling bbcore_reg in mod file)
        //and assign the correct offset in pdata (offset of vdata is in pdata[1])
        for (size_t i=0; i<instance.nodecount; i++)
        {
            //point pdata to the correct offset, and allocate vdata
            assert(dataOffset  <= dataCount);
            assert(vdataOffset <= vdataCount);
            floble_t * instanceData  = (floble_t*) &this->nt._data[dataOffset ];
            floble_t * instanceData2 = (floble_t*) &instance.data [i*mech->dataSize];
            offset_t *  instancePdata = (offset_t *) &instance.pdata[i*mech->pdataSize];
            assert (instanceData = instanceData2); //Make sure data offsets are good so far
            if (mech->pntMap>0)
            {
                assert( (mech->type == IClamp && mech->vdataSize == 1 && mech->pdataSize == 2)
                    || ((mech->type == ProbAMPANMDA_EMS || mech->type == ProbGABAAB_EMS)
                        && mech->vdataSize == 2 && mech->pdataSize == 3));

                void** instanceVdata = (void**) &this->nt._vdata[vdataOffset];
                Point_process * pp = (Point_process *) instanceVdata[0];
                (void) pp;
                (void) instanceVdata;
                if (mech->vdataSize==2)
                {
                    offset_t offsetPP  = instancePdata[1];
                    offset_t offsetRNG = instancePdata[2];
                    assert(offsetPP  < vdataCount);
                    assert(offsetRNG < vdataCount);
                    void * rng = instanceVdata[1];
                }
                //We will call bbcore_red on the correct vdata offset
                //this->vdata[vdataOffset];
            }
            else
            { assert (mech->vdataSize==0);}
            dataOffset  += mech->dataSize;
            pdataOffset += mech->pdataSize;
            vdataOffset += mech->vdataSize;
            assert(dataOffset  < 2^sizeof(offset_t));
            assert(pdataOffset < 2^sizeof(offset_t));
            assert(vdataOffset < 2^sizeof(offset_t));
            instancesOffset++;
        }
    }
    assert( dataOffset == dataCount);
    assert(pdataOffset == pdataCount);
    assert(vdataOffset == vdataCount);
    assert(instancesOffset == nodesIndicesCount);

    //Shadow arrays
    int shadowElemsCount = std::max(
                mechsInstances[mechanismsMap[ProbAMPANMDA_EMS]].nodecount,
                mechsInstances[mechanismsMap[ProbGABAAB_EMS]].nodecount
            );
    this->nt._shadow_rhs = new double[shadowElemsCount];
    this->nt._shadow_d   = new double[shadowElemsCount];

    //reconstructs parents or branching
    if (inputParams->multiSplix)
    {
      this->neuronTree = new Branch::NeuronTree;
      if (branchesCount>0)
      {
        this->neuronTree->branchesCount = branchesCount;
        this->neuronTree->branches = new hpx_t[branchesCount];
        memcpy(this->neuronTree->branches, branches, branchesCount*sizeof(hpx_t));
      }
      else
      {
        this->neuronTree->branchesCount = 0;
        this->neuronTree->branches = nullptr;
      }
    }
    else
    {
        this->neuronTree = nullptr;
    }

    //reconstructs netcons
    offset_t argsOffset=0;
    for (offset_t nc=0; nc<netconsCount; nc++)
    {
        this->netcons[ netConsPreId[nc] ].push_back(
                    new NetConX(netcons[nc].mechType, netcons[nc].mechInstance, netcons[nc].delay,
                    &netConsWeights[argsOffset], netcons[nc].argsCount, netcons[nc].active));
        argsOffset += netcons[nc].argsCount;
    }
}

Branch::~Branch()
{
    delete [] nt._v_parent_index;
    delete [] nt._actual_area;
    delete [] this->nt._data;

    if (neuronTree)
    {
      delete [] neuronTree->branches;
      delete [] neuronTree->branchesLCOs;
    }

    for (int m=0; m<mechanismsCount; m++)
    {
        delete [] mechsInstances[m].nodeindices;
        delete [] mechsInstances[m].pdata;
    }
    delete [] mechsInstances;
}

hpx_action_t Branch::init = 0;
int Branch::init_handler( const int nargs, const void *args[],
                          const size_t sizes[])
{
    neurox_hpx_pin(Branch);
    assert(nargs==14 || nargs==15);
    new(local) Branch(
        *(offset_t*) args[0], //number of compartments
        target, //current branch HPX addres
        (floble_t*) args[2], sizes[2]/sizeof(floble_t), //data (RHS, D, A, V, B, area, and mechs...)
        (offset_t*) args[3], sizes[3]/sizeof(offset_t), //pdata
        (offset_t*) args[4], sizes[4]/sizeof(offset_t), //instances count per mechanism
        (offset_t*) args[5], sizes[5]/sizeof(offset_t), //nodes indices
        (hpx_t*) args[6], sizes[6]/sizeof(hpx_t), //branches
        (offset_t*) args[7], sizes[7]/sizeof(offset_t), //parent index
        (floble_t*) args[8], sizes[8]/sizeof(floble_t), //vecplay T data
        (floble_t*) args[9], sizes[8]/sizeof(floble_t), //vecplay Y data
        (PointProcInfo*) args[10], sizes[10]/sizeof(PointProcInfo), //point processes info
        (NetConX*) args[11], sizes[11]/sizeof(NetConX), //netcons
        (neuron_id_t *) args[12], sizes[12]/sizeof(neuron_id_t), //netcons preneuron ids
        (floble_t *) args[13], sizes[13]/sizeof(floble_t), //netcons args
        (void**) args[14], sizes[14]/sizeof(void*));
    neurox_hpx_unpin;
}

void Branch::initMechanismsGraph(hpx_t target)
{
    //initializes mechanisms graphs (capacitance is excluded from graph)
    if (inputParams->multiMex)
    {
      this->mechsGraph = new MechanismsGraphLCO;
      this->mechsGraph->graphLCO = hpx_lco_and_new(mechanismsCount-1); //excludes 'capacitance'
      this->mechsGraph->mechsLCOs = new hpx_t[mechanismsCount];
      this->mechsGraph->mechsLCOs[mechanismsMap[CAP]] = HPX_NULL;
      size_t terminalMechanismsCount=0;
      for (size_t m=0; m<mechanismsCount; m++)
      {
        if (mechanisms[m]->type == CAP) continue; //exclude capacitance
        this->mechsGraph->mechsLCOs[m] = hpx_lco_reduce_new(
                    max((short) 1,mechanisms[m]->dependenciesCount),
                    sizeof(Mechanism::ModFunction), Mechanism::initModFunction,
                    Mechanism::reduceModFunction);
        if (mechanisms[m]->successorsCount==0) //bottom of mechs graph
            terminalMechanismsCount++;
        hpx_call(target,  Branch::MechanismsGraphLCO::nodeFunction,
                 this->mechsGraph->graphLCO, &mechanisms[m]->type, sizeof(int));
      }
      this->mechsGraph->endLCO = hpx_lco_and_new(terminalMechanismsCount);
    }
    else
    {
      this->mechsGraph = NULL;
    }
}

hpx_action_t Branch::initSoma = 0;
int Branch::initSoma_handler(const int nargs,
                             const void *args[], const size_t[])
{
    neurox_hpx_pin(Branch);
    assert(nargs==2);
    const neuron_id_t neuronId = *(const neuron_id_t*)    args[0];
    const floble_t APthreshold = *(const floble_t*) args[1];
    local->soma=new Neuron(neuronId, APthreshold);
    neurox_hpx_unpin;
}

hpx_action_t Branch::clear=0;
int Branch::clear_handler()
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::clear);
    //hpx_lco_wait(-local->mechsGraphLco);
    hpx_lco_delete_sync(local->mechsGraph->graphLCO);
    local->mechsGraph->graphLCO=HPX_NULL;
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

void Branch::initEventsQueue()
{
    hpx_lco_sema_p(this->eventsQueueMutex);
    for (size_t v=0; v<this->nt.n_vecplay; v++)
    {
        VecPlayContinuouX * vecplay = reinterpret_cast<VecPlayContinuouX*>(this->nt._vecplay[v]);
        eventsQueue.push(make_pair(vecplay->getFirstInstant(),
                                   (Event*) vecplay));
    }
    hpx_lco_sema_v_sync(this->eventsQueueMutex);
}

void Branch::callModFunction(const Mechanism::ModFunction functionId)
{
    if (functionId<BEFORE_AFTER_SIZE) return; //N/A

    //only for capacitance mechanism
    if (functionId == Mechanism::ModFunction::currentCapacitance
     || functionId == Mechanism::ModFunction::jacobCapacitance)
    {
      mechanisms[mechanismsMap[CAP]]->callModFunction(this, functionId);
    }
    //for all others except capacitance (mechanisms graph)
    else
    {
        if (this->mechsGraph!=NULL) //parallel
        {
          //launch execution on top nodes of the branch
          for (int m=0; m<mechanismsCount;m++)
          {
            if (mechanisms[m]->type == CAP)   continue; //not capacitance
            if (mechanisms[m]->dependenciesCount > 0) continue; //not a top branche
            hpx_lco_set(this->mechsGraph->mechsLCOs[m],
                        sizeof(functionId), &functionId, HPX_NULL, HPX_NULL);
          }
          //wait for the completion of the graph by waiting at 'end node' lco
          hpx_lco_wait_reset(this->mechsGraph->endLCO);
        }
        else //serial
        {
            for (int m=0; m<mechanismsCount; m++)
                if ( mechanisms[m]->type == CAP
                   && (  functionId == Mechanism::ModFunction::current
                      || functionId == Mechanism::ModFunction::jacob))
                    continue;
                else
                    mechanisms[m]->callModFunction(this, functionId);
        }
    }
}

hpx_action_t Branch::MechanismsGraphLCO::nodeFunction = 0;
int Branch::MechanismsGraphLCO::nodeFunction_handler(const int * mechType_ptr,
                                                     const size_t)
{
    neurox_hpx_pin(Branch);
    int type = *mechType_ptr;
    assert(type!=CAP); //capacitance should be outside mechanisms graph
    assert(local->mechsGraph->mechsLCOs[mechanismsMap[type]] != HPX_NULL);
    Mechanism * mech = getMechanismFromType(type);

    Mechanism::ModFunction functionId;
    while (local->mechsGraph->graphLCO != HPX_NULL)
    {
      //wait until all dependencies have completed, and get the argument
      //(function id) from the hpx_lco_set
      hpx_lco_get_reset(
              local->mechsGraph->mechsLCOs[mechanismsMap[type]],
              sizeof(Mechanism::ModFunction), &functionId);
      assert(functionId!=Mechanism::ModFunction::jacobCapacitance);
      assert(functionId!=Mechanism::ModFunction::currentCapacitance);
      mech->callModFunction(local, functionId);

      if (mech->successorsCount==0) //bottom mechanism
        hpx_lco_set(local->mechsGraph->endLCO, 0, NULL, HPX_NULL, HPX_NULL);
      else
        for (int c=0; c<mech->successorsCount; c++)
          hpx_lco_set(
              local->mechsGraph->mechsLCOs[mechanismsMap[mech->successors[c]]],
              sizeof(functionId), &functionId, HPX_NULL, HPX_NULL);
    }
    neurox_hpx_unpin;
}

//from eion.c
#define nparm 5
#define cur	3
#define dcurdv 4

void Branch::secondOrderCurrent()
{
    floble_t *& rhs = this->nt._actual_rhs;

    Memb_list * mechInstances= this->mechsInstances;
    for (int m=0; m<mechanismsCount; m++)
    {
        Mechanism * mech = mechanisms[m];
        if (!mech->isIon) continue;
        assert(nparm==mech->dataSize); //see '#define nparm 5' in eion.c

        for (int i=0; i<mechInstances[m].nodecount; i++)
        {
            floble_t * data = &mechInstances[m].data[i*mech->dataSize];
            offset_t & nodeIndex = mechInstances[m].nodeindices[i];
            data[cur] += data[dcurdv] * rhs[nodeIndex];
        }
    }
}

hpx_action_t Branch::addSpikeEvent = 0;
int Branch::addSpikeEvent_handler(
        const int nargs, const void *args[], const size_t[] )
{
    neurox_hpx_pin(Branch);
    assert(nargs==2);
    const int * preNeuronId = (const int *) args[0];
    const floble_t * deliveryTime = (const floble_t*) args[1];

    //netcvode::PreSyn::send()
    hpx_lco_sema_p(local->eventsQueueMutex);
    for (auto nc : local->netcons.at(*preNeuronId))
        local->eventsQueue.push( make_pair(*deliveryTime, (Event*) nc) );
    hpx_lco_sema_v_sync(local->eventsQueueMutex);
    neurox_hpx_unpin;
}

hpx_action_t Branch::initNeuronTreeLCO = 0;
int Branch::initNeuronTreeLCO_handler()
{
    neurox_hpx_pin(Branch);
    if (!inputParams->multiSplix)
    {
      local->neuronTree = nullptr;
    }
    else
    {
      assert(local->neuronTree);
      offset_t branchesCount = local->neuronTree->branchesCount;
      for (int i=0; i<NeuronTree::futuresSize; i++)
          local->neuronTree->localLCO[i] =
                  local->soma ? HPX_NULL : hpx_lco_future_new(sizeof(floble_t));
      local->neuronTree->branchesLCOs = branchesCount ?
                  new hpx_t[branchesCount][NeuronTree::futuresSize] : nullptr;

      //send my LCOs to children, and receive theirs
      if (branchesCount>0)
      {
        hpx_t * futures = branchesCount ? new hpx_t[branchesCount]  : nullptr;
        void ** addrs   = branchesCount ? new void*[branchesCount]  : nullptr;
        size_t* sizes   = branchesCount ? new size_t[branchesCount] : nullptr;
        for (offset_t c = 0; c < branchesCount; c++)
        {
          futures[c] = hpx_lco_future_new(sizeof (hpx_t)*NeuronTree::futuresSize);
          addrs[c]   = &local->neuronTree->branchesLCOs[c];
          sizes[c]   = sizeof(hpx_t)*NeuronTree::futuresSize;
          hpx_call(local->neuronTree->branches[c], Branch::initNeuronTreeLCO,
                   futures[c], local->neuronTree->localLCO,
                   sizeof(hpx_t)*NeuronTree::futuresSize); //pass my LCO down
        }
        hpx_lco_get_all(branchesCount, futures, sizes, addrs, NULL);
        hpx_lco_delete_all(branchesCount, futures, NULL);

        delete [] futures;
        delete [] addrs;
        delete [] sizes;
      }

      if (!local->soma) //send my LCO to parent
          neurox_hpx_unpin_continue(local->neuronTree->localLCO);
    }
    neurox_hpx_unpin;
}

void Branch::initialize()
{
    floble_t *& v = this->nt._actual_v;
    double t = this->nt._t;

    //set up by finitialize.c:nrn_finitialize(): if (setv)
    assert(inputParams->secondorder < sizeof(char));
    initEventsQueue();
    deliverEvents(t);

    //set up by finitialize.c:nrn_finitialize(): if (setv)
    for (int n=0; n<this->nt.end; n++)
        v[n]=inputParams->voltage;

    // the INITIAL blocks are ordered so that mechanisms that write
    // concentrations are after ions and before mechanisms that read
    // concentrations.
    callModFunction(Mechanism::ModFunction::before_initialize);
    callModFunction(Mechanism::ModFunction::initialize);
    callModFunction(Mechanism::ModFunction::after_initialize);

    //set up by finitialize.c:nrn_finitialize() -> fadvance_core.c:dt2thread()
    //local->cj = inputParams->secondorder ? 2.0/inputParams->dt : 1.0/inputParams->dt;
    //done when calling mechanisms //TODO have a copy per branch to speed-up?
    deliverEvents(t);
    setupTreeMatrixMinimal();
    deliverEvents(t);

    //part of nrn_fixed_step_group_minimal
    //1. multicore.c::nrn_thread_table_check()
    callModFunction(Mechanism::ModFunction::threadTableCheck);
}

void Branch::backwardEulerStep()
{
    floble_t *& v = this->nt._actual_v;
    floble_t *& rhs = this->nt._actual_rhs;
    double & t = this->nt._t;
    double & dt = this->nt._dt;

    //2. multicore.c::nrn_fixed_step_thread()
    //2.1 cvodestb::deliver_net_events(nth);
    //(send outgoing spikes netcvode.cpp::NetCvode::check_thresh() )
    bool reachedThresold = soma && v[0] >= soma->APthreshold;
    if (reachedThresold) soma->fireActionPotential(t);
    t += .5*dt;
    deliverEvents(t);
    fixedPlayContinuous();
    setupTreeMatrixMinimal();

    //eion.c : second_order_cur()
    //TODO when we use only CoreNeuron structs, replace this by original func;
    if (inputParams->secondorder)
        secondOrderCurrent();

    //fadvance_core.c : update() / Branch::updateV
    for (int i=0; i<this->nt.end; i++)
        v[i] += (inputParams->secondorder ? 2 : 1) * rhs[i];

    callModFunction(Mechanism::ModFunction::jacob);
    t += .5*dt;
    fixedPlayContinuous();
    callModFunction(Mechanism::ModFunction::state);
    callModFunction(Mechanism::ModFunction::after_solve);
    deliverEvents(t);

    //if we are at the output time instant output to file
    if (fmod(t, inputParams->dt_io) == 0)
    {
        //output
    }
    //make sure all synapses from N steps before were delivered
    //(thus other neurons wait for this neuron one before stepping)
    //>waitForSynapsesDelivery(stepsCountPerComm);
}

hpx_action_t Branch::finitialize = 0;
int Branch::finitialize_handler()
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::finitialize);
    local->initialize();
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t Branch::backwardEuler = 0;
int Branch::backwardEuler_handler()
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::backwardEuler);
    int i=0;
    double & t = local->nt._t;
    while (t <= inputParams->tstop)
    {
        if (local->soma) //neuron
        { //docs: last argument should be 0 for and-lco
          //hpx_t step = hpx_lco_array_at(timeMachine, i, 0);
          //hpx_lco_set(step, 0, NULL, HPX_NULL, HPX_NULL);
        }

        local->backwardEulerStep();
        if (local->soma && i>=4)
        {
            //TODO use lco_gencount maybe?
            //hpx_t step = hpx_lco_array_at(timeMachine, i-4, 0);
            //hpx_lco_wait(step);
        }
        i++;
    }
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

void Branch::setupTreeMatrixMinimal()
{
    floble_t *& d = this->nt._actual_d;
    floble_t *& rhs = this->nt._actual_rhs;

    for (int i=0; i<this->nt.end; i++)
    {
        rhs[i]=0;
        d[i]=0;
    }
    this->callModFunction(Mechanism::ModFunction::before_breakpoint);
    this->callModFunction(Mechanism::ModFunction::current);

    Solver::HinesSolver::forwardTriangulation(this);

    // calculate left hand side of
    //cm*dvm/dt = -i(vm) + is(vi) + ai_j*(vi_j - vi)
    //cx*dvx/dt - cm*dvm/dt = -gx*(vx - ex) + i(vm) + ax_j*(vx_j - vx)
    //with a matrix so that the solution is of the form [dvm+dvx,dvx] on the right
    //hand side after solving.
    //This is a common operation for fixed step, cvode, and daspk methods
    this->callModFunction(Mechanism::ModFunction::jacob);

    //finitialize.c:nrn_finitialize()->set_tree_matrix_minimal->nrn_rhs (treeset_core.c)
    //now the cap current can be computed because any change to cm
    //by another model has taken effect.
    this->callModFunction(Mechanism::ModFunction::jacobCapacitance);

    Solver::HinesSolver::backSubstitution(this);
}

void Branch::deliverEvents(floble_t t)
{
    hpx_lco_sema_p(this->eventsQueueMutex);
    while (!this->eventsQueue.empty() &&
           this->eventsQueue.top().first <= t)
    {
        Event * e = this->eventsQueue.top().second;
        e->deliver(t, this);
        this->eventsQueue.pop();
    }
    hpx_lco_sema_v_sync(this->eventsQueueMutex);
}

void Branch::fixedPlayContinuous()
{
    double t = this->nt._t;
    for (int v=0; v<this->nt.n_vecplay; v++)
    {
        void * vecplay_void = this->nt._vecplay[v];
        VecPlayContinuouX * vecplay = reinterpret_cast<VecPlayContinuouX*>(vecplay_void);
        vecplay->continuous(t);
    }
}

void Branch::registerHpxActions()
{
    neurox_hpx_register_action(2, Branch::init);
    neurox_hpx_register_action(2, Branch::initSoma);
    neurox_hpx_register_action(0, Branch::clear);
    neurox_hpx_register_action(2, Branch::addSpikeEvent);
    neurox_hpx_register_action(0, Branch::finitialize);
    neurox_hpx_register_action(0, Branch::backwardEuler);
    neurox_hpx_register_action(0, Branch::initNeuronTreeLCO);
    neurox_hpx_register_action(1, Branch::MechanismsGraphLCO::nodeFunction);
}
