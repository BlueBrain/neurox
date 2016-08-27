#include "neurox/neurox.h"
#include <cstring>
#include <algorithm>
#include <numeric>
#include <set>

using namespace neurox;

void* Branch::operator new(size_t bytes, void* addr) {
  return addr;
}

void Branch::operator delete(void* worker) {
  free(worker);
}

Branch::Branch (int n,
                hpx_t branchHpxAddr,
                double * data, size_t dataCount,
                int *pdata, size_t pdataCount,
                int * instancesCount, size_t recvMechanismsCount,
                int * nodesIndices, size_t nodesIndicesCount,
                hpx_t * branches, size_t branchesCount,
                int * p, size_t pCount,
                double * vecplayTdata, size_t vecplayTCount,
                double * vecplayYdata, size_t vecplayYCount,
                PointProcInfo * ppis, size_t vecplayCount,
                NetConX * netcons, size_t netconsCount,
                int * netConsPreId, size_t netConsPreIdsCount,
                double *netConsArgs, size_t netConsArgsCount,
                void** vdata, size_t vdataCount):
    soma(nullptr), n(n)
{
    this->t = inputParams->t;
    this->dt = inputParams->dt;
    this->secondOrder = inputParams->secondorder;
    this->eventsQueueMutex = hpx_lco_sema_new(1);
    this->initMechanismsGraph(branchHpxAddr);

    this->data = dataCount==0 ? nullptr : new double[dataCount];
    memcpy(this->data, data, dataCount*sizeof(double));
    this->vdata = vdataCount==0 ? nullptr :  new void*[vdataCount];
    memcpy(this->vdata, vdata, vdataCount*sizeof(void*));
    this->rhs = this->data + this->n*0;
    this->d = this->data + this->n*1;
    this->a = this->data + this->n*2;
    this->b = this->data + this->n*3;
    this->v = this->data + this->n*4;
    this->area = this->data + this->n*5;

    // reconstruct mechanisms
    assert (recvMechanismsCount == mechanismsCount);
    int dataOffset=6*this->n;
    int pdataOffset=0;
    int vdataOffset=0;
    int instancesOffset=0;
    this->mechsInstances = new MechanismInstance[mechanismsCount];

    for (int m=0; m<mechanismsCount; m++)
    {
        MechanismInstance & instance = this->mechsInstances[m];
        Mechanism * mech = mechanisms[m];
        instance.count = instancesCount[m];

        //data, pdata, and nodesIndices arrays
        instance.data = mech->dataSize ==0 ? nullptr : this->data+dataOffset;
        instance.pdata = mech->pdataSize==0 ? nullptr : new int[mech->pdataSize * instance.count];
        memcpy(instance.pdata, &pdata[pdataOffset], sizeof(int)*(mech->pdataSize * instance.count));
        instance.nodesIndices = instance.count>0 ? new int[instance.count] : nullptr;
        memcpy(instance.nodesIndices, &nodesIndices[instancesOffset], sizeof(int)*instance.count);

        //vdata: if is point process we need to allocate the vdata (by calling bbcore_reg in mod file)
        //and assign the correct offset in pdata (offset of vdata is in pdata[1])
        for (int i=0; i<instance.count; i++)
        {
            //point pdata to the correct offset, and allocate vdata
            assert(dataOffset  <= dataCount);
            assert(vdataOffset <= vdataCount);
            double * instanceData  = (double*) &this->data[dataOffset ];
            double * instanceData2 = (double*) &instance.data [i*mech->dataSize];
            int *    instancePdata = (int   *) &instance.pdata[i*mech->pdataSize];
            assert (instanceData = instanceData2); //Make sure data offsets are good so far
            if (mech->pntMap>0)
            {
                assert( (mech->type == IClamp && mech->vdataSize == 1 && mech->pdataSize == 2)
                    || ((mech->type == ProbAMPANMDA_EMS || mech->type == ProbGABAAB_EMS)
                        && mech->vdataSize == 2 && mech->pdataSize == 3));

                void** instanceVdata = (void**) &this->vdata[vdataOffset];
                Point_process * pp = (Point_process *) instanceVdata[0];
                (void) pp;
                (void) instanceVdata;
                if (mech->vdataSize==2)
                {
                    int offsetPP  = instancePdata[1];
                    int offsetRNG = instancePdata[2];
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
            instancesOffset++;
        }
    }
    assert( dataOffset == dataCount);
    assert(pdataOffset == pdataCount);
    assert(vdataOffset == vdataCount);
    assert(instancesOffset == nodesIndicesCount);

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

    if (pCount>0)
    {
        assert(pCount==n);
        this->p = new int[pCount];
        memcpy(this->p, p, n*sizeof(int));
    }
    else
    {
        this->p = nullptr;
    }

    //reconstructs vecplay
    this->vecplayCount = vecplayCount;
    this->vecplay = this->vecplayCount == 0 ? nullptr : new VecPlayContinuouX*[this->vecplayCount];

    int vOffset=0;
    for (int v=0; v<this->vecplayCount; v++)
    {
        int m = mechanismsMap[ppis[v].mechType];
        double *instancesData = this->mechsInstances[m].data;
        double *pd = &(instancesData[ppis->mechInstance*mechanisms[m]->dataSize+ppis->instanceDataOffset]);
        this->vecplay[v] = new VecPlayContinuouX(pd, &vecplayTdata[vOffset], &vecplayYdata[vOffset], ppis[v].size);
        vOffset += ppis[v].size;
    }

    //reconstructs netcons
    int argsOffset=0;
    //this->netcons = std::map<int, vector<NetConX*> > ();
    for (int nc=0; nc<netconsCount; nc++)
    {
        this->netcons[ netConsPreId[nc] ].push_back(
                    new NetConX(netcons[nc].mechType, netcons[nc].mechInstance, netcons[nc].delay,
                    &netConsArgs[argsOffset], netcons[nc].argsCount, netcons[nc].active));
        //(int mechType, int mechInstance, double delay, double * args, short int argsCount, bool active);
        argsOffset += netcons[nc].argsCount;
    }
}

Branch::~Branch()
{
    delete [] b;
    delete [] d;
    delete [] a;
    delete [] v;
    delete [] p;
    delete [] rhs;
    delete [] area;
    delete [] data;

    if (neuronTree)
    {
      delete [] neuronTree->branches;
      delete [] neuronTree->branchesLCOs;
    }

    for (int m=0; m<mechanismsCount; m++)
    {
        delete [] mechsInstances[m].nodesIndices;
        delete [] mechsInstances[m].pdata;
    }
    delete [] mechsInstances;
}

hpx_action_t Branch::init = 0;
int Branch::init_handler(const int nargs, const void *args[], const size_t sizes[])
{
    neurox_hpx_pin(Branch);
    assert(nargs==14 || nargs==15);
    new(local) Branch(
        *(int*) args[0], //number of compartments
        target, //current branch HPX addres
        (double*) args[2], sizes[2]/sizeof(double), //data (RHS, D, A, V, B, area, and mechs...)
        (int*) args[3], sizes[3]/sizeof(int), //pdata
        (int*) args[4], sizes[4]/sizeof(int), //instances count per mechanism
        (int*) args[5], sizes[5]/sizeof(int), //nodes indices
        (hpx_t*) args[6], sizes[6]/sizeof(hpx_t), //branches
        (int*) args[7], sizes[7]/sizeof(int), //parent index
        (double*) args[8], sizes[8]/sizeof(double), //vecplay T data
        (double*) args[9], sizes[8]/sizeof(double), //vecplay Y data
        (PointProcInfo*) args[10], sizes[10]/sizeof(PointProcInfo), //point processes info
        (NetConX*) args[11], sizes[11]/sizeof(NetConX), //netcons
        (int *) args[12], sizes[12]/sizeof(int), //netcons preneuron ids
        (double *) args[13], sizes[13]/sizeof(double), //netcons args
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
      this->mechsGraph->mechsLCOs[mechanismsMap[capacitance]] = HPX_NULL;
      int terminalMechanismsCount=0;
      for (int m=0; m<mechanismsCount; m++)
      {
        if (mechanisms[m]->type == capacitance) continue; //exclude capacitance
        this->mechsGraph->mechsLCOs[m] = hpx_lco_reduce_new(max((short) 1,mechanisms[m]->dependenciesCount),
                      sizeof(Mechanism::ModFunction), Mechanism::initModFunction, Mechanism::reduceModFunction);
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
int Branch::initSoma_handler(const int nargs, const void *args[], const size_t[])
{
    neurox_hpx_pin(Branch);
    assert(nargs==2);
    const int    neuronId    = *(const int*)    args[0];
    const double APthreshold = *(const double*) args[1];
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
    for (int v=0; v<this->vecplayCount; v++)
        eventsQueue.push(make_pair(this->vecplay[v]->getFirstInstant(),
                                   (Event*) this->vecplay[v]));
    hpx_lco_sema_v_sync(this->eventsQueueMutex);
}

void Branch::callModFunction(const Mechanism::ModFunction functionId)
{
    if (functionId<BEFORE_AFTER_SIZE) return; //N/A

    //only for capacitance mechanism
    if (functionId == Mechanism::ModFunction::currentCapacitance
     || functionId == Mechanism::ModFunction::jacobCapacitance)
    {
        mechanisms[mechanismsMap[capacitance]]->callModFunction(this, functionId);
    }
    //for all others except capacitance (mechanisms graph)
    else
    {
        if (this->mechsGraph!=NULL) //parallel
        {
          //launch execution on top nodes of the branch
          for (int m=0; m<mechanismsCount;m++)
          {
            if (mechanisms[m]->type == capacitance)   continue; //not capacitance
            if (mechanisms[m]->dependenciesCount > 0) continue; //not a top branche
            hpx_lco_set(this->mechsGraph->mechsLCOs[m], sizeof(functionId), &functionId, HPX_NULL, HPX_NULL);
          }
          //wait for the completion of the graph by waiting at the 'end node' lco
          hpx_lco_wait_reset(this->mechsGraph->endLCO);
        }
        else //serial
        {
            for (int m=0; m<mechanismsCount; m++)
                if ( mechanisms[m]->type == capacitance
                   && (  functionId == Mechanism::ModFunction::current
                      || functionId == Mechanism::ModFunction::jacob))
                    continue;
                else
                    mechanisms[m]->callModFunction(this, functionId);
        }
    }
}

hpx_action_t Branch::MechanismsGraphLCO::nodeFunction = 0;
int Branch::MechanismsGraphLCO::nodeFunction_handler(const int * mechType_ptr, const size_t)
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

void Branch::secondOrderCurrent()
{
    MechanismInstance * mechInstances= this->mechsInstances;
    for (int m=0; m<mechanismsCount; m++)
    {
        Mechanism * mech = mechanisms[m];
        if (!mech->isIon) continue;
        assert(nparm==mech->dataSize); //see '#define nparm 5' in eion.c

        for (int i=0; i<mechInstances[m].count; i++)
        {
            double * data = &mechInstances[m].data[i*mech->dataSize];
            int & nodeIndex = mechInstances[m].nodesIndices[i];
            data[cur] += data[dcurdv] * this->rhs[nodeIndex];
        }
    }
}

hpx_action_t Branch::addSpikeEvent = 0;
int Branch::addSpikeEvent_handler(const int nargs, const void *args[], const size_t [])
{
    neurox_hpx_pin(Branch);
    assert(nargs==2);
    const int * preNeuronId = (const int *) args[0];
    const double * deliveryTime = (const double*) args[1];

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
      int branchesCount = local->neuronTree->branchesCount;
      local->neuronTree->localLCO[0] = local->soma ? HPX_NULL : hpx_lco_future_new(sizeof(double));
      local->neuronTree->localLCO[1] = local->soma ? HPX_NULL : hpx_lco_future_new(sizeof(double));
      local->neuronTree->branchesLCOs = branchesCount ? new hpx_t[branchesCount*2] : nullptr;

      //send my LCOs to children, and receive theirs
      if (branchesCount>0)
      {
        hpx_t * futures = branchesCount ? new hpx_t[branchesCount]  : nullptr;
        void ** addrs   = branchesCount ? new void*[branchesCount]  : nullptr;
        size_t* sizes   = branchesCount ? new size_t[branchesCount] : nullptr;
        for (int c = 0; c < branchesCount; c++)
        {
          futures[c] = hpx_lco_future_new(sizeof (hpx_t)*2);
          addrs[c]   = &local->neuronTree->branchesLCOs[c*2];
          sizes[c]   = sizeof(hpx_t)*2;
          hpx_call(local->neuronTree->branches[c], Branch::initNeuronTreeLCO, futures[c],
                  local->neuronTree->localLCO, sizeof(hpx_t)*2); //pass my LCO down
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
    //set up by finitialize.c:nrn_finitialize(): if (setv)
    assert(inputParams->secondorder < sizeof(char));
    this->secondOrder = inputParams->secondorder;
    initEventsQueue();
    deliverEvents(this->t);

    //set up by finitialize.c:nrn_finitialize(): if (setv)
    for (int n=0; n<this->n; n++)
        this->v[n]=inputParams->voltage;

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
    if (this->secondOrder == 2)
        secondOrderCurrent();

    //fadvance_core.c : update() / Branch::updateV
    for (int i=0; i<n; i++)
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
    neurox_hpx_recursive_branch_async_call(Branch::backwardEuler);
    local->initialize();
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t Branch::backwardEuler = 0;
int Branch::backwardEuler_handler()
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::backwardEuler);
    while (local->t <= inputParams->tstop)
        local->backwardEulerStep();
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

void Branch::setupTreeMatrixMinimal()
{
    for (int i=0; i<this->n; i++)
    {
        this->rhs[i]=0;
        this->d[i]=0;
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

void Branch::deliverEvents(double t)
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
    for (int v=0; v<this->vecplayCount; v++)
        this->vecplay[v]->continuous(this->t);
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
