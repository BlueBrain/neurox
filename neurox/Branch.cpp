#include "neurox/neurox.h"
#include <cstring>
#include <algorithm>
#include <numeric>
#include <set>

#include "coreneuron/utils/randoms/nrnran123.h" //TODO

using namespace neurox;
using namespace neurox::Solver;
using namespace neurox::Tools;

void* Branch::operator new(size_t bytes, void* addr) {
  return addr;
}

void Branch::operator delete(void* worker) {}

Branch::Branch(offset_t n,
               int nrnThreadId,
               int thresholdVoffset,
               hpx_t branchHpxAddr,
               floble_t * data, size_t dataCount,
               offset_t *pdata, size_t pdataCount,
               offset_t * instancesCount, size_t recvMechanismsCount,
               offset_t * nodesIndices, size_t nodesIndicesCount,
               hpx_t topBranchAddr,
               hpx_t * branches, size_t branchesCount,
               offset_t * p, size_t pCount,
               floble_t * vecplayT, size_t vecplayTCount,
               floble_t * vecplayY, size_t vecplayYCount,
               PointProcInfo * vecplayPPI, size_t vecplayPPIcount,
               NetConX * netcons, size_t netconsCount,
               neuron_id_t * netConsPreId, size_t netConsPreIdsCount,
               floble_t *weights, size_t weightsCount,
               unsigned char* vdataSerialized, size_t vdataSerializedCount):
    soma(nullptr),nt(nullptr), thvar_ptr(nullptr)
{
    this->nt = (NrnThread*) malloc(sizeof(NrnThread));
    NrnThread * nt = this->nt;

    //all non usable values
    nt->_ml_list = NULL;
    nt->tml = NULL;
    nt->pntprocs = NULL;
    nt->presyns = NULL;
    nt->presyns_helper = NULL;
    nt->pnt2presyn_ix = NULL;
    nt->netcons = NULL;
    nt->n_pntproc=-1;
    nt->n_presyn=-1;
    nt->n_input_presyn = -1;
    nt->n_netcon = -1;
    nt->ncell=-1;
    nt->id = nrnThreadId;
    nt->_stop_stepping = -1;
    nt->_permute = NULL;
    nt->_sp13mat = NULL;
    nt->_ecell_memb_list = NULL;
    nt->_ctime = -1;
    for (int i=0; i<BEFORE_AFTER_SIZE; i++)
      nt->tbl[i] = NULL;
    nt->compute_gpu=0;
    nt->_net_send_buffer_size=-1;
    nt->_net_send_buffer_cnt=-1;
    nt->_net_send_buffer = NULL;
    nt->mapping = NULL;
    nt->_idata = NULL;
    nt->_nidata = -1;

    //assignemnts start here
    nt->_dt = inputParams->dt;
    nt->_t  = inputParams->tstart;
    nt->cj = (inputParams->secondorder ? 2.0 : 1.0 ) / inputParams->dt;
    nt->end = n;

    nt->_data = dataCount==0 ? nullptr : Vectorizer::new_<floble_t>(dataCount);
    memcpy(nt->_data, data, dataCount*sizeof(floble_t));
    nt->_ndata = dataCount;

    nt->weights = weightsCount==0 ? nullptr : new floble_t[weightsCount];
    memcpy(nt->weights, weights, sizeof(floble_t)*weightsCount);
    nt->n_weight = weightsCount;

    nt->_actual_rhs  = nt->_data + n*0;
    nt->_actual_d    = nt->_data + n*1;
    nt->_actual_a    = nt->_data + n*2;
    nt->_actual_b    = nt->_data + n*3;
    nt->_actual_v    = nt->_data + n*4;
    nt->_actual_area = nt->_data + n*5;

    //AP threshold offset
    this->thvar_ptr = thresholdVoffset==-1 ? NULL : &(this->nt->_actual_v[thresholdVoffset]);

    //events mutex
    this->eventsQueueMutex = hpx_lco_sema_new(1);

    //parent index
    if (pCount>0)
    {
        assert(pCount==n);
        this->nt->_v_parent_index = Vectorizer::new_<offset_t>(pCount);
        memcpy(this->nt->_v_parent_index, p, n*sizeof(offset_t));
    }
    else
    {
        this->nt->_v_parent_index = nullptr;
    }

    // reconstruct mechanisms
    assert (recvMechanismsCount == mechanismsCount);
    offset_t dataOffset=6*n;
    offset_t pdataOffset=0;
    offset_t instancesOffset=0;
    this->mechsInstances = new Memb_list[mechanismsCount];

    int maxMechId = 0;

    vector<void*> vdataPtrs;
    offset_t vdataOffset=0;

    for (offset_t m=0; m<mechanismsCount; m++)
    {
        Memb_list & instance = this->mechsInstances[m];
        Mechanism * mech = mechanisms[m];
        instance.nodecount = instancesCount[m];
        maxMechId = max(maxMechId, mech->type);

        //data, pdata, and nodesIndices arrays
        instance.data  = mech->dataSize ==0 || instance.nodecount==0 ? nullptr : this->nt->_data+dataOffset;
        instance.pdata = mech->pdataSize==0 || instance.nodecount==0 ? nullptr : Vectorizer::new_<offset_t>(mech->pdataSize * instance.nodecount);
        if (instance.pdata)
            memcpy(instance.pdata, &pdata[pdataOffset], sizeof(offset_t)*(mech->pdataSize * instance.nodecount));
        instance.nodeindices = instance.nodecount>0 ? Vectorizer::new_<offset_t>(instance.nodecount) : nullptr;
        if (instance.nodeindices)
            memcpy(instance.nodeindices, &nodesIndices[instancesOffset], sizeof(offset_t)*instance.nodecount);

        //init thread data :: nrn_setup.cpp->setup_ThreadData();
        if (mech->membFunc.thread_size_)
        {
            instance._thread = new ThreadDatum[mech->membFunc.thread_size_];
            if (mech->membFunc.thread_mem_init_)
                mech->membFunc.thread_mem_init_(instance._thread);
        }
        else
        {
            instance._thread = nullptr;
        }

        //vdata: if is point process we need to allocate the vdata (by calling bbcore_reg in mod file)
        //and assign the correct offset in pdata (offset of vdata is in pdata[1])
        for (size_t i=0; i<instance.nodecount; i++)
        {
            //point pdata to the correct offset, and allocate vdata
            assert(dataOffset  <= dataCount);
            assert(vdataOffset <= vdataSerializedCount);
            floble_t * instanceData  = (floble_t*) &this->nt->_data[dataOffset ];
            floble_t * instanceData2 = (floble_t*) &instance.data [i*mech->dataSize];
            offset_t * instancePdata = (offset_t *) &instance.pdata[i*mech->pdataSize];
            assert (instanceData = instanceData2); //Make sure data offsets are good so far

            if (mech->vdataSize>0 || mech->pntMap>0)
            {
                assert((mech->type == IClamp && mech->vdataSize == 1 && mech->pdataSize == 2 && mech->pntMap>0)
                    || (mech->type == StochKv && mech->vdataSize == 1 && mech->pdataSize == 5 && mech->pntMap==0)
                    || ((mech->type == ProbAMPANMDA_EMS || mech->type == ProbGABAAB_EMS)
                        && mech->vdataSize == 2 && mech->pdataSize == 3 && mech->pntMap>0 ));

                //ProbAMPANMDA_EMS, ProbAMPANMDA_EMS and IClamp:
                //pdata[0]: offset in data (area)
                //pdata[1]: offset for Point_process in vdata[0]
                //pdata[2]: offset for RNG in vdata[1]   (NOT for IClamp)

                //StochKv:
                //pdata[0]: offset in area (ion_ek)
                //pdata[1]: offset in area (ion_ik)
                //pdata[2]: offset in area (ion_dikdv)
                //pdata[3]: offset for RNG in vdata[0]
                //pdata[4]: offset in data (area)

                //copy Point_processes by replacing vdata pointers and pdata offset by the ones referring to a copy
                if (mech->type == IClamp || mech->type == ProbAMPANMDA_EMS || mech->type == ProbGABAAB_EMS)
                {
                    int pointProcOffsetInPdata = 1;
                    Point_process * pp = (Point_process *) (void*) &vdataSerialized[vdataOffset];
                    assert(pp->_i_instance>=0 && pp->_tid>=0 && pp->_type>=0);
                    Point_process * ppcopy = new Point_process;
                    memcpy(ppcopy, pp, sizeof(Point_process));
                    vdataOffset += sizeof(Point_process);
                    instancePdata[pointProcOffsetInPdata] = vdataPtrs.size();
                    vdataPtrs.push_back(ppcopy);
                }

                //copy RNG by replacing vdata pointers and pdata offset by the ones referring to a copy
                if (mech->type == StochKv || mech->type == ProbAMPANMDA_EMS || mech->type == ProbGABAAB_EMS)
                {
                    int rngOffsetInPdata = mech->type == StochKv ? 3 : 2;
                    nrnran123_State * rng = (nrnran123_State*) (void*) &vdataSerialized[vdataOffset];
                    nrnran123_State * rngcopy = new nrnran123_State;
                    memcpy(rngcopy, rng, sizeof(nrnran123_State));
                    vdataOffset += sizeof(nrnran123_State);
                    instancePdata[rngOffsetInPdata] = vdataPtrs.size();
                    vdataPtrs.push_back(rngcopy);
                }
            }
            dataOffset  += mech->dataSize;
            pdataOffset += mech->pdataSize;
            assert(dataOffset  < 2^sizeof(offset_t));
            assert(pdataOffset < 2^sizeof(offset_t));
            assert(vdataOffset < 2^sizeof(offset_t));
            instancesOffset++;
        }
    }
    assert( dataOffset ==  dataCount);
    assert(pdataOffset == pdataCount);
    assert(instancesOffset == nodesIndicesCount);

    //vdata pointers
    nt->_nvdata = vdataPtrs.size();
    nt->_vdata = vdataSerializedCount==0 ? nullptr :  new void*[vdataPtrs.size()];
    memcpy(nt->_vdata, vdataPtrs.data(), vdataPtrs.size()*sizeof(void*));
    vdataPtrs.clear();

    //nt->_ml_list
    nt->_ml_list = new Memb_list*[maxMechId+1];
    for (int i=0; i<=maxMechId; i++)
        nt->_ml_list[i] = NULL;


    int ionsCount=0;
    for (offset_t m=0; m<mechanismsCount; m++)
    {
        Mechanism * mech = mechanisms[m];
        Memb_list & instances = this->mechsInstances[m];
        this->nt->_ml_list[mech->type] = &instances;
        if (mech->isIon) ionsCount++;
    }
    assert(ionsCount==Mechanism::Ion::size_writeable_ions+1); //ttx excluded (no writes to ttx state)

    //vecplay
    nt->n_vecplay = vecplayPPIcount;
    nt->_vecplay = vecplayPPIcount == 0 ? nullptr : new void*[vecplayPPIcount];

    offset_t vOffset=0;
    for (size_t v=0; v < nt->n_vecplay; v++)
    {
        PointProcInfo & ppi = vecplayPPI[v];
        size_t size =  ppi.size;
        int m = mechanismsMap[ppi.mechType];
        floble_t *instancesData = this->mechsInstances[m].data;
        floble_t *pd = &(instancesData[ppi.mechInstance*mechanisms[m]->dataSize + ppi.instanceDataOffset]);
        floble_t * yvec = new floble_t[size];
        floble_t * tvec = new floble_t[size];
        for (size_t i=0; i<size; i++)
        {
            yvec[i] = vecplayY[vOffset+i];
            tvec[i] = vecplayT[vOffset+i];
        }
        nt->_vecplay[v] = new VecPlayContinuousX(pd,size, yvec,tvec, NULL);
        vOffset += size;
    }

    //Shadow arrays
    this->nt->shadow_rhs_cnt=0;
    this->nt->_shadow_d=NULL;
    this->nt->_shadow_rhs=NULL;

    for (int m=0; m<mechanismsCount; m++)
    {
        int shadowSize = 0;
        if (mechanisms[m]->membFunc.current && !mechanisms[m]->isIon) //ions have no updates
            shadowSize = this->mechsInstances[m].nodecount;

        Memb_list * ml = &mechsInstances[m];
        ml->_shadow_d   = shadowSize==0 ? nullptr : Vectorizer::new_<double>(shadowSize);
        ml->_shadow_rhs = shadowSize==0 ? nullptr : Vectorizer::new_<double>(shadowSize);

        for (int i=0; i<shadowSize; i++)
        {
            ml->_shadow_d[i] = 0;
            ml->_shadow_rhs[i] = 0;
        }

        if (mechanisms[m]->dependencyIonIndex >= Mechanism::Ion::size_writeable_ions)
            shadowSize = 0; //> only mechanisms with parent ions update I and DI/DV

        ml->_shadow_i            = shadowSize==0 ? nullptr : Vectorizer::new_<double>(shadowSize);
        ml->_shadow_didv         = shadowSize==0 ? nullptr : Vectorizer::new_<double>(shadowSize);
        ml->_shadow_i_offsets    = shadowSize==0 ? nullptr : Vectorizer::new_<int>(shadowSize);
        ml->_shadow_didv_offsets = shadowSize==0 ? nullptr : Vectorizer::new_<int>(shadowSize);
        for (int i=0; i<shadowSize; i++)
        {
            ml->_shadow_i[i] = 0;
            ml->_shadow_didv[i] = 0;
            ml->_shadow_i_offsets[i] = -1;
            ml->_shadow_didv_offsets[i] = -1;
        }
    }

    //reconstructs netcons
    offset_t weightsOffset=0;
    for (offset_t nc=0; nc<netconsCount; nc++)
    {
        this->netcons[ netConsPreId[nc] ].push_back(
                    new NetConX(netcons[nc].mechType, netcons[nc].mechInstance, netcons[nc].delay,
                    netcons[nc].weightIndex, netcons[nc].weightsCount, netcons[nc].active));
        assert(weightsOffset == netcons[nc].weightIndex);
        weightsOffset += netcons[nc].weightsCount;
    }

    //create branchTree and MechsGraph
    this->branchTree = inputParams->branchingDepth>0 ? new Branch::BranchTree(topBranchAddr, branches,branchesCount) : nullptr;
    this->mechsGraph = inputParams->multiMex         ? new Branch::MechanismsGraph(n)         : nullptr;
    if (this->mechsGraph) mechsGraph->initMechsGraph(branchHpxAddr);
    assert(weightsCount == weightsOffset);

//#if LAYOUT==1
    if (inputParams->vectorize)
        Tools::Vectorizer::vectorize(this);
//#endif
}

Branch::~Branch()
{
    Vectorizer::delete_(this->nt->_data);
    Vectorizer::delete_(this->nt->_v_parent_index);
    delete [] this->nt->weights;
    delete [] this->nt->_ml_list;

    hpx_lco_delete_sync(this->eventsQueueMutex);
    for (int m=0; m<mechanismsCount; m++)
    {
        Memb_list & instance = this->mechsInstances[m];
        if (mechanisms[m]->membFunc.thread_cleanup_)
            mechanisms[m]->membFunc.thread_cleanup_(instance._thread);

        Vectorizer::delete_(mechsInstances[m].nodeindices);
        Vectorizer::delete_(mechsInstances[m].pdata);
        delete [] mechsInstances[m]._thread;

        Vectorizer::delete_(mechsInstances[m]._shadow_d);
        Vectorizer::delete_(mechsInstances[m]._shadow_didv);
        Vectorizer::delete_(mechsInstances[m]._shadow_didv_offsets);
        Vectorizer::delete_(mechsInstances[m]._shadow_i);
        Vectorizer::delete_(mechsInstances[m]._shadow_rhs);
        Vectorizer::delete_(mechsInstances[m]._shadow_i_offsets);
    }
    delete [] mechsInstances;

    for (int i=0; i< this->nt->_nvdata; i++)
        delete nt->_vdata[i];
    delete [] this->nt->_vdata;

    for (int i=0; i< this->nt->n_vecplay; i++)
        delete (VecPlayContinuousX*) nt->_vecplay[i];
    delete [] this->nt->_vecplay;

    free(this->nt);

    for (auto & nc_pair : this->netcons)
        for (auto & nc : nc_pair.second)
            delete nc;

    delete soma;
    delete branchTree;
    delete mechsGraph;
}

hpx_action_t Branch::init = 0;
int Branch::init_handler( const int nargs, const void *args[],
                          const size_t sizes[])
{
    neurox_hpx_pin(Branch);
    assert(nargs==17 || nargs==18); //16 for normal init, 17 for benchmark (initializes, benchmarks, and clears memory)
    new(local) Branch(
        *(offset_t*) args[0], //number of compartments
        *(int*) args[1], //nrnThreadId (nt.id)
        *(int*) args[2], //offset  AP voltage threshold (-1 if none)
        target, //current branch HPX address
        (floble_t*) args[3], sizes[3]/sizeof(floble_t), //data (RHS, D, A, V, B, area, and mechs...)
        (offset_t*) args[4], sizes[4]/sizeof(offset_t), //pdata
        (offset_t*) args[5], sizes[5]/sizeof(offset_t), //instances count per mechanism
        (offset_t*) args[6], sizes[6]/sizeof(offset_t), //nodes indices
        *(hpx_t*) args[7], //top branchAddr
        (hpx_t*) args[8], sizes[8]/sizeof(hpx_t), //branches
        (offset_t*) args[9], sizes[9]/sizeof(offset_t), //parent index
        (floble_t*) args[10], sizes[10]/sizeof(floble_t), //vecplay T data
        (floble_t*) args[11], sizes[11]/sizeof(floble_t), //vecplay Y data
        (PointProcInfo*) args[12], sizes[12]/sizeof(PointProcInfo), //point processes info
        (NetConX*) args[13], sizes[13]/sizeof(NetConX), //netcons
        (neuron_id_t *) args[14], sizes[14]/sizeof(neuron_id_t), //netcons preneuron ids
        (floble_t *) args[15], sizes[15]/sizeof(floble_t), //netcons weights
        (unsigned char*) args[16], sizes[16]/sizeof(unsigned char)); //serialized vdata

    bool runBenchmarkAndClear=false;
    if (nargs==18)
        runBenchmarkAndClear = *(bool*) args[17];
    if (runBenchmarkAndClear)
    {
        //initialize soma (gid irrelevant, APthreshold very high (never spikes)
        local->soma=new Neuron( -1, 999);

        //initialize datatypes and offsets (required for mechs graph-parallelism shadow vecs)
        local->finitialize2();

        //benchmark execution time of a communication-step time-frame
        hpx_time_t now = hpx_time_now();
        for (int i=0; i< Neuron::CommunicationBarrier::commStepSize; i++)
            local->backwardEulerStep();
        double timeElapsed = hpx_time_elapsed_ms(now)/1e3;
        delete local;
        neurox_hpx_unpin_continue(timeElapsed);
    }
    neurox_hpx_unpin;
}

hpx_action_t Branch::initSoma=0;
int Branch::initSoma_handler(const int nargs,
                             const void *args[], const size_t[])
{
    neurox_hpx_pin(Branch);
    assert(nargs==2);
    const neuron_id_t neuronId = *(const neuron_id_t*) args[0];
    const floble_t APthreshold = *(const floble_t*) args[1];
    local->soma=new Neuron(neuronId, APthreshold);
    neurox_hpx_unpin;
}

hpx_action_t Branch::clear=0;
int Branch::clear_handler()
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::clear);
    delete local;
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

void Branch::initVecPlayContinous()
{
    //nrn_play_init
    for (size_t v=0; v<this->nt->n_vecplay; v++)
    {
        VecPlayContinuousX * vecplay = reinterpret_cast<VecPlayContinuousX*>(this->nt->_vecplay[v]);
        vecplay->play_init(this);
    }
}

void Branch::addEventToQueue(floble_t tt, Event * e)
{
    this->eventsQueue.push(make_pair(tt, e));
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
            if (mechanisms[m]->dependenciesCount > 0) continue; //not a top branch
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
                {
                    mechanisms[m]->callModFunction(this, functionId);
                }
        }
    }
}

//netcvode.cpp::PreSyn::send() --> NetCvode::bin_event.cpp
hpx_action_t Branch::addSpikeEvent = 0;
int Branch::addSpikeEvent_handler(
        const int nargs, const void *args[], const size_t[] )
{
    neurox_hpx_pin(Branch);
    assert(nargs == (inputParams->algorithm == Algorithm::BackwardEulerWithTimeDependencyLCO ? 3 : 2));

    //auto source = libhpx_parcel_get_source(p);
    const neuron_id_t preNeuronId = *(const neuron_id_t *) args[0];
    const spike_time_t spikeTime  = *(const spike_time_t*) args[1];

    assert(local->netcons.find(preNeuronId)!=local->netcons.end());
    auto & netcons = local->netcons.at(preNeuronId);
    hpx_lco_sema_p(local->eventsQueueMutex);
    for (auto nc : netcons)
    {
        floble_t deliveryTime = spikeTime+nc->delay;
        local->eventsQueue.push( make_pair(deliveryTime, (Event*) nc) );
    }
    hpx_lco_sema_v_sync(local->eventsQueueMutex);

    if (inputParams->algorithm == Algorithm::BackwardEulerWithTimeDependencyLCO)
    {
        //inform soma of this neuron of new time dependency update
        spike_time_t maxTime = *(const spike_time_t*) args[2];
        hpx_t topBranchAddr = local->soma ? target : local->branchTree->topBranchAddr;
        if (local->soma)
            local->soma->timeDependencies->updateTimeDependency(preNeuronId, maxTime);
        else
            hpx_call(topBranchAddr, Branch::updateTimeDependency, HPX_NULL,
                 &preNeuronId, sizeof(neuron_id_t), &maxTime, sizeof(spike_time_t));
    }
    neurox_hpx_unpin;
}

hpx_action_t Branch::updateTimeDependency = 0;
int Branch::updateTimeDependency_handler(
        const int nargs, const void *args[], const size_t[] )
{
    neurox_hpx_pin(Branch);
    assert(nargs==2 || nargs==3);

    //auto source = libhpx_parcel_get_source(p);
    const neuron_id_t preNeuronId = *(const neuron_id_t *) args[0];
    const spike_time_t maxTime  = *(const spike_time_t*) args[1];
    const bool initPhase = nargs==3 ? *(const bool*) args[2] : false;

    assert(local->soma);
    assert(local->soma->timeDependencies);
    local->soma->timeDependencies->updateTimeDependency(
                preNeuronId, (floble_t) maxTime, local->soma ? local->soma->gid : -1, initPhase);
    neurox_hpx_unpin;
}

void Branch::finitialize2()
{
    floble_t * v = this->nt->_actual_v;
    double t = this->nt->_t;

    //set up by finitialize.c:nrn_finitialize(): if (setv)
    assert(inputParams->secondorder < sizeof(char));
    callModFunction(Mechanism::ModFunction::threadTableCheck);
    initVecPlayContinous();
    deliverEvents(t);

    //set up by finitialize.c:nrn_finitialize(): if (setv)
    for (int n=0; n<this->nt->end; n++)
        v[n]=inputParams->voltage;

    // the INITIAL blocks are ordered so that mechanisms that write
    // concentrations are after ions and before mechanisms that read
    // concentrations.
    callModFunction(Mechanism::ModFunction::before_initialize);
    callModFunction(Mechanism::ModFunction::initialize);
    callModFunction(Mechanism::ModFunction::after_initialize);

    //initEvents(t); //not needed because we copy the status and weights of events
    deliverEvents(t);
    setupTreeMatrix();
    callModFunction(Mechanism::ModFunction::before_step);
    deliverEvents(t);
}

//fadvance_core.c::nrn_fixed_step_thread
void Branch::backwardEulerStep()
{
    double & t  = this->nt->_t;
    const double dt = this->nt->_dt;
    hpx_t spikesLco = HPX_NULL;

    if (soma && inputParams->algorithm==neurox::Algorithm::BackwardEulerWithTimeDependencyLCO)
    {
        //inform time dependants that must be notified in this step
        soma->timeDependencies->sendSteppingNotification(t, dt, this->soma->gid, this->soma->synapses);
        //wait until Im sure I can start and finalize this step at t+dt
        soma->timeDependencies->waitForTimeDependencyNeurons(t, dt, this->soma->gid);
    }

    //cvodestb.cpp::deliver_net_events()
    //netcvode.cpp::NetCvode::check_thresh(NrnThread*)
    if (this->soma)
    {
      //Soma waits for AIS to have threshold V value updated
      floble_t thresholdV;
      HinesSolver::synchronizeThresholdV(this, &thresholdV);
      if (soma->checkAPthresholdAndTransmissionFlag(thresholdV))
          spikesLco = soma->sendSpikes(nt->_t);
    }
    else if (this->thvar_ptr)
        //Axon Initial Segment send threshold  V to parent
        HinesSolver::synchronizeThresholdV(this);


    //netcvode.cpp::NetCvode::deliver_net_events()
    t += .5*dt;
    deliverEvents(t); //delivers events in the first HALF of the step
    fixedPlayContinuous();
    setupTreeMatrix();
    solveTreeMatrix();
    second_order_cur(this->nt, inputParams->secondorder );

    ////// fadvance_core.c : update()
    Solver::HinesSolver::updateV(this);

    callModFunction(Mechanism::ModFunction::currentCapacitance);

    ////// fadvance_core.::nrn_fixed_step_lastpart()
    //callModFunction(Mechanism::ModFunction::jacob);
    t += .5*dt;
    fixedPlayContinuous();
    callModFunction(Mechanism::ModFunction::state);
    callModFunction(Mechanism::ModFunction::after_solve);
    callModFunction(Mechanism::ModFunction::before_step);
    deliverEvents(t); //delivers events in the second HALF of the step

    //if we are at the output time instant output to file
    if (fmod(t, inputParams->dt_io) == 0) {}

    //wait for spikes sent 4 steps ago (queue has always size 3)
    if (soma)
    if (inputParams->algorithm == Algorithm::BackwardEulerWithSlidingTimeWindow
     || inputParams->algorithm == Algorithm::BackwardEulerWithAllReduceBarrier)
    {
      Neuron::SlidingTimeWindow * stw = this->soma->slidingTimeWindow;
      assert(stw->spikesLcoQueue.size() == Neuron::CommunicationBarrier::commStepSize-1);
      stw->spikesLcoQueue.push(spikesLco);
      hpx_t queuedSpikesLco = stw->spikesLcoQueue.front();
      stw->spikesLcoQueue.pop();
      if (queuedSpikesLco != HPX_NULL)
      {
          hpx_lco_wait(queuedSpikesLco);
          hpx_lco_delete_sync(queuedSpikesLco);
      }
    }

#if !defined(NDEBUG) 
    //fixed comm barrier and serial jobs can be compared at runtime
    if (inputParams->branchingDepth==0)
    if (inputParams->algorithm == Algorithm::BackwardEulerDebugWithCommBarrier
     || !inputParams->parallelDataLoading)
    {
        Input::Debugger::fixed_step_minimal(&nrn_threads[this->nt->id], inputParams->secondorder);
        Input::Debugger::compareBranch2(this);
    }
#endif
}

hpx_action_t Branch::backwardEuler = 0;
int Branch::backwardEuler_handler(const int * steps_ptr, const size_t size)
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::backwardEuler, steps_ptr, size);

    const int steps = *steps_ptr;
    if (local->soma)
    {
      //fixes crash for Algorithm::ALL when TimeDependency algorithm starts at t=inputParams->tend*2
      if (inputParams->algorithm == Algorithm::BackwardEulerWithTimeDependencyLCO)
        {
          //increase notification and dependencies time
          for (Neuron::Synapse *& s : local->soma->synapses)
              s->nextNotificationTime += local->nt->_t;
          local->soma->timeDependencies->increseDependenciesTime(local->nt->_t);
        }
    }

    if (inputParams->algorithm == Algorithm::BackwardEulerWithSlidingTimeWindow
    ||  inputParams->algorithm == Algorithm::BackwardEulerWithAllReduceBarrier)
    {
        const int reductionsPerCommStep = Neuron::SlidingTimeWindow::reductionsPerCommStep;
        const int stepsPerReduction = Neuron::CommunicationBarrier::commStepSize / Neuron::SlidingTimeWindow::reductionsPerCommStep;
        const int commStepSize = Neuron::CommunicationBarrier::commStepSize;
        const Neuron::SlidingTimeWindow * stw = local->soma ? local->soma->slidingTimeWindow : nullptr;

        for (int s=0; s<steps; s += commStepSize) //for every communication step
        {
          #ifdef NEUROX_TIME_STEPPING_VERBOSE
              if (hpx_get_my_rank()==0 && target == neurox::neurons->at(0))
              {
                  printf("-- t=%.4f ms\n", inputParams->dt*s);
                  fflush(stdout);
              }
          #endif
          for (int r=0; r<reductionsPerCommStep; r++) //for every reduction step
          {
              if (local->soma)
              {
                  if (s>= commStepSize) //first comm-window does not wait
                    hpx_lco_wait_reset(stw->allReduceFuture[r]);
                  else
                    //fixes crash for Algorithm::ALL when running two hpx-reduce -based algorithms in a row
                    hpx_lco_reset_sync(stw->allReduceFuture[r]);

                  hpx_process_collective_allreduce_join(stw->allReduceLco[r], stw->allReduceId[r], NULL, 0);
              }

              for (int n=0; n<stepsPerReduction; n++)
                  local->backwardEulerStep();
                  // Input::Coreneuron::Debugger::stepAfterStepBackwardEuler(local, &nrn_threads[this->nt->id], secondorder); //SMP ONLY
          }
        }
    }
    else
    {
        for (int step=0; step<steps; step++)
            local->backwardEulerStep();
            // Input::Coreneuron::Debugger::stepAfterStepBackwardEuler(local, &nrn_threads[this->nt->id], secondorder); //SMP ONLY

        if (local->soma)
        if (inputParams->algorithm == Algorithm::BackwardEulerDebugWithCommBarrier) //end of comm-step
            if (local->soma->commBarrier->allSpikesLco != HPX_NULL) //was set/used once
                hpx_lco_wait(local->soma->commBarrier->allSpikesLco); //wait if needed

        #ifdef NEUROX_TIME_STEPPING_VERBOSE
            if (local->soma && inputParams->algorithm == Algorithm::BackwardEulerWithTimeDependencyLCO)
                printf("-- neuron %d finished\n", local->soma->gid);
        #endif
    }
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t Branch::backwardEulerOnLocality = 0;
int Branch::backwardEulerOnLocality_handler(const int * steps_ptr, const size_t size)
{
    neurox_hpx_pin(uint64_t);
    assert(inputParams->allReduceAtLocality);
    assert(Neuron::SlidingTimeWindow::AllReduceLocality::localityNeurons);
    assert(inputParams->algorithm == Algorithm::BackwardEulerWithSlidingTimeWindow
        || inputParams->algorithm == Algorithm::BackwardEulerWithAllReduceBarrier);

    const int localityNeuronsCount = Neuron::SlidingTimeWindow::AllReduceLocality::localityNeurons->size();
    const hpx_t localityNeuronsLCO = hpx_lco_and_new(localityNeuronsCount);
    const int commStepSize = Neuron::CommunicationBarrier::commStepSize;
    const int reductionsPerCommStep = Neuron::SlidingTimeWindow::reductionsPerCommStep;
    const int stepsPerReduction = Neuron::CommunicationBarrier::commStepSize / Neuron::SlidingTimeWindow::reductionsPerCommStep;
    const int steps = *steps_ptr;

    for (int s=0; s<steps; s+= commStepSize)
    {
        for (int r=0; r< reductionsPerCommStep; r++)
        {
            if (s>= commStepSize) //first comm-window does not wait
                hpx_lco_wait_reset(Neuron::SlidingTimeWindow::AllReduceLocality::allReduceFuture[r]);
            else
                //fixes crash for Algorithm::ALL when running two hpx-reduce -based algorithms in a row
                hpx_lco_reset_sync(Neuron::SlidingTimeWindow::AllReduceLocality::allReduceFuture[r]);

            hpx_process_collective_allreduce_join(
                        Neuron::SlidingTimeWindow::AllReduceLocality::allReduceLco[r],
                        Neuron::SlidingTimeWindow::AllReduceLocality::allReduceId[r], NULL, 0);

            for (int i=0; i<localityNeuronsCount; i++)
               hpx_call(Neuron::SlidingTimeWindow::AllReduceLocality::localityNeurons->at(i), Branch::backwardEuler,
                        localityNeuronsLCO, &stepsPerReduction, sizeof(int));
            hpx_lco_wait_reset(localityNeuronsLCO);
        }
    }
    hpx_lco_delete_sync(localityNeuronsLCO);
    neurox_hpx_unpin;
}

hpx_action_t Branch::threadTableCheck =0;
int Branch::threadTableCheck_handler()
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::threadTableCheck);
    local->callModFunction(Mechanism::ModFunction::threadTableCheck);
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t Branch::finitialize = 0;
int Branch::finitialize_handler()
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::finitialize);
    local->finitialize2();
#if !defined(NDEBUG)
    //Input::Debugger::stepAfterStepFinitialize(local, &nrn_threads[local->nt->id]);
#endif
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

void Branch::setupTreeMatrix()
{
    //treeset_core.c::nrn_rhs: Set up Right-Hand-Side of Matrix-Vector multiplication
    Solver::HinesSolver::resetMatrixRHSandD(this);

    this->callModFunction(Mechanism::ModFunction::before_breakpoint);
    this->callModFunction(Mechanism::ModFunction::current);

    Solver::HinesSolver::setupMatrixRHS(this);

    //treeset_core.c::nrn_lhs: Set up Left-Hand-Side of Matrix-Vector multiplication
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

    Solver::HinesSolver::setupMatrixDiagonal(this);
}

void Branch::solveTreeMatrix()
{
    Solver::HinesSolver::backwardTriangulation(this);
    Solver::HinesSolver::forwardSubstituion(this);
}

void Branch::deliverEvents(floble_t til) //til=t+0.5*dt
{
    //delivers events in the preivous half-step
    floble_t tsav = this->nt->_t; //copying cvodestb.cpp logic
    hpx_lco_sema_p(this->eventsQueueMutex);
    while (!this->eventsQueue.empty()
         && this->eventsQueue.top().first <= til)
     {
        auto event_it = this->eventsQueue.top();
        floble_t & tt = event_it.first;
        Event *& e = event_it.second;

        //TODO if (tt>=0) //until CoreNeuron issue of negative delivery times is fixed
        //must be delivered in previous half step (not earlier, or it's been missed)
        //{   assert(tt >= til-0.5*nt->_dt-Neuron::teps && tt<=til+Neuron::teps); }
        //{   assert(tt >= til-0.5*nt->_dt && tt<=til); }

        e->deliver(tt, this);
        this->eventsQueue.pop();
    }
    hpx_lco_sema_v_sync(this->eventsQueueMutex);
    this->nt->_t = tsav;
}

void Branch::fixedPlayContinuous()
{
    double t = this->nt->_t;
    for (int v=0; v<this->nt->n_vecplay; v++)
    {
        void * vecplay_void = this->nt->_vecplay[v];
        VecPlayContinuousX * vecplay = reinterpret_cast<VecPlayContinuousX*>(vecplay_void);
        vecplay->continuous(t);
    }
}

//////////////////// Branch::NeuronTree ///////////////////////

Branch::BranchTree::BranchTree(
        hpx_t topBranchAddr, hpx_t * branches, size_t branchesCount)
    : topBranchAddr(topBranchAddr), branches(nullptr), branchesCount(branchesCount), withChildrenLCOs(nullptr)
{
    if (branchesCount>0)
    {
      this->branches = new hpx_t[branchesCount];
      memcpy(this->branches, branches, branchesCount*sizeof(hpx_t));
    }
}

Branch::BranchTree::~BranchTree()
{
    delete [] branches;
    delete [] withChildrenLCOs;
}

hpx_action_t Branch::BranchTree::initLCOs = 0;
int Branch::BranchTree::initLCOs_handler()
{
    neurox_hpx_pin(Branch);
    BranchTree * branchTree = local->branchTree;
    if (branchTree)
    {
      offset_t branchesCount = branchTree->branchesCount;
      for (int i=0; i<BranchTree::futuresSize; i++)
          branchTree->withParentLCO[i] = local->soma ? HPX_NULL : hpx_lco_future_new(sizeof(floble_t));
      branchTree->withChildrenLCOs =  branchesCount ? new hpx_t[branchesCount][BranchTree::futuresSize] : nullptr;

      //send my LCOs to children, and receive theirs
      if (branchesCount>0)
      {
        hpx_t * futures = branchesCount ? new hpx_t[branchesCount]  : nullptr;
        void ** addrs   = branchesCount ? new void*[branchesCount]  : nullptr;
        size_t* sizes   = branchesCount ? new size_t[branchesCount] : nullptr;
        for (offset_t c = 0; c < branchesCount; c++)
        {
          futures[c] = hpx_lco_future_new(sizeof (hpx_t)*BranchTree::futuresSize);
          addrs[c]   = &branchTree->withChildrenLCOs[c];
          sizes[c]   = sizeof(hpx_t)*BranchTree::futuresSize;
          hpx_call(branchTree->branches[c], Branch::BranchTree::initLCOs,
                   futures[c], branchTree->withParentLCO,
                   sizeof(hpx_t)*BranchTree::futuresSize); //pass my LCO down
        }
        hpx_lco_get_all(branchesCount, futures, sizes, addrs, NULL);
        hpx_lco_delete_all(branchesCount, futures, NULL);

        delete [] futures;
        delete [] addrs;
        delete [] sizes;
      }

      if (!local->soma) //send my LCO to parent
          neurox_hpx_unpin_continue(branchTree->withParentLCO);
    }
    neurox_hpx_unpin;
}

//////////////////// Branch::MechanismsGraph ///////////////////////

Branch::MechanismsGraph::MechanismsGraph(int compartmentsCount)
{
    //initializes mechanisms graphs (capacitance is excluded from graph)
    this->graphLCO  = hpx_lco_and_new(mechanismsCount-1); //excludes 'capacitance'
    this->mechsLCOs = new hpx_t[mechanismsCount];
    this->mechsLCOs[mechanismsMap[CAP]] = HPX_NULL;
    size_t terminalMechanismsCount=0;
    for (size_t m=0; m<mechanismsCount; m++)
    {
      if (mechanisms[m]->type == CAP) continue; //exclude capacitance
      this->mechsLCOs[m] = hpx_lco_reduce_new(
                  max((short) 1,mechanisms[m]->dependenciesCount),
                  sizeof(Mechanism::ModFunction),
                  Branch::MechanismsGraph::init,
                  Branch::MechanismsGraph::reduce);
      if (mechanisms[m]->successorsCount==0) //bottom of mechs graph
          terminalMechanismsCount++;
    }
    this->endLCO = hpx_lco_and_new(terminalMechanismsCount);

    this->rhs_d_mutex = hpx_lco_sema_new(1);
    for (int i=0; i<Mechanism::Ion::size_writeable_ions;i++)
        this->i_didv_mutex[i] = hpx_lco_sema_new(1);
}

void Branch::MechanismsGraph::initMechsGraph(hpx_t branchHpxAddr)
{
    for (size_t m=0; m<mechanismsCount; m++)
      if (mechanisms[m]->type != CAP) //exclude capacitance
        hpx_call(branchHpxAddr, Branch::MechanismsGraph::mechFunction,
               this->graphLCO, &mechanisms[m]->type, sizeof(int));
}

Branch::MechanismsGraph::~MechanismsGraph()
{
    hpx_lco_delete_sync(endLCO);
    hpx_lco_delete_sync(graphLCO);

    for (int i=0; i<mechanismsCount; i++)
        if (i != mechanismsMap[CAP]) //HPX_NULL
            hpx_lco_delete_sync(mechsLCOs[i]);
    delete [] mechsLCOs;

    hpx_lco_delete_sync(rhs_d_mutex);
    for (int i=0; i<Mechanism::Ion::size_writeable_ions; i++)
        hpx_lco_delete_sync(i_didv_mutex[i]);
}

hpx_action_t Branch::MechanismsGraph::mechFunction = 0;
int Branch::MechanismsGraph::mechFunction_handler(
        const int * mechType_ptr, const size_t)
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

void Branch::MechanismsGraph::accumulate_rhs_d( NrnThread* nt, Memb_list * ml, int, void* args)
{
    MechanismsGraph * mg = (MechanismsGraph*) args;
    hpx_lco_sema_p(mg->rhs_d_mutex);
    for (int n=0; n<ml->nodecount; n++)
    {
        int & i = ml->nodeindices[n];
        nt->_actual_rhs[i] += ml->_shadow_rhs[n];
        nt->_actual_d[i] += ml->_shadow_d[n];
    }
    hpx_lco_sema_v_sync(mg->rhs_d_mutex);
}

void Branch::MechanismsGraph::accumulate_i_didv(NrnThread* nt,  Memb_list * ml, int type, void* args)
{
    MechanismsGraph * mg = (MechanismsGraph*) args;
    Mechanism * mech = getMechanismFromType(type);
    assert(mech->dependencyIonIndex<Mechanism::Ion::size_writeable_ions);
    hpx_lco_sema_p(mg->i_didv_mutex[mech->dependencyIonIndex]); //TODO should mutex be at branch and not mech level?
    for (int n=0; n<ml->nodecount; n++)
    {
        //TODO does it need *_STRIDE to vectorize?
        int & i_offset = ml->_shadow_i_offsets[n];
        int & didv_offset = ml->_shadow_didv_offsets[n];
        assert(i_offset>=0 && didv_offset>=0);
        nt->_data[i_offset] += ml->_shadow_i[n];
        nt->_data[didv_offset] += ml->_shadow_didv[n];
    }
    hpx_lco_sema_v_sync(mg->i_didv_mutex[mech->dependencyIonIndex]);
}

hpx_action_t Branch::MechanismsGraph::init = 0;
void Branch::MechanismsGraph::init_handler(Mechanism::ModFunction *, const size_t)
{}

hpx_action_t Branch::MechanismsGraph::reduce = 0;
void Branch::MechanismsGraph::reduce_handler
    (Mechanism::ModFunction * lhs, const Mechanism::ModFunction *rhs, const size_t)
{ *lhs = *rhs; }

void Branch::registerHpxActions()
{
    neurox_hpx_register_action(neurox_zero_var_action,     Branch::clear);
    neurox_hpx_register_action(neurox_zero_var_action,     Branch::finitialize);
    neurox_hpx_register_action(neurox_zero_var_action,     Branch::threadTableCheck);
    neurox_hpx_register_action(neurox_zero_var_action,     Branch::BranchTree::initLCOs);
    neurox_hpx_register_action(neurox_single_var_action,   Branch::backwardEuler);
    neurox_hpx_register_action(neurox_single_var_action,   Branch::backwardEulerOnLocality);
    neurox_hpx_register_action(neurox_single_var_action,   Branch::MechanismsGraph::mechFunction);
    neurox_hpx_register_action(neurox_several_vars_action, Branch::init);
    neurox_hpx_register_action(neurox_several_vars_action, Branch::initSoma);
    neurox_hpx_register_action(neurox_several_vars_action, Branch::addSpikeEvent);
    neurox_hpx_register_action(neurox_several_vars_action, Branch::updateTimeDependency);
    neurox_hpx_register_action(neurox_reduce_op_action,    Branch::MechanismsGraph::init);
    neurox_hpx_register_action(neurox_reduce_op_action,    Branch::MechanismsGraph::reduce);
}
