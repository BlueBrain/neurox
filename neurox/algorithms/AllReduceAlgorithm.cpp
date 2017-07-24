#include "neurox/algorithms/AllReduceAlgorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

hpx_t* DERIVED_CLASS_NAME::allReduces = nullptr;

DERIVED_CLASS_NAME::DERIVED_CLASS_NAME()
{
    AllReduceAlgorithm::AllReducesInfo::reductionsPerCommStep = DERIVED_CLASS_NAME::allReducesCount;
}

DERIVED_CLASS_NAME::~DERIVED_CLASS_NAME() {}

const AlgorithmType DERIVED_CLASS_NAME::getType()
{
    return AlgorithmType::BackwardEulerAllReduce;
}

const char* DERIVED_CLASS_NAME::getTypeString()
{
    return "BackwardEulerAllReduce";
}

void DERIVED_CLASS_NAME::Init() {
    SubscribeAllReduces(DERIVED_CLASS_NAME::allReduces,
                        DERIVED_CLASS_NAME::allReducesCount);
}

void DERIVED_CLASS_NAME::Clear() {
    UnsubscribeAllReduces(DERIVED_CLASS_NAME::allReduces,
                          DERIVED_CLASS_NAME::allReducesCount);
}

double DERIVED_CLASS_NAME::Launch()
{
    int totalSteps = Algorithm::getTotalStepsCount();
    hpx_time_t now = hpx_time_now();
    if (inputParams->allReduceAtLocality)
        hpx_bcast_rsync(Branch::BackwardEulerOnLocality, &totalSteps, sizeof(int));
    else
        neurox_hpx_call_neurons_lco(Branch::BackwardEuler, &totalSteps, sizeof(int));
    double elapsedTime = hpx_time_elapsed_ms(now)/1e3;
    input::Debugger::RunCoreneuronAndCompareAllBranches();
    return elapsedTime;
}

void DERIVED_CLASS_NAME::StepBegin(Branch*) {}

void DERIVED_CLASS_NAME::StepEnd(Branch* b, hpx_t spikesLco)
{
    WaitForSpikesDelivery(b, spikesLco);
    input::Debugger::SingleNeuronStepAndCompare(&nrn_threads[b->nt->id], b, inputParams->secondorder);
}

void DERIVED_CLASS_NAME::Run(Branch*b, const void* args)
{
    Run2(b, args);
}

hpx_t DERIVED_CLASS_NAME::SendSpikes(Neuron* b, double tt, double)
{
    return SendSpikes2(b,tt);
}

void DERIVED_CLASS_NAME::SubscribeAllReduces(hpx_t *& allReduces, size_t allReducesCount)
{
    assert(allReduces==nullptr);
    allReduces = new hpx_t[allReducesCount];

    hpx_bcast_rsync(AllReduceAlgorithm::AllReducesInfo::SetReductionsPerCommStep, &allReducesCount, sizeof(int));

    for (int i=0; i<allReducesCount; i++)
        allReduces[i] = hpx_process_collective_allreduce_new(
                    0, AllReduceAlgorithm::AllReducesInfo::Init, AllReduceAlgorithm::AllReducesInfo::Reduce);

    if (inputParams->allReduceAtLocality)
        hpx_bcast_rsync(AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::SubscribeAllReduce,
                        allReduces, sizeof(hpx_t)*allReducesCount);
    else
        neurox_hpx_call_neurons_lco(AllReduceAlgorithm::AllReducesInfo::SubscribeAllReduce,
                        allReduces, sizeof(hpx_t)*allReducesCount);

    for (int i=0; i<allReducesCount; i++)
        hpx_process_collective_allreduce_subscribe_finalize(allReduces[i]);
}

void DERIVED_CLASS_NAME::UnsubscribeAllReduces(hpx_t *& allReduces, size_t allReducesCount)
{
    assert(allReduces!=nullptr);
    if (inputParams->allReduceAtLocality)
        hpx_bcast_rsync(AllReduceAlgorithm::AllReducesInfo::AllReduceLocality::UnsubscribeAllReduce,
                        allReduces, sizeof(hpx_t)*allReducesCount);
    else
        neurox_hpx_call_neurons_lco(AllReduceAlgorithm::AllReducesInfo::UnsubscribeAllReduce,
                        allReduces, sizeof(hpx_t)*allReducesCount);

    for (int i=0; i<allReducesCount; i++)
        hpx_process_collective_allreduce_delete(allReduces[i]);

    delete [] allReduces; allReduces=nullptr;
}

void DERIVED_CLASS_NAME::WaitForSpikesDelivery(Branch *b, hpx_t spikesLco)
{
    //wait for spikes sent 4 steps ago (queue has always size 3)
    if (b->soma)
    {
      AllReducesInfo * stw = (AllReducesInfo*) b->soma->algorithmMetaData;
      std::queue<hpx_t> q = stw->spikesLcoQueue;
      assert(stw->spikesLcoQueue.size() == CoreneuronDebugAlgorithm::CommunicationBarrier::commStepSize-1);
      stw->spikesLcoQueue.push(spikesLco);
      hpx_t queuedSpikesLco = stw->spikesLcoQueue.front();
      stw->spikesLcoQueue.pop();
      if (queuedSpikesLco != HPX_NULL)
      {
          hpx_lco_wait(queuedSpikesLco);
          hpx_lco_delete_sync(queuedSpikesLco);
      }
    }
}

hpx_t DERIVED_CLASS_NAME::SendSpikes2(Neuron *neuron, double tt)
{
    hpx_t newSynapsesLco = hpx_lco_and_new(neuron->synapses.size());
    for (Neuron::Synapse *& s : neuron->synapses)
        hpx_call(s->branchAddr, Branch::AddSpikeEvent, newSynapsesLco,
            &neuron->gid, sizeof(neuron_id_t), &tt, sizeof(spike_time_t));
    return newSynapsesLco;
}

void DERIVED_CLASS_NAME::Run2(Branch *b, const void *args)
{
    int steps = *(int*)args;
    const int reductionsPerCommStep = AllReduceAlgorithm::AllReducesInfo::reductionsPerCommStep;
    const int commStepSize = CoreneuronDebugAlgorithm::CommunicationBarrier::commStepSize;
    const int stepsPerReduction = commStepSize / reductionsPerCommStep;
    const AllReducesInfo * stw = b->soma ? (AllReducesInfo*) b->soma->algorithmMetaData : nullptr;

    for (int s=0; s<steps; s += commStepSize) //for every communication step
    {
      #ifndef NDEBUG
          if (hpx_get_my_rank()==0 && b->nt->id == 0)
          {
              printf("-- t=%.4f ms\n", inputParams->dt*s);
              fflush(stdout);
          }
      #endif
      for (int r=0; r<reductionsPerCommStep; r++) //for every reduction step
      {
          if (b->soma)
          {
              if (s>= commStepSize) //first comm-window does not wait
                hpx_lco_wait_reset(stw->allReduceFuture[r]);
              else
                //fixes crash for Algorithm::ALL when running two hpx-reduce -based algorithms in a row
                hpx_lco_reset_sync(stw->allReduceFuture[r]);

              hpx_process_collective_allreduce_join(stw->allReduceLco[r], stw->allReduceId[r], NULL, 0);
          }

          for (int n=0; n<stepsPerReduction; n++)
              b->BackwardEulerStep();
              // Input::Coreneuron::Debugger::stepAfterStepBackwardEuler(local, &nrn_threads[this->nt->id], secondorder); //SMP ONLY
      }
    }
}


DERIVED_CLASS_NAME::AllReducesInfo::AllReducesInfo()
{
    for (int s=0; s< CoreneuronDebugAlgorithm::CommunicationBarrier::commStepSize-1; s++)
        this->spikesLcoQueue.push(HPX_NULL);
}

DERIVED_CLASS_NAME::AllReducesInfo::~AllReducesInfo()
{
    for (int i=0; i<spikesLcoQueue.size(); i++)
    {
        hpx_t queuedSpikesLco = spikesLcoQueue.front();
        if (queuedSpikesLco != HPX_NULL)
            hpx_lco_delete_sync(queuedSpikesLco);
        spikesLcoQueue.pop();
    }
}

hpx_action_t DERIVED_CLASS_NAME::AllReducesInfo::SubscribeAllReduce = 0;
int DERIVED_CLASS_NAME::AllReducesInfo::SubscribeAllReduce_handler(const hpx_t * allreduces, const size_t size)
{
    neurox_hpx_pin(Branch);
    AllReducesInfo * stw = (AllReducesInfo*) local->soma->algorithmMetaData;
    stw->allReduceFuture = new hpx_t[AllReducesInfo::reductionsPerCommStep];
    stw->allReduceLco = new hpx_t[AllReducesInfo::reductionsPerCommStep];
    stw->allReduceId = new int[AllReducesInfo::reductionsPerCommStep];
    for (int i=0; i<size/sizeof(hpx_t); i++)
    {
        stw->allReduceLco[i] = allreduces[i];
        stw->allReduceFuture[i] = hpx_lco_future_new(0); //no value to be reduced
        stw->allReduceId[i] = hpx_process_collective_allreduce_subscribe(
                allreduces[i], hpx_lco_set_action, stw->allReduceFuture[i]);
    }
    neurox_hpx_unpin;
}

hpx_action_t DERIVED_CLASS_NAME::AllReducesInfo::UnsubscribeAllReduce = 0;
int DERIVED_CLASS_NAME::AllReducesInfo::UnsubscribeAllReduce_handler(const hpx_t * allreduces, const size_t size)
{
    neurox_hpx_pin(Branch);
    AllReducesInfo * stw = (AllReducesInfo*) local->soma->algorithmMetaData;
    for (int i=0; i<size/sizeof(hpx_t); i++)
    {
      hpx_process_collective_allreduce_unsubscribe(allreduces[i], stw->allReduceId[i]);
      if (stw->allReduceFuture[i]!= HPX_NULL)
          hpx_lco_delete_sync(stw->allReduceFuture[i]);
    }
    delete [] stw->allReduceLco; stw->allReduceLco=nullptr;
    delete [] stw->allReduceFuture; stw->allReduceFuture=nullptr;
    delete [] stw->allReduceId; stw->allReduceId=nullptr;
    neurox_hpx_unpin;
}

int DERIVED_CLASS_NAME::AllReducesInfo::reductionsPerCommStep = -1;
std::vector<hpx_t>* DERIVED_CLASS_NAME::AllReducesInfo::AllReduceLocality::localityNeurons = nullptr;
hpx_t* DERIVED_CLASS_NAME::AllReducesInfo::AllReduceLocality::allReduceFuture = nullptr;
hpx_t* DERIVED_CLASS_NAME::AllReducesInfo::AllReduceLocality::allReduceLco = nullptr;
int* DERIVED_CLASS_NAME::AllReducesInfo::AllReduceLocality::allReduceId = nullptr;

hpx_action_t DERIVED_CLASS_NAME::AllReducesInfo::SetReductionsPerCommStep = 0;
int DERIVED_CLASS_NAME::AllReducesInfo::SetReductionsPerCommStep_handler(const int* val, const size_t)
{
    neurox_hpx_pin(uint64_t);
    reductionsPerCommStep = *val;
    neurox_hpx_unpin;
}

hpx_action_t DERIVED_CLASS_NAME::AllReducesInfo::AllReduceLocality::SubscribeAllReduce = 0;
int DERIVED_CLASS_NAME::AllReducesInfo::AllReduceLocality::SubscribeAllReduce_handler(const hpx_t * allreduces, const size_t size)
{
    neurox_hpx_pin(uint64_t);
    assert(inputParams->allReduceAtLocality);
    AllReduceLocality::allReduceLco = new hpx_t[AllReducesInfo::reductionsPerCommStep];
    AllReduceLocality::allReduceFuture = new hpx_t[AllReducesInfo::reductionsPerCommStep];
    AllReduceLocality::allReduceId = new int[AllReducesInfo::reductionsPerCommStep];
    for (int i=0; i<size/sizeof(hpx_t); i++)
    {
        allReduceLco[i] = allreduces[i];
        allReduceFuture[i] = hpx_lco_future_new(0); //no value to be reduced
        allReduceId[i] = hpx_process_collective_allreduce_subscribe(
                allreduces[i], hpx_lco_set_action, allReduceFuture[i]);
    }
    neurox_hpx_unpin;
}

hpx_action_t DERIVED_CLASS_NAME::AllReducesInfo::AllReduceLocality::UnsubscribeAllReduce = 0;
int DERIVED_CLASS_NAME::AllReducesInfo::AllReduceLocality::UnsubscribeAllReduce_handler(const hpx_t * allreduces, const size_t size)
{
    neurox_hpx_pin(uint64_t);
    assert(inputParams->allReduceAtLocality);
    for (int i=0; i<size/sizeof(hpx_t); i++)
    {
      hpx_process_collective_allreduce_unsubscribe(allreduces[i], allReduceId[i]);
      hpx_lco_delete_sync(allReduceFuture[i]);
    }
    delete [] allReduceLco; allReduceLco=nullptr;
    delete [] allReduceFuture; allReduceFuture=nullptr;
    delete [] allReduceId; allReduceId=nullptr;
    neurox_hpx_unpin;
}

hpx_action_t DERIVED_CLASS_NAME::AllReducesInfo::Init = 0;
void DERIVED_CLASS_NAME::AllReducesInfo::Init_handler
    (const void*, const size_t) {}

hpx_action_t DERIVED_CLASS_NAME::AllReducesInfo::Reduce = 0;
void DERIVED_CLASS_NAME::AllReducesInfo::Reduce_handler
    (void* , const void* , const size_t) {}


void DERIVED_CLASS_NAME::AllReducesInfo::registerHpxActions()
{
    neurox_hpx_register_action(neurox_single_var_action, AllReducesInfo::SubscribeAllReduce);
    neurox_hpx_register_action(neurox_single_var_action, AllReducesInfo::UnsubscribeAllReduce);
    neurox_hpx_register_action(neurox_single_var_action, AllReducesInfo::AllReduceLocality::SubscribeAllReduce);
    neurox_hpx_register_action(neurox_single_var_action, AllReducesInfo::AllReduceLocality::UnsubscribeAllReduce);
    neurox_hpx_register_action(neurox_single_var_action, AllReducesInfo::SetReductionsPerCommStep);
    neurox_hpx_register_action(neurox_reduce_op_action,  AllReducesInfo::Init);
    neurox_hpx_register_action(neurox_reduce_op_action,  AllReducesInfo::Reduce);
}
