#include "neurox/algorithms/BackwardEulerAllReduceAlgorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

hpx_t* DERIVED_CLASS_NAME::allReduces = nullptr;

DERIVED_CLASS_NAME::DERIVED_CLASS_NAME()
{
    Neuron::SlidingTimeWindow::reductionsPerCommStep = DERIVED_CLASS_NAME::allReducesCount;
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

void DERIVED_CLASS_NAME::Finalize() {
    UnsubscribeAllReduces(DERIVED_CLASS_NAME::allReduces,
                          DERIVED_CLASS_NAME::allReducesCount);
}

double DERIVED_CLASS_NAME::Run()
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

void DERIVED_CLASS_NAME::CommStepBegin(Branch*) {}

void DERIVED_CLASS_NAME::CommStepEnd(Branch*) {}

hpx_t DERIVED_CLASS_NAME::SendSpikes(Neuron* b, double tt, double)
{
    return SendSpikes2(b,tt);
}

void DERIVED_CLASS_NAME::SubscribeAllReduces(hpx_t *& allReduces, size_t allReducesCount)
{
    assert(allReduces==nullptr);
    allReduces = new hpx_t[allReducesCount];

    hpx_bcast_rsync(neurox::Neuron::SlidingTimeWindow::SetReductionsPerCommStep, &allReducesCount, sizeof(int));

    for (int i=0; i<allReducesCount; i++)
        allReduces[i] = hpx_process_collective_allreduce_new(0, Neuron::SlidingTimeWindow::Init, Neuron::SlidingTimeWindow::Reduce);

    if (inputParams->allReduceAtLocality)
        hpx_bcast_rsync(neurox::Neuron::SlidingTimeWindow::AllReduceLocality::SubscribeAllReduce,
                        allReduces, sizeof(hpx_t)*allReducesCount);
    else
        neurox_hpx_call_neurons_lco(Neuron::SlidingTimeWindow::SubscribeAllReduce,
                        allReduces, sizeof(hpx_t)*allReducesCount);

    for (int i=0; i<allReducesCount; i++)
        hpx_process_collective_allreduce_subscribe_finalize(allReduces[i]);
}

void DERIVED_CLASS_NAME::UnsubscribeAllReduces(hpx_t *& allReduces, size_t allReducesCount)
{
    assert(allReduces!=nullptr);
    if (inputParams->allReduceAtLocality)
        hpx_bcast_rsync(neurox::Neuron::SlidingTimeWindow::AllReduceLocality::UnsubscribeAllReduce,
                        allReduces, sizeof(hpx_t)*allReducesCount);
    else
        neurox_hpx_call_neurons_lco(Neuron::SlidingTimeWindow::UnsubscribeAllReduce,
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
      Neuron::SlidingTimeWindow * stw = b->soma->slidingTimeWindow;
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
}

hpx_t DERIVED_CLASS_NAME::SendSpikes2(Neuron *neuron, double tt)
{
    hpx_t newSynapsesLco = hpx_lco_and_new(neuron->synapses.size());
    for (Neuron::Synapse *& s : neuron->synapses)
        hpx_call(s->branchAddr, Branch::AddSpikeEvent, newSynapsesLco,
            &neuron->gid, sizeof(neuron_id_t), &tt, sizeof(spike_time_t));
    return newSynapsesLco;
}
