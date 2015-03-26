#include "neurox/neurox.h"

#include <algorithm>
#include <cstring>
#include <utility>

using namespace neurox;
using namespace neurox::solver;
using namespace neurox::algorithms;

Neuron::Neuron(neuron_id_t neuronId, floble_t APthreshold)
    : gid(neuronId), threshold(APthreshold), algorithmMetaData(nullptr) {
  this->synapsesTransmissionFlag = false;
  this->synapsesMutex = hpx_lco_sema_new(1);
  this->refractoryPeriod = 0;
  this->algorithmMetaData = AlgorithmMetaData::New(input_params->algorithm);
  assert(this->algorithmMetaData != nullptr);
  assert(
      TimeDependencyLCOAlgorithm::TimeDependencies::notificationIntervalRatio >
          0 &&
      TimeDependencyLCOAlgorithm::TimeDependencies::notificationIntervalRatio <=
          1);
  assert(CoreneuronAlgorithm::CommunicationBarrier::commStepSize %
             AllReduceAlgorithm::AllReducesInfo::reductionsPerCommStep ==
         0);
}

Neuron::~Neuron() {
  if (synapsesMutex != HPX_NULL) hpx_lco_delete_sync(synapsesMutex);
  for (Synapse*& s : synapses) delete s;
  delete algorithmMetaData;
}

Neuron::Synapse::Synapse(hpx_t branchAddr, floble_t minDelay,
                         hpx_t topBranchAddr, int destinationGid)
    : branchAddr(branchAddr), minDelay(minDelay), topBranchAddr(topBranchAddr) {
  const double& teps = TimeDependencyLCOAlgorithm::TimeDependencies::teps;
  const double& notifRatio =
      TimeDependencyLCOAlgorithm::TimeDependencies::notificationIntervalRatio;
  this->nextNotificationTime =
      input_params->tstart + teps + this->minDelay * notifRatio;
  this->previousSpikeLco = hpx_lco_future_new(0);
  hpx_lco_set_rsync(
      this->previousSpikeLco, 0,
      NULL);  // starts as set and will be reset when synapses happen
#ifndef NDEBUG
  this->destinationGid = destinationGid;
#endif
}

Neuron::Synapse::~Synapse() {
  if (previousSpikeLco != HPX_NULL) hpx_lco_delete_sync(previousSpikeLco);
}

size_t Neuron::GetSynapsesCount() {
  size_t size = -1;
  hpx_lco_sema_p(synapsesMutex);
  size = synapses.size();
  hpx_lco_sema_v_sync(synapsesMutex);
  return size;
}

void Neuron::AddSynapse(Synapse* syn) {
  hpx_lco_sema_p(synapsesMutex);
  synapses.push_back(syn);
  synapses.shrink_to_fit();
  hpx_lco_sema_v_sync(synapsesMutex);
}

// netcvode.cpp::static bool pscheck(...)
bool Neuron::CheckAPthresholdAndTransmissionFlag(floble_t v) {
  // can only spike if AP threshold has been reach and spikes havent already
  // been transmitted
  if (v > threshold) {
    if (synapsesTransmissionFlag == false) {
      synapsesTransmissionFlag = true;
      return true;
    }
  } else {
    synapsesTransmissionFlag = false;
  }
  return false;
}

hpx_t Neuron::SendSpikes(floble_t t)  // netcvode.cpp::PreSyn::send()
{
  const spike_time_t tt =
      (spike_time_t)t + 1e-10;  // Coreneuron logic, do not change!
#if !defined(NDEBUG)
  printf("== Neuron gid %d spiked at %.3f ms\n", this->gid, tt);
#endif

  if (synapses.size() == 0) return HPX_NULL;
  return algorithm->SendSpikes(this, tt, t);
}
