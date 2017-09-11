#include "neurox/algorithms/cvodes_algorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

CvodesAlgorithm::CvodesAlgorithm() {
}

CvodesAlgorithm::~CvodesAlgorithm() {}

const AlgorithmType CvodesAlgorithm::GetType() {
  return AlgorithmType::kCvodes;
}

const char* CvodesAlgorithm::GetTypeString() {
  return "CVODES";
}

void CvodesAlgorithm::Init() {
}

void CvodesAlgorithm::Clear() {
}

double CvodesAlgorithm::Launch() {
}

void CvodesAlgorithm::StepBegin(Branch*) {}

void CvodesAlgorithm::StepEnd(Branch* b, hpx_t spikesLco) {
}

void CvodesAlgorithm::Run(Branch* b, const void* args)
{
    
}

hpx_t CvodesAlgorithm::SendSpikes(Neuron* b, double tt, double) {
}
