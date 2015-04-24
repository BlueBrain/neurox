#include "neurox/algorithms/BackwardEulerSlidingTimeWindowAlgorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

hpx_t* DERIVED_CLASS_NAME::allReduces = nullptr;

DERIVED_CLASS_NAME::DERIVED_CLASS_NAME() {}

DERIVED_CLASS_NAME::~DERIVED_CLASS_NAME() {}

AlgorithmType DERIVED_CLASS_NAME::getType()
{
    return AlgorithmType::BackwardEulerSlidingTimeWindow;
}

void DERIVED_CLASS_NAME::Init() {
    BackwardEulerAllReduceAlgorithm::SubscribeAllReduces(
                DERIVED_CLASS_NAME::allReduces,
                DERIVED_CLASS_NAME::allReducesCount);
}

void DERIVED_CLASS_NAME::Clear() {
    BackwardEulerAllReduceAlgorithm::UnsubscribeAllReduces(
                DERIVED_CLASS_NAME::allReduces,
                DERIVED_CLASS_NAME::allReducesCount);
}

void DERIVED_CLASS_NAME::StepBegin(Branch*) {}

void DERIVED_CLASS_NAME::StepEnd(Branch*) {}

void DERIVED_CLASS_NAME::CommStepBegin(Branch*) {}

void DERIVED_CLASS_NAME::CommStepEnd(Branch*) {}

void DERIVED_CLASS_NAME::AfterSpike(Branch*) {}
