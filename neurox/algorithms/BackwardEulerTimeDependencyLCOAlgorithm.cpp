#include "neurox/algorithms/BackwardEulerTimeDependencyLCOAlgorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

DERIVED_CLASS_NAME::DERIVED_CLASS_NAME() {}

DERIVED_CLASS_NAME::~DERIVED_CLASS_NAME() {}

AlgorithmType DERIVED_CLASS_NAME::getType()
{
    return AlgorithmType::BackwardEulerTimeDependencyLCO;
}

void DERIVED_CLASS_NAME::Init() {}

void DERIVED_CLASS_NAME::Clear() {}

void DERIVED_CLASS_NAME::StepBegin(Branch*) {}

void DERIVED_CLASS_NAME::StepEnd(Branch*) {}

void DERIVED_CLASS_NAME::CommStepBegin(Branch*) {}

void DERIVED_CLASS_NAME::CommStepEnd(Branch*) {}

void DERIVED_CLASS_NAME::AfterSpike(Branch*) {}
