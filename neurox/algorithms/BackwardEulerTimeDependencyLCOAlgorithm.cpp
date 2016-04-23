#include "neurox/algorithms/BackwardEulerTimeDependencyLCOAlgorithm.h"

using namespace neurox;
using namespace neurox::algorithms;

DERIVED_CLASS_NAME::DERIVED_CLASS_NAME() {}

DERIVED_CLASS_NAME::~DERIVED_CLASS_NAME() {}

AlgorithmType DERIVED_CLASS_NAME::getType()
{
    return AlgorithmType::BackwardEulerTimeDependencyLCO;
}

void DERIVED_CLASS_NAME::init(Branch*) {}

void DERIVED_CLASS_NAME::clear(Branch*) {}

void DERIVED_CLASS_NAME::stepBegin(Branch*) {}

void DERIVED_CLASS_NAME::stepEnd(Branch*) {}

void DERIVED_CLASS_NAME::commStepBegin(Branch*) {}

void DERIVED_CLASS_NAME::commStepEnd(Branch*) {}

void DERIVED_CLASS_NAME::afterSpike(Branch*) {}
