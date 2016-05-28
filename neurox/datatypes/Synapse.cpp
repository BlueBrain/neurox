#include "neurox/neurox.h"
#include <cstring>

using namespace std;

Synapse::Synapse(const double weight, const double delay, const int mechType, const int mechInstance)
    :weight(weight), delay(delay), mechType(mechType), mechInstance(mechInstance) {};

Synapse::~Synapse(){};
