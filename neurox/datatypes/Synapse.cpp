#include "neurox/neurox.h"
#include <cstring>

using namespace std;

Synapse::Synapse(const double weight, const double delay, const hpx_t target, const int mechOffset, const int mechInstance)
   :weight(weight), delay(delay), target(target), mechOffset(mechOffset), mechInstance(mechInstance) {};

Synapse::~Synapse(){};
