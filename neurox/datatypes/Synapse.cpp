#include "neurox/neurox.h"
#include <cstring>

using namespace std;

Synapse::Synapse(const double weight, const double delay, const hpx_t target)
   :weight(weight), delay(delay), target(target) {};

Synapse::~Synapse(){};
