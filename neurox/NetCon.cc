#include <cstring>
#include "neurox/neurox.h"

using namespace neurox;

NetConX::NetConX() {}

NetConX::NetConX(int mechType, offset_t mechInstance, floble_t delay,
                 offset_t weightIndex, unsigned short weightsCount, bool active)
    : mechType(mechType),
      weightsCount(weightsCount),
      mechInstance(mechInstance),
      delay(delay),
      weightIndex(weightIndex),
      active(active) {}

NetConX::~NetConX() {}

void NetConX::Deliver(floble_t tt, Branch* branch) {
  if (this->active)
    GetMechanismFromType(mechType)->CallModFunction(
        branch, Mechanism::ModFunctions::kNetReceive, this, tt);
}
