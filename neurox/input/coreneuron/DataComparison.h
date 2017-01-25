#pragma once

#include "neurox/neurox.h"

#include <map>
#include <vector>
#include <memory>

using namespace std;
namespace neurox {
namespace Input {
namespace Coreneuron {

/**
 * @brief The DataComparison class
 * Compares CoreNeuron and HPX-based data structures
 */
class DataComparison
{
  public:
    DataComparison()=delete;
    ~DataComparison()=delete;

    static void coreNeuronFinitialize();
    static void coreNeuronFinitialize2();
    static void compareDataStructuresWithCoreNeuron(Branch * branch);

  private:
};

}; //Coreneuron
}; //Input
}; //NeuroX
