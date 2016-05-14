#pragma once

#include <vector>

using namespace std;

/**
 * @brief The MechanismsDependencies class
 * Compressed Row Storage of dependencies between mechanisms
 */
class MechanismsDependencies
{
  public:
    MechanismsDependencies(){};
    ~MechanismsDependencies(){};

    vector<short int> count;
    vector<int> offsets;
    vector<int> ids;


  private:
};
