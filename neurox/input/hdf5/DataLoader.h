#pragma once

#include "neurox/neurox.h"

using namespace std;

//TODO should be part of Neurox::Input::Hdf5 namespace

class Hdf5DataLoader //: IDataLoader
{
  public:
    Hdf5DataLoader()=delete;
    ~Hdf5DataLoader()=delete;

    static void loadData(int argc, char ** argv){};

  private:
};
