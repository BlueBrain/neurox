#pragma once

#include "neurox/neurox.h"

using namespace std;

//TODO should be part of Neurox::Input::Hdf5 namespace

class BlueConfigDataLoader //: IDataLoader
{
  public:
    BlueConfigDataLoader()=delete;
    ~BlueConfigDataLoader()=delete;

    static void loadData(int argc, char ** argv){};

  private:
};
