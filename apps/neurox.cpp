#include "neurox/neurox.h"
#include "coreneuron/nrnmpi/nrnmpi.h"

/// if  ADDITIONAL_MECHS is defined, this is automatically generated  by CMake in build/coreneuron/mod_func.c
extern "C" {extern void modl_reg(void);}

/// Declare an empty function if Neurodamus mechanisms are not used, otherwise register them in mechs/cfile/mod_func.c
#ifndef ADDITIONAL_MECHS
void modl_reg() {}
#endif

int main(int argc, char** argv)
{
    neurox::RegisterHpxActions();
    neurox::Branch::RegisterHpxActions();
    neurox::tools::Statistics::RegisterHpxActions();
    neurox::tools::LoadBalancing::RegisterHpxActions();
    neurox::input::DataLoader::RegisterHpxActions();
    neurox::algorithms::AllReduceAlgorithm::AllReducesInfo::RegisterHpxActions();
#if !defined(NDEBUG)
    neurox::input::Debugger::RegisterHpxActions();
#endif

    if (hpx_init(&argc, &argv) != 0)
    {
        printf("HPX failed to initialize!\n");
        return 1;
    }

    //parse command line arguments
    neurox::inputParams = new tools::CmdLineParser(argc, argv);

    ///all compute nodes load the data (mechs info is accessible to all)
    neurox::input::DataLoader::InitAndLoadCoreneuronData(argc, argv, false, false);

    int e = hpx_run(&neurox::Main, NULL);
    hpx_finalize();
    return e;
}
