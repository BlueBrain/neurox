#include "neurox/neurox.h"

#ifdef ADDITIONAL_MECHS
extern "C" {extern void modl_reg(void);} ///> generated by CMake in build/coreneuron/mod_func.c
#else
void modl_reg() {} ///No additional mechs, dont register any external
#endif

int Main(int argc, char** argv)
{
    neurox::RegisterHpxActions();
    neurox::Branch::registerHpxActions();
    neurox::Neuron::registerHpxActions();
    neurox::tools::Statistics::RegisterHpxActions();
    neurox::tools::LoadBalancing::RegisterHpxActions();
    neurox::input::DataLoader::registerHpxActions();
#if !defined(NDEBUG)
    neurox::input::Debugger::RegisterHpxActions();
#endif

    if (hpx_init(&argc, &argv) != 0)
    {
        printf("HPX failed to initialize!\n");
        return 1;
    }

    neurox::inputParams = new tools::CmdLineParser(argc, argv);

    if (inputParams->parallelDataLoading) //coreneuron data loading
        neurox::input::DataLoader::InitAndLoadCoreneuronData(argc, argv);
    else if (hpx_get_my_rank()==0) //one loads all neurons and spreads them
        neurox::input::DataLoader::InitAndLoadCoreneuronData(argc, argv);

    int e = hpx_run(&neurox::Main, NULL);
    hpx_finalize();
    return e;
}
