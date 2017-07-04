#include "neurox/neurox.h"

#ifdef ADDITIONAL_MECHS
extern "C" {extern void modl_reg(void);} ///> generated by CMake in build/coreneuron/mod_func.c
#else
void modl_reg() {} ///No additional mechs, dont register any external
#endif

int main(int argc, char** argv)
{
    neurox::registerHpxActions();
    neurox::Branch::registerHpxActions();
    neurox::Neuron::registerHpxActions();
    neurox::Tools::Statistics::registerHpxActions();
    neurox::Tools::LoadBalancing::registerHpxActions();
    neurox::Input::DataLoader::registerHpxActions();
#if !defined(NDEBUG)
    neurox::Input::Debugger::registerHpxActions();
#endif

    if (hpx_init(&argc, &argv) != 0)
    {
        printf("HPX failed to initialize!\n");
        return 1;
    }

    neurox::inputParams = new Tools::CmdLineParser(argc, argv);

    if (inputParams->parallelDataLoading) //coreneuron data loading
        neurox::Input::DataLoader::initAndLoadCoreneuronData(argc, argv);
    else if (hpx_get_my_rank()==0) //one loads all neurons and spreads them
        neurox::Input::DataLoader::initAndLoadCoreneuronData(argc, argv);

    int e = hpx_run(&neurox::main, NULL);
    hpx_finalize();
    return e;
}
