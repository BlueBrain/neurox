#include "neurox/neurox.h"

#ifdef ADDITIONAL_MECHS
extern "C" {extern void modl_reg(void);} ///> generated by CMake in build/coreneuron/mod_func.c
#else
void modl_reg() {} ///No additional mechs, dont register any external
#endif

int main(int argc, char** argv)
{
    printf("neurox::registerHpxActions...\n");
    neurox::registerHpxActions();
    neurox::Branch::registerHpxActions();
    neurox::Mechanism::registerHpxActions();
    neurox::Misc::Statistics::registerHpxActions();
#if defined(CORENEURON_H)
    neurox::Input::Coreneuron::DataLoader::registerHpxActions();
#if !defined(NDEBUG)
    neurox::Input::Coreneuron::Debugger::registerHpxActions();
#endif
#endif

    printf("neurox::hpx_init...\n");
    if (hpx_init(&argc, &argv) != 0)
    {
        printf("HPX failed to initialize!\n");
        return 1;
    }

    //MPI parallel data loading
    neurox::Input::Coreneuron::DataLoader::initAndLoadCoreneuronData(argc, argv);

    printf("neurox::main...\n");
    int e = hpx_run(&neurox::main, NULL, &argv, argc);
    hpx_finalize();
    return e;
}
