#include "neurox/neurox.h"

#ifdef ADDITIONAL_MECHS
extern "C" {extern void modl_reg(void);} ///> generated by CMake in build/coreneuron/mod_func.c
#else
void modl_reg() {} ///No additional mechs, dont register any external
#endif

int main(int argc, char** argv)
{
    //register HPX methods
    neurox::registerHpxActions();
    neurox::Branch::registerHpxActions();
    neurox::Mechanism::registerHpxActions();
    neurox::Misc::Statistics::registerHpxActions();
    neurox::Solver::HinesSolver::registerHpxActions();
    neurox::Input::Coreneuron::DataLoader::registerHpxActions();

    //Init HPX
    if (hpx_init(&argc, &argv) != 0)
    {
        printf("HPX failed to initialize!\n");
        return 1;
    }

    int e = hpx_run(&neurox::main, argv, argc);
    hpx_finalize();
    return e;
}
