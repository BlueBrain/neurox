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
    neurox::Solver::HinesSolver::registerHpxActions();
    neurox::Misc::Statistics::registerHpxActions();

    //Init HPX
    if (hpx_init(&argc, &argv) != 0)
    {
        printf("HPX failed to initialize!\n");
        return 1;
    }
    printf("neurox started (localities: %d, threads/locality: %d.)\n", HPX_LOCALITIES, HPX_THREADS);

    //Run main
    int e = hpx_run(&neurox::main, argv, argc);

    //clean up
    hpx_finalize();
    return e;
}
