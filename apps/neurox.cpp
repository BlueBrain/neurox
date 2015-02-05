#include "neurox/neurox.h"
#include "coreneuron/nrnmpi/nrnmpi.h"

#ifdef ADDITIONAL_MECHS
extern "C" {extern void modl_reg(void);} ///> generated by CMake in build/coreneuron/mod_func.c
#else
void modl_reg() {} ///No additional mechs, dont register any external
#endif

int main(int argc, char** argv)
{
    neurox::RegisterHpxActions();
    neurox::Branch::registerHpxActions();
    neurox::tools::Statistics::RegisterHpxActions();
    neurox::tools::LoadBalancing::RegisterHpxActions();
    neurox::input::DataLoader::registerHpxActions();
    neurox::algorithms::AllReduceAlgorithm::AllReducesInfo::registerHpxActions();
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

    if (inputParams->branchingDepth>0 && hpx_get_num_ranks()>1)
    {
        //Branches will be split across several localities.
        //To load and instantiate all mechs in all nodes:
        //1. load all morphologies on all nodes
        //2. get and initialize the mechanisms for those neurons
        //3. clean neurons data and keep mechs
        //4. perform regular parallel loading of neurons
        //(This is because mechs have static variables that can't be accessed
        //by or communicated to other cpus, so mechs info can't be serialized)

        neurox::input::DataLoader::InitAndLoadCoreneuronData(argc, argv, true);
        hpx_run(&neurox::InitMechanismsAndQuit, NULL);
    }

    neurox::input::DataLoader::InitAndLoadCoreneuronData(argc, argv, false);
    int e = hpx_run(&neurox::Main, NULL);
    hpx_finalize();
    return e;
}
