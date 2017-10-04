#include "neurox/neurox.h"
#include "coreneuron/nrnmpi/nrnmpi.h"

int main(int argc, char** argv) {
  neurox::RegisterHpxActions();
  neurox::Branch::RegisterHpxActions();
  neurox::tools::Statistics::RegisterHpxActions();
  neurox::tools::LoadBalancing::RegisterHpxActions();
  neurox::input::DataLoader::RegisterHpxActions();
  neurox::algorithms::AllreduceAlgorithm::AllReducesInfo::RegisterHpxActions();
  neurox::BranchCvodes::RegisterHpxActions();
#if !defined(NDEBUG)
  neurox::input::Debugger::RegisterHpxActions();
#endif

  if (hpx_init(&argc, &argv) != 0) {
    printf("HPX failed to initialize!\n");
    return 1;
  }

  // parse command line arguments
  neurox::input_params_ = new tools::CmdLineParser(argc, argv);

  /// all compute nodes load the data (mechs info is accessible to all)
  neurox::input::DataLoader::InitAndLoadCoreneuronData(argc, argv, false,
                                                       false);

  int e = hpx_run(&neurox::Main, NULL);
  hpx_finalize();
  return e;
}
