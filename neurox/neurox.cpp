#include "neurox/neurox.h"
#include <cstring>
#include <map>

#include "nrniv/nrniv_decl.h"
#include "nrniv/nrn_stats.h"

namespace neurox
{

int neuronsCount=-1;
hpx_t neuronsAddr = HPX_NULL;

int mechanismsCount=-1;
extern int * mechanismsMap = nullptr;
neurox::Mechanism ** mechanisms = nullptr;
Input::InputParams * inputParams = nullptr;

hpx_t getNeuronAddr(int i) {
    return hpx_addr_add(neuronsAddr, sizeof(Branch)*i, sizeof(Branch));
}

Mechanism * getMechanismFromType(int type) {
    assert(mechanismsMap[type]!=-1);
    return mechanisms[mechanismsMap[type]];
}

hpx_action_t setNeurons = 0;
int setNeurons_handler(const int nargs, const void *args[], const size_t[])
{
    /** nargs=2 where
     * args[0] = neuronsCount
     * args[1] = neuronsAddr
     */

    neurox_hpx_pin(uint64_t);
    assert(nargs==2);
    neuronsCount = *(int*)args[0];
    neuronsAddr = *(hpx_t*)args[1];
    neurox_hpx_unpin;
}

hpx_action_t setInputParams = 0;
int setInputParams_handler(const Input::InputParams * data, const size_t)
{
    neurox_hpx_pin(uint64_t);
    if (inputParams!=nullptr)
        delete [] neurox::inputParams;

    inputParams = new neurox::Input::InputParams();
    memcpy(neurox::inputParams, data, sizeof(neurox::Input::InputParams));
    neurox_hpx_unpin;
}

hpx_action_t setMechanisms = 0;
int setMechanisms_handler(const int nargs, const void *args[], const size_t sizes[])
{
    /**
     * nargs=3 where:
     * args[0] = array of all mechanisms info
     * args[1] = array of all mechanisms successors (children in mechanisms tree)
     * args[2] = array of all mechanisms names (sym)
     */

    neurox_hpx_pin(uint64_t);
    assert(nargs==3);
    mechanismsCount = sizes[0]/sizeof(Mechanism);
    mechanisms = new Mechanism*[mechanismsCount];

    int offsetSuccessors=0;
    int offsetSym=0;
    int maxMechType=-1;
    for (int m=0; m<mechanismsCount; m++)
    {
        Mechanism & mech = ((Mechanism*) args[0])[m];
        int * successorsIds = mech.successorsCount == 0 ? nullptr : &((int*) args[1])[offsetSuccessors];
        char * sym = mech.symLength == 0 ? nullptr : &((char*) args[2])[offsetSym];
        mechanisms[m] = new Mechanism(
                    mech.type, mech.dataSize, mech.pdataSize,
                    mech.isArtificial, mech.pntMap, mech.isIon,
                    mech.symLength, sym,
                    mech.dependenciesCount, mech.successorsCount, successorsIds);
        offsetSuccessors +=  mech.successorsCount;
        offsetSym += mech.symLength;
        if (mech.type > maxMechType)
            maxMechType = mech.type;
    }

    //initializes map of mechanisms ids to offset
    mechanismsMap = new int[maxMechType+1];
    for (int i=0; i<maxMechType+1; i++)
        mechanismsMap[i]=-1;
    for (int m=0; m<mechanismsCount; m++)
        mechanismsMap[mechanisms[m]->type]=m;
    neurox_hpx_unpin;
}

hpx_action_t main = 0;
static int main_handler(char **argv, size_t argc)
{
    //parse command line arguments and broadcasts it to other localities
    Input::InputParams inputParams(argc, argv);
    printf("neurox started (localities: %d, threads/locality: %d)\n", HPX_LOCALITIES, HPX_THREADS);
    hpx_bcast_rsync(neurox::setInputParams, &inputParams, sizeof (Input::InputParams));

    // many steps with large dt so that cells start at their resting potential
    assert(neurox::inputParams->forwardSkip == 0); //not supported yet

    //reads morphology data
    printf("neurox::Input::Coreneuron::DataLoader::loadData...\n");
    neurox::Input::Coreneuron::DataLoader::loadData(argc, argv);
    printf("neurox::Misc::Statistics::printSimulationSize...\n", neuronsCount);
    Misc::Statistics::printSimulationSize();
    printf("neurox::Misc::Statistics::printMechanismsDistribution...\n", neuronsCount);
    Misc::Statistics::printMechanismsDistribution();

    printf("neurox::BackwardEuler::initNeuronTreeLCO...\n");
    hpx_par_for_sync( [&] (int i, void*) -> int
     {  return hpx_call_sync(getNeuronAddr(i), Branch::initNeuronTreeLCO, NULL, 0, HPX_NULL, 0);
     }, 0, neuronsCount, NULL);

    printf("neurox::BackwardEuler::finitialize...\n");
    hpx_par_for_sync( [&] (int i, void*) -> int
    {  return hpx_call_sync(getNeuronAddr(i), Branch::finitialize, NULL, 0, HPX_NULL, 0);
    }, 0, neuronsCount, NULL);

    printf("neurox::BackwardEuler::solver...\n");
    hpx_t mainLCO = hpx_lco_and_new(neuronsCount);
    hpx_time_t now = hpx_time_now();
    for (int i=0; i<neuronsCount;i++)
        hpx_call(getNeuronAddr(i), Branch::backwardEuler, mainLCO);
    hpx_lco_wait(mainLCO);
    double elapsed = hpx_time_elapsed_ms(now)/1e3;

    printf("neurox finished (solver time: %.2f secs).\n", elapsed);
    hpx_exit(HPX_SUCCESS);
}

hpx_action_t clear = 0;
int clear_handler()
{
    delete [] mechanisms;
}

void registerHpxActions()
{
    neurox_hpx_register_action(1,neurox::main);
    neurox_hpx_register_action(0,neurox::clear);
    neurox_hpx_register_action(2,neurox::setNeurons);
    neurox_hpx_register_action(1,neurox::setInputParams);
    neurox_hpx_register_action(2,neurox::setMechanisms);
}
};
