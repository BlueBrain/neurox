#include "neurox/neurox.h"
#include <cstring>

#include "coreneuron/nrnoc/multicore.h" //nrn_threads

Brain * brain = nullptr; //global variable (defined in neurox.h)

Brain::~Brain()
{
    delete [] mechsTypes;
}

Brain::Brain(const int neuronsCount,
             const hpx_t neuronsAddr, const Mechanism * mechanisms,
             const size_t mechsTypesCount, const int * mechDependencies)
    : neuronsCount(neuronsCount), neuronsAddr(neuronsAddr), mechsTypesCount(mechsTypesCount)
{
    //add mechanisms information
    int offset=0;
    this->mechsTypes = new Mechanism[mechsTypesCount];
    for (int m=0; m<mechsTypesCount; m++)
    {
        this->mechsTypes[m]=Mechanism(m, mechanisms[m].dataSize, mechanisms[m].pdataSize,
                                      mechanisms[m].dependenciesCount, mechanisms[m].pntMap,
                                      mechanisms[m].isArtificial, &mechDependencies[offset],
                                      mechanisms[m].isIon, mechanisms[m].conci, mechanisms[m].conco, mechanisms[m].charge //for ions
                                      );
        offset += mechanisms[m].dependenciesCount;
    }

    //TODO will work on single node only
    //Copy function pointers for BA calls
    beforeAfterFunction = new (mod_f_t[BEFORE_AFTER_SIZE])(mechsTypesCount);
    for (int m=0; m<mechsTypesCount; m++)
        for (int ba=0; ba<BEFORE_AFTER_SIZE; ba++)
            beforeAfterFunction[m][ba] = 0 (mod_f_t);

    for (int nt=0; nt<nrn_nthread; nt++)
        for (NrnThreadBAList *tbl = nt->tbl[bat]; tbl; tbl = tbl->next)
        {
            mod_f_t f = tbl->bam->f;
            int type = tbl->bam->type;
            beforeAfterFunction[type][bat]=f;
        }
};


int Brain::callBeforeAfterMethod(const int functionId)
{
    for (int m=0; m<mechsTypesCount; m++)
        if (beforeAfterFunction[m][functionId])
            beforeAfterFunction[m][functionId](NULL, NULL, m);
}

hpx_action_t Brain::clear = 0;
int Brain::clear_handler()
{
    neurox_hpx_pin(Brain);
    delete brain;
    neurox_hpx_unpin ;
}

hpx_action_t Brain::finitialize = 0;
int Brain::finitialize_handler()
{
    neurox_hpx_pin(Brain)

    //set up by finitialize.c:nrn_finitialize() -> fadvance_core.c:dt2thread()
    double cj = inputParams->secondorder ? 2.0/inputParams->dt : 1.0/inputParams->dt;
    hpx_par_for_sync(
       [&] (int i, void*) { hpx_call_sync(local->getNeuronAddr(i), Neuron::setCj);},
       0, local->neuronsCount, NULL, cj);

    //TODO detect min delay (must be divisable by dt)
    //in netpar.c: usable_mindelay_ = floor(mindelay_ * dt1_ + 1e-9) * dt;
    //probably should be solver-depentend eg Solver::init(...)
    //TODO missing nrn_thread_table_check()

    //set up by finitialize.c:nrn_finitialize(): if (setv)
    double v = inputParams->voltage;
    hpx_par_for_sync(
       [&] (int i, void*) { hpx_call_sync(local->getNeuronAddr(i), Neuron::setV);},
       0, local->neuronsCount, NULL, v);

    callBeforeAfterMethod(BEFORE_INITIAL);

    // the INITIAL blocks are ordered so that mechanisms that write
    // concentrations are after ions and before mechanisms that read
    // concentrations.
    for (int m=0; m<local->mechsTypesCount; m++)
        if (local->mechsTypes[m].initialize)
            local->mechsTypes[m].initialize(NULL, NULL, m);

    callBeforeAfterMethod(AFTER_INITIAL);

    //now add the axial currents

    //finitialize.c:nrn_finitialize()->set_tree_matrix_minimal->nrn_rhs
    callBeforeAfterMethod(BEFORE_BREAKPOINT);

    //note that CAP has no current
    for (int m=0; m<local->mechsTypesCount; m++)
        if (local->mechsTypes[m].current)
            local->mechsTypes[m].current(NULL, NULL, m);

    //finitialize.c:nrn_finitialize()->set_tree_matrix_minimal->nrn_lhs (treeset_core.c)
    // now the internal axial currents.
    //The extracellular mechanism contribution is already done.
    //	rhs += ai_j*(vi_j - vi)
    hpx_par_for_sync(
       [&] (int i, void*) { hpx_call_sync(local->getNeuronAddr(i), Neuron::setupMatrixRHS);},
       0, local->neuronsCount, NULL);

    // calculate left hand side of
    //cm*dvm/dt = -i(vm) + is(vi) + ai_j*(vi_j - vi)
    //cx*dvx/dt - cm*dvm/dt = -gx*(vx - ex) + i(vm) + ax_j*(vx_j - vx)
    //with a matrix so that the solution is of the form [dvm+dvx,dvx] on the right
    //hand side after solving.
    //This is a common operation for fixed step, cvode, and daspk methods
    // note that CAP has no jacob
    for (int m=0; m<local->mechsTypesCount; m++)
        if (local->mechsTypes[m].jacob)
            local->mechsTypes[m].jacob(NULL, NULL, m);

    //finitialize.c:nrn_finitialize()->set_tree_matrix_minimal->nrn_rhs (treeset_core.c)
    //now the cap current can be computed because any change to cm
    //by another model has taken effect. note, the first is CAP
    local->mechsTypes[CAP].nrn_cap_jacob(NULL, NULL);

    //now add the axial currents
    hpx_par_for_sync(
       [&] (int i, void*) { hpx_call_sync(local->getNeuronAddr(i), Neuron::setupMatrixLHS);},
       0, local->neuronsCount, NULL);

    //TODO: missing nrn_spike_exchange step
    neurox_hpx_unpin;
}

hpx_action_t Brain::init = 0;
int Brain::init_handler(const int neuronsCount,
           const hpx_t neuronsAddr, const Mechanism * mechanisms,
           const size_t mechanismsCount, const int * mechDependencies)
{
    neurox_hpx_pin(Brain);

    //initialize global variable brain
    brain = new Brain(neuronsCount, neuronsAddr, mechanisms, mechanismsCount, mechDependencies);

    //clean up
    delete [] mechanisms;
    delete [] mechDependencies;

    neurox_hpx_unpin;
}

hpx_action_t Brain::solve = 0;
int Brain::solve()
{
    neurox_hpx_pin(Brain);

    neurox_hpx_unpin;
}

void Brain::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  init, init_handler, HPX_INT, HPX_ADDR, HPX_POINTER, HPX_SIZE_T, HPX_POINTER);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  finitialize, finitialize_handler);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  setV, setV_handler);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  clear, clear_handler);
}

