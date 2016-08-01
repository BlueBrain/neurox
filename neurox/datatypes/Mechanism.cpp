#include "neurox/Neurox.h"
#include <cstring>

#include "coreneuron/nrnoc/membdef.h" //Memb_func, BAMech

using namespace std;
using namespace Neurox;

Mechanism::Mechanism():
 type(-1), dataSize(-1), pdataSize(-1), dependenciesCount(-1), pntMap(0),
 isArtificial(0), dependencies(nullptr), isIon(0), conci(-1), conco(-1),
 charge(-1), sym(nullptr)
{};

Mechanism::Mechanism(const short int type, const short int dataSize, const short int pdataSize,
                     const char isArtificial, const char pntMap, const char isIon,
                     const short int symLength, const char * sym,
                     const short int dependenciesCount, const short int * dependencies):
    type(type), dataSize(dataSize), pdataSize(pdataSize),
    dependenciesCount(dependenciesCount),
    pntMap(pntMap), isArtificial(isArtificial),
    symLength(symLength), isIon(isIon),
    dependencies(nullptr), sym(nullptr),
    conci(-1), conco(-1), charge(-1)
{
    if (dependencies != nullptr && dependenciesCount>0)
    {
        this->dependencies = new  short int[dependenciesCount];
        std::memcpy(this->dependencies, dependencies, dependenciesCount*sizeof(short int));
    }

    if (sym != nullptr && symLength>0)
    {
        this->sym = new char[symLength+1];
        std::memcpy(this->sym, sym, symLength);
        this->sym[symLength] = '\0';
    }

    registerMechFunctions();
    //registerBAFunctions();
    if (isIon)
        registerIonicCharges();

#ifdef DEBUG
    printf("DEBUG Mechanism: type %d, dataSize %d, pdataSize %d, isArtificial %d,\n"
           "      pntMap %d, isIon %d, symLength %d, sym %s, dependenciesCount %d\n"
           "      conci %.2f, conco %.2f, charge %.2f\n",
           type, dataSize, pdataSize, isArtificial, pntMap, isIon, symLength,
           sym!=nullptr ? sym : "", dependenciesCount, conci, conco, charge);
#endif
};

void Mechanism::registerBAFunctions()
{
    //Copy Before-After functions
    //register_mech.c::hoc_reg_ba()
    for (int i=0; i< BEFORE_AFTER_SIZE; i++)
    {
        //BUG not implemented by CoreNeuron, undefined bam leads to SEGFAULT:
        this->BAfunctions[i] = nrn_threads[0].tbl[i]->bam->f;
    }
}

void Mechanism::registerMechFunctions()
{
    int type = this->type;
    //register functions //TODO will not work in more than 1 compute node
    this->membFunc.alloc = memb_func[type].alloc;
    this->membFunc.setdata_ = memb_func[type].setdata_;
    this->membFunc.destructor = memb_func[type].destructor;
    this->membFunc.current =  memb_func[type].current;
    this->membFunc.jacob = memb_func[type].jacob;
    this->membFunc.state = memb_func[type].state;
    this->membFunc.initialize = memb_func[type].initialize;
    this->membFunc.thread_cleanup_ = memb_func[type].thread_cleanup_;
    this->membFunc.thread_mem_init_ = memb_func[type].thread_mem_init_;
    this->membFunc.thread_size_ = memb_func[type].thread_size_;
    //look up tables are created and destroyed inside mod files, not accessible via coreneuron
    this->membFunc.thread_table_check_ = memb_func[type].thread_table_check_;
}

void Mechanism::registerIonicCharges()
{
    double ** ion_global_map = get_ion_global_map(); // added to membfunc.h and eion.c;
    conci = ion_global_map[type][0];
    conco = ion_global_map[type][1];
    charge = ion_global_map[type][2];
}

Mechanism::~Mechanism(){
    delete [] dependencies;
    delete [] sym;
};

hpx_action_t Mechanism::callModFunction = 0;
int Mechanism::callModFunction_handler(const int nargs, const void *args[], const size_t sizes[])
{
    neurox_hpx_pin(Branch); //we are in a Branch's memory space
    assert(nargs==2);

    /** nargs=2 where:
     * args[0] = Id of the mechanism to be called
     * args[1] = Id of the function to be called
     */
    short int mechType = *(short int*)args[0];
    Mechanism::ModFunction functionId = *(ModFunction*)args[1];
    Mechanism & mech = mechanisms[mechType];

    //call this function in all mechanisms that should run first
    //(ie the children on the tree of mechanisms dependencies)
    neurox_hpx_recursive_mechanism_sync(mechType, Mechanism::callModFunction,
                                        args[1], sizes[1]);

    Memb_list membList;
    NrnThread nrnThread;

    //Note:The Jacob updates D and nrn_cur updates RHS, so we need a mutex for compartments
    //The state function does not write to compartment, only reads, so no mutex needed

    Input::Coreneuron::DataLoader::fromHpxToCoreneuronDataStructs(local, membList, nrnThread, mechType);
    if (functionId<BEFORE_AFTER_SIZE)
        mech.BAfunctions[functionId](&nrnThread, &membList, mechType);
    else
    {
        switch(functionId)
        {
            case Mechanism::ModFunction::alloc:
                if (mech.membFunc.alloc)
                    mech.membFunc.alloc(membList.data, membList.pdata, mechType);
                break;
            case Mechanism::ModFunction::current:
                if (mech.membFunc.current)
                    mech.membFunc.current(&nrnThread, &membList, mechType);
                break;
            case Mechanism::ModFunction::state:
                if (mech.membFunc.state)
                    mech.membFunc.state(&nrnThread, &membList, mechType);
                break;
            case Mechanism::ModFunction::jacob:
                if (mech.membFunc.jacob)
                    mech.membFunc.jacob(&nrnThread, &membList, mechType);
                break;
            case Mechanism::ModFunction::initialize:
                if (mech.membFunc.initialize)
                    mech.membFunc.initialize(&nrnThread, &membList, mechType);
                break;
            case Mechanism::ModFunction::destructor:
                if (mech.membFunc.destructor)
                    mech.membFunc.destructor();
                break;
            case Mechanism::ModFunction::threadMemInit:
                if (mech.membFunc.thread_mem_init_)
                    mech.membFunc.thread_mem_init_(membList._thread);
                break;
            case Mechanism::ModFunction::threadTableCheck:
                if (mech.membFunc.thread_table_check_)
                    mech.membFunc.thread_table_check_
                        (0, membList.nodecount, membList.data, membList.pdata, membList._thread, &nrnThread, mechType);
                break;
            case Mechanism::ModFunction::threadCleanup:
                if (mech.membFunc.thread_cleanup_)
                    mech.membFunc.thread_cleanup_(membList._thread);
                break;
            default:
                printf("ERROR! Unknown function %d!!\n", functionId);
              }
            }

    neurox_hpx_unpin;
}

hpx_action_t Mechanism::callNetReceiveFunction = 0;
int Mechanism::callNetReceiveFunction_handler(const int nargs, const void *args[], const size_t sizes[])
{
    neurox_hpx_pin(Branch); //we are in a Branch's memory space
    assert(nargs==4);
    /** nargs=4 where:
     * args[0] = Id of the mechanism to be called
     * args[1] = function flag: NetReceiveInit (1) or NetReceive(0)
     * args[2] = actual time
     * args[3] = timestep size
     */
    short int mechType = *(short int*)args[0];
    const char * initFunction = (const char*) args[1];
    const double * t = (const double *) args[2];
    const double * dt = (const double *) args[3];
    Mechanism & mech = mechanisms[mechType];

    //call this function in all mechanisms that should run first
    //(ie the children on the tree of mechanisms dependencies)
    neurox_hpx_recursive_mechanism_sync(mechType, Mechanism::callModFunction,
                                        args[1], sizes[1], args[2], sizes[2], args[3], sizes[3]);

    Memb_list membList;
    NrnThread nrnThread;
    Point_process pp;
    while (!local->spikesQueue.empty() &&
           local->spikesQueue.top().deliveryTime < *t+*dt)
    {
        Spike spike = local->spikesQueue.top();
        if (spike.netcon->active)
        {
            short int mechType = spike.netcon->mechType;
            Input::Coreneuron::DataLoader::fromHpxToCoreneuronDataStructs(local, membList, nrnThread, mechType);
            pp._i_instance = spike.netcon->mechInstance;
            pp._presyn = NULL;
            pp._tid = -1;
            pp._type = spike.netcon->mechType;
            nrnThread._t = spike.deliveryTime;
            //see netcvode.cpp:433 (NetCon::deliver)
            //We have to pass NrnThread, MembList, and deliveryTime instead
            if (*initFunction)
            {
                    //mechanisms[spike.netcon->mechType].pnt_receive_init(&nrnThread, &membList, &pp, spike.netcon->args, 0);
            }
            else
            {
                    //mechanisms[spike.netcon->mechType].pnt_receive(&nrnThread, &membList, &pp, spike.netcon->args, 0);
            }
        }
        local->spikesQueue.pop();
    }
    neurox_hpx_unpin;
}

void Mechanism::registerHpxActions()
{
    neurox_hpx_register_action(2, Mechanism::callModFunction);
    neurox_hpx_register_action(2, Mechanism::callNetReceiveFunction);
}
