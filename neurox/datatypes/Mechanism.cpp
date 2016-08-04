#include "neurox/Neurox.h"
#include <cstring>

#include "coreneuron/nrnoc/membdef.h" //Memb_func, BAMech
#include "coreneuron/nrnoc/multicore.h" //NrnThread
#include "coreneuron/nrnoc/nrnoc_ml.h" //Memb_list
#include "coreneuron/nrniv/netcon.h"

using namespace std;
using namespace Neurox;

Mechanism::Mechanism(const int type, const short int dataSize, const short int pdataSize,
                     const char isArtificial, const char pntMap, const char isIon,
                     const short int symLength, const char * sym,
                     const char isTopMechanism,
                     const short int childrenCount, const int * children):
    type(type), dataSize(dataSize), pdataSize(pdataSize),
    childrenCount(childrenCount),
    pntMap(pntMap), isArtificial(isArtificial),
    isTopMechanism(isTopMechanism),
    symLength(symLength), isIon(isIon),
    children(nullptr), sym(nullptr),
    conci(-1), conco(-1), charge(-1) //registered in registerIonicCharges()
{
    if (children != nullptr)
    {
        assert(childrenCount>0);
        this->children = new int[childrenCount];
        std::memcpy(this->children, children, childrenCount*sizeof(int));
    }

    if (sym != nullptr)
    {
        assert(symLength>0);
        this->sym = new char[symLength+1];
        std::memcpy(this->sym, sym, symLength);
        this->sym[symLength] = '\0';
    }

    registerMechFunctions();
    registerBAFunctions();
    if (isIon)
        registerIonicCharges();

#ifdef DEBUG
    printf("DEBUG Mechanism: type %d, dataSize %d, pdataSize %d, isArtificial %d,\n"
           "      pntMap %d, isIon %d, symLength %d, sym %s, childrenCount %d\n"
           "      conci %.2f, conco %.2f, charge %.2f\n",
           this->type, this->dataSize, this->pdataSize, this->isArtificial,
           this->pntMap, this->isIon, this->symLength, this->sym, this->childrenCount,
           this->conci, this->conco, this->charge);
#endif
};

void Mechanism::registerBAFunctions()
{
    //Copy Before-After functions
    //register_mech.c::hoc_reg_ba()
    for (int i=0; i< BEFORE_AFTER_SIZE; i++)
    {
        //BUG not implemented by CoreNeuron, undefined bam leads to SEGFAULT:
        //this->BAfunctions[i] = nrn_threads[0].tbl[i]->bam->f;
        this->BAfunctions[i] = NULL;
    }
}

void Mechanism::registerMechFunctions()
{
    if (this->sym && strcmp("capacitance", this->sym))
    {
        return;
    }

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

Mechanism::Mechanism(){}

Mechanism::~Mechanism(){
    delete [] sym;
    delete [] children;
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
    int mechType = *(int*)args[0];
    Mechanism::ModFunction functionId = *(ModFunction*)args[1];
    Mechanism & mech = getMechanismFromType(mechType);

    Memb_list membList;
    NrnThread nrnThread;

    //Note:The Jacob updates D and nrn_cur updates RHS, so we need a mutex for compartments
    //The state function does not write to compartment, only reads, so no mutex needed

    Input::Coreneuron::DataLoader::fromHpxToCoreneuronDataStructs(local, membList, nrnThread, mechType);
    if (functionId<BEFORE_AFTER_SIZE)
        if (mech.BAfunctions[functionId])
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

    //call this function in all mechanisms that depend on this one
    //(ie the children on the tree of mechanisms dependencies)
    neurox_hpx_recursive_mechanism_sync(mechType, Mechanism::callModFunction,
                                        args[1], sizes[1]);

    neurox_hpx_unpin;
}

void Mechanism::callNetReceiveFunction(
        const void * branch, const Spike * spike,
        const char isInitFunction, const double t, const double dt)
{
    Memb_list membList;
    NrnThread nrnThread;
    Input::Coreneuron::DataLoader::fromHpxToCoreneuronDataStructs
            ((Branch*) branch, membList, nrnThread, type);

    Point_process pp;
    pp._i_instance = spike->netcon->mechInstance;
    pp._presyn = NULL;
    pp._tid = -1;
    pp._type = spike->netcon->mechType;
    nrnThread._t = spike->deliveryTime;
    //see netcvode.cpp:433 (NetCon::deliver)
    //We have to pass NrnThread, MembList, and deliveryTime instead
    if (isInitFunction)
    {
        //mechanisms[spike.netcon->mechType].pnt_receive_init(&nrnThread, &membList, &pp, spike.netcon->args, 0);
    }
    else
    {
        //mechanisms[spike.netcon->mechType].pnt_receive(&nrnThread, &membList, &pp, spike.netcon->args, 0);
    }
}

void Mechanism::registerHpxActions()
{
    neurox_hpx_register_action(2, Mechanism::callModFunction);
}
