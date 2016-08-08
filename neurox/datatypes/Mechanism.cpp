#include "neurox/Neurox.h"
#include <cstring>

#include "coreneuron/nrnoc/membfunc.h" //mod_f_t
#include "coreneuron/nrnoc/membdef.h" //Memb_func, BAMech
#include "coreneuron/nrnoc/multicore.h" //NrnThread
#include "coreneuron/nrnoc/nrnoc_ml.h" //Memb_list
#include "coreneuron/nrniv/netcon.h"

using namespace std;
using namespace NeuroX;

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
    conci(-1), conco(-1), charge(-1)
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

    if (sym==nullptr && children == nullptr) //not constructing, just serialized for transmission
        return;

    disableMechFunctions();
    if (this->type == CAP) //capacitance: capac.c
        registerCapacitance();
    else if (this->isIon)  //ion: eion.c
        registerIon();
    else //general mechanism: mod file
        registerModMechanism();

    //registerBeforeAfterFunctions();

#ifdef DEBUG
    printf("DEBUG Mechanism: type %d, dataSize %d, pdataSize %d, isArtificial %d,\n"
           "      pntMap %d, isIon %d, symLength %d, sym %s, childrenCount %d\n"
           "      conci %.2f, conco %.2f, charge %.2f\n",
           this->type, this->dataSize, this->pdataSize, this->isArtificial,
           this->pntMap, this->isIon, this->symLength, this->sym, this->childrenCount,
           this->conci, this->conco, this->charge);
#endif
};

void Mechanism::registerBeforeAfterFunctions()
{
    //Copy Before-After functions
    //register_mech.c::hoc_reg_ba()
    for (int i=0; i< BEFORE_AFTER_SIZE; i++)
        this->BAfunctions[i] = nrn_threads[0].tbl[i]->bam->f;
}

void Mechanism::disableMechFunctions()
{
    for (int i=0; i< BEFORE_AFTER_SIZE; i++)
        this->BAfunctions[i] = NULL;

    this->membFunc.alloc = NULL;
    this->membFunc.current = NULL;
    this->membFunc.state = NULL;
    this->membFunc.jacob = NULL;
    this->membFunc.initialize = NULL;
    this->membFunc.destructor = NULL;
    this->membFunc.thread_mem_init_ = NULL;
    this->membFunc.thread_cleanup_ = NULL;
    this->membFunc.thread_table_check_ = NULL;
    this->membFunc.setdata_ = NULL;

    this->membFunc.thread_size_ = -1;
    this->membFunc.is_point = -1;
}

void Mechanism::registerModMechanism()
{
    int type = this->type;
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

//from coreneuron/nrnoc/capac.c
extern void cap_alloc(double*, int*, int type);
extern void cap_init(struct NrnThread*, Memb_list*, int);
extern void nrn_capacity_current(struct NrnThread*, Memb_list*, int);
extern void nrn_cap_jacob(struct NrnThread*, Memb_list*, int);

void Mechanism::registerCapacitance()
{
    assert(this->sym && strcmp("capacitance", this->sym)==0);
    this->membFunc.alloc = cap_alloc;
    this->membFunc.initialize = cap_init;
    this->membFunc.current = nrn_capacity_current;
    this->membFunc.jacob = nrn_cap_jacob;
    //this->membFunc.alloc = Capacitance::cap_alloc;
    //this->membFunc.initialize = Capacitance::cap_init;
    //this->membFunc.current = Capacitance::nrn_capacity_current;
    //this->membFunc.jacob = Capacitance::nrn_cap_jacob;
}

//from eion.c
extern void ion_init(NrnThread* nt, Memb_list* ml, int type);
extern void ion_cur(NrnThread* nt, Memb_list* ml, int type);
extern void ion_alloc(double*, int*, int type);

void Mechanism::registerIon()
{
    assert(this->isIon);

    double ** ion_global_map = get_ion_global_map(); // added to membfunc.h and eion.c;
    conci = ion_global_map[type][0];
    conco = ion_global_map[type][1];
    charge = ion_global_map[type][2];

    this->membFunc.alloc = ion_alloc; //assert(0) should never be called
    this->membFunc.initialize = ion_init;
    this->membFunc.current = ion_cur;
}

Mechanism::~Mechanism(){
    delete [] sym;
    delete [] children;
};

hpx_action_t Mechanism::callModFunction = 0;
int Mechanism::callModFunction_handler(const int nargs, const void *args[], const size_t sizes[])
{
    assert(nargs==3);
    /** nargs=3 where:
     * args[0] = Memory Location of the branch
     * args[1] = Id of the function to be called
     * args[2] = Id of the mechanism to be called
     */
    Branch * branch = *(Branch**) args[0];
    Mechanism::ModFunction functionId = *(ModFunction*)args[1];
    int mechType = *(int*)args[2];
    Mechanism * mech = getMechanismFromType(mechType);

    if (branch->mechsInstances[mech->type].instancesCount > 0)
    {
        Memb_list membList;
        NrnThread nrnThread;

        //Note:The Jacob updates D and nrn_cur updates RHS, so we need a mutex for compartments
        //The state function does not write to compartment, only reads, so no mutex needed (TODO)
        Input::Coreneuron::DataLoader::fromHpxToCoreneuronDataStructs(branch, membList, nrnThread, mechType);
        switch(functionId)
        {
            case Mechanism::before_initialize:
            case Mechanism::after_initialize:
            case Mechanism::before_breakpoint:
            case Mechanism::after_solve:
            case Mechanism::before_step:
                   if (mech->BAfunctions[(int) functionId])
                       mech->BAfunctions[(int) functionId](&nrnThread, &membList, mechType);
                break;
            case Mechanism::ModFunction::alloc:
                if (mech->membFunc.alloc)
                    mech->membFunc.alloc(membList.data, membList.pdata, mechType);
                break;
            case Mechanism::ModFunction::currentCapacitance:
                assert(mechType == CAP);
                assert(mech->membFunc.current != NULL);
                mech->membFunc.current(&nrnThread, &membList, mechType);
                break;
            case Mechanism::ModFunction::current:
                assert(mechType != CAP);
                if (mech->membFunc.current)
                    mech->membFunc.current(&nrnThread, &membList, mechType);
                break;
            case Mechanism::ModFunction::state:
                if (mech->membFunc.state)
                    mech->membFunc.state(&nrnThread, &membList, mechType);
                break;
            case Mechanism::ModFunction::jacobCapacitance:
                assert(mechType == CAP);
                assert(mech->membFunc.jacob != NULL);
                mech->membFunc.jacob(&nrnThread, &membList, mechType);
                break;
            case Mechanism::ModFunction::jacob:
                assert(mechType != CAP);
                if (mech->membFunc.jacob)
                    mech->membFunc.jacob(&nrnThread, &membList, mechType);
                break;
            case Mechanism::ModFunction::initialize:
                if (mech->membFunc.initialize)
                    mech->membFunc.initialize(&nrnThread, &membList, mechType);
                break;
            case Mechanism::ModFunction::destructor:
                if (mech->membFunc.destructor)
                    mech->membFunc.destructor();
                break;
            case Mechanism::ModFunction::threadMemInit:
                if (mech->membFunc.thread_mem_init_)
                    mech->membFunc.thread_mem_init_(membList._thread);
                break;
            case Mechanism::ModFunction::threadTableCheck:
                if (mech->membFunc.thread_table_check_)
                    mech->membFunc.thread_table_check_
                        (0, membList.nodecount, membList.data, membList.pdata, membList._thread, &nrnThread, mechType);
                break;
            case Mechanism::ModFunction::threadCleanup:
                if (mech->membFunc.thread_cleanup_)
                    mech->membFunc.thread_cleanup_(membList._thread);
                break;
            default:
                printf("ERROR: Unknown ModFunction with id %d.\n", functionId);
                exit(1);
        }
    }

    //recursive mechanisms graph call (not for "capacitance")
    if (functionId  != Mechanism::ModFunction::jacobCapacitance
     && functionId  != Mechanism::ModFunction::jacob)
    {
    //call this function in all mechanisms that depend on this one
    //(ie the children on the tree of mechanisms dependencies)
      short int childrenCount = mech->childrenCount;
      hpx_addr_t lco = childrenCount > 0 ? hpx_lco_and_new(childrenCount) : HPX_NULL;

      for (short int c=0; c<childrenCount; c++)
          hpx_call(HPX_HERE, Mechanism::callModFunction, lco,
                   args[0], sizes[0],
                   args[1], sizes[1],
                   &(mech->children[c]), sizes[2]);
      if (childrenCount>0)
      {
        hpx_lco_wait(lco);
        hpx_lco_delete(lco, HPX_NULL);
      }
    }
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
