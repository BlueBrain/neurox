#include "neurox/neurox.h"
#include "coreneuron/nrnoc/membfunc.h"
#include <cstring>

using namespace std;
using namespace neurox;

Mechanism::Mechanism(const int type, const short int dataSize,
                     const short int pdataSize, const char isArtificial,
                     const char pntMap, const char isIon,
                     const short int symLength, const char * sym,
                     const short int dependenciesCount,
                     const short int successorsCount, const int * successors):
    type(type), dataSize(dataSize), pdataSize(pdataSize), vdataSize(0),
    successorsCount(successorsCount), dependenciesCount(dependenciesCount),
    symLength(symLength), pntMap(pntMap), isArtificial(isArtificial),
    isIon(isIon), successors(nullptr), sym(nullptr),
    conci(-1), conco(-1), charge(-1)
{
    if (pntMap>0)
    {
        switch (type)
        {
            case IClamp           : vdataSize=1; break;
            case ProbAMPANMDA_EMS : vdataSize=2; break;
            case ProbGABAAB_EMS   : vdataSize=2; break;
            default  : throw std::invalid_argument
                    ("Unknown vdataSize for mech type (FLAG2).");
        }
        assert(vdataSize == pdataSize -1); //always?
    }

    if (successors != nullptr){
        assert(successorsCount>0);
        this->successors = new int[successorsCount];
        std::memcpy(this->successors, successors, successorsCount*sizeof(int));
    }

    if (sym != nullptr)
    {
        assert(symLength>0);
        this->sym = new char[symLength+1];
        std::memcpy(this->sym, sym, symLength);
        this->sym[symLength] = '\0';
    }

    if (sym==nullptr && successors == nullptr) //not constructing, just serialized for transmission
        return;

    disableMechFunctions();
    if (this->type == CAP) //capacitance: capac.c
        registerCapacitance();
    else if (this->isIon)  //ion: eion.c
        registerIon();
    else //general mechanism: mod file
        registerModFunctions(this->type);

    //registerBeforeAfterFunctions();

#ifndef NDEBUG
    if (HPX_LOCALITY_ID ==0)
        printf("DEBUG Mechanism: type %d, dataSize %d, pdataSize %d, isArtificial %d,\n"
           "      pntMap %d, isIon %d, symLength %d, sym %s, successorsCount %d\n"
           "      dependenciesCount %d, conci %.2f, conco %.2f, charge %.2f\n",
           this->type, this->dataSize, this->pdataSize, this->isArtificial,
           this->pntMap, this->isIon, this->symLength, this->sym, this->successorsCount,
           this->dependenciesCount, this->conci, this->conco, this->charge);
#endif
};

hpx_action_t Mechanism::initModFunction = 0;
void Mechanism::initModFunction_handler(ModFunction *function_ptr, const size_t)
{}

hpx_action_t Mechanism::reduceModFunction = 0;
void Mechanism::reduceModFunction_handler
    (ModFunction * lhs, const ModFunction *rhs, const size_t)
{ *lhs = *rhs; }

void Mechanism::registerBeforeAfterFunctions()
{
    //Copy Before-After functions
    //register_mech.c::hoc_reg_ba()
    for (int i=0; i< BEFORE_AFTER_SIZE; i++)
        this->BAfunctions[i] = get_BA_function(this->sym, i);
        //this->BAfunctions[i] = nrn_threads[0].tbl[i]->bam->f;
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
    this->nrn_bbcore_read = NULL;
    this->membFunc.thread_size_ = -1;
    this->membFunc.is_point = -1;
    this->pnt_receive = NULL;
    this->pnt_receive_init = NULL;
}

void Mechanism::registerModFunctions(int type)
{
    this->membFunc.alloc = NULL; // memb_func[type].alloc;
    this->membFunc.setdata_ = NULL; // memb_func[type].setdata_;
    this->membFunc.destructor = NULL; // memb_func[type].destructor;
    this->membFunc.current = get_cur_function(this->sym); //memb_func[type].current;
    this->membFunc.jacob = NULL; //get_jacob_function(this->sym); //memb_func[type].jacob;
    this->membFunc.state = get_state_function(this->sym); //memb_func[type].state;
    this->membFunc.initialize = get_init_function(this->sym); //memb_func[type].initialize;
    this->membFunc.thread_cleanup_ = NULL; //memb_func[type].thread_cleanup_;
    this->membFunc.thread_mem_init_ = NULL; //TODO memb_func[type].thread_mem_init_;
    this->membFunc.thread_size_ = 0; //TODO memb_func[type].thread_size_;
    this->membFunc.thread_table_check_ = NULL; //TODO memb_func[type].thread_table_check_;
    this->nrn_bbcore_read = NULL; //TODO nrn_bbcore_read_[type];
}

void Mechanism::registerCapacitance()
{
    assert(this->sym && strcmp("capacitance", this->sym)==0);
    this->membFunc.alloc = nrn_alloc_capacitance;
    this->membFunc.initialize = nrn_init_capacitance;
    this->membFunc.current = nrn_cur_capacitance;
    this->membFunc.jacob = nrn_jacob_capacitance;
    //this->membFunc.alloc = Capacitance::cap_alloc;
    //this->membFunc.initialize = Capacitance::cap_init;
    //this->membFunc.current = Capacitance::nrn_capacity_current;
    //this->membFunc.jacob = Capacitance::nrn_cap_jacob;
}

void Mechanism::registerIon()
{
    assert(this->isIon);

    double ** ion_global_map = nrn_ion_global_map; //get_ion_global_map();
    conci  = (floble_t) ion_global_map[type][0];
    conco  = (floble_t) ion_global_map[type][1];
    charge = (floble_t) ion_global_map[type][2];

    //this->membFunc.alloc = ion_alloc; //assert(0) should never be called
    this->membFunc.initialize = nrn_init_ion;
    this->membFunc.current = nrn_cur_ion;
}

Mechanism::~Mechanism(){
    delete [] sym;
    delete [] successors;
};

void Mechanism::callModFunction(const void * branch_ptr,
                                const Mechanism::ModFunction functionId)
{
    Branch * branch = (Branch*) branch_ptr;
    assert(branch);
    NrnThread * nrnThread = branch->nt;
    assert(nrnThread);
    Memb_list * membList = &branch->mechsInstances[mechanismsMap[this->type]];
    assert(membList);

    if (membList->nodecount>0)
    switch(functionId)
    {
        case Mechanism::before_initialize:
        case Mechanism::after_initialize:
        case Mechanism::before_breakpoint:
        case Mechanism::after_solve:
        case Mechanism::before_step:
               if (BAfunctions[(int) functionId])
                   BAfunctions[(int) functionId](nrnThread, membList, type);
            break;
        case Mechanism::ModFunction::alloc:
            if (membFunc.alloc)
                membFunc.alloc(membList->data, membList->pdata, type);
            break;
        case Mechanism::ModFunction::currentCapacitance:
            assert(type == CAP);
            assert(membFunc.current != NULL);
            membFunc.current(nrnThread, membList, type);
            break;
        case Mechanism::ModFunction::current:
            assert(type != CAP);
            if (membFunc.current)
                membFunc.current(nrnThread, membList, type);
            break;
        case Mechanism::ModFunction::state:
            if (membFunc.state)
                membFunc.state(nrnThread, membList, type);
            break;
        case Mechanism::ModFunction::jacobCapacitance:
            assert(type == CAP);
            assert(membFunc.jacob != NULL);
            membFunc.jacob(nrnThread, membList, type);
            break;
        case Mechanism::ModFunction::jacob:
            assert(type != CAP);
            if (membFunc.jacob)
                membFunc.jacob(nrnThread, membList, type);
            break;
        case Mechanism::ModFunction::initialize:
            if (membFunc.initialize)
                membFunc.initialize(nrnThread, membList, type);
            break;
        case Mechanism::ModFunction::destructor:
            if (membFunc.destructor)
                membFunc.destructor();
            break;
        case Mechanism::ModFunction::threadMemInit:
            if (membFunc.thread_mem_init_)
                membFunc.thread_mem_init_(membList->_thread);
            break;
        case Mechanism::ModFunction::threadTableCheck:
            if (membFunc.thread_table_check_)
                membFunc.thread_table_check_
                    (0, membList->nodecount, membList->data, membList->pdata,
                     membList->_thread, nrnThread, type);
            break;
        case Mechanism::ModFunction::threadCleanup:
            if (membFunc.thread_cleanup_)
                membFunc.thread_cleanup_(membList->_thread);
            break;
        default:
            printf("ERROR: Unknown ModFunction with id %d.\n", functionId);
            exit(1);
    }
}

void Mechanism::callNetReceiveFunction(
        const void * branch_ptr, const NetConX * netcon,
        const floble_t t, const char callNetReceiveInit)
{
    Branch * branch = (Branch*) branch_ptr;
    //Memb_list & membList = branch->mechsInstances[mechanismsMap];
    NrnThread * nt = branch->nt;

    Point_process pp;
    pp._i_instance = netcon->mechInstance;
    pp._tid = -1;

    pp._type = netcon->mechType;
    nt->_t = nt->_dt + netcon->delay; //delivery time
    //see netcvode.cpp:433 (NetCon::deliver)
    //We have to pass NrnThread, MembList, and deliveryTime instead
    if (callNetReceiveInit)
    {   //Not possible to pass all arguments, will fail!
        //getMechanismFromType(netcon->mechType)->pnt_receive_init
                //(&nrnThread, &membList, &pp, netcon->args, 0);
    }
    else
    {
        //getMechanismFromType(netcon->mechType)->pnt_receive
               // (&nrnThread, &membList, &pp, netcon->args, 0);
    }
}

void Mechanism::registerHpxActions()
{
    neurox_hpx_register_action(5, Mechanism::initModFunction);
    neurox_hpx_register_action(5, Mechanism::reduceModFunction);
}
