#include "neurox/neurox.h"
#include "coreneuron/nrnoc/membfunc.h"
#include <cstring>

using namespace std;
using namespace neurox;

Mechanism::Mechanism(const int type, const short int dataSize,
                     const short int pdataSize, const char isArtificial,
                     const char pntMap, const char isIon,
                     const short int symLength, const char * sym,
                     const short int dependenciesCount, const int *dependencies,
                     const short int successorsCount,   const int *successors):
    type(type), dataSize(dataSize), pdataSize(pdataSize), vdataSize(0),
    successorsCount(successorsCount), dependenciesCount(dependenciesCount),
    symLength(symLength), pntMap(pntMap), isArtificial(isArtificial),
    isIon(isIon), dependencies(nullptr), successors(nullptr), sym(nullptr)
{

    switch (type)
    {
        case IClamp           : vdataSize=1; break;
        case ProbAMPANMDA_EMS : vdataSize=2; break;
        case ProbGABAAB_EMS   : vdataSize=2; break;
        case StochKv          :	vdataSize=1; break;
        default               : vdataSize = 0;
    }
    this->membFunc.is_point = pntMap > 0 ? 1 : 0; //TODO never used

    if (dependencies != nullptr){
        assert(dependenciesCount>0);
        this->dependencies = new int[dependenciesCount];
        std::memcpy(this->dependencies, dependencies, dependenciesCount*sizeof(int));
    }
    //ion index will be set later when all mechanisms are created
    dependencyIonIndex = Mechanism::Ion::no_ion;

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

    DisableMechFunctions();
    if (this->type == CAP) //capacitance: capac.c
        RegisterCapacitance();
    else if (this->isIon)  //ion: eion.c
        RegisterIon();

    //general mechanism: mod file registration and thread tables
    RegisterModFunctions();

    RegisterBeforeAfterFunctions();

#ifndef NDEBUG
    if (HPX_LOCALITY_ID ==0)
        printf("- %s (%d), dataSize %d, pdataSize %d, isArtificial %d, pntMap %d, "
               "isIon %d, symLength %d, %d successors, %d dependencies\n",
           this->sym, this->type, this->dataSize, this->pdataSize, this->isArtificial, this->pntMap,
           this->isIon, this->symLength, this->successorsCount, this->dependenciesCount);
#endif
};

Mechanism::Ion Mechanism::GetIonIndex()
{
    assert(this->sym);
    if (strcmp("na_ion",  this->sym)==0) return Mechanism::Ion::na;
    if (strcmp("k_ion",   this->sym)==0) return Mechanism::Ion::k;
    if (strcmp("ttx_ion", this->sym)==0) return Mechanism::Ion::ttx;
    if (strcmp("ca_ion",  this->sym)==0) return Mechanism::Ion::ca;
    return Mechanism::Ion::no_ion;
}

void Mechanism::RegisterBeforeAfterFunctions()
{
    //Copy Before-After functions
    //register_mech.c::hoc_reg_ba()
    for (int i=0; i< BEFORE_AFTER_SIZE; i++)
        this->BAfunctions[i] = get_BA_function(this->sym, i);
        //this->BAfunctions[i] = nrn_threads[0].tbl[i]->bam->f;
}

void Mechanism::DisableMechFunctions()
{
    for (int i=0; i< BEFORE_AFTER_SIZE; i++)
        this->BAfunctions[i] = NULL;

    this->membFunc.alloc = NULL;
    this->membFunc.current = NULL;
    this->membFunc.current_parallel = NULL;
    this->membFunc.state = NULL;
    this->membFunc.jacob = NULL;
    this->membFunc.initialize = NULL;
    this->membFunc.destructor = NULL;
    this->membFunc.thread_mem_init_ = NULL;
    this->membFunc.thread_cleanup_ = NULL;
    this->membFunc.thread_table_check_ = NULL;
    this->nrn_bbcore_read = NULL;
    this->membFunc.thread_size_ = -1;
    this->membFunc.is_point = -1;
    this->pnt_receive = NULL;
    this->pnt_receive_init = NULL;
}

void Mechanism::RegisterModFunctions()
{
    int type = this->type;
    //TODO this should be done by exposing pointers
    //Even with cach-table disable, mechs (eg StockKv) use some of these functions
    this->membFunc.is_point = memb_func[type].is_point;
    this->membFunc.setdata_ = memb_func[type].setdata_; //TODO never used?
    this->membFunc.thread_cleanup_ = memb_func[type].thread_cleanup_;
    this->membFunc.thread_mem_init_ = memb_func[type].thread_mem_init_;
    this->membFunc.thread_size_ = memb_func[type].thread_size_;
    this->membFunc.thread_table_check_ = memb_func[type].thread_table_check_;

    if (this->type==CAP || this->isIon) return;

    this->membFunc.alloc = NULL;  //N/A memb_func[type].alloc; //TODO never used?
    this->membFunc.destructor = NULL; //N/A memb_func[type].destructor;
    this->nrn_bbcore_read = NULL; //N/A nrn_bbcore_read_[type];

    this->membFunc.current = get_cur_function(this->sym); //memb_func[type].current;
    this->membFunc.current_parallel = get_cur_parallel_function(this->sym);
    this->membFunc.jacob =  NULL; //get_jacob_function(this->sym); //memb_func[type].jacob;
    this->membFunc.state = get_state_function(this->sym); //memb_func[type].state;
    this->membFunc.initialize = get_init_function(this->sym); //memb_func[type].initialize;
    this->pnt_receive = get_net_receive_function(this->sym);
}

void Mechanism::RegisterCapacitance()
{
    assert(this->sym && strcmp("capacitance", this->sym)==0);
    this->membFunc.alloc = nrn_alloc_capacitance;
    this->membFunc.initialize = nrn_init_capacitance;
    this->membFunc.current = nrn_cur_capacitance;
    this->membFunc.current_parallel = nrn_cur_parallel_capacitance;
    this->membFunc.jacob = nrn_jacob_capacitance;
}

void Mechanism::RegisterIon()
{
    assert(this->isIon);
    assert(nrn_ion_global_map);
    //this->membFunc.alloc = ion_alloc; //assert(0) should never be called
    this->membFunc.initialize = nrn_init_ion;
    this->membFunc.current = nrn_cur_ion;
    this->membFunc.current_parallel = nrn_cur_parallel_ion;
}

Mechanism::~Mechanism(){
    delete [] sym;
    delete [] successors;
};

void Mechanism::CallModFunction(const void * branch_ptr,
                                const Mechanism::ModFunction functionId,
                                const NetConX * netcon , //for net_receive only
                                const floble_t tt        //for net_receive only
        )
{
    Branch * branch = (Branch*) branch_ptr;
    assert(branch);
    NrnThread * nrnThread = branch->nt;
    assert(nrnThread);

    if (functionId == Mechanism::ModFunction::netReceive
     || functionId == Mechanism::ModFunction::netReceiveInit)
    {
        assert(functionId != Mechanism::ModFunction::netReceiveInit); //N/A yet
        assert(this->pnt_receive);

        Memb_list * membList = &branch->mechsInstances[mechanismsMap[netcon->mechType]];
        assert(membList);
        int iml = netcon->mechInstance;
        int weightIndex = netcon->weightIndex;
        nrnThread->_t = tt; //as seen in netcvode.cpp:479 (NetCon::deliver)
        this->pnt_receive(nrnThread, membList, iml, weightIndex, 0);
        return;
    }

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
            if (membFunc.current) //has a current function
            {
                if (branch->mechsGraph //parallel execution
                 && strcmp(this->sym, "CaDynamics_E2")!=0 //not CaDynamics_E2 (no updates in cur function)
                 && !this->isIon) //not ion (updates in nrn_cur_ion function)
                {
                    if (this->dependenciesCount>0)
                        membFunc.current_parallel(
                                    nrnThread, membList, type,
                                    Branch::MechanismsGraph::AccumulateRHSandD,
                                    Branch::MechanismsGraph::AccumulateIandDIDV,
                                    branch->mechsGraph);
                    else
                        membFunc.current_parallel(
                                    nrnThread, membList, type,
                                    Branch::MechanismsGraph::AccumulateRHSandD,
                                    NULL, //no accummulation of i and di/dv
                                    branch->mechsGraph);
                }
                else //regular version
                    membFunc.current(nrnThread, membList, type);
            }
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
            {
                assert(0); //No jacob function pointers yet (get_jacob_function(xxx))
                membFunc.jacob(nrnThread, membList, type);
            }
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
            assert(0); //TODO should be called only by constructor Branch(...)
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
            assert(0); //TODO should only be called by destructor ~Branch(...)
            if (membFunc.thread_cleanup_)
                membFunc.thread_cleanup_(membList->_thread);
        break;
        default:
            printf("ERROR: Unknown ModFunction with id %d.\n", functionId);
            exit(1);
    }
}
