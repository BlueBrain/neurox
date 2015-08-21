#include "neurox/neurox.h"
#include "coreneuron/nrnoc/membfunc.h"
#include <cstring>

using namespace std;
using namespace neurox;

Mechanism::Mechanism(const int type, const short int dataSize,
                     const short int pdataSize, const char isArtificial,
                     char pntMap, const char isIon,
                     const short int symLength, const char * sym,
                     Memb_func & memb_func,
                     const short int dependenciesCount, const int *dependencies,
                     const short int successorsCount,   const int *successors):
    type(type), dataSize(dataSize), pdataSize(pdataSize), vdataSize(0),
    successorsCount(successorsCount), dependenciesCount(dependenciesCount),
    symLength(symLength), pntMap(pntMap), isArtificial(isArtificial),
    isIon(isIon), dependencies(nullptr), successors(nullptr)
{

    //to be set by neuronx::UpdateMechanismsDependencies
    this->dependencyIonIndex = Mechanism::Ion::no_ion;

    //set function pointers
    memcpy(&this->membFunc, &memb_func, sizeof(Memb_func));
    assert(symLength>0);

    //set name (overwrites existing, pointing to CN data structures)
    this->membFunc.sym = new char[symLength+1];
    std::memcpy(this->membFunc.sym, sym, symLength);
    this->membFunc.sym[symLength] = '\0';

    //copy dparam_semantics (overwrites existing, pointing to CN data structures)
    this->membFunc.dparam_semantics = new int[pdataSize];
    memcpy(this->membFunc.dparam_semantics, memb_func.dparam_semantics, sizeof(int)*pdataSize);

    //non-coreneuron functions
    if (this->type!=CAP && !this->isIon)
    {
        this->membFunc.current_parallel = get_cur_parallel_function(this->membFunc.sym);
        this->pnt_receive = get_net_receive_function(this->membFunc.sym);
    }
    else if (this->type == CAP) //capacitance: capac.c
    {
        //these are not registered by capac.c
        this->membFunc.current = nrn_cur_capacitance;
        this->membFunc.current_parallel = nrn_cur_parallel_capacitance;
        this->membFunc.jacob = nrn_jacob_capacitance;
    }
    else if (this->isIon)  //ion: eion.c
    {
        //these are not registered by eion.c
        this->membFunc.current = nrn_cur_ion;
        this->membFunc.current_parallel = nrn_cur_parallel_ion;
    }

    switch (type)
    {
        case IClamp           : vdataSize=1; break;
        case ProbAMPANMDA_EMS : vdataSize=2; break;
        case ProbGABAAB_EMS   : vdataSize=2; break;
        case StochKv          :	vdataSize=1; break;
        default               : vdataSize = 0;
    }
    this->membFunc.is_point = pntMap > 0 ? 1 : 0;

    if (dependencies != nullptr){
        assert(dependenciesCount>0);
        this->dependencies = new int[dependenciesCount];
        std::memcpy(this->dependencies, dependencies, dependenciesCount*sizeof(int));
    }

    if (successors != nullptr){
        assert(successorsCount>0);
        this->successors = new int[successorsCount];
        std::memcpy(this->successors, successors, successorsCount*sizeof(int));
    }
};

Mechanism::Ion Mechanism::GetIonIndex()
{
    assert(this->membFunc.sym);
    if (strcmp("na_ion",  this->membFunc.sym)==0) return Mechanism::Ion::na;
    if (strcmp("k_ion",   this->membFunc.sym)==0) return Mechanism::Ion::k;
    if (strcmp("ttx_ion", this->membFunc.sym)==0) return Mechanism::Ion::ttx;
    if (strcmp("ca_ion",  this->membFunc.sym)==0) return Mechanism::Ion::ca;
    return Mechanism::Ion::no_ion;
}

void Mechanism::RegisterBeforeAfterFunctions()
{
    //Copy Before-After functions
    //register_mech.c::hoc_reg_ba()
    for (int i=0; i< BEFORE_AFTER_SIZE; i++)
        this->BAfunctions[i] = get_BA_function(this->membFunc.sym, i); //NULL at the moment
        //this->BAfunctions[i] = nrn_threads[0].tbl[i]->bam->f;
}

Mechanism::~Mechanism(){
    delete [] membFunc.sym;
    delete [] successors;
    delete [] dependencies;
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

        Memb_list * membList = &branch->mechsInstances[mechanisms_map[netcon->mechType]];
        assert(membList);
        int iml = netcon->mechInstance;
        int weightIndex = netcon->weightIndex;
        nrnThread->_t = tt; //as seen in netcvode.cpp:479 (NetCon::deliver)
        this->pnt_receive(nrnThread, membList, iml, weightIndex, 0);
        return;
    }

    Memb_list * membList = &branch->mechsInstances[mechanisms_map[this->type]];
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
                if (input_params->multiMex //parallel execution
                 && strcmp(this->membFunc.sym, "CaDynamics_E2")!=0 //not CaDynamics_E2 (no updates in cur function)
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
            assert(0); //should be called only by constructor Branch(...)
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
            assert(0); //should only be called by destructor ~Branch(...)
            if (membFunc.thread_cleanup_)
                membFunc.thread_cleanup_(membList->_thread);
        break;
        default:
            printf("ERROR: Unknown ModFunction with id %d.\n", functionId);
            exit(1);
    }
}
