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
                     const short int dependenciesCount, const int * dependencies):
    type(type), dataSize(dataSize), pdataSize(pdataSize),
    pntMap(pntMap), isArtificial(isArtificial),
    symLength(symLength), isIon(isIon),
    dependencies(nullptr), sym(nullptr),
    conci(-1), conco(-1), charge(-1)
{
    if (dependencies != nullptr && dependenciesCount>0)
    {
        this->dependencies = new int[dependenciesCount];
        std::memcpy(this->dependencies, dependencies, dependenciesCount*sizeof(int));
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

