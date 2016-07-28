#include "neurox/Neurox.h"
#include <cstring>
//#include "coreneuron/nrnoc/membfunc.h" //Memb_func, BAMech

using namespace std;
using namespace Neurox;

//from capac.c
extern void cap_alloc(double*, Datum*, int type);
extern void cap_init(NrnThread*, Memb_list*, int);

//from eion.c
extern void ion_alloc();
extern void ion_cur(NrnThread*, Memb_list*, int);
extern void ion_init(NrnThread*, Memb_list*, int);

Mechanism::Mechanism():
 type(-1), dataSize(-1), pdataSize(-1), dependenciesCount(-1), pntMap(0),
 isArtificial(0), dependencies(nullptr), isIon(0), conci(-1), conco(-1),
 charge(-1), sym(nullptr)
{};

Mechanism::Mechanism(const short int type, const short int dataSize, const short int pdataSize,
                     const char isArtificial, const char pntMap, const char isIon,
                     const double conci, const double conco, const double charge,
                     const short int symLength, const char * sym,
                     const short int dependenciesCount, const int * dependencies):
    type(type), dataSize(dataSize), pdataSize(pdataSize), dependenciesCount(dependenciesCount),
    pntMap(pntMap), isArtificial(isArtificial),
    isIon(isIon), conci(conci), conco(conco), charge(charge)
{
    if (dependencies != nullptr && dependenciesCount>0)
    {
        this->dependencies = new int[dependenciesCount];
        std::memcpy(this->dependencies, dependencies, dependenciesCount*sizeof(int));
    }
    else this->dependencies=nullptr;

    if (sym != nullptr && symLength>0)
    {
        this->sym = new char[symLength+1];
        std::memcpy(this->sym, sym, symLength);
        this->sym[symLength] = '\0';
    }
    else this->sym=nullptr;

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

    //Copy Before-After functions
    //register_mech.c::hoc_reg_ba()
    for (int i=0; i< BEFORE_AFTER_SIZE; i++)
    {
        this->BAfunctions[i] = nullptr;
        //not implemented by CoreNeuron, undefined bam leads to SEGFAULT:
        //this->BAfunctions[i] = nrn_threads[0].tbl[i]->bam->f;
    }
};

Mechanism::~Mechanism(){
    delete [] dependencies;
    delete [] sym;
};

