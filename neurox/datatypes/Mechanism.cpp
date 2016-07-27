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

Mechanism::Mechanism(const short int type, const short int dataSize, const short int pdataSize,
                     const short int dependenciesCount, const char pntMap,
                     const char isArtificial, const int * dependencies,
                     const char isIon, const double conci, const double conco, const double charge):
    type(type), dataSize(dataSize), pdataSize(pdataSize),
    dependenciesCount(dependenciesCount), pntMap(pntMap), isArtificial(isArtificial),
    isIon(isIon), conci(conci), conco(conco), charge(charge)
{
    dependencies = new int[dependenciesCount];
    std::memcpy(this->dependencies, dependencies, dependenciesCount*sizeof(int));

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
    memcpy(this->name, memb_func[type].sym, 64);

    //Copy Before-After functions
    //register_mech.c::hoc_reg_ba()
    for (int i=0; i< BEFORE_AFTER_SIZE; i++)
        this->BAfunctions[i] = nrn_threads[0].tbl[i]->bam->f;
};

Mechanism::~Mechanism(){
    delete [] dependencies;
    //TODO anything missing
};

hpx_action_t Mechanism::setMechanisms = 0;
int Mechanism::setMechanisms_handler(const Mechanism * mechanisms, const size_t mechanisms_size)
{
    neurox_hpx_pin(uint64_t);
    if (Neurox::mechanisms==NULL)
        delete [] Neurox::mechanisms;

    Neurox::mechanisms = new Mechanism[mechanisms_size*sizeof(Mechanism)];
    memcpy(Neurox::mechanisms, mechanisms, mechanisms_size);
    neurox_hpx_unpin;
}

void Mechanism::registerHpxActions()
{
    neurox_hpx_register_action(1, Mechanism::setMechanisms);
}
