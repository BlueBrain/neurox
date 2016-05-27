#include "neurox/neurox.h"
#include <cstring>
#include "coreneuron/nrnoc/membfunc.h" //Memb_func, BAMech

//BAMech = Before and After running mechanism

using namespace std;

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
    this->alloc = memb_func[type].alloc;
    this->thread_mem_init = memb_func[type].thread_mem_init_;
    this->thread_cleanup = memb_func[type].thread_cleanup_;
    this->thread_table_check = memb_func[type].thread_table_check_;
    this->setdata = memb_func[type].setdata_;
    this->destructor = memb_func[type].destructor;
    memcpy(this->name, memb_func[type].sym, 64);

    //finitialize.c->nrn_finitialize()->nrn_ba()
    functions[Function::current] = memb_func[type].current;
    functions[Function::jacob] = memb_func[type].jacob;
    functions[Function::state] = memb_func[type].state;
    functions[Function::initialize] = memb_func[type].initialize;

    //Before After Functions will be set after (by Brain)
    //register_mech.c::hoc_reg_ba()
    for (int i=Function::initialize; i< Function::functionsCount; i++)
        functions[i] = (mod_f_t) 0; //TODO

     //look up tables: multicore.c::nrn_thread_table_check()
    //where's the look up table data
    if (*thread_table_check)
    {
        //TODO
        //(*thread_table_check)(0, ml->nodecount, ml->data, ml->pdata, ml->_thread, nt, tml->index);
        //look up tables are created and destroyed inside mod files, not accessible via coreneuron
    }
};

Mechanism::~Mechanism(){
    delete [] dependencies;
    //TODO anything missing
};
