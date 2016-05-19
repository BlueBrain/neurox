#include "neurox/neurox.h"
#include <cstring>
#include <algorithm>

Branch::Branch(const int n, const double *a, const double *b, const double *d,
               const double *v, const double *rhs, const double *area,
               const int mechsCount, const int * mechsIds, const double *data,
               const Datum *pdata, const int childrenCount, const hpx_t * children)
    :n(n), childrenCount(childrenCount)
{
    this->a = new double[n];
    this->b = new double[n];
    this->d = new double[n];
    this->v = new double[n];
    this->rhs = new double[n];
    this->area = new double[n];
    this->mechsIds = new int[mechsCount];
    this->children = new hpx_t[childrenCount];

    memcpy(this->a,a,n*sizeof(double));
    memcpy(this->b,b,n*sizeof(double));
    memcpy(this->d,d,n*sizeof(double));
    memcpy(this->v,v,n*sizeof(double));
    memcpy(this->rhs,rhs,n*sizeof(double));
    memcpy(this->area,area,n*sizeof(double));
    memcpy(this->mechsIds,mechsIds,n*sizeof(int));
    memcpy(this->children, children, childrenCount*sizeof(hpx_t));

    int dataSize=0, pdataSize=0;
    for_each(brain->mechanisms, brain->mechanisms+brain->mechanismsCount, [&] (Mechanism & m){
        dataSize += m.dataSize; pdataSize += m.pdataSize;
    });

    this->data = new double[dataSize];
    this->pdata = new Datum[pdataSize];
    memcpy(this->data, data, dataSize*sizeof(double));
    memcpy(this->pdata, pdata, pdataSize*sizeof(Datum));
}

Branch::~Branch()
{
    delete [] b;
    delete [] d;
    delete [] a;
    delete [] v;
    delete [] rhs;
    delete [] area;
    delete [] mechsIds;
    delete [] data;
    delete [] pdata;
    delete [] children;
}

hpx_action_t Branch::initialize = 0;
int Branch::initialize_handler(const int n, const double *a, const double *b, const double *d,
                               const double *v, const double *rhs, const double *area,
                               const int mechsCount, const int * mechsIds, const double *data,
                               const Datum *pdata, const int childrenCount, const hpx_t * children)
{
    //Make sure message arrived correctly, and pin memory
    hpx_t branch_addr = hpx_thread_current_target();
    Branch * branch = NULL;
    if (!hpx_gas_try_pin(branch_addr, (void**) &branch))
        return HPX_RESEND;

    //do the work
    branch = new Branch(n,a,b,d,v,rhs,area,mechsCount,mechsIds, data, pdata, childrenCount, children);

    //unpin and return success
    hpx_gas_unpin(branch_addr);
    return HPX_SUCCESS;
}

void Branch::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  initialize, initialize_handler, HPX_INT, HPX_POINTER,
                        HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_INT, HPX_POINTER,
                        HPX_POINTER, HPX_POINTER, HPX_INT, HPX_POINTER);
}
