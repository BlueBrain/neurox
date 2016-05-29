#include "neurox/neurox.h"
#include <cstring>
#include <algorithm>

Branch::Branch(const int n, const double *a, const double *b, const double *d,
               const double *v, const double *rhs, const double *area,
               const int m, const int * mechsCounts, const double *data,
               const Datum *pdata, const int childrenCount, const hpx_t * children)
    :n(n), m(m), childrenCount(childrenCount)
{
    mutex = hpx_lco_sema_new(1);

    this->a = new double[n];
    this->b = new double[n];
    this->d = new double[n];
    this->v = new double[n];
    this->rhs = new double[n];
    this->area = new double[n];
    this->mechsOffsets = new int[brain->mechsTypesCount];
    this->children = new hpx_t[childrenCount];

    memcpy(this->a,a,n*sizeof(double));
    memcpy(this->b,b,n*sizeof(double));
    memcpy(this->d,d,n*sizeof(double));
    memcpy(this->v,v,n*sizeof(double));
    memcpy(this->rhs,rhs,n*sizeof(double));
    memcpy(this->area,area,n*sizeof(double));
    memcpy(this->children, children, childrenCount*sizeof(hpx_t));

    //calculate offsets based on count
    int dataSize=0, pdataSize=0;
    for (int i=0; i<brain->mechsTypesCount; i++)
    {
        dataSize += brain->mechsTypes[i].dataSize * mechsCounts[i];
        pdataSize += brain->mechsTypes[i].pdataSize * mechsCounts[i];
        mechsOffsets[i] = i==0 ? 0 : mechsOffsets[i-1]+ mechsCounts[i];
    }

    this->data = new double[dataSize];
    this->pdata = new Datum[pdataSize];
    memcpy(this->data, data, dataSize*sizeof(double));
    memcpy(this->pdata, pdata, pdataSize*sizeof(Datum));

    //for recursive calls
    this->futures = childrenCount ? new hpx_t[childrenCount] : nullptr;
    this->futuresAddrs = childrenCount ? new (void*)[childrenCount] : nullptr;
    this->futuresSizes = childrenCount ? new int[childrenCount] : nullptr;
    this->futuresData = childrenCount ? new FwSubFutureData[childrenCount] : nullptr;
}

Branch::~Branch()
{
    delete [] b;
    delete [] d;
    delete [] a;
    delete [] v;
    delete [] rhs;
    delete [] area;
    delete [] mechsOffsets;
    delete [] data;
    delete [] pdata;
    delete [] children;

    delete [] futures;
    delete [] futuresAddrs;
    delete [] futuresSizes;
    delete [] futuresData;
}

hpx_action_t Branch::init = 0;
int Branch::init_handler(const int n, const double *a, const double *b, const double *d,
                               const double *v, const double *rhs, const double *area,
                               const int m, const int * mechsCount, const double *data,
                               const Datum *pdata, const int childrenCount, const hpx_t * children)
{
    neurox_hpx_pin(Branch);;
    local = new Branch(n,a,b,d,v,rhs,area,m,mechsCount, data, pdata, childrenCount, children);
    neurox_hpx_unpin;;
}

hpx_action_t Branch::setV = 0;
int Branch::setV_handler(const double v)
{
    neurox_hpx_pin(Branch);
    //do the work
    for (int n=0; n<local->n; n++)
    {
        local->v[n]=v;
        local->rhs[n]=0;
        local->d[n]=0;
    }
    neurox_hpx_unpin;
}

hpx_action_t Branch::setupMatrixRHS = 0;
int Branch::setupMatrixRHS_handler(const char isSoma, const double parentV)
{
    neurox_hpx_pin(Branch);

    double *a   = local->a;
    double *b   = local->b;
    double *d   = local->d;
    double *rhs = local->rhs;
    int childrenCount = local->childrenCount;

    double returnValue = -1; //contribution to upper branch
    double dv=-1;

    //for (i = i2; i < i3; ++i))
    if (!isSoma)
    {
        dv = parentV-v[0];
        rhs[0] -= b[0]*dv;
        returnValue = a[0]*dv;
    }
    for (int i=1; i<local->n; i++)
    {
        dv = v[i-1]-v[i];
        rhs[i] -= b[i]*dv;
        rhs[i-1] += a[i]*dv;
    }

    //send/receive contribution to/from children
    hpx_t * futures = new hpx_t[childrenCount];
    void  ** addrs  = new (void*)[childrenCount];
    size_t * sizes  = new size_t[childrenCount];
    double * values = new double[childrenCount];
    char notSoma=0;
    for (int c = 0; c < local->childrenCount; c++)
    {
        futures[c] = hpx_lco_future_new(sizeof(double));
        addrs[c]   = &values[c];
        sizes[c]   = sizeof(double);
        hpx_call(local->children[c], setupMatrixRHS_handler, futures[c], &local->v[n-1], &notSoma);
    }

    if (childrenCount > 0) //required or fails
        hpx_lco_get_all(childrenCount, futures, sizes, addrs, NULL);

    //received contributions, can now update value
    for (int c = 0; c < childrenCount; c++)
    {
        rhs[n-1] += values[c];
        hpx_lco_delete(futures[c], HPX_NULL);
    }

    delete [] futures;
    delete [] addrs;
    delete [] sizes;
    delete [] values;

    if (!isSoma)
        neurox_hpx_unpin_continue(returnValue);
    else
        neurox_hpx_unpin;
}

hpx_action_t Branch::setupMatrixLHS = 0;
int Branch::setupMatrixLHS_handler(const char isSoma)
{
    neurox_hpx_pin(Branch);

    double *a   = local->a;
    double *b   = local->b;
    double *d   = local->d;
    int childrenCount = local->childrenCount;

    double returnValue = -1; //contribution to upper branch

    //for (i = i2; i < i3; ++i))
    if (!isSoma)
    {
        d[0] -= b[0];
        returnValue = a[0];
    }
    for (int i=1; i<local->n; i++)
    {
        d[i]   -= b[i];
        d[i-1] -= a[i];
    }

    //send/receive contribution to/from children
    hpx_t * futures = new hpx_t[childrenCount];
    void  ** addrs  = new (void*)[childrenCount];
    size_t * sizes  = new size_t[childrenCount];
    double * values = new double[childrenCount];
    char notSoma=0;
    for (int c = 0; c<local->childrenCount; c++)
    {
        futures[c] = hpx_lco_future_new(sizeof(double));
        addrs[c]   = &values[c];
        sizes[c]   = sizeof(double);
        hpx_call(local->children[c], setupMatrixRHS_handler, futures[c], &local->v[n-1], &notSoma);
    }

    if (childrenCount > 0) //required or fails
        hpx_lco_get_all(childrenCount, futures, sizes, addrs, NULL);

    //received contributions, can now update value
    for (int c = 0; c < childrenCount; c++)
    {
        d[n-1] -= values[c];
        hpx_lco_delete(futures[c], HPX_NULL);
    }

    delete [] futures;
    delete [] addrs;
    delete [] sizes;
    delete [] values;

    if (!isSoma)
        neurox_hpx_unpin_continue(returnValue);
    else
        neurox_hpx_unpin;
}

hpx_action_t Branch::callMechsFunction = 0;
int Branch::callMechsFunction_handler(const Mechanism::Function functionId)
{
    neurox_hpx_pin(Branch);
    for (int m=0; m<brain->mechsTypesCount; m++)
        if (brain->mechsTypes[m].functions[functionId])
            brain->mechsTypes[m].functions[functionId](NULL, NULL, m);
    neurox_hpx_unpin;
}

void Branch::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  setupMatrixRHS, setupMatrixRHS_handler, HPX_CHAR, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  setupMatrixLHS, setupMatrixLHS_handler, HPX_CHAR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  setV, setV_handler, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  callMechsFunction, callMechsFunction_handler, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  init, init_handler, HPX_INT, HPX_POINTER,
                        HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_INT, HPX_POINTER,
                        HPX_POINTER, HPX_POINTER, HPX_INT, HPX_POINTER);
}
