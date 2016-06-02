#include "neurox/Neurox.h"
#include <cstring>
#include <algorithm>

using namespace Neurox;

Branch::Branch(const int n, const double *a, const double *b, const double *d,
               const double *v, const double *rhs, const double *area,
               const int mechsCount, const int * mechsCounts, const double *data,
               const Datum *pdata, const int branchesCount, const hpx_t * branches)
    :n(n), mechsCount(mechsCount), branchesCount(branchesCount)
{
    mutex = hpx_lco_sema_new(1);

    this->a = new double[n];
    this->b = new double[n];
    this->d = new double[n];
    this->v = new double[n];
    this->rhs = new double[n];
    this->area = new double[n];
    this->mechsDataOffsets = new int[brain->mechanismsCount];
    this->branches = new hpx_t[branchesCount];

    memcpy(this->a,a,n*sizeof(double));
    memcpy(this->b,b,n*sizeof(double));
    memcpy(this->d,d,n*sizeof(double));
    memcpy(this->v,v,n*sizeof(double));
    memcpy(this->rhs,rhs,n*sizeof(double));
    memcpy(this->area,area,n*sizeof(double));
    memcpy(this->branches, branches, branchesCount*sizeof(hpx_t));

    //calculate offsets based on count
    int dataSize=0, pdataSize=0;
    for (int i=0; i<brain->mechanismsCount; i++)
    {
        dataSize += brain->mechanisms[i].dataSize * mechsCounts[i];
        pdataSize += brain->mechanisms[i].pdataSize * mechsCounts[i];
        mechsDataOffsets[i] = i==0 ? 0 : mechsDataOffsets[i-1]+ mechsCounts[i];
    }

    this->data = new double[dataSize];
    this->pdata = new Datum[pdataSize];
    memcpy(this->data, data, dataSize*sizeof(double));
    memcpy(this->pdata, pdata, pdataSize*sizeof(Datum));

    //for recursive calls
    this->futures = branchesCount ? new hpx_t[branchesCount] : nullptr;
    this->futuresAddrs = branchesCount ? new (void*)[branchesCount] : nullptr;
    this->futuresSizes = branchesCount ? new int[branchesCount] : nullptr;
    this->futuresData = branchesCount ? new BackTriangFuture[branchesCount] : nullptr;
}

Branch::~Branch()
{
    delete [] b;
    delete [] d;
    delete [] a;
    delete [] v;
    delete [] rhs;
    delete [] area;
    delete [] mechsDataOffsets;
    delete [] data;
    delete [] pdata;
    delete [] branches;

    delete [] futures;
    delete [] futuresAddrs;
    delete [] futuresSizes;
    delete [] futuresData;
}

hpx_action_t Branch::init = 0;
int Branch::init_handler(const int n, const double *a, const double *b, const double *d,
                               const double *v, const double *rhs, const double *area,
                               const int m, const int * mechsCount, const double *data,
                               const Datum *pdata, const int branchesCount, const hpx_t * branches)
{
    neurox_hpx_pin(Branch);;
    local = new Branch(n,a,b,d,v,rhs,area,m,mechsCount, data, pdata, branchesCount, branches);
    neurox_hpx_unpin;;
}

hpx_action_t Branch::setV = 0;
int Branch::setV_handler(const double v)
{
    neurox_hpx_pin(Branch);
    for (int n=0; n<local->n; n++)
        local->v[n]=v;
    neurox_hpx_recursive_branch_call(Branch::setV, v);
    neurox_hpx_unpin;
}


hpx_action_t Branch::updateV = 0;
static int updateV_handler(const int secondOrder)
{
    neurox_hpx_pin(Branch);
    for (int n=0; n<local->n; n++)
        local->v[i] += (secondOrder ? 2 : 1) * local->rhs[i];
    neurox_hpx_recursive_branch_call(Branch::updateV, v);
    neurox_hpx_unpin;
}

hpx_action_t Branch::setupMatrixInitValues = 0;
int Branch::setupMatrixInitValues_handler()
{
    neurox_hpx_pin(Branch);
    for (int n=0; n<local->n; n++)
    {
        local->rhs[n]=0;
        local->d[n]=0;
    }
    neurox_hpx_recursive_branch_call(Branch::setupMatrixInitValues);
    neurox_hpx_unpin;
}


hpx_action_t Branch::setupMatrixRHS = 0;
int Branch::setupMatrixRHS_handler(const char isSoma, const double parentV)
{
    neurox_hpx_pin(Branch);

    double *a   = local->a;
    double *b   = local->b;
    double *rhs = local->rhs;
    int branchesCount = local->branchesCount;

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

    //send/receive contribution to/from branches
    hpx_t * futures = new hpx_t[branchesCount];
    void  ** addrs  = new (void*)[branchesCount];
    size_t * sizes  = new size_t[branchesCount];
    double * values = new double[branchesCount];
    char isSoma=0;
    for (int c = 0; c < local->branchesCount; c++)
    {
        futures[c] = hpx_lco_future_new(sizeof(double));
        addrs[c]   = &values[c];
        sizes[c]   = sizeof(double);
        hpx_call(local->branches[c], setupMatrixRHS_handler, futures[c], &local->v[n-1], &isSoma);
    }

    if (branchesCount > 0) //required or fails
        hpx_lco_get_all(branchesCount, futures, sizes, addrs, NULL);

    //received contributions, can now update value
    for (int c = 0; c < branchesCount; c++)
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

struct BackTriangFuture
{
    double rhs;
    double b;
}; ///> future value of the back-triangulation method

hpx_action_t Branch::gaussianBackTriangulation = 0;
int Branch::gaussianBackTriangulation_handler(const char isSoma)
{
    neurox_hpx_pin(Branch);
    double *a   = local->a;
    double *b   = local->b;
    double *d   = local->d;
    double *rhs = local->rhs;
    int branchesCount = local->branchesCount;

    hpx_t * futures = branchesCount ? new hpx_t[branchesCount]   : nullptr;
    void  ** addrs  = branchesCount ? new (void*)[branchesCount] : nullptr;
    size_t * sizes  = branchesCount ? new size_t[branchesCount]  : nullptr;
    BackTriangFuture * values = branchesCount ? new BackTriangFuture[branchesCount] : nullptr;

    char isSoma=0;
    for (int c = 0; c < branchesCount; c++)
    {
        futures[c] = hpx_lco_future_new(sizeof (BackTriangFuture));
        addrs[c]   = &values[c];
        sizes[c]   = sizeof(BackTriangFuture);
        hpx_call(local->branches[c], Branch:::gaussianBackTriangulation, futures[c], isSoma);
    }

    if (branchesCount > 0) //required or fails
        hpx_lco_get_all(branchesCount, futures, sizes, addrs, NULL);

    //bottom compartment can now be updated with children's contribution
    for (int c = 0; c < branchesCount; c++)
    {
        d[n-1]   -= values[c].b;
        rhs[n-1] -= values[c].rhs;
        hpx_lco_delete(futures[c], HPX_NULL);
    }

    double q;
    for (int i = local->n - 1; i >= 1; --i)
    {
        q = a[i]/d[i];
        d[i-1]   -= q * b[i];
        rhs[i-1] -= q * rhs[i];
    }

    delete [] futures;
    delete [] addrs;
    delete [] sizes;
    delete [] values;

    //value to be decremented will be sent to parent branch (except soma)
    if (!isSoma)
    {
        BackTriangFuture futureData;
        q = a[0] / d[0];
        futureData.d   = q * b[0];
        futureData.rhs = q * rhs[0];
        neurox_hpx_unpin_continue(futureData);
    }
    else
        neurox_hpx_unpin;
}


hpx_action_t Branch::gaussianFwdSubstitution = 0;
int Branch::gaussianFwdSubstitution_handler(const char isSoma, const double parentRHS)
{
    neurox_hpx_pin(Branch);
    double *b   = local->b;
    double *d   = local->d;
    double *rhs = local->rhs;
    int n = local->n;
    int branchesCount = local->branchesCount;

    if(isSoma)
    {
        rhs[0] /= d[0];
    }
    else for (int i=0; i<n; i++)
    {
        rhs[i] -= b[i] * (i==0 ? parentRHS : rhs[i-1]);
        rhs[i] /= d[i];
    }

    char isSoma=0;
    double childrenRHS=rhs[n-1];
    neurox_hpx_recursive_branch_call(Branch::gaussianFwdSubstitution, isSoma, childrenRHS);
    neurox_hpx_unpin;
}


hpx_action_t Branch::setupMatrixLHS = 0;
int Branch::setupMatrixLHS_handler(const char isSoma)
{
    neurox_hpx_pin(Branch);
    double *a   = local->a;
    double *b   = local->b;
    double *d   = local->d;
    int branchesCount = local->branchesCount;

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

    //send/receive contribution to/from branches
    hpx_t * futures = new hpx_t[branchesCount];
    void  ** addrs  = new (void*)[branchesCount];
    size_t * sizes  = new size_t[branchesCount];
    double * values = new double[branchesCount];
    char notSoma=0;
    for (int c = 0; c<local->branchesCount; c++)
    {
        futures[c] = hpx_lco_future_new(sizeof(double));
        addrs[c]   = &values[c];
        sizes[c]   = sizeof(double);
        hpx_call(local->branches[c], setupMatrixRHS_handler, futures[c], &local->v[n-1], &notSoma);
    }

    if (branchesCount > 0) //required or fails
        hpx_lco_get_all(branchesCount, futures, sizes, addrs, NULL);

    //received contributions, can now update value
    for (int c = 0; c < branchesCount; c++)
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
int Branch::callMechsFunction_handler(const Mechanism::modFunctionId functionId)
{
    neurox_hpx_pin(Branch);
    for (int m=0; m<mechanismsCount; m++)
        if (mechanisms[m].functions[functionId])
            mechanisms[m].functions[functionId](
                    (mechsDataOffsets[m+1]-mechsDataOffsets[m])/dataSize, //instances count
                    mechanisms[m].dataSize,  &mechsDataOffsets[m],
                    mechanisms[m].pdataSize, &mechsPDataOffsets[m],
                    mechsNodesIndices[m]);
    neurox_hpx_recursive_branch_call(Branch::callMechsFunction, functionId);
    neurox_hpx_unpin;
}

//eion.c:
#define	nparm 5
#define cur 3
#define dcurdv 4
static void ion_alloc() { assert(0); } //used in secondOrderCurrent (below)

hpx_action_t Branch::secondOrderCurrent = 0;
int Branch::secondOrderCurrent_handler()
{
    neurox_hpx_pin(Branch);
    for (int m=0; m<mechanismsCount; m++)
        if (mechanisms[m].functions[Mechanism::modFunctionId::alloc] == ion_alloc())
        {
            int * nodeIndices = &local->mechsNodesIndices[m];
            int indicesCount = local->mechsNodesIndices[m+1] - local->mechsNodesIndices[m];
            double * mechData = &local->data[local->mechsDataOffsets[m]];
            for (int i=0; i<indicesCount; i++)
            {
                double * data = &mechData[i*nparm];
                data[cur] += data[dcurdv] * local->rhs[nodeIndices[i]]; // cur += dcurdv * rhs(ni[i])
            }
        }
    neurox_hpx_recursive_branch_call(Branch::secondOrderCurrent);
    neurox_hpx_unpin;
}


hpx_action_t Branch::queueSpike = 0;
int Branch::queueSpike_handler(const Synapse * syn, size_t)
{
    neurox_hpx_pin(Branch);
    //netcvode::PreSyn::send()
    hpx_lco_sema_p(local->mutex);
    local->queuedSynapses.push(*syn);
    hpx_lco_sema_v_sync(local->mutex);
    neurox_hpx_unpin;
}


hpx_action_t Branch::deliverSpikes = 0;
int Branch::deliverSpikes_handler()
{
    neurox_hpx_pin(Branch);
    //netcvode.cpp::NetCvode::deliver_net_events()
    //            ->NetCvode::deliver_events()
    //            ->NetCvode::deliver_event()
    //            ->NetCon::deliver()
    // (runs NET_RECEIVE on mod files)

    //Launch in all sub branches
    hpx_addr_t lco = hpx_lco_and_new(local->branchesCount);
    for (int c=0; c<local->branchesCount; c++)
        hpx_call(local->branches[c], Branch::deliverSpikes,lco);

    //Deliver all spikes
    while (local->queuedSynapses.front().deliveryTime < t+dt)
    {
        hpx_lco_sema_p(local->mutex);
        Synapse & syn = local->queuedSynapses.front();
        (*brain->mechanisms[syn.mechType].pnt_receive_t)((void *) syn.mechInstance, (void *) syn.weight, syn.deliveryTime);
        local->queuedSynapses.pop();
        hpx_lco_sema_v_sync(local->mutex);
    }

    hpx_lco_wait(lco);
    hpx_lco_delete(lco, HPX_NULL);
    neurox_hpx_unpin;
}

void Branch::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  setupMatrixRHS, setupMatrixRHS_handler, HPX_CHAR, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  setupMatrixLHS, setupMatrixLHS_handler, HPX_CHAR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  updateV, updateV_handler, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  gaussianFwdSubstitution, gaussianFwdSubstitution_handler, HPX_CHAR, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  gaussianBackTriangulation, gaussianBackTriangulation_handler, HPX_CHAR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  setV, setV_handler, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  callMechsFunction, callMechsFunction_handler, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_ATTR_NONE,   setupMatrixInitValues, setupMatrixInitValues_handler);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_MARSHALLED,  init, init_handler, HPX_INT, HPX_POINTER,
                        HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_INT, HPX_POINTER,
                        HPX_POINTER, HPX_POINTER, HPX_INT, HPX_POINTER);
}
