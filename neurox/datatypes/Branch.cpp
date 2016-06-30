#include "neurox/Neurox.h"
#include <cstring>
#include <algorithm>
#include <numeric>

using namespace Neurox;

Branch::Branch(const int n, const double *a, const double *b, const double *d,
               const double *v, const double *rhs, const double *area,
               const int * mechsCounts, const double *data,
               const int *pdata, const int branchesCount, const hpx_t * branches)
    :n(n), branchesCount(branchesCount)
{
    synapsesQueueMutex = hpx_lco_sema_new(1);

    this->a = new double[n];
    this->b = new double[n];
    this->d = new double[n];
    this->v = new double[n];
    this->rhs = new double[n];
    this->area = new double[n];
    mechsInstances.dataOffsets = new int[Neurox::mechanismsCount+1];
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
    for (int i=0; i<Neurox::mechanismsCount; i++)
    {
        dataSize  += mechanisms[i].dataSize * mechsCounts[i];
        pdataSize += mechanisms[i].pdataSize * mechsCounts[i];
        mechsInstances.dataOffsets[i] = i==0 ? 0 : mechsInstances.dataOffsets[i-1]+mechsCounts[i];
    }
    //copy children addresses
    for (int b=0; b<branchesCount; b++)
        this->branches[b]=branches[b];

    mechsInstances.data = new double[dataSize];
    mechsInstances.pdata = new int[pdataSize];
    memcpy(mechsInstances.data, data,  dataSize*sizeof(double));
    memcpy(mechsInstances.pdata, pdata, pdataSize*sizeof(int));
}

Branch::~Branch()
{
    delete [] b;
    delete [] d;
    delete [] a;
    delete [] v;
    delete [] rhs;
    delete [] area; //TODO is this ever used?
    delete [] branches;
    delete [] mechsInstances.data;
    delete [] mechsInstances.dataOffsets;
    delete [] mechsInstances.pdata;
    delete [] mechsInstances.pdataOffsets;
    delete [] mechsInstances.nodesIndices;
    delete [] mechsInstances.nodesIndicesOffsets;

}

hpx_action_t Branch::init = 0;
int Branch::init_handler(const int n, const double *a, const double *b, const double *d,
                               const double *v, const double *rhs, const double *area,
                               const int * mechsCount, const double *data,
                               const Datum *pdata, const int branchesCount, const hpx_t * branches)
{
    neurox_hpx_pin(Branch);
    local = new Branch(n,a,b,d,v,rhs,area,mechsCount, data, pdata, branchesCount, branches);
    neurox_hpx_unpin;;
}

hpx_action_t Branch::setV = 0;
int Branch::setV_handler(const double v)
{
    neurox_hpx_pin(Branch);
    for (int n=0; n<local->n; n++)
        local->v[n]=v;
    neurox_hpx_recursive_branch_call_sync(Branch::setV, v);
    neurox_hpx_unpin;
}


hpx_action_t Branch::updateV = 0;
int Branch::updateV_handler(const int secondOrder)
{
    neurox_hpx_pin(Branch);
    for (int i=0; i<local->n; i++)
        local->v[i] += (secondOrder ? 2 : 1) * local->rhs[i];
    neurox_hpx_recursive_branch_call_sync(Branch::updateV, secondOrder);
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
    neurox_hpx_recursive_branch_call_sync(Branch::setupMatrixInitValues);
    neurox_hpx_unpin;
}


hpx_action_t Branch::setupMatrixRHS = 0;
int Branch::setupMatrixRHS_handler(const char isSoma, const double parentV)
{
    neurox_hpx_pin(Branch);

    int n = local->n;
    double *a   = local->a;
    double *b   = local->b;
    double *v   = local->v;
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
    void  ** addrs  = new void*[branchesCount];
    size_t * sizes  = new size_t[branchesCount];
    double * values = new double[branchesCount];
    char isSomaFlag=0;
    for (int c = 0; c < local->branchesCount; c++)
    {
        futures[c] = hpx_lco_future_new(sizeof(double));
        addrs[c]   = &values[c];
        sizes[c]   = sizeof(double);
        hpx_call(local->branches[c], Branch::setupMatrixRHS, futures[c], local->v[n-1], isSomaFlag);
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
    int n = local->n;
    double *a   = local->a;
    double *b   = local->b;
    double *d   = local->d;
    double *rhs = local->rhs;
    int branchesCount = local->branchesCount;

    hpx_t * futures = branchesCount ? new hpx_t[branchesCount]  : nullptr;
    void  ** addrs  = branchesCount ? new void*[branchesCount]  : nullptr;
    size_t * sizes  = branchesCount ? new size_t[branchesCount] : nullptr;
    BackTriangFuture * values = branchesCount ? new BackTriangFuture[branchesCount] : nullptr;

    char isSomaFlag=0;
    for (int c = 0; c < branchesCount; c++)
    {
        futures[c] = hpx_lco_future_new(sizeof (BackTriangFuture));
        addrs[c]   = &values[c];
        sizes[c]   = sizeof(BackTriangFuture);
        hpx_call(local->branches[c], Branch::gaussianBackTriangulation, futures[c], isSomaFlag);
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
    for (int i = n - 1; i >= 1; --i)
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
        futureData.b   = q * b[0];
        futureData.rhs = q * rhs[0];
        neurox_hpx_unpin_continue(futureData);
    }
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

    if(isSoma)
    {
        rhs[0] /= d[0];
    }
    else for (int i=0; i<n; i++)
    {
        rhs[i] -= b[i] * (i==0 ? parentRHS : rhs[i-1]);
        rhs[i] /= d[i];
    }

    char isSomaFlag=0;
    double childrenRHS=rhs[n-1];
    neurox_hpx_recursive_branch_call_sync(Branch::gaussianFwdSubstitution, isSomaFlag, childrenRHS);
    neurox_hpx_unpin;
}


hpx_action_t Branch::setupMatrixLHS = 0;
int Branch::setupMatrixLHS_handler(const char isSoma)
{
    neurox_hpx_pin(Branch);
    int n = local->n;
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
    void  ** addrs  = new void*[branchesCount];
    size_t * sizes  = new size_t[branchesCount];
    double * values = new double[branchesCount];
    char isSomaFlag=0;
    for (int c = 0; c<local->branchesCount; c++)
    {
        futures[c] = hpx_lco_future_new(sizeof(double));
        addrs[c]   = &values[c];
        sizes[c]   = sizeof(double);
        hpx_call(local->branches[c], Branch::setupMatrixLHS, futures[c], &isSomaFlag);
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
    neurox_hpx_unpin;
}

hpx_action_t Branch::callMechsFunction = 0;
int Branch::callMechsFunction_handler(const Mechanism::Functions functionId, const double t, const double dt)
{
    neurox_hpx_pin(Branch);

    Memb_list membList;
    NrnThread nrnThread;

    Branch::MechanismInstances & mechs = local->mechsInstances;
    //TODO parallelize this loop
    //PS what if different mechanisms overwrite the same nodes voltage?
    for (int m=0; m<mechanismsCount; m++)
        if (mechanisms[m].BAfunctions[functionId])
        {
            membList.data  = &mechs.data[mechs.dataOffsets[m]];
            membList.pdata = &mechs.pdata[mechs.pdataOffsets[m]];
            membList.nodecount = (mechs.dataOffsets[m+1] - mechs.dataOffsets[m])/mechanisms[m].dataSize; //instances count
            membList.nodeindices = &mechs.nodesIndices[mechs.nodesIndicesOffsets[m]];
            membList._thread = NULL; //TODO: ThreadDatum never used ?
            nrnThread._actual_d = local->d;
            nrnThread._actual_rhs = local->rhs;
            nrnThread._actual_v = local->v;
            nrnThread._actual_area = local->area;
            nrnThread._dt = dt;
            nrnThread._t = t;
            if (functionId<BEFORE_AFTER_SIZE)
                mechanisms[m].BAfunctions[functionId](&nrnThread, &membList, m);
            else
                switch(functionId)
                {
                case Mechanism::Functions::alloc:
                   mechanisms[m].membFunc.alloc(membList.data, membList.pdata, m); break;
                case Mechanism::Functions::current:
                    mechanisms[m].membFunc.current(&nrnThread, &membList, m); break;
                case Mechanism::Functions::state:
                    mechanisms[m].membFunc.state(&nrnThread, &membList, m); break;
                case Mechanism::Functions::jacob:
                    mechanisms[m].membFunc.jacob(&nrnThread, &membList, m); break;
                case Mechanism::Functions::initialize:
                    mechanisms[m].membFunc.initialize(&nrnThread, &membList, m); break;
                case Mechanism::Functions::destructor:
                    mechanisms[m].membFunc.destructor(); break;
                case Mechanism::Functions::threadMemInit:
                    mechanisms[m].membFunc.thread_mem_init_(membList._thread); break;
                case Mechanism::Functions::threadTableCheck:
                    mechanisms[m].membFunc.thread_table_check_
                        (0, membList.nodecount, membList.data, membList.pdata, membList._thread, &nrnThread, m);
                    break;
                case Mechanism::Functions::threadCleanup:
                    mechanisms[m].membFunc.thread_cleanup_(membList._thread); break;
                case Mechanism::Functions::setData:
                    mechanisms[m].membFunc.setdata_(membList.data, membList.pdata); break;
                default:
                    printf("ERROR! Unknown function!!\n"); //TODO
                }
        }
    //TODO should be called at the beginning of this function, not the end (maybe we'll need mutex)
    neurox_hpx_recursive_branch_call_sync(Branch::callMechsFunction, functionId, t, dt);
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
    MechanismInstances & mechs = local->mechsInstances;
    for (int m=0; m<mechanismsCount; m++)
        if (mechanisms[m].BAfunctions[Mechanism::Functions::alloc] ) //TODO used to be if == ion_alloc()
        {
            int * nodeIndices = &mechs.nodesIndices[m];
            int indicesCount = mechs.nodesIndicesOffsets[m+1] - mechs.nodesIndicesOffsets[m];
            double * mechData = &mechs.data[mechs.dataOffsets[m]];
            for (int i=0; i<indicesCount; i++)
            {
                double * data = &mechData[i*nparm];
                data[cur] += data[dcurdv] * local->rhs[nodeIndices[i]]; // cur += dcurdv * rhs(ni[i])
            }
        }
    neurox_hpx_recursive_branch_call_sync(Branch::secondOrderCurrent);
    neurox_hpx_unpin;
}


hpx_action_t Branch::queueSpike = 0;
int Branch::queueSpike_handler(const int preNeuronId, const double deliveryTime)
{
    neurox_hpx_pin(Branch);
    //netcvode::PreSyn::send()
    hpx_lco_sema_p(local->synapsesQueueMutex);
    local->synapsesQueue.push(*syn);
    hpx_lco_sema_v_sync(local->synapsesQueueMutex);
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
    if (local->synapsesQueue.size()>0)
    while (local->synapsesQueue.top().deliveryTime < t+dt)
    {
        hpx_lco_sema_p(local->synapsesQueueMutex);
        Synapse syn = local->synapsesQueue.top();
        //(mechanisms[syn.mechType].functions[Mechanism::Functions::NetReceive])((void *) syn.mechInstance, (void *) syn.weight, syn.deliveryTime);
        (mechanisms[syn.mechType].BAfunctions[Mechanism::Functions::NetReceive])(NULL,NULL,0);
        //( short instancesCount, short dataSize, double * data, short pdataSize, int * pdata, int * nodesIndices)
        local->synapsesQueue.pop();
        hpx_lco_sema_v_sync(local->synapsesQueueMutex);
    }

    hpx_lco_wait(lco);
    hpx_lco_delete(lco, HPX_NULL);
    neurox_hpx_unpin;
}

void Branch::registerHpxActions()
{
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,  setupMatrixRHS, setupMatrixRHS_handler, HPX_CHAR, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,  setupMatrixLHS, setupMatrixLHS_handler, HPX_CHAR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,  updateV, updateV_handler, HPX_INT);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,  gaussianFwdSubstitution, gaussianFwdSubstitution_handler, HPX_CHAR, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,  gaussianBackTriangulation, gaussianBackTriangulation_handler, HPX_CHAR);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,  setV, setV_handler, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,  callMechsFunction, callMechsFunction_handler, HPX_INT, HPX_DOUBLE, HPX_DOUBLE);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,  setupMatrixInitValues, setupMatrixInitValues_handler);
    HPX_REGISTER_ACTION(HPX_DEFAULT, HPX_PINNED,  init, init_handler, HPX_INT, HPX_POINTER,
                        HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_POINTER, HPX_INT, HPX_POINTER,
                        HPX_POINTER, HPX_POINTER, HPX_INT, HPX_POINTER);
}
