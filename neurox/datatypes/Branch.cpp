#include "neurox/Neurox.h"
#include <cstring>
#include <algorithm>
#include <numeric>

//Serialization data structures
//#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

using namespace Neurox;

Branch::~Branch()
{
    delete [] b;
    delete [] d;
    delete [] a;
    delete [] v;
    delete [] rhs;
    delete [] area;
    delete [] branches;
    for (int m=0; m<mechanismsCount; m++)
    {
        delete [] mechsInstances[m].data;
        delete [] mechsInstances[m].pdata;
        delete [] mechsInstances[m].nodesIndices;
    }
    delete [] mechsInstances;
}

hpx_action_t Branch::init = 0;
int Branch::init_handler(const int n, std::vector<double> * a, const std::vector<double> *b,
                         const std::vector<double> * d, const std::vector<double> * v,
                         const std::vector<double> * rhs, const std::vector<double> * area,
                         const std::vector<int> * instancesCount, const std::vector<std::vector<double> > * data,
                         const std::vector<std::vector<double> > * pdata, std::vector<std::vector<int>> * nodesIndices,
                         const std::vector<hpx_t> * branches, std::map<int, std::vector<NetConX> > * netcons)
{
    neurox_hpx_pin(Branch);
    local->n = n;
    local->branchesCount = branches->size();
    local->spikesQueueMutex = hpx_lco_sema_new(1);
    local->netcons.insert(netcons->begin(), netcons->end());

    local->a = new double[n];
    local->b = new double[n];
    local->d = new double[n];
    local->v = new double[n];
    local->rhs = new double[n];
    local->area = new double[n];
    local->branches = new hpx_t[local->branchesCount];

    std::copy(a->begin(), a->end(), local->a);
    std::copy(b->begin(), b->end(), local->b);
    std::copy(d->begin(), d->end(), local->d);
    std::copy(v->begin(), v->end(), local->v);
    std::copy(rhs->begin(), rhs->end(), local->rhs);
    std::copy(area->begin(), area->end(), local->area);
    std::copy(branches->begin(), branches->end(), local->branches);

    //copy mechanisms instances
    local->mechsInstances = new MechanismInstances[mechanismsCount];
    for (int m=0; m<mechanismsCount; m++)
    {
        local->mechsInstances[m].instancesCount = instancesCount->at(m);
        local->mechsInstances[m].data = new double[instancesCount->at(m)*mechanisms[m].dataSize];
        local->mechsInstances[m].pdata = new int[instancesCount->at(m)*mechanisms[m].pdataSize];
        local->mechsInstances[m].nodesIndices = new int[instancesCount->at(m)];

        std::copy(data->at(m).begin(),  data->at(m).end(),  local->mechsInstances[m].data);
        std::copy(pdata->at(m).begin(), pdata->at(m).end(), local->mechsInstances[m].pdata);
        std::copy(nodesIndices->at(m).begin(), nodesIndices->at(m).end(), local->mechsInstances[m].nodesIndices);
    }

    //inform pre-synaptic neurons that we connect (my hpx address is stored in variable "target")
    hpx_addr_t lco =  local->netcons.size() ?  local->netcons.size() : HPX_NULL;
    for (auto nc = local->netcons.begin(); nc != local->netcons.end(); nc++)
    {
        int preNeuronId = getNeuronAddr(nc->first);
        int e = hpx_call(preNeuronId, Neuron::addSynapseTarget, lco, &target, sizeof(target)) ;
        assert(e==HPX_SUCCESS);
    }

    neurox_hpx_unpin;;
}

hpx_action_t Branch::setV = 0;
int Branch::setV_handler(const double * v, const size_t v_size)
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::setV, &v, v_size);
    for (int n=0; n<local->n; n++)
        local->v[n]=*v;
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t Branch::updateV = 0;
int Branch::updateV_handler(const int * secondOrder, size_t secondOrder_size)
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::updateV, &secondOrder, secondOrder_size);
    for (int i=0; i<local->n; i++)
        local->v[i] += (*secondOrder ? 2 : 1) * local->rhs[i];
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t Branch::setupMatrixInitValues = 0;
int Branch::setupMatrixInitValues_handler()
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::setupMatrixInitValues);
    for (int n=0; n<local->n; n++)
    {
        local->rhs[n]=0;
        local->d[n]=0;
    }
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t Branch::setupMatrixRHS = 0;
int Branch::setupMatrixRHS_handler(const int nargs, const void *args[], const size_t[])
{
    neurox_hpx_pin(Branch);
    assert(nargs==2);
    const char * isSoma = (const char*) args[0];
    const double * parentV = (const double*) args[1];

    int n = local->n;
    double *a   = local->a;
    double *b   = local->b;
    double *v   = local->v;
    double *rhs = local->rhs;
    int branchesCount = local->branchesCount;

    double returnValue = -1; //contribution to upper branch
    double dv=-1;

    //for (i = i2; i < i3; ++i))
    if (!*isSoma)
    {
        dv = *parentV-v[0];
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
        hpx_call(local->branches[c], Branch::setupMatrixRHS, futures[c],
                 &isSomaFlag, sizeof(isSomaFlag), &v[n-1], sizeof(v[n-1]));
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
int Branch::gaussianBackTriangulation_handler(const char * isSoma, const size_t isSoma_size)
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
        hpx_call(local->branches[c], Branch::gaussianBackTriangulation, futures[c], &isSomaFlag, isSoma_size);
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
    if (!*isSoma)
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
int Branch::gaussianFwdSubstitution_handler(const int nargs, const void *args[], const size_t sizes[])
{
    neurox_hpx_pin(Branch);
    assert(nargs==2);
    const char * isSoma = (const char *) args[0];
    const double *parentRHS = (const double *) args[1];

    double *b   = local->b;
    double *d   = local->d;
    double *rhs = local->rhs;
    int n = local->n;

    if(*isSoma)
    {
        rhs[0] /= d[0];
    }
    else for (int i=0; i<n; i++)
    {
        rhs[i] -= b[i] * (i==0 ? *parentRHS : rhs[i-1]);
        rhs[i] /= d[i];
    }

    char isSomaFlag=0;
    double childrenRHS=rhs[n-1];
    neurox_hpx_recursive_branch_sync(Branch::gaussianFwdSubstitution, &isSomaFlag, sizeof(isSomaFlag), &childrenRHS, sizeof(childrenRHS));
    neurox_hpx_unpin;
}

hpx_action_t Branch::setupMatrixLHS = 0;
int Branch::setupMatrixLHS_handler(const char * isSoma, const size_t isSoma_size)
{
    neurox_hpx_pin(Branch);
    int n = local->n;
    double *a   = local->a;
    double *b   = local->b;
    double *d   = local->d;
    int branchesCount = local->branchesCount;

    double returnValue = -1; //contribution to upper branch

    //for (i = i2; i < i3; ++i))
    if (!*isSoma)
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
        hpx_call(local->branches[c], Branch::setupMatrixLHS, futures[c], &isSomaFlag, isSoma_size);
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

hpx_action_t Branch::callNetReceiveFunction = 0;
int Branch::callNetReceiveFunction_handler(const int nargs, const void *args[], const size_t sizes[])
{
    neurox_hpx_pin(Branch);

    assert(nargs==3);
    /** nargs=3 where:
     * args[0] = function flag: NetReceiveInit (1) or NetReceive(0)
     * args[1] = actual time
     * args[2] = timestep size
     */
    neurox_hpx_recursive_branch_sync(Branch::callNetReceiveFunction,
        args[0], sizes[0], args[1], sizes[1], args[2], sizes[2]);

    const char isInitFunction = *(const char*) args[0];
    const double t = *(const double *) args[1];
    const double dt = *(const double *) args[2];

    //*sequential* delivery of received spikes
    while (!local->spikesQueue.empty() &&
           local->spikesQueue.top().deliveryTime < t+dt)
    {
        Spike spike = local->spikesQueue.top();
        if (spike.netcon->active)
        {
             int mechType = spike.netcon->mechType;
             getMechanismFromType(mechType).callNetReceiveFunction
                     (local, &spike, isInitFunction, t, dt);
        }
        local->spikesQueue.pop();
    }
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t Branch::callModFunction = 0;
int Branch::callModFunction_handler(const Mechanism::ModFunction * functionId_ptr, const size_t functionId_size)
{
    neurox_hpx_pin(Branch);
    neurox_hpx_recursive_branch_async_call(Branch::callModFunction, functionId_ptr, functionId_size);

    int topDependenciesCount = 0;
    for (int m=0; m<mechanismsCount; m++)
        if (mechanisms[m].isTopMechanism)
            topDependenciesCount++;
    assert(topDependenciesCount>0);

    //*parallel* execution of independent mechanisms
    hpx_t lco_mechs = hpx_lco_and_new(topDependenciesCount);
    for (int m=0; m<mechanismsCount; m++)
    {
        Mechanism & mech = mechanisms[m];
        if (mech.isTopMechanism)
        {
            int e = hpx_call(target, Mechanism::callModFunction, lco_mechs,
                              &mech.type, sizeof(mech.type),
                              functionId_ptr, functionId_size);
            assert(e==HPX_SUCCESS);
        }
    }
    hpx_lco_wait(lco_mechs);
    hpx_lco_delete(lco_mechs, HPX_NULL);

    neurox_hpx_recursive_branch_async_wait;
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
    neurox_hpx_recursive_branch_async_call(Branch::secondOrderCurrent);
    MechanismInstances * mechInstances= local->mechsInstances;
    for (int m=0; m<mechanismsCount; m++)
    {
        Mechanism & mech = mechanisms[m];
        if (mech.BAfunctions[Mechanism::ModFunction::alloc] ) //TODO used to be if == ion_alloc()
        {
            for (int i=0; i<mechInstances[m].instancesCount; i++)
            {
                double * data = &mechInstances[m].data[i*nparm];
                int & nodeIndex = mechInstances[m].nodesIndices[i];
                data[cur] += data[dcurdv] * local->rhs[nodeIndex]; // cur += dcurdv * rhs(ni[i])
            }
        }
    }
    neurox_hpx_recursive_branch_async_wait;
    neurox_hpx_unpin;
}

hpx_action_t Branch::queueSpikes = 0;
int Branch::queueSpikes_handler(const int nargs, const void *args[], const size_t sizes[])
{
    neurox_hpx_pin(Branch);
    assert(nargs==2);
    const int * preNeuronId = (const int *) args[0];
    const double * deliveryTime = (const double*) args[1];

    //netcvode::PreSyn::send()
    for (auto nc = local->netcons.at(*preNeuronId).begin();
         nc!=local->netcons.at(*preNeuronId).end(); nc++)
    {
      hpx_lco_sema_p(local->spikesQueueMutex);
      local->spikesQueue.push( Spike(*deliveryTime, &(*nc)) );
      hpx_lco_sema_v_sync(local->spikesQueueMutex);
    }
    neurox_hpx_unpin;
}

void Branch::registerHpxActions()
{
    neurox_hpx_register_action(2, Branch::setupMatrixRHS);
    neurox_hpx_register_action(1, Branch::setupMatrixLHS);
    neurox_hpx_register_action(1, Branch::updateV);
    neurox_hpx_register_action(1, Branch::gaussianBackTriangulation);
    neurox_hpx_register_action(2, Branch::gaussianFwdSubstitution);
    neurox_hpx_register_action(1, Branch::setV);
    neurox_hpx_register_action(1, Branch::callModFunction);
    neurox_hpx_register_action(2, Branch::callNetReceiveFunction);
    neurox_hpx_register_action(0, Branch::setupMatrixInitValues);
    neurox_hpx_register_action(2, Branch::init);
    neurox_hpx_register_action(2, Branch::queueSpikes);
}
