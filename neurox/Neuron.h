#pragma once

#include "neurox/neurox.h"
#include "libhpx/libhpx.h"
#include <deque>

namespace neurox
{

/**
 * @brief The Neuron class
 * Represents the soma structure and outgoing spikes network of a neuron
 */
class Neuron
{
  public:
    Neuron() = delete;
    ~Neuron();

    Neuron(neuron_id_t neuronId, floble_t APthreshold, floble_t * thvar_ptr);

    neuron_id_t gid;     ///> neuron global id
    floble_t threshold;  ///> Action Potential threshold (PreSyn->_threshold)
    floble_t* thvar_ptr;       ///> pointer to variable holding threashold var (typically RHS[0] is AIS's first compartment)
    floble_t refractoryPeriod; ///> refractory period

    void setupTreeMatrixMinimal(); ///>set_tree_matrix_minimal
    bool checkAPthresholdAndTransmissionFlag (floble_t v); ///> checks if AP threshold has been reached and whether spikes can be transmitted  (PreSynHelper->flag)

    // the outgoing neurons:
    class Synapse
    {
      public:
        Synapse()=delete;
        Synapse(hpx_t branchAddr, floble_t minDelay, hpx_t topBranchAddr = HPX_NULL, int destinationGid=-1);
        ~Synapse();
        hpx_t branchAddr;          ///> address of destination
        hpx_t topBranchAddr; ///> addres of top-branch (soma) of destination neuron
#ifndef NDEBUG
        int destinationGid;
#endif
        floble_t nextNotificationTime; ///> next time this post-syn neuron needs to be informed of my actual time
        floble_t minDelay; ///>interval of notification in case of no spykes
                           ///(fastest Netcon from current neuron to dependant-neuron)
        hpx_t previousSpikeLco; ///>lco controlling spikes delivery
    };

    hpx_t sendSpikes(floble_t t); ///> fires AP, returns LCO for sent synapses
    void addSynapse(Synapse * target);///> add hpx address of post-synaptic branch
    size_t getSynapseCount(); ///> get size of vector synapse

    class CommunicationBarrier
    {
      public:
        CommunicationBarrier();
        ~CommunicationBarrier();

        hpx_t allSpikesLco; ///> LCO for all spikes of previous Comm Step (for fixed step methods and debug)
        static constexpr int commStepSize = 4; ///> Fixed communication step size
    } * commBarrier;

    class TimeDependencies  ///from incoming neuron connections
    {
      public:
        TimeDependencies();
        ~TimeDependencies();

        void waitForTimeDependencyNeurons(floble_t t, floble_t dt, int gid);
        void sendSteppingNotification(floble_t t, floble_t dt, int gid, std::vector<Synapse*> & synapses); ///> inform my outgoing-connection neurons that I stepped
        void updateTimeDependency(neuron_id_t srcGid, floble_t dependencyNotificationTime, neuron_id_t myGid = -1, bool initialization = false);
        floble_t getDependenciesMinTime();
        size_t getDependenciesCount();
        void increseDependenciesTime(floble_t t);

        static constexpr floble_t notificationIntervalRatio = 1; ///> ration of notification interval (0,1]
        static constexpr double teps =1e-8; ///>time-epsilon to correct wrong delivery of events due to floating point rounding

      private:
        std::map<neuron_id_t, floble_t> dependenciesMap; ///> map to previous structure (pointing to vector values)
        libhpx_cond_t dependenciesWaitCondition;
        libhpx_mutex_t dependenciesLock;
        floble_t dependenciesTimeNeuronWaitsFor;

    } * timeDependencies;

    class SlidingTimeWindow
    {
      public:
        SlidingTimeWindow();
        ~SlidingTimeWindow();

        //set by initNodeLevelInformaion
        static int reductionsPerCommStep;

        //for node level reduction only (initialized by initNodeLevelInformaion)
        class AllReduceLocality
        {
          public:
            static std::vector<hpx_t> * localityNeurons;
            static hpx_t * allReduceFuture;
            static hpx_t * allReduceLco;
            static int * allReduceId;

            static hpx_action_t subscribeAllReduce;
            static hpx_action_t unsubscribeAllReduce;

            static int subscribeAllReduce_handler(const hpx_t*, const size_t);
            static int unsubscribeAllReduce_handler(const hpx_t*, const size_t);
        };

        //initiated by constructor (one per neuron)
        std::queue<hpx_t> spikesLcoQueue;
        hpx_t *allReduceFuture;
        hpx_t *allReduceLco;
        int *allReduceId;

        //all-reduce functions
        static hpx_action_t init;
        static hpx_action_t reduce;
        static hpx_action_t subscribeAllReduce;
        static hpx_action_t unsubscribeAllReduce;
        static hpx_action_t setReductionsPerCommStep;

        static void init_handler(const void*, const size_t);
        static void reduce_handler(void* rhs, const void* lhs, const size_t);
        static int subscribeAllReduce_handler(const hpx_t*, const size_t);
        static int unsubscribeAllReduce_handler(const hpx_t*, const size_t);
        static int setReductionsPerCommStep_handler(const int*, const size_t);
    } * slidingTimeWindow;

    std::vector<Synapse*> synapses; ///> out-going synapse information

    static void registerHpxActions(); ///> Register all HPX actions

  private:
    //the outgoing neuron connections:
    hpx_t synapsesMutex; ///> protects synapses
    bool synapsesTransmissionFlag; ///> PreSynHelper* psh -> flag (to flag whether spikes for a given AP have been sent or not
}; //Neuron
}; //namespace neurox
