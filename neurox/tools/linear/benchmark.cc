#include <stdlib.h>

//#define LINEAR

#ifdef LINEAR
  #include "map.h"
  #include "priority_queue.h"
  #include "vector.h"
  using namespace neurox::tools;
#else
  #include <map>
  #include <queue>
  #include <vector>
#endif

typedef double Synapse;
typedef double Event;
typedef unsigned NetconX; //offset to element in data array
typedef std::pair<double, Event*> TimedEvent;
typedef double Time;
typedef int neuron_id_t;

class Neuron
{
  public:

    unsigned char *buffer;
    size_t buffer_size;

#ifdef LINEAR
    linear::Vector<Synapse> *v; //out synapses
    linear::Map<neuron_id_t, NetconX> *m; //in synapses
    linear::Map<neuron_id_t, Time> *n; //4 instances
    linear::PriorityQueue<neuron_id_t, Time, TimedEvent>* q;
#else
    std::vector<Synapse*> v;
    std::map<neuron_id_t, std::vector<NetconX*> > m;
    std::map<neuron_id_t, Time > n;
    std::priority_queue<TimedEvent, std::vector<TimedEvent>, std::greater_equal<TimedEvent> > q;
#endif

    /// get size of data structure with padding
    size_t SizeOf(size_t size) {
      const int kSOAPadding=4;
      const int layout=0;
      return this->soa_padded_size<kSOAPadding>(size, layout);
    }

    Neuron() = delete;

    ~Neuron() {
        delete [] buffer;
#ifdef LINEAR
        delete v;
        delete m;
        delete n;
        delete q;
#else
        v.clear();
        m.clear();
        n.clear();
        q=std::priority_queue<TimedEvent, std::vector<TimedEvent>, std::greater_equal<TimedEvent>>();
#endif
    }


#ifndef LINEAR
    Neuron(size_t buffer_size)
    {
        buffer = new unsigned char[buffer_size]();
        //remaining vars are set externally
    }
#else LINEAR
    Neuron(size_t buffer_size,
           std::vector<Synapse*> vv,
           std::map<neuron_id_t, std::vector<NetconX*> >& mm,
           std::map<neuron_id_t, Time>& nn,
           size_t q_keys_count,
           neuron_id_t *q_keys,
           size_t *q_max_vals_per_key
           ) {

        size_t v_size = SizeOf(linear::Vector<Synapse>::Size(vv.size()));

        size_t keys_count = mm.size();
        size_t * count_per_key = new size_t[keys_count];
        int i=0;
        for (auto key_val : mm)
             count_per_key[i++]=key_val.second.size();
        size_t m_size = SizeOf(linear::Map<neuron_id_t, NetconX>::Size(
                          keys_count, count_per_key));

        size_t n_size = SizeOf(linear::Map<neuron_id_t, Time>::Size(
                        nn.size()));

        size_t q_size = SizeOf(linear::PriorityQueue<neuron_id_t, Time, TimedEvent>::Size(
                        q_keys_count, q_max_vals_per_key));;

        unsigned int buffer_it=buffer_size;
        buffer_size += q_size + v_size + m_size;
        buffer = new unsigned char[buffer_size]();

        v = (linear::Vector<Synapse>*) &buffer[buffer_it];
        new (v) linear::Vector<Synapse>(vv, &buffer[buffer_it]);
        buffer_it += v_size;

        m = (linear::Map<neuron_id_t, NetconX> *)&buffer[buffer_it];
        new (m) linear::Map<neuron_id_t, NetconX>(mm, &buffer[buffer_it]);
        buffer_it += m_size;

        n = (linear::Map<neuron_id_t, Time> *)&buffer[buffer_it];
        new (n) linear::Map<neuron_id_t, Time>(nn, &buffer[buffer_it]);
        buffer_it += n_size;

        q = (linear::PriorityQueue<neuron_id_t, Time, TimedEvent> *)&buffer[buffer_it];
        new (q) linear::PriorityQueue<neuron_id_t, Time, TimedEvent>(
            q_keys_count, q_keys, q_max_vals_per_key, &buffer[buffer_it]);
        buffer_it += q_size;
        assert(buffer_it==buffer_size);
    }
#endif

    /// Random integer number between min and max
    static int irand(int min, int max)
    {
        return min+ (std::rand() % (max-min + 1));
    }

  private:
    template <int chunk>
    inline int soa_padded_size(int cnt, int layout) {
        int imod = cnt % chunk;
        if (layout == 1)
            return cnt;
        if (imod) {
            int idiv = cnt / chunk;
            return (idiv + 1) * chunk;
        }
        return cnt;
    }
};

void benchmark(const Neuron *neurons,
               const size_t neuron_count,
               const Time sim_time)
{

    // For a fair comparison of structs, we will
    //benchmark 4 steps per iteration interval.
    double dumb;
    const Time dt = 0.025;
    int comm_it=0;
    for (Time time=0; time<sim_time; time+=0.1, comm_it++)
    {
      for (int n=0; n<neuron_count; n++)
      {
        const Neuron * neuron = &neurons[n];

        for (Time t=time; t<time+0.1; t+=dt)
        {
          if (comm_it % 100 == 0)
          {
            // V: if 10 ms have past (outgoing spike)
#ifdef LINEAR
            for (auto & v_it : &(neuron->v))
                dumb += v_it;
#else
            for (auto & v_it : neuron->v)
                dumb += *v_it;
#endif

            // M: f 10 ms have past (incoming netcon)
#ifdef LINEAR
            for (auto & id_nc_it : neuron->m)
                for (auto & nc : id_nc_it.second)
                  dumb += *nc;
#else
            while (!neuron->q.empty() &&
                   neuron->q.top().first <= t+dt) {
              auto q_it = neuron->q.top();
              dumb += q_it.first  + *q_it.second;
            }
#endif
          }


          // N: time or pre-synaptic id (4 arrays)
          for (int x=0; x<4; x++)
#ifdef LINEAR
            for (auto & n_it : neuron->n)
                 dumb += n_it.first + n_it.second;
#else
            for (auto & n_it : neuron->n)
                 dumb += n_it.first + n_it.second;
#endif

          // Q: deliver events for next step
          for (int x=0; x<4; x++)
#ifdef LINEAR
            for (auto & n_it : neuron->n)
                 dumb += n_it.first + n_it.second;
#else
            for (auto & n_it : neuron->n)
                 dumb += n_it.first + n_it.second;
#endif
        }
      }
    }
}

int main()
{
#ifdef LINEAR
    printf("Benchark starting (LINEAR data structs)\n");
#else
    printf("Benchark starting (std data structs)\n");
#endif

    const size_t buffer_size = 3*1024*1024; //3MB
    const Time sim_time=10; //ms
    const float netcons_per_syn=5;

    for (int scale=0; scale<=1; scale*=2)
    {
      const size_t neuron_count = 10*scale;
      const size_t synapse_count = std::min(0.8*neuron_count, 10000./netcons_per_syn);

      std::vector<Neuron> neurons(neurons);

      for (int n=0; n<neuron_count; n++)
      {
        //synapses to outgoing neurons at random times
        std::vector<Synapse*> vv(synapse_count);
        for (int i=0; i<vv.size(); i++)
          vv.at(i) = new Synapse;

#ifndef LINEAR
        neurons.at(n)=Neuron(buffer_size);
        std::copy (vv.begin(), vv.end(), neurons.at(n).v.begin());
      }

      for (int n=0; n<neuron_count; n++){
#endif
        //map of incoming netcons per neurons and dependencies time
        std::map<neuron_id_t, std::vector<NetconX*> > mm;

        for (int i=0; i<synapse_count; i++)
        {
          size_t pre_id = Neuron::irand(0,neuron_count);
          while (mm.find(pre_id) != mm.end())
            pre_id = Neuron::irand(0,neuron_count);

          //map of netcons per incomming synapse
          mm.at(pre_id) = std::vector<NetconX*>(netcons_per_syn);
          for (int i=0; i<netcons_per_syn; i++)
          {
              int offset = Neuron::irand(0,buffer_size);
              mm.at(pre_id).at(i) = new NetconX(offset);
          }
        }
#ifndef LINEAR
        neurons.at(n).m.insert(mm.begin(), mm.end());
      }

      for (int n=0; n<neuron_count; n++)
      {
#endif
        //map of dependencies time
        std::map<neuron_id_t, Time> nn;
        for (int i=0; i<synapse_count; i++)
        {
          size_t pre_id = Neuron::irand(0,neuron_count);
          nn.at(pre_id) = -1;
        }

#ifndef LINEAR
        neurons.at(n).n.insert(nn.begin(), nn.end());
      }

      for (int n=0; n<neuron_count; n++)
      {
#endif

        std::vector<std::pair<neuron_id_t, Time>> evs;
#ifdef LINEAR
        std::vector<neuron_id_t> q_keys;
        std::vector<size_t> q_max_vals_per_key;
#endif
        for (neuron_id_t pre_id=0; pre_id<neuron_count; pre_id++)
        {
#ifdef LINEAR
          q_keys.push_back(pre_id);
          q_max_vals_per_key.push_back(100);
#endif
          for (int i=0; i<synapse_count; i++)
            for (Time t=0; t<sim_time; t+=1)
              evs.push_back(std::make_pair(pre_id, t));
        }

        for (auto ev : evs)
        {
            Event *e = new Event;
#ifdef LINEAR
            neurons.at(n).q->Push(ev.first, std::make_pair(ev.second, e));
#else
            neurons.at(n).q.push(std::make_pair(ev.second, e));
#endif
        }

#ifdef LINEAR
        neurons.at(n)=Neuron(buffer_size,
                             vv, mm, nn,
                             synapse_count, //q_keys_count
                             q_keys.data(),
                             q_max_vals_per_key.data()
                             );
#endif
      }
      printf("Running scale %d\n", scale);
      benchmark(neurons.data(), neurons.size(), sim_time);
    }

    return 0;
}
