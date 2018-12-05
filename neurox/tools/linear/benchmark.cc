#include <stdlib.h>

#define LINEAR

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

typedef std::pair<int, double> Synapse;
typedef unsigned NetconX; //offset to element in data array
typedef std::pair<double, double> TimedEvent;
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

    Neuron(size_t buffer_size,
           std::vector<Synapse*> vv,
           std::map<neuron_id_t, std::vector<NetconX*> >& mm,
           std::map<neuron_id_t, Time>& nn,
           size_t q_keys_count,
           neuron_id_t *q_keys,
           size_t *q_max_vals_per_key
           ) {

#ifdef LINEAR
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
#endif
        buffer = new unsigned char[buffer_size];

#ifdef LINEAR
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
#else
        v = vv;
        n = nn;
        m = mm;
        //q = qq; //WRONG
#endif

    }

    /// Random number between min and max
    static int rand(int min, int max)
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
               const float sim_time)
{

}

int main()
{
#ifdef LINEAR
    printf("Benchark starting (LINEAR data structs)\n");
#else
    printf("Benchark starting (std data structs)\n");
#endif

    const size_t buffer_size = 3*1024*1024; //3MB
    const float sim_time=10; //ms
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
          vv.at(i) = new Synapse(Neuron::rand(0,neuron_count), -1);

        //map of incoming netcons per neurons
        std::map<neuron_id_t, std::vector<NetconX*> > mm(synapse_count);

        //map of dependencies time
        std::map<neuron_id_t, Time> nn;
        for (int i=0; i<synapse_count; i++)
        {
          size_t pre_id = Neuron::rand(0,neuron_count);
          while (mm.find(pre_id) != mm.end())
            pre_id = new Synapse(Neuron::rand(0,neuron_count));

          //map of netcons per incomming synapse
          mm.at(pre_id) = new std::vector<NetconX*>(netcons_per_syn);
          for (int i=0; i<netcons_per_syn; i++)
          {
              int offset = Neuron::rand(0,buffer_size);
              mm.at(pre_id).at(i) = new NetconX(offset);
          }

          //map of dependencies time
          nn.at(pre_id) = -1;
        }

        //events queue
        size_t q_keys_count = synapse_count;
        std::vector<neuron_id_t> q_keys;
        std::vector<size_t> q_max_vals_per_key;

        neurons.at(n)=Neuron(buffer_size,
                             vv, mm, nn,
                             q_keys_count,
                             q_keys.data(),
                             q_max_vals_per_key.data()
                             );
      }
      benchmark(neurons.data(), neurons.size(), sim_time);
    }

    return 0;
}
