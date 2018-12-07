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

typedef double Synapse;
typedef double Event;
typedef double NetconX;
typedef double Time;
typedef std::pair<Time, Event*> TimedEvent;
typedef int neuron_id_t;

const double sample = -1;

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

    Neuron() {};

    ~Neuron() {
        delete [] buffer; //for Linear case it deletes all structures
#ifndef LINEAR
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
#else 
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
	delete [] count_per_key;

        size_t n_size = SizeOf(linear::Map<neuron_id_t, Time>::Size(
                        nn.size()));

        size_t q_size = SizeOf(linear::PriorityQueue<neuron_id_t, Time, TimedEvent>::Size(
                        q_keys_count, q_max_vals_per_key));; 

        unsigned int buffer_it=buffer_size;
        buffer_size += q_size + v_size + m_size + n_size;
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

double benchmark(Neuron **neurons,
               const size_t neuron_count,
               const Time sim_time,
               const size_t synapse_count,
               const size_t buffer_size)
{

#ifdef LINEAR
    neuron_id_t *m_keys, *n_keys;
    size_t m_key_count, n_key_count;
    size_t m_count;
    NetconX *m_val;
    std::vector<TimedEvent> q_events;

    linear::Map<neuron_id_t, NetconX> * map_m;
    linear::Map<neuron_id_t, NetconX> * map_n;
#endif

    // For a fair comparison of structs, we will
    //benchmark 4 steps per iteration interval.
    double dumb=0;
    const Time dt = 0.025;
    int comm_it=0;
    Neuron * neuron;
    for (Time time=0; time<sim_time; time+=0.1, comm_it++)
    {
      for (int n=0; n<neuron_count; n++)
      {
        neuron = neurons[n];

        //run through the structure
        //for (size_t b=0; b<buffer_size; b+=256)
        //    dumb+= neuron->buffer[b];

#ifdef LINEAR
        q_events.clear();
        m_key_count = neuron->m->KeysCount();
        m_keys = neuron->m->Keys();
        n_key_count = neuron->n->KeysCount();
        n_keys = neuron->n->Keys();
        map_m = neuron->m;
        map_n = neuron->n;
#endif

        for (Time t=time; t<time+0.1; t+=dt)
        {
          if (comm_it % 100 == 0)
          {
            // V: if 10 ms have past (outgoing spike)
#ifdef LINEAR
            for (int i=0; i< neuron->v->Count(); i++)
                dumb += *neuron->v->At(i);
#else
            for (int i=0; i< neuron->v.size(); i++)
                dumb += *neuron->v.at(i);
#endif

            // M: f 10 ms have past (incoming netcon)
#ifdef LINEAR
            for (int k=0; k<m_key_count; k++)
            {
              map_m->At(m_keys[k], m_count, m_val);
              for (int i=0; i<m_count; i++)
                  dumb += m_val[i];
            }
#else
            for (std::pair<neuron_id_t, std::vector<NetconX*> > id_nc_it : neuron->m)
                for (NetconX *& nc : id_nc_it.second)
                  dumb += *nc;
#endif
          }



          // N: time or pre-synaptic id (4 arrays)
          for (int x=0; x<4; x++)
	  {
#ifdef LINEAR
            for (int k=0; k<n_key_count; k++)
                dumb += n_keys[k] + *map_n->At(n_keys[k]);
#else
            for (auto & n_it : neuron->n)
                 dumb += n_it.first + n_it.second;
#endif
	  }


          // Q: at every 3 ms, create events
          if (comm_it + n % 30 == 0)
          {
            for (size_t i=0; i<synapse_count; i++)
            {
              //create events at "random" positions in queue
              Time til = t + dt + i*0.000001*(i%2==0 ? -1 : 1);
#ifdef LINEAR
              neurons[n]->q->Push(i, std::make_pair(til, (Event*)&sample));
#else
              neurons[n]->q.push(std::make_pair(til, (Event*) &sample));
#endif
            }
          }


          // Q: deliver events for next step
#ifdef LINEAR
            neuron->q->PopAllBeforeTime(t+dt, q_events);
            for (auto q_it : q_events)
              dumb += q_it.first + *q_it.second;
#else
            while (!neuron->q.empty() &&
                   neuron->q.top().first <= t+dt) {
              auto q_it = neuron->q.top();
              dumb += q_it.first  + *q_it.second;
              neuron->q.pop();
	    }
#endif
        }
      }
    }
    return dumb;
}

int main(int argc, char** argv)
{
    if (argc!=2)
    {
      printf("Usage: %s <neuron-count>\n", argv[0]);
      exit(1);
    }
    const size_t neuron_count = atoi(argv[1]);
#ifdef LINEAR
    printf("Benchark starting (LINEAR data structs) with %d neurons\n", (int) neuron_count);
#else
    printf("Benchark starting (std data structs) with %d neurons\n", (int) neuron_count);
#endif

    const size_t buffer_size = 4*1024*1024; //4MB
    const Time sim_time=10; //ms
    const float netcons_per_syn=5;

      const size_t synapse_count = std::min(0.8*(double)neuron_count, 10000./(double)netcons_per_syn);

      std::vector<Neuron*> neurons(neuron_count);

      for (int n=0; n<neuron_count; n++)
      {
        //synapses to outgoing neurons at random times
        std::vector<Synapse*> vv(synapse_count);
        for (int i=0; i<synapse_count; i++)
          vv.at(i) = (Synapse*)&sample;

#ifndef LINEAR
        neurons.at(n)=new Neuron(buffer_size);
        neurons.at(n)->v = vv;
        //std::copy (vv.begin(), vv.end(), neurons.at(n)->v.begin());
      }

      for (int n=0; n<neuron_count; n++){
#endif
        //map of incoming netcons per neurons and dependencies time
        std::map<neuron_id_t, std::vector<NetconX*> > mm;

        for (int i=0; i<synapse_count; i++)
        {
          mm[i] = std::vector<NetconX*>(netcons_per_syn);
          for (int j=0; j<netcons_per_syn; j++)
            mm.at(i).at(j)=(NetconX*)&sample;
        }
#ifndef LINEAR
        neurons.at(n)->m.insert(mm.begin(), mm.end());
      }

      for (int n=0; n<neuron_count; n++)
      {
#endif
        //map of dependencies time
        std::map<neuron_id_t, Time> nn;
        for (int i=0; i<synapse_count; i++)
          nn[i] = -1;

#ifndef LINEAR
        neurons.at(n)->n.insert(nn.begin(), nn.end());
#else
        std::vector<neuron_id_t> q_keys;
        std::vector<size_t> q_max_vals_per_key;
        for (neuron_id_t pre_id=0; pre_id<synapse_count; pre_id++)
        {
          q_keys.push_back(pre_id);
          q_max_vals_per_key.push_back(5);
        }

        neurons.at(n)=new Neuron(buffer_size,
                             vv, mm, nn,
                             synapse_count, //q_keys_count
                             q_keys.data(),
                             q_max_vals_per_key.data()
                             );
#endif
      }
      printf("Running %f msecs...\n", sim_time);
      clock_t begin = clock();
      double dumb = benchmark(neurons.data(), neurons.size(), sim_time, synapse_count, buffer_size);
      clock_t end = clock();
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
      printf("Finished: %f secs. Sum check=%.4f\n", time_spent, dumb);

      for (auto & neuron : neurons)
          delete neuron;

    return 0;
}
