/*
# =============================================================================
# Copyright (c) 2015 - 2021-2022 Blue Brain Project/EPFL
#
# See top-level LICENSE file for details.
# =============================================================================.
*/
#include <math.h>  // fabs
#include <stdlib.h>
#include <algorithm>  // std::random_shuffle
#include <cstdlib>    // std::rand, std::srand
#include <ctime>      // std::time
#include <iterator>   // std::distance
#include <vector>     // std::vector

// This block enables compilation of the code with and without LIKWID in place
#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

#define LINEAR
#define SCHEDULER

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
typedef std::pair<Time, Event *> TimedEvent;
typedef int neuron_id_t;

const double sample = 1;

class Neuron {
 public:
  unsigned char *buffer;
  size_t buffer_size;

#ifdef LINEAR
  linear::Vector<Synapse> *v;            // out synapses
  linear::Map<neuron_id_t, NetconX> *m;  // in synapses
  linear::Map<neuron_id_t, Time> *n;     // 4 instances
  linear::PriorityQueue<neuron_id_t, Time, TimedEvent> *q;
#else
  std::vector<Synapse *> v;
  std::map<neuron_id_t, std::vector<NetconX *>> m;
  std::map<neuron_id_t, Time> n;
  std::priority_queue<TimedEvent, std::vector<TimedEvent>,
                      std::greater_equal<TimedEvent>>
      q;
#endif

  /// get size of data structure with padding
  size_t SizeOf(size_t size) {
    const int kSOAPadding = 4;
    const int layout = 0;
    return this->soa_padded_size<kSOAPadding>(size, layout);
  }

  Neuron() {}

  ~Neuron() {
    delete[] buffer;  // for Linear case it deletes all structures
#ifndef LINEAR
    v.clear();
    m.clear();
    n.clear();
    q = std::priority_queue<TimedEvent, std::vector<TimedEvent>,
                            std::greater_equal<TimedEvent>>();
#endif
  }

#ifdef LINEAR
  inline static linear::Map<neuron_id_t, NetconX> *get_m(Neuron *neuron) {
    return neuron->m;
  }
#else
  inline static std::map<neuron_id_t, std::vector<NetconX *>> *get_m(
      Neuron *neuron) {
    return &neuron->m;
  }
#endif

#ifdef LINEAR
  inline static linear::Map<neuron_id_t, Time> *get_n(Neuron *neuron) {
    return neuron->n;
  }
#else
  inline static std::map<neuron_id_t, Time> *get_n(Neuron *neuron) {
    return &neuron->n;
  }
#endif

#ifndef LINEAR
  Neuron(size_t buffer_size) {
    buffer = new unsigned char[buffer_size]();
    // remaining vars are set externally
  }
#else
  Neuron(size_t buffer_size, std::vector<Synapse *> vv,
         std::map<neuron_id_t, std::vector<NetconX *>> &mm,
         std::map<neuron_id_t, Time> &nn, size_t q_keys_count,
         neuron_id_t *q_keys, size_t *q_max_vals_per_key) {
    size_t v_size = SizeOf(linear::Vector<Synapse>::Size(vv.size()));

    size_t keys_count = mm.size();
    size_t *count_per_key = new size_t[keys_count];
    int i = 0;
    for (auto key_val : mm) count_per_key[i++] = key_val.second.size();
    size_t m_size = SizeOf(
        linear::Map<neuron_id_t, NetconX>::Size(keys_count, count_per_key));
    delete[] count_per_key;

    size_t n_size = SizeOf(linear::Map<neuron_id_t, Time>::Size(nn.size()));

    size_t q_size =
        SizeOf(linear::PriorityQueue<neuron_id_t, Time, TimedEvent>::Size(
            q_keys_count, q_max_vals_per_key));
    ;

    unsigned int buffer_it = buffer_size;
    buffer_size += q_size + v_size + m_size + n_size;
    buffer = new unsigned char[buffer_size]();

    v = (linear::Vector<Synapse> *)&buffer[buffer_it];
    new (v) linear::Vector<Synapse>(vv, &buffer[buffer_it]);
    buffer_it += v_size;

    m = (linear::Map<neuron_id_t, NetconX> *)&buffer[buffer_it];
    new (m) linear::Map<neuron_id_t, NetconX>(mm, &buffer[buffer_it]);
    buffer_it += m_size;

    n = (linear::Map<neuron_id_t, Time> *)&buffer[buffer_it];
    new (n) linear::Map<neuron_id_t, Time>(nn, &buffer[buffer_it]);
    buffer_it += n_size;

    q = (linear::PriorityQueue<neuron_id_t, Time, TimedEvent>
             *)&buffer[buffer_it];
    new (q) linear::PriorityQueue<neuron_id_t, Time, TimedEvent>(
        q_keys_count, q_keys, q_max_vals_per_key, &buffer[buffer_it]);
    buffer_it += q_size;
    assert(buffer_it == buffer_size);
  }
#endif

 private:
  template <int chunk>
  inline int soa_padded_size(int cnt, int layout) {
    int imod = cnt % chunk;
    if (layout == 1) return cnt;
    if (imod) {
      int idiv = cnt / chunk;
      return (idiv + 1) * chunk;
    }
    return cnt;
  }
};

int irand(int userBeg, int userEnd) {
  return rand() % ((userEnd - userBeg) + 1) + userBeg;
}

double benchmark(Neuron **neurons, const size_t neuron_count,
                 const Time sim_time, const size_t buffer_size,
                 const size_t *random_ids) {
#ifdef LINEAR
  size_t m_count;
  NetconX *m_val;
  std::vector<TimedEvent> q_events;

  linear::Map<neuron_id_t, NetconX> *map_m;
  linear::Map<neuron_id_t, NetconX> *map_n;
#else
  std::map<neuron_id_t, std::vector<NetconX *>> *map_m;
  std::map<neuron_id_t, Time> *map_n;
  std::vector<NetconX *> *m_val;
#endif

  // For a fair comparison of structs, we will
  // benchmark 4 steps per iteration interval.
  double dumb = 0;
  const Time dt = 0.025;
  Neuron *neuron;
  size_t *neuron_it = new size_t[neuron_count]();
  neuron_id_t neuron_id = -1;
  size_t i, j;

#ifndef SCHEDULER
  for (Time time = 0; time < sim_time; time += 0.1)
    for (size_t n = 0; n < neuron_count; n++) {
      neuron = neurons[n];

      map_m = Neuron::get_m(neuron);
      map_n = Neuron::get_n(neuron);

      // run through dataset to make sure it wont stay in cache
      for (size_t b = 0; b < buffer_size; b += 256) dumb += neuron->buffer[b];

      LIKWID_MARKER_START("benchmark");
      for (Time t = time; t < time + 0.1; t += dt, neuron_it[n]++) {
#else

  std::vector<double> neuron_times(neuron_count);
  for (size_t n = 0; n < neuron_count; n++) neuron_times.at(n) = 0;

  while (true) {
    int min_e = std::distance(
        neuron_times.begin(),
        std::min_element(neuron_times.begin(), neuron_times.end()));
    int max_e = std::distance(
        neuron_times.begin(),
        std::max_element(neuron_times.begin(), neuron_times.end()));
    Time min_time = neuron_times.at(min_e);
    Time max_time = neuron_times.at(max_e);

    if (min_time == max_time && fabs(sim_time - min_time) < 0.001)
      break;  // end of execution

    size_t n = min_e;
    Time t = min_time;
    neuron = neurons[n];
    double step_count = (double)irand(4, 20);
    if (neuron_times[n] + step_count * dt > sim_time)
      step_count = (sim_time - neuron_times[n]) / dt;
    double step_end_time = neuron_times[n] + step_count * dt;
    // printf("---- neuron %d does %.1f steps from %f to %f\n", n, step_count,
    // neuron_times.at(n), step_end_time);
    neuron_times[n] = step_end_time;

    map_m = Neuron::get_m(neuron);
    map_n = Neuron::get_n(neuron);

    // run through dataset to make sure it wont stay in cache
    for (size_t b = 0; b < buffer_size; b += 256) dumb += neuron->buffer[b];

    LIKWID_MARKER_START("benchmark");
    for (; t < step_end_time; t += dt, neuron_it[n]++) {
#endif

        // M: handle one random incoming netcon per neuron per step
        for (i = 0; i < neuron_count; i++) {
          neuron_id = random_ids[i];
#ifdef LINEAR
          map_m->At(neuron_id, m_count, m_val);
          for (j = 0; j < m_count; j++) dumb += m_val[j];
#else
        m_val = &map_m->at(neuron_id);
        for (j = 0; j < m_val->size(); j++) dumb += *m_val->at(j);
#endif
        }

        // one spike check at every 1 ms
        if (neuron_it[n] % 40 == 0) {
#ifdef LINEAR
          for (i = 0; i < neuron->v->Count(); i++) dumb += *neuron->v->At(i);
#else
        for (i = 0; i < neuron->v.size(); i++) dumb += *neuron->v.at(i);
#endif
        }

        // N: time of pre-synaptic pre-synaptic random id (4 arrays)
        for (i = 0; i < 4; i++) {
          for (j = 0; j < neuron_count; j++) {
            neuron_id = random_ids[j];
#ifdef LINEAR
            dumb += *map_n->At(neuron_id);
#else
          dumb += map_n->at(neuron_id);
#endif
          }
        }

        // Q: at every 1ms, create events
        if (neuron_it[n] % 40 == 0) {
          for (i = 0; i < neuron_count; i++) {
            // create events at "random" positions in queue
            Time til = t + dt + i * 0.01 * (i % 2 == 0 ? -1 : 1);
            neuron_id_t id = random_ids[i];
#ifdef LINEAR
            neurons[id]->q->Push(n, std::make_pair(til, (Event *)&sample));
#else
          neurons[id]->q.push(std::make_pair(til, (Event *)&sample));
#endif
            dumb++;
          }
        }

        // Q: deliver events for next step
#ifdef LINEAR
        neuron->q->PopAllBeforeTime(t + dt, q_events);
        for (auto q_it : q_events) {
          dumb += q_it.first + *q_it.second;
        }
#else
      while (!neuron->q.empty() && neuron->q.top().first <= t + dt) {
        auto q_it = neuron->q.top();
        dumb += q_it.first + *q_it.second;
        neuron->q.pop();
      }
#endif
      }  // end of steping
      LIKWID_MARKER_STOP("benchmark");
    }
  delete[] neuron_it;
  return dumb;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    printf("Usage: %s <neuron-count>\n", argv[0]);
    exit(1);
  }
  LIKWID_MARKER_INIT;
  const size_t neuron_count = atoi(argv[1]);
#ifdef LINEAR
  printf("Benchmark starting (LINEAR data structs) with %d neurons\n",
         (int)neuron_count);
#else
  printf("Benchmark starting (std data structs) with %d neurons\n",
         (int)neuron_count);
#endif

  // shuffled list of elements;
  std::srand(unsigned(12345));
  // std::srand ( unsigned ( std::time(0) ) );
  std::vector<size_t> random_ids(neuron_count);
  for (int id = 0; id < neuron_count; id++) random_ids[id] = id;
  std::random_shuffle(random_ids.begin(), random_ids.end());

  const size_t buffer_size = 32 * 1024 * 1024;  // 4MB
  const Time sim_time = 10;                     // ms
  const float netcons_per_syn = 5;

  std::vector<Neuron *> neurons(neuron_count);

  for (int n = 0; n < neuron_count; n++) {
    // synapses to outgoing (all) neurons at random times
    std::vector<Synapse *> vv(neuron_count);
    for (int i = 0; i < neuron_count; i++) vv.at(i) = (Synapse *)&sample;

#ifndef LINEAR
    neurons.at(n) = new Neuron(buffer_size);
    neurons.at(n)->v = vv;
    // std::copy (vv.begin(), vv.end(), neurons.at(n)->v.begin());
  }

  for (int n = 0; n < neuron_count; n++) {
#endif
    // map of incoming netcons per neurons and dependencies time
    std::map<neuron_id_t, std::vector<NetconX *>> mm;

    for (int i = 0; i < neuron_count; i++) {
      mm[i] = std::vector<NetconX *>(netcons_per_syn);
      for (int j = 0; j < netcons_per_syn; j++)
        mm.at(i).at(j) = (NetconX *)&sample;
    }
#ifndef LINEAR
    neurons.at(n)->m.insert(mm.begin(), mm.end());
  }

  for (int n = 0; n < neuron_count; n++) {
#endif
    // map of dependencies time
    std::map<neuron_id_t, Time> nn;
    for (int i = 0; i < neuron_count; i++) nn[i] = -1;

#ifndef LINEAR
    neurons.at(n)->n.insert(nn.begin(), nn.end());
#else
    std::vector<neuron_id_t> q_keys;
    std::vector<size_t> q_max_vals_per_key;
    for (int i = 0; i < neuron_count; i++) {
      q_keys.push_back(i);
      q_max_vals_per_key.push_back(neuron_count);
    }

    neurons.at(n) = new Neuron(buffer_size, vv, mm, nn,
                               neuron_count,  // q_keys_count
                               q_keys.data(), q_max_vals_per_key.data());
#endif
  }
  printf("Running %f msecs...\n", sim_time);
  clock_t begin = clock();
  double dumb = benchmark(neurons.data(), neurons.size(), sim_time, buffer_size,
                          random_ids.data());
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Finished: %f secs. Sum check=%.4f\n", time_spent, dumb);

  for (auto &neuron : neurons) delete neuron;

  LIKWID_MARKER_CLOSE;
  return 0;
}
