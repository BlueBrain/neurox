#include <stdlib.h>

#define LINEAR

#ifdef LINEAR
  #include "map.h"
  #include "priority_queue.h"
  #include "vector.h"
  using namespace neurox::tools;
#else
  #include <map>
  #include <priority_queue>
  #include <vector>
#endif

typedef std::pair<int, double> Synapse;
typedef double NetconX;
typedef std::pair<double, double> TimedEvent;
typedef double Time;
typedef int neuron_id_t;

class Neuron
{
  public:

    unsigned char *dumb;

#ifdef LINEAR
    linear::Vector<Synapse> *v;
    linear::Map<neuron_id_t, NetconX> *m;
    linear::PriorityQueue<neuron_id_t, Time, TimedEvent>* q;
#else
    std::vector<Synapse> *v;
    std::map<neuron_id_t, std::vector<NetconX*> > *m;
    std::priority_queue<TimedEvent, std::vector<TimedEvent>, std::greater_equal<TimedEvent> > *q;
#endif

    Neuron(): dumb(nullptr) {}
    ~Neuron() { delete dumb; }


    Neuron(size_t size_MB)
    {
    }

    static void main()
    {

    }
};
