#include <stdio.h>
#include <string>
#include <set>
#include <map>
#include <list>
#include <tuple>
#include <algorithm>
#include <numeric>

#include "coreneuron/nrnconf.h"
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrnoc/nrnoc_decl.h"
#include "coreneuron/nrniv/netcon.h"
#include "coreneuron/nrniv/netcvode.h"
#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrniv/nrn_assert.h"
#include "coreneuron/nrniv/nrn_setup.h"
#include "coreneuron/utils/memory_utils.h"
#include "coreneuron/nrnoc/nrnoc_decl.h" //nrn_is_ion()
#include "coreneuron/nrniv/vrecitem.h"

#include "neurox/neurox.h"
#include "neurox/datatypes/Mechanism.h"
#include "neurox/datatypes/Branch.h"

using namespace std;
using namespace neurox::Input;
using namespace neurox::Input::Coreneuron;

void DataComparison::compareDataStructuresWithCoreNeuron(Branch * branch)
{
    int nrnThreadId = branch->soma->nrnThreadId;
    assert(sizeof(floble_t) == sizeof(double)); //only works with doubles!
    assert(branch->soma); //only non-branched neurons
    NrnThread & nt = nrn_threads[nrnThreadId];

    //make sure all morphology and mechs data is correct
    for (int i=0; i<nt._ndata; i++)
    {   assert(nt._data[i] == branch->nt._data[i]); }

    for (int i=0; i<6*branch->nt.end; i++)
    {   assert(nt._data[i] == branch->nt._data[i]); }

    for (offset_t i=0; i<branch->nt.end; i++)
    {
        assert(nt._actual_a[i] == branch->nt._actual_a[i]);
        assert(nt._actual_b[i] == branch->nt._actual_b[i]);
        assert(nt._actual_d[i] == branch->nt._actual_d[i]);
        assert(nt._actual_v[i] == branch->nt._actual_v[i]);
        assert(nt._actual_rhs[i] == branch->nt._actual_rhs[i]);
        assert(nt._actual_area[i] == branch->nt._actual_area[i]);
        if (branch->nt._v_parent_index)
        {  assert(nt._v_parent_index[i] == branch->nt._v_parent_index[i]); }
    }

    int mechCount=0;
    int vdataOffset=0;
    for (NrnThreadMembList* tml = nt.tml; tml!=NULL; tml = tml->next) //For every mechanism
    {
        int type = tml->index;
        int m = mechanismsMap[type];
        Memb_list * ml = tml->ml; //Mechanisms application to each compartment
        Memb_list & instances = branch->mechsInstances[m];
        assert(ml->nodecount == instances.nodecount);
        //assert(ml->_nodecount_padded == instance.instancesCount);
        short dataSize  = mechanisms[m]->dataSize;
        short pdataSize = mechanisms[m]->pdataSize;
        for (int n=0; n<ml->nodecount; n++) //for every mech instance
        {
            assert(ml->nodeindices[n]==instances.nodeindices[n]);
            for (int i=0; i<dataSize; i++)
            {   assert(ml->data[n*dataSize+i]==instances.data[n*dataSize+i]); }

            for (int i=0; i<pdataSize; i++)
            {
                int ptype = memb_func[type].dparam_semantics[i];
                assert(ml->pdata[n*pdataSize+i] == instances.pdata[n*pdataSize+i]);
            }

            /* We comment this because it runs for NULL presyn
            if (mechanisms[m]->pntMap)
            {
                //compare point_processes (index 1)
                Point_process * pp = (Point_process *) nt._vdata[vdataOffset+1];
                Point_process * pp2 = (Point_process *) branch->vdata[vdataOffset+1];
                assert(pp->_type == pp2->_type );
                assert(pp->_i_instance == pp2->_i_instance );
                assert(pp->_tid == pp2->_tid );
                assert(pp->_presyn == pp2->_presyn );
                vdataOffset+= mechanisms[m]->vdataSize;
            }
            */
        }
        mechCount++;
    }
    assert(mechCount==mechanismsCount);
}

void DataComparison::coreNeuronFinitialize2()
{
    nrn_finitialize(inputParams->voltage != 1000., inputParams->voltage);
}

void DataComparison::coreNeuronFinitialize() //can be deleted
{
    //nrn_finitialize
    int i;
    NrnThread* _nt;
    dt2thread(-1.);
    nrn_thread_table_check();
    clear_event_queue();
    nrn_spike_exchange_init();
    nrn_play_init(); /* Vector.play */
        ///Play events should be executed before initializing events
    for (i=0; i < nrn_nthread; ++i) {
        nrn_deliver_events(nrn_threads + i); /* The play events at t=0 */
    }
          for (_nt = nrn_threads; _nt < nrn_threads + nrn_nthread; ++_nt) {
            for (i=0; i < _nt->end; ++i) {
            (_nt->_actual_v[(i)]) = inputParams->voltage;
            }
          }
    for (i=0; i < nrn_nthread; ++i) {
        nrn_ba(nrn_threads + i, BEFORE_INITIAL);
    }

    /* the memblist list in NrnThread is already so ordered */
    for (i=0; i < nrn_nthread; ++i) {
        NrnThread* nt = nrn_threads + i;
        NrnThreadMembList* tml;
        for (tml = nt->tml; tml; tml = tml->next) {
            mod_f_t s = memb_func[tml->index].initialize;
            if (s) {
                // TODO COMMENTED (*s)(nt, tml->ml, tml->index);
            }
        }
    }
}

