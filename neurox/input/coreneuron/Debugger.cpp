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

using namespace std;
using namespace neurox::Input;
using namespace neurox::Input::Coreneuron;

void Debugger::compareMechanismsFunctionPointers( std::list<NrnThreadMembList*> & uniqueMechs)
{
    printf("NDEBUG::comparing Mechanisms functions...\n");
    for (auto & tml : uniqueMechs)
    {
        Memb_func & mf_cn = memb_func[tml->index]; //coreneuron
        Memb_func & mf_nx = mechanisms[mechanismsMap[tml->index]]->membFunc; //neurox
        if (tml->index != CAP)
        { //we call MechFunctions::jacobCapacitance and currentCapacitance
          //so value is set. Coreneuron sets to null and calls nrn_cur_capacitance instead
          assert (mf_cn.jacob == mf_nx.jacob);
          assert (mf_cn.current == mf_nx.current);
        }
        //assert (mf_cn.alloc == mf_nx.alloc); //TODO never used?
        assert (mf_cn.destructor == mf_nx.destructor);
        assert (mf_cn.initialize == mf_nx.initialize);
        assert (mf_cn.state == mf_nx.state);
        assert (mf_cn.thread_cleanup_ == mf_nx.thread_cleanup_);
        assert (mf_cn.thread_mem_init_ == mf_nx.thread_mem_init_);
        assert (mf_cn.thread_table_check_ == mf_nx.thread_table_check_);
    }

}

hpx_action_t Debugger::compareBranch = 0;
int Debugger::compareBranch_handler()
{
    neurox_hpx_pin(Branch);
    Branch * branch = local;
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
    neurox_hpx_unpin;
}

void Debugger::coreNeuronFinitialize()
{
    nrn_finitialize(inputParams->voltage != 1000., inputParams->voltage);
}

void Debugger::registerHpxActions()
{
    neurox_hpx_register_action(0, Debugger::compareBranch);
}
