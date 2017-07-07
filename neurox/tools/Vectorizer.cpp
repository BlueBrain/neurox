#include "neurox/neurox.h"

using namespace std;
using namespace neurox;

//TODO replace these functions by Coreneuron's "coreneuron/nrniv/memory.h"?

size_t Tools::Vectorizer::sizeof_(size_t size)
{
    size_t size_padded = size % NEUROX_MEM_ALIGNMENT ==0 ? size : size + (NEUROX_MEM_ALIGNMENT - size % NEUROX_MEM_ALIGNMENT);
    assert(size_padded % NEUROX_MEM_ALIGNMENT == 0);
    return size_padded;
}

void Tools::Vectorizer::vectorize(Branch * b)
{
   //NOTE: arrays with memory-aligned allocation:
   //ml->pdata, nt->data, ml->nodeindices, v->parent_index, and ml->shadow_*;
   //data needs to add gaps, pdata needs new offset values

   assert(memb_func);

   //get total counts
   int N = b->nt->end;
   size_t totalDataSize  = 6*sizeof_(N);
   assert(totalDataSize %NEUROX_MEM_ALIGNMENT==0);

   std::vector<int> dataOffsets;
   for (int m=0; m<mechanismsCount; m++)
   {
       b->mechsInstances[m]._nodecount_padded = sizeof_(b->mechsInstances[m].nodecount);
       totalDataSize  += b->mechsInstances[m]._nodecount_padded * mechanisms[m]->dataSize;
   }

   //old thvar offset
   int thvar_idx = b->thvar_ptr ? b->thvar_ptr - &b->nt->_actual_v[0] : -1; //compartment id (typically 2)
   assert(thvar_idx==-1 || thvar_idx>=0 && thvar_idx < N);

   double* dataNew = new_<double>(totalDataSize);
   size_t dataNewOffset =0;

   //add padding to data for RHS, D, A, B, V and area
   for (int i=0; i<6; i++)
   {
        switch(i)
        {
            case 0: b->nt->_actual_rhs  = &dataNew[dataNewOffset]; break;
            case 1: b->nt->_actual_d    = &dataNew[dataNewOffset]; break;
            case 2: b->nt->_actual_a    = &dataNew[dataNewOffset]; break;
            case 3: b->nt->_actual_b    = &dataNew[dataNewOffset]; break;
            case 4: b->nt->_actual_v    = &dataNew[dataNewOffset]; break;
            case 5: b->nt->_actual_area = &dataNew[dataNewOffset]; break;
        }

        for (size_t j=0; j<sizeof_(b->nt->end); j++, dataNewOffset++)
            if (j < b->nt->end)
            {
                dataNew[dataNewOffset] = b->nt->_data[i*N+j];
                dataOffsets.push_back(dataNewOffset);
            }
   }

   //convert AP-threshold pointer
   b->thvar_ptr = b->thvar_ptr ? &b->nt->_actual_v[thvar_idx] : nullptr;

   //add padding and converting AoS->SoA in nt->data and update ml->data for mechs instances
   for (int m=0; m<neurox::mechanismsCount; m++)
   {
       Memb_list * instances = &b->mechsInstances[m];
       double* instanceDataNew = &dataNew[dataNewOffset];
       for (size_t d=0; d<mechanisms[m]->dataSize; d++) //for every variable
           for (int i=0; i<instances->_nodecount_padded; i++, dataNewOffset++) //for every node
               if (i < instances->nodecount)
               {
                   int oldOffset = mechanisms[m]->dataSize*i + d; //d-th variable in i-th instance
                   dataNew[dataNewOffset] = instances->data[oldOffset];
                   dataOffsets.push_back(dataNewOffset);
               }
       instances->data = instanceDataNew;
   }
   dataOffsets.resize(dataOffsets.size());

   //convert VecPlay continuous pointers
   for (int v=0; v<b->nt->n_vecplay; v++)
   {
       VecPlayContinuousX * vc = (VecPlayContinuousX*) b->nt->_vecplay[v];
       int vc_offset = vc->pd_ - &b->nt->_data[0];
       assert(vc_offset>=0 && vc_offset <= b->nt->_ndata);
       int vc_offset_new = dataOffsets.at(vc_offset);
       vc->pd_ = &dataNew[ vc_offset_new ];
   }

   delete_(b->nt->_data);
   b->nt->_data  = dataNew;

   //add padding and update ml->pdata for mechs instances
   for (int m=0; m<mechanismsCount; m++)
   {
       Memb_list * instances = &b->mechsInstances[m];
       int totalPDataSize = instances->_nodecount_padded * mechanisms[m]->pdataSize;
       int* pdataNew = new_<int>(totalPDataSize);
       size_t pdataNewOffset =0;
       for (size_t p=0; p<mechanisms[m]->pdataSize; p++) //for every pointer
           for (int i=0; i<instances->_nodecount_padded; i++, pdataNewOffset++) //for every node
           {
               if (i < instances->nodecount)
               {
                   int oldOffset = mechanisms[m]->pdataSize*i + p; //p-th pointer in i-th instance
                   pdataNew[pdataNewOffset] = instances->pdata[oldOffset];

                   int ptype = memb_func[mechanisms[m]->type].dparam_semantics[p];
                   bool isPointer = ptype==-1 || (ptype>0 && ptype<1000);
                   if (isPointer) //true data pointer to area or ion
                   {
                     int & pd = pdataNew[pdataNewOffset];
                     pd = dataOffsets.at(pd); //point to new offset (with gaps)
                   }
               }
               else
               {
                   pdataNew[pdataNewOffset]=-1;
               }
           }
        delete_(instances->pdata);
        instances->pdata = pdataNew;
    }
}
