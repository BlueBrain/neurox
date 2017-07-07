#include "neurox/neurox.h"

using namespace std;
using namespace neurox;

size_t Tools::Vectorizer::sizeof_(size_t size)
{
    size_t size_padded = size==0 ? 0 : size + (NEUROX_MEM_ALIGNMENT - size % NEUROX_MEM_ALIGNMENT);
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
   size_t totalDataSize = 6*sizeof_(N);
   assert(totalDataSize %NEUROX_MEM_ALIGNMENT==0);
   //int *  dataNewOffset = new int[b->nt->_ndata];
   std::vector<int> dataNewOffset;
   for (int m=0; m<mechanismsCount; m++)
       totalDataSize += b->mechsInstances[m].nodecount * sizeof_(mechanisms[m]->dataSize);

   //old thvar offset
   int thvar_idx = b->thvar_ptr ? b->thvar_ptr - &b->nt->_actual_v[0] : -1; //compartment id (typically 2)
   assert(thvar_idx==-1 || thvar_idx>=0 && thvar_idx < N);

   double* data_new = new_<double>(totalDataSize);
   size_t dataOffset =0;

   //add padding to data for RHS, D, A, B, V and area
   for (int i=0; i<6; i++)
   {
        switch(i)
        {
            case 0: b->nt->_actual_rhs  = &data_new[dataOffset]; break;
            case 1: b->nt->_actual_d    = &data_new[dataOffset]; break;
            case 2: b->nt->_actual_a    = &data_new[dataOffset]; break;
            case 3: b->nt->_actual_b    = &data_new[dataOffset]; break;
            case 4: b->nt->_actual_v    = &data_new[dataOffset]; break;
            case 5: b->nt->_actual_area = &data_new[dataOffset]; break;
        }

        for (size_t j=0; j<sizeof_(b->nt->end); j++, dataOffset++)
            if (j < b->nt->end)
            {
                data_new[dataOffset] = b->nt->_data[i*N+j];
                //dataNewOffset[i*N+j] = dataOffset;
                dataNewOffset.push_back(dataOffset);
            }
   }

   //add padding to nt->data and update ml->data for mechs instances
   for (int m=0; m<neurox::mechanismsCount; m++)
   {
       short dataSize = mechanisms[m]->dataSize;
       Memb_list * instances = &b->mechsInstances[m];
       double* instance_data_new = &data_new[dataOffset];
       for (int i=0; i<instances->nodecount; i++)
           for (size_t d=0; d<sizeof_(dataSize); d++, dataOffset++)
               if (d < dataSize)
               {
                   data_new[dataOffset] = instances->data[dataSize*i+d];
                   dataNewOffset.push_back(dataOffset);
               }
       instances->data = instance_data_new;
   }
   dataNewOffset.resize(dataNewOffset.size());

   //convert AP-threshold pointer
   b->thvar_ptr = b->thvar_ptr ? &b->nt->_actual_v[thvar_idx] : nullptr;

   //convert VecPlay continuous pointers
   for (int v=0; v<b->nt->n_vecplay; v++)
   {
       VecPlayContinuousX * vc = (VecPlayContinuousX*) b->nt->_vecplay[v];
       int vc_offset = vc->pd_ - &b->nt->_data[0];
       assert(vc_offset>=0 && vc_offset <= b->nt->_ndata);
       int vc_offset_new = dataNewOffset.at(vc_offset);
       vc->pd_ = &data_new[ vc_offset_new ];
   }

   delete_(b->nt->_data);
   b->nt->_data  = data_new;

   //convert ml->pdata offsets for data vector
   for (int m=0; m<mechanismsCount; m++)
   {
       short pdataSize = mechanisms[m]->pdataSize;
       Memb_list * instances = &b->mechsInstances[m];
       for (int i=0; i<instances->nodecount; i++)
           for (int p=0; p<pdataSize; p++)
           {
               int ptype = memb_func[mechanisms[m]->type].dparam_semantics[p];
               //true data pointer to area or ion, false for vdata pointers or other types
               bool isPointer = ptype==-1 || (ptype>0 && ptype<1000);
               if (isPointer)
               {
                   int & pd = instances->pdata[pdataSize*i+p];
                   pd = dataNewOffset.at(pd); //point to new offset (with gaps)
               }
           }
    }

//swap rows and columns
   //convert AoS to AoS
#if LAYOUT ==1
   //TODO
#endif
}
