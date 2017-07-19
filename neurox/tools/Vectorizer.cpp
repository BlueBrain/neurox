#include "neurox/neurox.h"

#include <iostream>

using namespace std;
using namespace neurox;

void tools::Vectorizer::ConvertToSOA(Branch * b)
{
   //NOTE: arrays with memory-aligned allocation:
   //ml->pdata, nt->data, ml->nodeindices, v->parent_index, and ml->shadow_*;
   //data needs to add gaps, pdata needs new offset values

   assert(memb_func);
   assert(LAYOUT==0);

   //get total counts
   int N = b->nt->end;
   size_t oldDataSize  = 6*N;
   size_t newDataSize  = 6*SizeOf(N);
   assert(newDataSize % NEUROX_SOA_PADDING==0);

   fflush(stdout);
   for (int i=0; i<b->nt->_ndata; i++)
       std::cout << "## dataOld[" << i << "]=" << b->nt->_data[i] << std::endl;

   for (int m=0; m<mechanismsCount; m++)
   {
       b->mechsInstances[m]._nodecount_padded = SizeOf(b->mechsInstances[m].nodecount);
       newDataSize  += b->mechsInstances[m]._nodecount_padded * mechanisms[m]->dataSize;
       oldDataSize  += b->mechsInstances[m].nodecount * mechanisms[m]->dataSize;
   }
   std::vector<int> dataOffsets(oldDataSize);

   //old thvar offset
   int thvar_idx = b->thvar_ptr ? b->thvar_ptr - &b->nt->_actual_v[0] : -1; //compartment id (typically 2)
   assert(thvar_idx==-1 || (thvar_idx>=0 && thvar_idx < N));

   double* dataNew = New<double>(newDataSize);

   //add padding to data for RHS, D, A, B, V and area
   size_t newOffset =0;
   for (int i=0; i<6; i++)
        for (size_t j=0; j<SizeOf(N); j++, newOffset++)
            if (j < N)
            {
                int oldOffset = N*i + j;
                dataNew[newOffset] = b->nt->_data[oldOffset];
                dataOffsets[oldOffset] = newOffset;
            }

   assert(newOffset == SizeOf(N)*6);
   b->nt->_actual_rhs  = &dataNew[SizeOf(N)*0];
   b->nt->_actual_d    = &dataNew[SizeOf(N)*1];
   b->nt->_actual_a    = &dataNew[SizeOf(N)*2];
   b->nt->_actual_b    = &dataNew[SizeOf(N)*3];
   b->nt->_actual_v    = &dataNew[SizeOf(N)*4];
   b->nt->_actual_area = &dataNew[SizeOf(N)*5];

   //convert AP-threshold pointer
   b->thvar_ptr = b->thvar_ptr ? &b->nt->_actual_v[thvar_idx] : nullptr;

   //add padding and converting AoS->SoA in nt->data and update ml->data for mechs instances
   unsigned oldOffsetAcc = N*6;
   for (int m=0; m<neurox::mechanismsCount; m++)
   {
       Memb_list * instances = &b->mechsInstances[m];
       double* instanceDataNew = &dataNew[newOffset];
       for (size_t i=0; i<mechanisms[m]->dataSize; i++) //for every variable
           for (int n=0; n<instances->_nodecount_padded; n++, newOffset++) //for every node
               if (n < instances->nodecount)
               {
                   int oldOffset = mechanisms[m]->dataSize*n + i; //d-th variable in i-th instance
                   assert(b->nt->_data[oldOffsetAcc+oldOffset] == instances->data[oldOffset]);
                   dataNew[newOffset] = instances->data[oldOffset];
                   dataOffsets[oldOffsetAcc+oldOffset]=newOffset;
               }

       oldOffsetAcc += mechanisms[m]->dataSize*instances->nodecount;

       //all instances processed, replace pointer by new padded data
       instances->data = instanceDataNew;
   }

   for (int i=0; i<newDataSize; i++)
       std::cout << "## dataNew[" << i << "]=" << dataNew[i] << std::endl;
   for (int i=0; i<dataOffsets.size(); i++)
       std::cout << "## dataOffsets[" << i << "]=" << dataOffsets[i] << std::endl;

   assert(newOffset == newDataSize);
   dataOffsets.shrink_to_fit();
   b->nt->_ndata = newDataSize;

   //convert VecPlay continuous pointers
   for (int v=0; v<b->nt->n_vecplay; v++)
   {
       VecPlayContinuousX * vc = (VecPlayContinuousX*) b->nt->_vecplay[v];
       int pd_offset = vc->pd_ - &b->nt->_data[0];
       assert(pd_offset>=0 && pd_offset <= b->nt->_ndata);
       int pd_offset_new = dataOffsets.at(pd_offset);
       vc->pd_ = &dataNew[ pd_offset_new ];
       printf("### post-vectorize: pd %f (offset %lld -> %lld)\n", *vc->pd_, pd_offset, pd_offset_new);
   }

   Delete(b->nt->_data);
   b->nt->_data  = dataNew;

   //add padding and update ml->pdata for mechs instances
   oldOffsetAcc = N*6;
   for (int m=0; m<mechanismsCount; m++)
   {
       Memb_list * instances = &b->mechsInstances[m];
       int totalPDataSize = instances->_nodecount_padded * mechanisms[m]->pdataSize;
       int* pdataNew = New<int>(totalPDataSize);
       size_t pdataNewOffset =0;
       for (size_t i=0; i<mechanisms[m]->pdataSize; i++) //for every pointer
           for (int n=0; n<instances->_nodecount_padded; n++, pdataNewOffset++) //for every node
           {
               if (n < instances->nodecount)
               {
                   //padding:
                   int oldOffset = mechanisms[m]->pdataSize*n + i; //p-th pointer in i-th instance
                   pdataNew[pdataNewOffset] = instances->pdata[oldOffset];

                   //new SoA pointer values:
                   int ptype = memb_func[mechanisms[m]->type].dparam_semantics[i];
                   bool isPointer = ptype==-1 || (ptype>0 && ptype<1000);
                   if (isPointer) //true for pointer to area in nt->data, or ion instance data
                     pdataNew[pdataNewOffset] = dataOffsets.at(oldOffsetAcc + oldOffset); //point to new offset (with gaps)
               }
               else
               {
                   pdataNew[pdataNewOffset]=-1;
               }
           }
        oldOffsetAcc += mechanisms[m]->dataSize*instances->nodecount;
        Delete(instances->pdata);
        instances->pdata = pdataNew;
    }
}
