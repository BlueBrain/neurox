#include "neurox/neurox.h"

using namespace std;
using namespace neurox;

size_t Tools::Vectorizer::sizeof_(size_t size)
{
    return size % NEUROX_MEM_ALIGNMENT + 1;
}

void Tools::Vectorizer::vectorize(Branch * b)
{
   //Vectors to expand:
   //pdata, data, nodeindices, parent index
   //ml->shadow_*

   //Vectors of offsets to be updateed
   //nt->pdata(ml->pdata)

   //get total counts
   size_t totalDataSize = 6*sizeof_(b->nt->end);
   int *  dataNewOffset = new int[b->nt->_ndata];
   for (int m=0; m<mechanismsCount; m++)
       totalDataSize += b->mechsInstances[m].nodecount * sizeof_(mechanisms[m]->dataSize);

   double* data_new = new_<double>(totalDataSize);

   //copy RHS, D, A, B, V and area to new data
   size_t dataOffset =0;
   for (int i=0; i<6; i++) //RHS, D, A, B, V and area
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

        for (int j=0; j<sizeof_(b->nt->end); j++, dataOffset++)
            if (j<b->nt->end)
            {
                data_new[dataOffset] = b->nt->_data[i*6+j];
                dataNewOffset[i*6+j] = dataOffset;
            }
   }

   //create data and pdata for mechs instances
   for (int m=0; m<mechanismsCount; m++)
   {
       Memb_list * instances = &b->mechsInstances[m];
       short dataSize = mechanisms[m]->dataSize;
       short pdataSize = mechanisms[m]->pdataSize;

       for (int i=0; i<instances->nodecount; i++)
       {
           Memb_list & instance = instances[m];

           //instances data are part of nt->data
           double * instance_data_new = &data_new[dataOffset];
           for (int d=0; d<sizeof_(dataSize); d++, dataOffset++)
               if (d<sizeof_(dataSize))
                   data_new[dataOffset] = instance.data[d];
           instance.data = instance_data_new;

           //instances pdata are single arrays
           size_t instance_pdata_size = sizeof_(pdataSize);
           int * instance_pdata_new = new_<int>(instance_pdata_size);
           size_t pdataOffset=0;
           for (int p=0; p<sizeof_(pdataSize); p++, pdataOffset++)
               if (p<sizeof_(pdataSize))
                   instance_pdata_new[pdataOffset] = instance.pdata[p];
           delete_(instance.pdata);
           instance.pdata = instance_pdata_new;
       }
   }
   delete_(b->nt->_data);
   b->nt->_data  = data_new;

   //TODO swap rows and columns!
   //todo update pdata values
}
