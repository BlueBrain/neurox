#include "neurox/neurox.h"

#include <iostream>

using namespace std;
using namespace neurox;

void tools::Vectorizer::ConvertToSOA(Branch* b) {
  // NOTE: arrays with memory-aligned allocation:
  // ml->pdata, nt->data, ml->nodeindices, v->parent_index, and ml->shadow_*;
  // data needs to add gaps, pdata needs new offset values

  assert(memb_func);
  assert(LAYOUT == 0);

  // get total counts
  int N = b->nt->end;
  size_t oldDataSize = 6 * N;
  size_t newDataSize = 6 * SizeOf(N);
  assert(newDataSize % NEUROX_SOA_PADDING == 0);

  for (int m = 0; m < mechanisms_count; m++) {
    b->mechsInstances[m]._nodecount_padded =
        SizeOf(b->mechsInstances[m].nodecount);
    newDataSize +=
        b->mechsInstances[m]._nodecount_padded * mechanisms[m]->dataSize;
    oldDataSize += b->mechsInstances[m].nodecount * mechanisms[m]->dataSize;
  }

  // old-to-new offset map
  std::vector<int> dataOffsets(oldDataSize, -99999);

  // old thvar offset
  int thvar_idx = b->thvar_ptr ? b->thvar_ptr - &b->nt->_actual_v[0]
                               : -1;  // compartment id (typically 2)
  assert(thvar_idx == -1 || (thvar_idx >= 0 && thvar_idx < N));

  double* dataNew = New<double>(newDataSize);
  double* dataOld = b->nt->_data;

  // add padding to data for RHS, D, A, B, V and area
  for (int i = 0; i < 6; i++)
    for (size_t j = 0; j < N; j++) {
      int oldOffset = N * i + j;
      int newOffset = SizeOf(N) * i + j;
      dataNew[newOffset] = b->nt->_data[oldOffset];
      dataOffsets.at(oldOffset) = newOffset;
    }

  b->nt->_actual_rhs = &dataNew[SizeOf(N) * 0];
  b->nt->_actual_d = &dataNew[SizeOf(N) * 1];
  b->nt->_actual_a = &dataNew[SizeOf(N) * 2];
  b->nt->_actual_b = &dataNew[SizeOf(N) * 3];
  b->nt->_actual_v = &dataNew[SizeOf(N) * 4];
  b->nt->_actual_area = &dataNew[SizeOf(N) * 5];

  // convert AP-threshold pointer
  b->thvar_ptr = b->thvar_ptr ? &b->nt->_actual_v[thvar_idx] : nullptr;

  // add padding and converting AoS->SoA in nt->data and update ml->data for
  // mechs instances
  unsigned oldOffsetAcc = N * 6;
  unsigned newOffsetAcc = SizeOf(N) * 6;

  for (int m = 0; m < neurox::mechanisms_count; m++) {
    Memb_list* instances = &b->mechsInstances[m];
    double* instanceDataNew = &dataNew[newOffsetAcc];

    int totalPDataSize =
        instances->_nodecount_padded * mechanisms[m]->pdataSize;
    int* pdataNew = New<int>(totalPDataSize);
    int* pdataOld = instances->pdata;

    for (int n = 0; n < instances->nodecount; n++)  // for every node
    {
      for (size_t i = 0; i < mechanisms[m]->dataSize;
           i++)  // for every variable
      {
        int oldOffset = mechanisms[m]->dataSize * n + i;
        int newOffset = SizeOf(instances->nodecount) * i + n;
        dataNew[newOffsetAcc + newOffset] = instances->data[oldOffset];
        assert(b->nt->_data[oldOffsetAcc + oldOffset] ==
               instances->data[oldOffset]);
        dataOffsets.at(oldOffsetAcc + oldOffset) = newOffsetAcc + newOffset;
      }

      for (size_t i = 0; i < mechanisms[m]->pdataSize;
           i++)  // for every pointer
      {
        // padding
        int oldOffset = mechanisms[m]->pdataSize * n + i;      // SoA
        int newOffset = instances->_nodecount_padded * i + n;  // AoS
        pdataNew[newOffset] = pdataOld[oldOffset];

        // get correct pdata offset: without branching, offsets are already
        // correct for both LAYOUTs and padding
        if (input_params->branchingDepth > 0) {
          int ptype = memb_func[mechanisms[m]->type].dparam_semantics[i];
          bool isPointer = ptype == -1 || (ptype > 0 && ptype < 1000);
          if (isPointer)  // true for pointer to area in nt->data, or ion
                          // instance data
          {
            pdataNew[newOffset] =
                dataOffsets.at(pdataOld[oldOffset]);  // point to new offset
            assert(dataNew[pdataNew[newOffset]] ==
                   dataOld[pdataOld[oldOffset]]);
          }
        }
        assert(pdataNew[newOffset] != -99999);
      }
    }

    oldOffsetAcc += mechanisms[m]->dataSize * instances->nodecount;
    newOffsetAcc += mechanisms[m]->dataSize * SizeOf(instances->nodecount);

    // all instances processed, replace pointer by new padded data
    Delete(instances->pdata);
    instances->pdata = pdataNew;
    instances->data = instanceDataNew;
  }

  // convert VecPlay continuous pointers
  for (int v = 0; v < b->nt->n_vecplay; v++) {
    VecPlayContinuousX* vc = (VecPlayContinuousX*)b->nt->_vecplay[v];
    int pd_offset = vc->pd_ - &b->nt->_data[0];
    assert(pd_offset >= 0 && pd_offset <= b->nt->_ndata);
    int pd_offset_new = dataOffsets.at(pd_offset);
    vc->pd_ = &dataNew[pd_offset_new];
  }

  Delete(b->nt->_data);
  b->nt->_data = dataNew;
  b->nt->_ndata = newDataSize;
}
