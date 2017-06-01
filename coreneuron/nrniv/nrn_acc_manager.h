#ifndef _nrn_device_manager_
#define _nrn_device_manager_

#if defined(_OPENACC)
#include <openacc.h>
#endif

#include "coreneuron/nrnoc/multicore.h"

void setup_nrnthreads_on_device(NrnThread* threads, int nthreads);
void update_nrnthreads_on_host(NrnThread* threads, int nthreads);
void update_nrnthreads_on_device(NrnThread* threads, int nthreads);
void modify_data_on_device(NrnThread* threads, int nthreads);
void dump_nt_to_file(char* filename, NrnThread* threads, int nthreads);
void finalize_data_on_device();

#ifdef __cplusplus
extern "C" {
#endif

void update_matrix_from_gpu(NrnThread* _nt);
void update_matrix_to_gpu(NrnThread* _nt);
void update_net_receive_buffer(NrnThread* _nt);
void realloc_net_receive_buffer(NrnThread* nt, Memb_list* ml);
void update_net_send_buffer_on_host(NrnThread* nt, NetSendBuffer_t* nsb);

#ifdef __cplusplus
}
#endif

#endif  // _nrn_device_manager_
