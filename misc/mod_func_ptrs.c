#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "coreneuron/coreneuron.h"

extern void 
  _nrn_init__Ca(NrnThread*, Memb_list*, int),
  _nrn_init__CaDynamics_E2(NrnThread*, Memb_list*, int),
  _nrn_init__Ca_HVA(NrnThread*, Memb_list*, int),
  _nrn_init__Ca_LVAst(NrnThread*, Memb_list*, int),
  _nrn_init__Ih(NrnThread*, Memb_list*, int),
  _nrn_init__Im(NrnThread*, Memb_list*, int),
  _nrn_init__K_Pst(NrnThread*, Memb_list*, int),
  _nrn_init__K_Tst(NrnThread*, Memb_list*, int),
  _nrn_init__KdShu2007(NrnThread*, Memb_list*, int),
  _nrn_init__NaTa_t(NrnThread*, Memb_list*, int),
  _nrn_init__NaTs2_t(NrnThread*, Memb_list*, int),
  _nrn_init__Nap_Et2(NrnThread*, Memb_list*, int),
  _nrn_init__ProbAMPANMDA_EMS(NrnThread*, Memb_list*, int),
  _nrn_init__ProbGABAAB_EMS(NrnThread*, Memb_list*, int),
  _nrn_init__SK_E2(NrnThread*, Memb_list*, int),
  _nrn_init__SKv3_1(NrnThread*, Memb_list*, int),
  _nrn_init__StochKv(NrnThread*, Memb_list*, int),
  _nrn_init__InhPoissonStim(NrnThread*, Memb_list*, int),
  _nrn_init__cc(NrnThread*, Memb_list*, int),
  _nrn_init__ExpSyn(NrnThread*, Memb_list*, int),
  _nrn_init__hh(NrnThread*, Memb_list*, int),
  _nrn_init__NetStim(NrnThread*, Memb_list*, int),
  _nrn_init__pas(NrnThread*, Memb_list*, int),
  _nrn_init__PatternStim(NrnThread*, Memb_list*, int),
  _nrn_init__IClamp(NrnThread*, Memb_list*, int);

mod_f_t get_init_function(const char * sym)
{
  if (strcmp(sym, "Ca") == 0)  return _nrn_init__Ca;
  if (strcmp(sym, "CaDynamics_E2") == 0)  return _nrn_init__CaDynamics_E2;
  if (strcmp(sym, "Ca_HVA") == 0)  return _nrn_init__Ca_HVA;
  if (strcmp(sym, "Ca_LVAst") == 0)  return _nrn_init__Ca_LVAst;
  if (strcmp(sym, "Ih") == 0)  return _nrn_init__Ih;
  if (strcmp(sym, "Im") == 0)  return _nrn_init__Im;
  if (strcmp(sym, "K_Pst") == 0)  return _nrn_init__K_Pst;
  if (strcmp(sym, "K_Tst") == 0)  return _nrn_init__K_Tst;
  if (strcmp(sym, "KdShu2007") == 0)  return _nrn_init__KdShu2007;
  if (strcmp(sym, "NaTa_t") == 0)  return _nrn_init__NaTa_t;
  if (strcmp(sym, "NaTs2_t") == 0)  return _nrn_init__NaTs2_t;
  if (strcmp(sym, "Nap_Et2") == 0)  return _nrn_init__Nap_Et2;
  if (strcmp(sym, "ProbAMPANMDA_EMS") == 0)  return _nrn_init__ProbAMPANMDA_EMS;
  if (strcmp(sym, "ProbGABAAB_EMS") == 0)  return _nrn_init__ProbGABAAB_EMS;
  if (strcmp(sym, "SK_E2") == 0)  return _nrn_init__SK_E2;
  if (strcmp(sym, "SKv3_1") == 0)  return _nrn_init__SKv3_1;
  if (strcmp(sym, "StochKv") == 0)  return _nrn_init__StochKv;
  if (strcmp(sym, "InhPoissonStim") == 0)  return _nrn_init__InhPoissonStim;
  if (strcmp(sym, "cc") == 0)  return _nrn_init__cc;
  if (strcmp(sym, "ExpSyn") == 0)  return _nrn_init__ExpSyn;
  if (strcmp(sym, "hh") == 0)  return _nrn_init__hh;
  if (strcmp(sym, "NetStim") == 0)  return _nrn_init__NetStim;
  if (strcmp(sym, "pas") == 0)  return _nrn_init__pas;
  if (strcmp(sym, "PatternStim") == 0)  return _nrn_init__PatternStim;
  if (strcmp(sym, "IClamp") == 0)  return _nrn_init__IClamp;
  return NULL;
}


extern void 
  _nrn_cur__Ca(NrnThread*, Memb_list*, int),
  _nrn_cur__CaDynamics_E2(NrnThread*, Memb_list*, int),
  _nrn_cur__Ca_HVA(NrnThread*, Memb_list*, int),
  _nrn_cur__Ca_LVAst(NrnThread*, Memb_list*, int),
  _nrn_cur__Ih(NrnThread*, Memb_list*, int),
  _nrn_cur__Im(NrnThread*, Memb_list*, int),
  _nrn_cur__K_Pst(NrnThread*, Memb_list*, int),
  _nrn_cur__K_Tst(NrnThread*, Memb_list*, int),
  _nrn_cur__KdShu2007(NrnThread*, Memb_list*, int),
  _nrn_cur__NaTa_t(NrnThread*, Memb_list*, int),
  _nrn_cur__NaTs2_t(NrnThread*, Memb_list*, int),
  _nrn_cur__Nap_Et2(NrnThread*, Memb_list*, int),
  _nrn_cur__ProbAMPANMDA_EMS(NrnThread*, Memb_list*, int),
  _nrn_cur__ProbGABAAB_EMS(NrnThread*, Memb_list*, int),
  _nrn_cur__SK_E2(NrnThread*, Memb_list*, int),
  _nrn_cur__SKv3_1(NrnThread*, Memb_list*, int),
  _nrn_cur__StochKv(NrnThread*, Memb_list*, int),
  _nrn_cur__cc(NrnThread*, Memb_list*, int),
  _nrn_cur__ExpSyn(NrnThread*, Memb_list*, int),
  _nrn_cur__hh(NrnThread*, Memb_list*, int),
  _nrn_cur__pas(NrnThread*, Memb_list*, int),
  _nrn_cur__IClamp(NrnThread*, Memb_list*, int);

mod_f_t get_cur_function(const char * sym)
{
  if (strcmp(sym, "Ca") == 0)  return _nrn_cur__Ca;
  if (strcmp(sym, "CaDynamics_E2") == 0)  return _nrn_cur__CaDynamics_E2;
  if (strcmp(sym, "Ca_HVA") == 0)  return _nrn_cur__Ca_HVA;
  if (strcmp(sym, "Ca_LVAst") == 0)  return _nrn_cur__Ca_LVAst;
  if (strcmp(sym, "Ih") == 0)  return _nrn_cur__Ih;
  if (strcmp(sym, "Im") == 0)  return _nrn_cur__Im;
  if (strcmp(sym, "K_Pst") == 0)  return _nrn_cur__K_Pst;
  if (strcmp(sym, "K_Tst") == 0)  return _nrn_cur__K_Tst;
  if (strcmp(sym, "KdShu2007") == 0)  return _nrn_cur__KdShu2007;
  if (strcmp(sym, "NaTa_t") == 0)  return _nrn_cur__NaTa_t;
  if (strcmp(sym, "NaTs2_t") == 0)  return _nrn_cur__NaTs2_t;
  if (strcmp(sym, "Nap_Et2") == 0)  return _nrn_cur__Nap_Et2;
  if (strcmp(sym, "ProbAMPANMDA_EMS") == 0)  return _nrn_cur__ProbAMPANMDA_EMS;
  if (strcmp(sym, "ProbGABAAB_EMS") == 0)  return _nrn_cur__ProbGABAAB_EMS;
  if (strcmp(sym, "SK_E2") == 0)  return _nrn_cur__SK_E2;
  if (strcmp(sym, "SKv3_1") == 0)  return _nrn_cur__SKv3_1;
  if (strcmp(sym, "StochKv") == 0)  return _nrn_cur__StochKv;
  if (strcmp(sym, "cc") == 0)  return _nrn_cur__cc;
  if (strcmp(sym, "ExpSyn") == 0)  return _nrn_cur__ExpSyn;
  if (strcmp(sym, "hh") == 0)  return _nrn_cur__hh;
  if (strcmp(sym, "pas") == 0)  return _nrn_cur__pas;
  if (strcmp(sym, "IClamp") == 0)  return _nrn_cur__IClamp;
  return NULL;
}


extern void 
  _nrn_state__Ca(NrnThread*, Memb_list*, int),
  _nrn_state__CaDynamics_E2(NrnThread*, Memb_list*, int),
  _nrn_state__Ca_HVA(NrnThread*, Memb_list*, int),
  _nrn_state__Ca_LVAst(NrnThread*, Memb_list*, int),
  _nrn_state__Ih(NrnThread*, Memb_list*, int),
  _nrn_state__Im(NrnThread*, Memb_list*, int),
  _nrn_state__K_Pst(NrnThread*, Memb_list*, int),
  _nrn_state__K_Tst(NrnThread*, Memb_list*, int),
  _nrn_state__KdShu2007(NrnThread*, Memb_list*, int),
  _nrn_state__NaTa_t(NrnThread*, Memb_list*, int),
  _nrn_state__NaTs2_t(NrnThread*, Memb_list*, int),
  _nrn_state__Nap_Et2(NrnThread*, Memb_list*, int),
  _nrn_state__ProbAMPANMDA_EMS(NrnThread*, Memb_list*, int),
  _nrn_state__ProbGABAAB_EMS(NrnThread*, Memb_list*, int),
  _nrn_state__SK_E2(NrnThread*, Memb_list*, int),
  _nrn_state__SKv3_1(NrnThread*, Memb_list*, int),
  _nrn_state__StochKv(NrnThread*, Memb_list*, int),
  _nrn_state__InhPoissonStim(NrnThread*, Memb_list*, int),
  _nrn_state__cc(NrnThread*, Memb_list*, int),
  _nrn_state__ExpSyn(NrnThread*, Memb_list*, int),
  _nrn_state__hh(NrnThread*, Memb_list*, int),
  _nrn_state__NetStim(NrnThread*, Memb_list*, int),
  _nrn_state__pas(NrnThread*, Memb_list*, int),
  _nrn_state__PatternStim(NrnThread*, Memb_list*, int),
  _nrn_state__IClamp(NrnThread*, Memb_list*, int);

mod_f_t get_state_function(const char * sym)
{
  if (strcmp(sym, "Ca") == 0)  return _nrn_state__Ca;
  if (strcmp(sym, "CaDynamics_E2") == 0)  return _nrn_state__CaDynamics_E2;
  if (strcmp(sym, "Ca_HVA") == 0)  return _nrn_state__Ca_HVA;
  if (strcmp(sym, "Ca_LVAst") == 0)  return _nrn_state__Ca_LVAst;
  if (strcmp(sym, "Ih") == 0)  return _nrn_state__Ih;
  if (strcmp(sym, "Im") == 0)  return _nrn_state__Im;
  if (strcmp(sym, "K_Pst") == 0)  return _nrn_state__K_Pst;
  if (strcmp(sym, "K_Tst") == 0)  return _nrn_state__K_Tst;
  if (strcmp(sym, "KdShu2007") == 0)  return _nrn_state__KdShu2007;
  if (strcmp(sym, "NaTa_t") == 0)  return _nrn_state__NaTa_t;
  if (strcmp(sym, "NaTs2_t") == 0)  return _nrn_state__NaTs2_t;
  if (strcmp(sym, "Nap_Et2") == 0)  return _nrn_state__Nap_Et2;
  if (strcmp(sym, "ProbAMPANMDA_EMS") == 0)  return _nrn_state__ProbAMPANMDA_EMS;
  if (strcmp(sym, "ProbGABAAB_EMS") == 0)  return _nrn_state__ProbGABAAB_EMS;
  if (strcmp(sym, "SK_E2") == 0)  return _nrn_state__SK_E2;
  if (strcmp(sym, "SKv3_1") == 0)  return _nrn_state__SKv3_1;
  if (strcmp(sym, "StochKv") == 0)  return _nrn_state__StochKv;
  if (strcmp(sym, "InhPoissonStim") == 0)  return _nrn_state__InhPoissonStim;
  if (strcmp(sym, "cc") == 0)  return _nrn_state__cc;
  if (strcmp(sym, "ExpSyn") == 0)  return _nrn_state__ExpSyn;
  if (strcmp(sym, "hh") == 0)  return _nrn_state__hh;
  if (strcmp(sym, "NetStim") == 0)  return _nrn_state__NetStim;
  if (strcmp(sym, "pas") == 0)  return _nrn_state__pas;
  if (strcmp(sym, "PatternStim") == 0)  return _nrn_state__PatternStim;
  if (strcmp(sym, "IClamp") == 0)  return _nrn_state__IClamp;
  return NULL;
}


extern void 
  _nrn_cur_lock__Ca(NrnThread*, Memb_list*, int),
  _nrn_cur_lock__CaDynamics_E2(NrnThread*, Memb_list*, int),
  _nrn_cur_lock__Ca_HVA(NrnThread*, Memb_list*, int),
  _nrn_cur_lock__Ca_LVAst(NrnThread*, Memb_list*, int),
  _nrn_cur_lock__Ih(NrnThread*, Memb_list*, int),
  _nrn_cur_lock__Im(NrnThread*, Memb_list*, int),
  _nrn_cur_lock__K_Pst(NrnThread*, Memb_list*, int),
  _nrn_cur_lock__K_Tst(NrnThread*, Memb_list*, int),
/*
  _nrn_cur_lock__KdShu2007(NrnThread*, Memb_list*, int),
  */
  _nrn_cur_lock__NaTa_t(NrnThread*, Memb_list*, int),
  _nrn_cur_lock__NaTs2_t(NrnThread*, Memb_list*, int),
  _nrn_cur_lock__Nap_Et2(NrnThread*, Memb_list*, int),
  _nrn_cur_lock__ProbAMPANMDA_EMS(NrnThread*, Memb_list*, int),
  _nrn_cur_lock__ProbGABAAB_EMS(NrnThread*, Memb_list*, int),
  _nrn_cur_lock__SK_E2(NrnThread*, Memb_list*, int),
  _nrn_cur_lock__SKv3_1(NrnThread*, Memb_list*, int),
  /*
  _nrn_cur_lock__StochKv(NrnThread*, Memb_list*, int),
  _nrn_cur_lock__InhPoissonStim(NrnThread*, Memb_list*, int),
  _nrn_cur_lock__cc(NrnThread*, Memb_list*, int),
  _nrn_cur_lock__ExpSyn(NrnThread*, Memb_list*, int),
  _nrn_cur_lock__hh(NrnThread*, Memb_list*, int),
  _nrn_cur_lock__NetStim(NrnThread*, Memb_list*, int),
  */
  _nrn_cur_lock__pas(NrnThread*, Memb_list*, int),
/*
  _nrn_cur_lock__PatternStim(NrnThread*, Memb_list*, int),
  */
  _nrn_cur_lock__IClamp(NrnThread*, Memb_list*, int)
;

mod_lock_f_t get_cur_lock_function(const char * sym)
{
  if (strcmp(sym, "Ca") == 0)  return _nrn_cur_lock__Ca;
  if (strcmp(sym, "CaDynamics_E2") == 0)  return _nrn_cur_lock__CaDynamics_E2;
  if (strcmp(sym, "Ca_HVA") == 0)  return _nrn_cur_lock__Ca_HVA;
  if (strcmp(sym, "Ca_LVAst") == 0)  return _nrn_cur_lock__Ca_LVAst;
  if (strcmp(sym, "Ih") == 0)  return _nrn_cur_lock__Ih;
  if (strcmp(sym, "Im") == 0)  return _nrn_cur_lock__Im;
  if (strcmp(sym, "K_Pst") == 0)  return _nrn_cur_lock__K_Pst;
  if (strcmp(sym, "K_Tst") == 0)  return _nrn_cur_lock__K_Tst;
  /*
  if (strcmp(sym, "KdShu2007") == 0)  return _nrn_cur_lock__KdShu2007;
  */
  if (strcmp(sym, "NaTa_t") == 0)  return _nrn_cur_lock__NaTa_t;
  if (strcmp(sym, "NaTs2_t") == 0)  return _nrn_cur_lock__NaTs2_t;
  if (strcmp(sym, "Nap_Et2") == 0)  return _nrn_cur_lock__Nap_Et2;
  if (strcmp(sym, "ProbAMPANMDA_EMS") == 0)  return _nrn_cur_lock__ProbAMPANMDA_EMS;
  if (strcmp(sym, "ProbGABAAB_EMS") == 0)  return _nrn_cur_lock__ProbGABAAB_EMS;
  if (strcmp(sym, "SK_E2") == 0)  return _nrn_cur_lock__SK_E2;
  if (strcmp(sym, "SKv3_1") == 0)  return _nrn_cur_lock__SKv3_1;
  /*
  if (strcmp(sym, "StochKv") == 0)  return _nrn_cur_lock__StochKv;
  if (strcmp(sym, "InhPoissonStim") == 0)  return _nrn_cur_lock__InhPoissonStim;
  if (strcmp(sym, "cc") == 0)  return _nrn_cur_lock__cc;
  if (strcmp(sym, "ExpSyn") == 0)  return _nrn_cur_lock__ExpSyn;
  if (strcmp(sym, "hh") == 0)  return _nrn_cur_lock__hh;
  if (strcmp(sym, "NetStim") == 0)  return _nrn_cur_lock__NetStim;
  */
  if (strcmp(sym, "pas") == 0)  return _nrn_cur_lock__pas;
  /*
  if (strcmp(sym, "PatternStim") == 0)  return _nrn_cur_lock__PatternStim;
  */
  if (strcmp(sym, "IClamp") == 0)  return _nrn_cur_lock__IClamp;
  return NULL;
}


extern void 
  _net_receive2__ProbAMPANMDA_EMS(NrnThread*, Memb_list*, int, int, double),
  _net_receive2__ProbGABAAB_EMS(NrnThread*, Memb_list*, int, int, double),
  _net_receive2__NetStim(NrnThread*, Memb_list*, int, int, double),
  _net_receive2__ExpSyn(NrnThread*, Memb_list*, int, int, double),
  _net_receive2__PatternStim(NrnThread*, Memb_list*, int, int, double),
  _net_receive2__InhPoissonStim(NrnThread*, Memb_list*, int, int, double);

pnt_receive2_t get_net_receive_function(const char * sym)
{
  if (strcmp(sym, "ProbAMPANMDA_EMS") == 0)  return _net_receive2__ProbAMPANMDA_EMS;
  if (strcmp(sym, "ProbGABAAB_EMS") == 0)  return _net_receive2__ProbGABAAB_EMS;
  if (strcmp(sym, "NetStim") == 0)  return _net_receive2__NetStim;
  if (strcmp(sym, "ExpSyn") == 0)  return _net_receive2__ExpSyn;
  if (strcmp(sym, "PatternStim") == 0)  return _net_receive2__PatternStim;
  if (strcmp(sym, "InhPoissonStim") == 0)  return _net_receive2__InhPoissonStim;
  return NULL;
}

mod_f_t get_BA_function(const char * sym, int BA_func_id)
{
  return NULL;
}
