/*       ----------------------------------------------------

	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.

	 Name of Program: Local Control Model
	 Version: Documented Version, C
	 Date: November 2003

	 --------------------------------------------------

	 mpi_master - routines that control computations
	              in all processes

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <cassert>
#include <omp.h>
#include "StochModel.h"

#if USE_MPI
#include <mpi.h>
#endif


void StochModel::send_calc_fru_avg()
{
  double CaSS, CaJSR, CaTOT_SS, CaTOT_JSR;
  double JRyR, NRyR_Open, NNorm_Mode, NVinact, NLType_Open, NIto2_Open, NRyR_ready;

  CaSS = 0.0;
  CaTOT_SS = 0.0;
  CaJSR = 0.0;
  CaTOT_JSR = 0.0;
  JRyR = 0.0;
  NRyR_Open = 0.0;
  NNorm_Mode = 0.0;
  NRyR_ready = 0.0;
  NVinact = 0.0;
  NLType_Open = 0.0;
  NIto2_Open = 0.0;

  double CaSLtot = 0;
  double Ca_avg = 0;

  for (int i = 0; i < NFRU; i++) {
    CaSS += FRUs[i].Get_CaSS_Avg();
    CaTOT_SS += FRUs[i].Get_CaTOT_SS();
    CaJSR += FRUs[i].FRU_states[index_frustates_CaJSR];
    CaTOT_JSR += FRUs[i].Get_CaTOT_JSR();
    JRyR += FRUs[i].Get_JRyR();
    NRyR_Open += FRUs[i].Get_RyR_Open();
    NNorm_Mode += FRUs[i].Get_LCC_NormMode();
    NRyR_ready += FRUs[i].Get_RyR_Ready();
    NVinact += FRUs[i].Get_LCC_Vinact();
    NLType_Open += FRUs[i].Get_LCC_Open();
    NIto2_Open += FRUs[i].Get_Ito2_Open();
    //printf("CaSS %d %g\n",i,FRUs[i].Get_CaSS_Avg()*1e3);
  }
  for (int i = 0; i < NFRU; i++) {
    CaSLtot += state[N_states+i] * (1.0 + CMDNtot / (KmCMDN + state[N_states+i]) + EGTAtot / (KmEGTA + state[N_states+i]) + BSLtot_SL / (KBSL + state[N_states+i]));
    Ca_avg += state[N_states+i];
    //printf("CaSL %d %g\n",i,state[N_states+i]*1e3);
  }

  double f_lcc_act = (beta_flag) ? beta_frac_active : lcc_frac_active;

  otherstates[index_CaSSavg]  = CaSS / (double)(NFRU);
  otherstates[index_CaJSRavg] = CaJSR / (double)(NFRU);
  otherstates[index_JRyRtot] = JRyR * NFRU_scale * VSS / Vmyo;
  otherstates[index_PRyR_Open]  = NRyR_Open / (double)(NRyRs_per_cleft * NFRU);
  otherstates[index_PRyR_ready]  = NRyR_ready / (double)(NRyRs_per_cleft * NFRU);
  otherstates[index_CaTOT2] = 1.e6 * ((CaTOT_SS * VSS + CaTOT_JSR * VJSR ) * NFRU_scale
				      + CaSLtot * NFRU_scale * VSL
                                      + Vmyo * (state[index_LTRPNCa] * LTRPNtot + state[index_HTRPNCa] * HTRPNtot
                                          + state[index_Cai] * (1.0 + CMDNtot / (KmCMDN + state[index_Cai]) + EGTAtot / (KmEGTA + state[index_Cai])))
                                      + VNSR * state[index_CaNSR]);

  otherstates[index_PNorm_Mode] = (NNorm_Mode / (double)(NFRU * Max_LCCs_per_cleft * f_lcc_act)); // / (1-f_lcc_inact)) - (f_lcc_inact/(1-f_lcc_inact));
  otherstates[index_PnotVinact] = (NVinact / (double)(NFRU * Max_LCCs_per_cleft * f_lcc_act)); // / (1-f_lcc_inact)) - (f_lcc_inact/(1-f_lcc_inact)); 
  otherstates[index_PLType_Open] = (NLType_Open / (double)(NFRU * Max_LCCs_per_cleft * f_lcc_act)); // / (1-f_lcc_inact)); // - (f_lcc_inact/(1-f_lcc_inact));
  otherstates[index_PIto2_Open] = NIto2_Open / (double)(NFRU);
  otherstates[index_CaJSRtot] = 1.e6 * CaTOT_JSR * VJSR * NFRU_scale;
  otherstates[index_CaSStot] = 1.e6 * CaTOT_SS * VSS * NFRU_scale;
  otherstates[index_CaSRtot] = (state[index_CaNSR] * VNSR + CaTOT_JSR * VJSR * NFRU_scale)  / Vmyo;
  otherstates[index_CaSLtot] = 1e6 * CaSLtot * NFRU_scale * VSL;
  state[index_CaSL] = Ca_avg/NFRU;
}


/*
   void distrib_simFRU

   Have slave processes perform computations ie distributed simFRU

*/
void StochModel::distrib_simFRU(double st_time, double end_time,
                                double FRUdep_states[Nstates_FRUdep],
				double FRUdep_statesf[Nstates_FRUdep],
				std::vector<double> &states0,
				std::vector<double> &statesf,
                                double *Jxfer,
                                double *Jtr,
                                double *ICa,
                                double *Ito2)
{
  double num_stat[5], total_num[5];
  double V, sum_Jtr, sum_Jxfer, sum_CaSS, sum_CaJSR;
  double ICa_numerator, sum_ICa, sum_Ito2, NIto2_Open;
  double VF_over_RT, VFsq_over_RT, exp_2VFRT, exp_VFRT;

  V = FRUdep_states[index_frudep_V];
  VF_over_RT = V / RT_over_F;
  VFsq_over_RT = (1000.0 * Faraday) * VF_over_RT;
  exp_VFRT = exp(VF_over_RT);
  exp_2VFRT = exp_VFRT * exp_VFRT;
  //double Cai = FRUdep_states[index_frudep_Cai];
  double CaNSR = FRUdep_states[index_frudep_CaNSR];

  for(int i = 0; i < 5; i++) {
    num_stat[i] = 0;
  }

  if (end_time > st_time) {
//#pragma omp parallel for //USE OMP ONLY FOR STOCH MODEL TESTING
    for (int i = 0; i < NFRU; i++) {

      double FRUdep_i_0[Nstates_FRUdep];
      double FRUdep_i_f[Nstates_FRUdep];
      memcpy(FRUdep_i_0, FRUdep_states, sizeof(double)*Nstates_FRUdep);
      memcpy(FRUdep_i_f, FRUdep_statesf, sizeof(double)*Nstates_FRUdep); 
      
      FRUdep_i_0[index_frudep_CaSL] = states0[N_states+i];
      FRUdep_i_f[index_frudep_CaSL] = statesf[N_states+i];

      int mythread = omp_get_thread_num();
      
      FRUs[i].simfru(st_time,
                     end_time,
                     FRUdep_i_0,
		     FRUdep_i_f,
		     mt[mythread], 
		     mti[mythread],
                     (int)beta_flag,
                     Cao,
		     ryr_kplus_scale);
    }
  }

  num_stat[1] = 0;
  for (int i = 0; i < NFRU; i++) {
    FRUs[i].calc_fru_flux(num_stat);
  }
  sum_Jtr = (CaNSR*NFRU - num_stat[1])/tautr;

  sum_CaSS = num_stat[0];
  //sum_CaJSR = num_stat[1];
  total_num[2] = num_stat[2];
  total_num[3] = num_stat[3];
  NIto2_Open = num_stat[4];

  // Compute actual currents

  if (fabs(V) < 1.e-6) { // First order Taylor expansion
    ICa_numerator = total_num[2] - total_num[3] * Cao * 0.341;
    sum_ICa = PCa * 2.0 * 1000.0 * Faraday * ICa_numerator;
    sum_ICa = sum_ICa / Acap;		// divide by uF(Acap) to get current normalized to surface area
    sum_Ito2 = ((double)NIto2_Open) * PCl * 1000.0 * Faraday * (Clo - Cli);
    sum_Ito2 = sum_Ito2 / Acap;	// divide by uF(Acap) to get current normalized to surface area
  } else {
    ICa_numerator = total_num[2] * exp_2VFRT - total_num[3] * Cao * 0.341;
    sum_ICa = PCa * 4.0 * VFsq_over_RT * ICa_numerator / (exp_2VFRT - 1.0);
    sum_ICa = sum_ICa / Acap;		// divide by uF(Acap) to get current normalized to surface area
    sum_Ito2 = ((double)NIto2_Open) * PCl * VFsq_over_RT * (Cli - Clo * exp_VFRT) / (1.0 - exp_VFRT);
    sum_Ito2 = sum_Ito2 / Acap;	// divide by uF(Acap) to get current normalized to surface area
  }

  //sum_Jtr = (((double)NFRU) * CaNSR - sum_CaJSR) / tautr;
  sum_Jxfer = 0; //(sum_CaSS - ((double)NFRU * (double)Nclefts_FRU) * Cai) / tauxfer;

  *Jxfer = NFRU_scale * sum_Jxfer;
  *Jtr = NFRU_scale * sum_Jtr;
  *ICa = NFRU_scale * sum_ICa;
  *Ito2 = NFRU_scale * sum_Ito2;
}

/*
   void send_save_state(void)

   save current state

*/
void StochModel::send_save_state(void)
{
  FRUs_hold = FRUs;
  for (int i = 0; i < OMP_THREAD_MAX; i++) {
    memcpy(mt_hold[i], mt[i], sizeof(unsigned long)*(mtN+1));
  }
  memcpy(mti_hold, mti, sizeof(int)*OMP_THREAD_MAX);
}

/*
   void send_resume_state(void)

   resume saved state

*/
void StochModel::send_resume_state(void)
{
  FRUs = FRUs_hold;
  for (int i = 0; i < OMP_THREAD_MAX; i++) {
    memcpy(mt[i], mt_hold[i], sizeof(unsigned long)*(mtN+1));
  }
  memcpy(mti, mti_hold, sizeof(int)*OMP_THREAD_MAX);
}
