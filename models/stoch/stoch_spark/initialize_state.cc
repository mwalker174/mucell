/*       ----------------------------------------------------

				 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
				 Copyright 2003, The Johns Hopkins University
				 School of Medicine. All rights reserved.

				 Name of Program: Local Control Model
				 Version: Documented Version, C
				 Date: November 2003

				 --------------------------------------------------

				 initialize_state.c - Initialize states globally and for FRUs

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <limits.h>
#include <omp.h>

#include "StochModel.h"

void StochModel::Initialize_Default_State()
{
  int iFRU;

  state = std::vector<double>(N_states+NFRU);

  state[index_V] = -93.4;
  state[index_mNa] = 0.2654309E-03;
  state[index_hNa] = 0.9985519E+00;
  state[index_jNa] = 0.9987711E+00;
  state[index_Nai] = Default_Nai;
  state[index_Ki] = Default_Ki;
  state[index_xKs] = 0.1501720E-03; 
  state[index_C0Kv43] = 0.590958E+00;
  state[index_C1Kv43] = 0.155217E+00;
  state[index_C2Kv43] = 0.152881E-01;
  state[index_C3Kv43] = 0.669242E-03;
  state[index_OKv43] = 0.109861E-04;
  state[index_CI0Kv43] = 0.220999E-00;
  state[index_CI1Kv43] = 0.143513E-01;
  state[index_CI2Kv43] = 0.596808E-03;
  state[index_CI3Kv43] = 0.642013E-04;
  state[index_OIKv43] = 0.184528E-02;
  state[index_C0Kv14] = 0.7646312E+00;
  state[index_C1Kv14] = 0.7393727E-01;
  state[index_C2Kv14] = 0.2681076E-02;
  state[index_C3Kv14] = 0.4321218E-04;
  state[index_OKv14] = 0.2620394E-06;
  state[index_CI0Kv14] = 0.1493312E+00;
  state[index_CI1Kv14] = 0.6441895E-02;
  state[index_CI2Kv14] = 0.2012634E-02;
  state[index_CI3Kv14] = 0.6128624E-03;
  state[index_OIKv14] = 0.3084292E-03;
  state[index_CaTOT] = 0;
  state[index_C1Herg] = 0.99e+00;
  state[index_C2Herg] = 0.8e-02;
  state[index_C3Herg] = 0.2e-02;
  state[index_OHerg] = 0.0e+00;
  state[index_IHerg] = 0.0e+00;

  state[index_Cai] = Default_Cai;
  state[index_CaSL] = Default_Cai;
  state[index_CaNSR] = Default_CaSR;
  state[index_LTRPNCa] = Default_Cai / (Default_Cai + (kltrpn_minus/kltrpn_plus));
  state[index_HTRPNCa] = Default_Cai / (Default_Cai + (khtrpn_minus/khtrpn_plus)); 
  state[index_xKs2] = 0;

  for (int i = 0; i < NFRU; i++) {
    state[N_states+i] = Default_Cai;
  }

  ReleaseUnit r;
  r.RyR_state = 0;
  r.FRU_states[1] = Default_Cai;
  r.Ito2_state = 0;
  r.FRU_states[index_frustates_CaJSR] = Default_CaSR;   //  local JSR

  FRUs = std::vector<ReleaseUnit>(NFRU, r);
  int mythread = omp_get_thread_num();

  //Orphaned release sites
  for (int i = 0; i < NFRU; i++) {
    double u1 = Get_Rand_Unif(mt[mythread],mti[mythread]);
    if (u1 < frac_orphan) {
      FRUs[i].Set_Orphan(1);
      Num_Orphans++;
    }
  }
  
  FRU_neighbs = std::vector<std::vector<int> >(NFRU);
  FRU_neighbs_distance = std::vector<std::vector<double> >(NFRU);
  int cell_idx = 0;

  for (int i = 0; i < NFRU_X; i++) {
    for (int j = 0; j < NFRU_Y; j++) {
      for (int k = 0; k < NFRU_Z; k++) {
	std::vector<int> neighbors_init;

	double vol = 1; //1.0 + generateGaussian(mt[mythread],mti[mythread],0, 0.25);
	vol = (vol < 0.2) ? 0.2 : vol;
	vol = (vol > 1.8) ? 1.8 : vol;

	int n = cell_idx-1;
	if (k > 0 && n < NFRU) {
	  neighbors_init.push_back(n);
	}

	n = cell_idx+1;
	if (k < NFRU_Z-1 && n < NFRU) {
	  neighbors_init.push_back(n);
	}

	n = cell_idx-NFRU_Z;
	if (j > 0 && n < NFRU) {
	  neighbors_init.push_back(n);
	}

	n = cell_idx+NFRU_Z;
	if (j < NFRU_Y-1 && n < NFRU) {
	  neighbors_init.push_back(n);
	}

	n = cell_idx-NFRU_Y*NFRU_Z;
	if (i > 0 && n < NFRU) {
	  neighbors_init.push_back(n);
	}

	n = cell_idx+NFRU_Y*NFRU_Z;
	if (i < NFRU_X-1 && n < NFRU) {
	  neighbors_init.push_back(n);
	}

	FRU_neighbs[cell_idx] = neighbors_init;
	FRU_neighbs_distance[cell_idx] = std::vector<double>(neighbors_init.size());
	cell_idx++;
      }
    }
  }
  
  //Heterogeneous release site spacing (diffusion rates)
  
  //double dist0 = 0.248;
  double dist_avg = 1.0; //1.0-dist0; //0.139;
  double dist_std = 0.25;
  double dist_trans = 2.0;
  
  //  DisjointSets components(NFRU);

  for (int i = 0; i < FRU_neighbs.size(); i++) {
    for (int j = 0; j < FRU_neighbs[i].size(); j++) {
      int n = FRU_neighbs[i][j];
      double u1 = Get_Rand_Unif(mt[mythread],mti[mythread]);
      double dist;
      if (abs(i-n) == 1) {
	dist = dist_trans;
      } else {
	//dist = -dist_avg*log(u1) + dist0;
	dist = 1; //dist_avg + generateGaussian(mt[mythread],mti[mythread],0, dist_std);
	//dist = (dist_avg - 6 * 1.5*log(u1))/(dist_avg+6*1.5);
	if (dist < 0.05) dist = 0.05;
      }
      FRU_neighbs_distance[i][j] = dist;
      for (int k = 0; k < FRU_neighbs[n].size(); k++) {
	if (FRU_neighbs[n][k] == i) {
	  FRU_neighbs_distance[n][k] = dist;
	}
      }
      /*if (dist < dist_avg-dist_std) {
	components.Union(i,n);
	}*/
    }
  }

  //fprintf(stdout,"Number of superclusters:%d\n", components.NumSets());
  //fprintf(stdout,"Mean subclusters per supercluster: %g\n", (double)NFRU/((double)components.NumSets()));

  //Spark tests
  //FRUs[NFRU/2].RyR_state = 5;
  //FRUs[0].RyR_state = 5;
  /*
  for(iFRU = 0; iFRU < 100; iFRU++) {
  	//state[idx_Cai + iFRU] = 1e-3;
  	//state[idx_CaNSR + iFRU] = 1;
  }
  for(iFRU = 0; iFRU < 4; iFRU++) {
  	for (k = 0; k < 5; k++) {
  		//RyR_state[iFRU][0][k] = 3;
  	}
  	//state[idx_Cai+iFRU] = 0.1;
  	}
  
  
  for(icleft = 0; icleft < Nclefts_FRU; icleft++) {
    FRUs[NFRU/2].RyR_state[icleft] = NRyRs_per_cleft;
  }
  */
  /*for (int i = 0; i < NFRU; i+=NFRU_Z) {
    FRUs[i].RyR_state = NRyRs_per_cleft;
    }*/
  for (int i = 0; i < NFRU; i++) {
    //FRUs[i].RyR_state = 10;
  }
  /*for (int i = NFRU_Z-1; i < NFRU; i+=NFRU_Z) {
    FRUs[i].RyR_state = NRyRs_per_cleft;
    }*/
  /*for (int i = NFRU_Z/2; i < NFRU; i+=NFRU_Z) {
    FRUs[i].RyR_state = NRyRs_per_cleft;
    }*/
  //FRUs[0].RyR_state = NRyRs_per_cleft;

  //Randomly set LCC phosphorylation
  //Note there are much better algorithms for this but only need to do it once...

  //Assign exact percentage to random channels
  int N_LCC = NFRU*Max_LCCs_per_cleft;
  int N_Phosph = (int)(N_LCC*(beta_flag ? beta_frac_active : lcc_frac_active));
  int N_Assigned = 0;
  while (N_Assigned < N_Phosph) {
    bool bDone = 0;
    while (!bDone) {
      double runif = Get_Rand_Unif(mt[mythread],mti[mythread]);
      int FRU_attempt = (int)(NFRU*runif);
      if (FRUs[FRU_attempt].N_LCC_Active < Max_LCCs_per_cleft) {
				int n_act = FRUs[FRU_attempt].N_LCC_Active;
				FRUs[FRU_attempt].LCC_States[n_act] = 1;;
				FRUs[FRU_attempt].LCC_Vdep[n_act] = Oy_LType;
				FRUs[FRU_attempt].LCC_Mode2[n_act] = 0;
				FRUs[FRU_attempt].N_LCC_Active++;
				bDone = 1;
				N_Assigned++;
      }
    }
  }

  //Assign exact percentage to random channels
  int N_Active = N_Assigned;
  N_Phosph = (int)(N_Active*mode2_frac);
  N_Assigned = 0;
  while (N_Assigned < N_Phosph) {
    bool bDone = 0;
    while (!bDone) {
      double runif = Get_Rand_Unif(mt[mythread],mti[mythread]);
      int FRU_attempt = (int)(NFRU*runif);
      for (int i = 0; i < FRUs[FRU_attempt].N_LCC_Active; i++) {
	if (FRUs[FRU_attempt].LCC_Mode2[i] == 0) {
	  FRUs[FRU_attempt].LCC_Mode2[i] = 1;
	  bDone = 1;
	  N_Assigned++;
	  break;
	}
      }
    }
  }
	
	/*
  for(iFRU = 0; iFRU < NFRU; iFRU++) {
    for(icleft = 0; icleft < Nclefts_FRU; icleft++) {
      for (int j = 0; j < NLCCs_per_cleft; j++) {
				double runif;
				MersenneTwister_fast(&mti, mt, 1, &runif);
				if ( runif < (beta_flag ? beta_mode2_frac : lcc_mode2_frac)) {
					FRUs[iFRU].LType_state[icleft][j][index_LCC_Mode2] = 2;
				} else {
					FRUs[iFRU].LType_state[icleft][j][index_LCC_Mode2] = 1;
				}
      }
    }
  }
  */
  //Initialize residuals
  for (int i = 0; i < NFRU; i++) {
    double runif = Get_Rand_Unif(mt[mythread],mti[mythread]);
    FRUs[i].Ri = log(runif);
  }

  // This section initalizes the calculation of total call Ca,
  // which is used as one of the model algorithm verification tests

  double CaTOT_Trpn = state[index_LTRPNCa] * LTRPNtot + state[index_HTRPNCa] * HTRPNtot;
  double CaTOT_cyto = state[index_Cai] * (1.0 + CMDNtot / (KmCMDN + state[index_Cai]) + EGTAtot / (KmEGTA + state[index_Cai]));
  double CaTOT_NSR = state[index_CaNSR];
  double CaTOT_SL = 0; //state[index_CaSL] * (1.0 + CMDNtot / (KmCMDN + state[index_CaSL]) + EGTAtot / (KmEGTA + state[index_CaSL]) + BSLtot_SL / (KBSL + state[index_CaSL]));

  for (int i = 0; i < NFRU; i++) {
    CaTOT_SL += VSL * state[N_states+i] * (1.0 + CMDNtot / (KmCMDN + state[N_states+i]) + EGTAtot / (KmEGTA + state[N_states+i]) + BSLtot_SL / (KBSL + state[N_states+i]));
  }

  double CaTOT_SS = 0.0;
  double CaTOT_JSR = 0.0;

  for(iFRU = 0; iFRU < NFRU; iFRU++) {
    CaTOT_SS += FRUs[iFRU].Get_CaTOT_SS();
    CaTOT_JSR += FRUs[iFRU].Get_CaTOT_JSR();
  }

  state[index_CaTOT] = 1.e6 * (NFRU_scale * (CaTOT_JSR * VJSR + CaTOT_SS * VSS ) + CaTOT_SL*NFRU_scale + ((CaTOT_cyto + CaTOT_Trpn) * Vmyo  + CaTOT_NSR * VNSR)) ; //picomoles

  errweight = std::vector<double>(state.size());
  errweight[index_V] = 0; //1.e-2; // not independent
  errweight[index_mNa] = 1.0;
  errweight[index_hNa] = 1.0;
  errweight[index_jNa] = 1.;
  errweight[index_Nai] = 0.1;
  errweight[index_Ki] = 1 / 140.0;
  errweight[index_xKs] = 1.0;
  errweight[index_C0Kv43] = 1.0;
  errweight[index_C1Kv43] = 1.0;
  errweight[index_C2Kv43] = 1.0;
  errweight[index_C3Kv43] = 1.0;
  errweight[index_OKv43] = 1.0;
  errweight[index_CI0Kv43] = 1.0;
  errweight[index_CI1Kv43] = 1.0;
  errweight[index_CI2Kv43] = 1.0;
  errweight[index_CI3Kv43] = 1.0;
  errweight[index_OIKv43] = 0.0; // 1.0, not independent
  errweight[index_C0Kv14] = 1.0;
  errweight[index_C1Kv14] = 1.0;
  errweight[index_C2Kv14] = 1.0;
  errweight[index_C3Kv14] = 1.0;
  errweight[index_OKv14] = 1.0;
  errweight[index_CI0Kv14] = 1.0;
  errweight[index_CI1Kv14] = 1.0;
  errweight[index_CI2Kv14] = 1.0;
  errweight[index_CI3Kv14] = 1.0;
  errweight[index_OIKv14] = 0.0; // 1.0, not independent
  errweight[index_CaTOT] = 0.1;
  errweight[index_C1Herg] = 1.0;
  errweight[index_C2Herg] = 1.0;
  errweight[index_C3Herg] = 1.0;
  errweight[index_OHerg] = 1.0;
  errweight[index_IHerg] = 0.0; // 1.0, not independent
  errweight[index_Cai] = 1.0;
  errweight[index_CaNSR] = 1.0;
  errweight[index_CaSL] = 1.0;
  errweight[index_LTRPNCa] = 1.0;
  errweight[index_HTRPNCa] = 1.0;
  errweight[index_xKs2] = 1.0;
  for (int i = 0; i < NFRU; i++) {
    errweight[N_states+i] = 1.0/(double)NFRU;
  }
}









