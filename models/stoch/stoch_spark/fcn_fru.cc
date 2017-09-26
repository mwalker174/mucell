/*       ----------------------------------------------------

	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.

	 Name of Program: Local Control Model
	 Version: Documented Version, C
	 Date: September 2003

	 --------------------------------------------------*/

// NB: Indices fixed

#include <math.h>

#include "StochModel.h"


// Computes and returns the velocity field for CaSS and CaJSR
// within the FRU

void ReleaseUnit::fcn_fru(const double FRU_states[Nstates_FRU],
                          const double FRUdep_states[Nstates_FRUdep],
                          double dFRU_states1[Nstates_FRU],
                          int beta_flag,
			  double Cao)
{
  double V, CaJSR;
  double CaSS_1, CaSS_2, CaSS_3, CaSS_4;
  double VF_over_RT, VFsq_over_RT, exp_VFRT;
  double a1, a2;

  double OScale = bOrphan ? ORPHAN_SCALE : 1.0;

  const double JDHconstant = 1.0 / (2.0 * VSS * Faraday * 1000.0);

  double Jtr, beta_JSR; //, JRyRtot;
  double Jxfer_1;
  double JRyR_1;
  double JDHPR_1;
  double beta_SS_1;
  double LType_open; 
  LType_open = Get_LCC_Open_Cleft();
  double NRyR_open;
  NRyR_open = Get_RyR_Open_Cleft();

  V = FRUdep_states[index_frudep_V];
  double CaSL = FRUdep_states[index_frudep_CaSL];
  double CaNSR = FRUdep_states[index_frudep_CaNSR];

  CaJSR = FRU_states[index_frustates_CaJSR];
  CaSS_1 = FRU_states[1];

  Jtr = (CaNSR - CaJSR) / tautr;

  JRyR_1 = JRyRmax * NRyR_open * (CaJSR - CaSS_1);

  Jxfer_1 = (CaSS_1 - CaSL)/tauxfer;

  VF_over_RT = V / RT_over_F;
  VFsq_over_RT = (1000.0 * Faraday) * VF_over_RT;
  exp_VFRT = FRUdep_states[index_frudep_exp_VFRT]; //exp(2.0 * VF_over_RT);

  if (fabs(V) < 1.e-6) {
    JDHPR_1 = -PCa * 2.0 * 1000.0 * Faraday * (CaSS_1 * LCC_GAMMA - Cao * 0.341) * LType_open * JDHconstant;
  } else {
    // The same unrolled
    a2 = PCa * 4.0 * VFsq_over_RT / (exp_VFRT - 1.0) * JDHconstant;
    JDHPR_1 = -(CaSS_1 * exp_VFRT * LCC_GAMMA - Cao * 0.341) * a2 * LType_open;
  }

  // The same unrolled
  a1 = BSRtot * KBSR / (OScale*(CaSS_1 + KBSR) * (CaSS_1 + KBSR)); // + (165*1.1 / (OScale*(CaSS_1+1.1)*(CaSS_1+1.1)));
  a2 = BSLtot * KBSL / (OScale*(CaSS_1 + KBSL) * (CaSS_1 + KBSL));
  beta_SS_1 = 1.0 / (1.0 + a1 + a2) ;

  a1 = CSQNtot * KmCSQN / ((CaJSR + KmCSQN) * (CaJSR + KmCSQN));
  beta_JSR = 1.0 / (1.0 + a1);

  dFRU_states1[index_frustates_CaJSR] = beta_JSR * (Jtr - (VSS / VJSR) * JRyR_1); // dCaJSR

  dFRU_states1[1] = beta_SS_1 * ( JDHPR_1 + JRyR_1 - Jxfer_1 ) / OScale ; // dCaSS1
}
