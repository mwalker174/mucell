/*       ----------------------------------------------------

	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.

	 Name of Program: Local Control Model
	 Version: Documented Version, C
	 Date: November 2003

	 --------------------------------------------------

	 calc_fru_avg_local.c - compute average FRU properties

*/

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "StochModel.h"

// Calculate additional output quantities that are algebraically
// related to the state variables and stored in 'otherstate'

void ReleaseUnit::calc_fru_avg(double num_stat[Nstat])
{
  double CaSS, CaJSR, CaTOT_SS, CaTOT_JSR, CaTOT_Cyto, CaTOT_NSR, CaTOT_TRPN;
  double JRyR, PRyR_Open, PNorm_Mode, PnotVinact, PLType_Open, PIto2_Open, NRyR_ready;

  CaJSR = 0.0;
  CaTOT_JSR = 0.0;
  JRyR = 0.0;
  PRyR_Open = 0.0;
  PNorm_Mode = 0.0;
  PnotVinact = 0.0;
  PLType_Open = 0.0;
  PIto2_Open = 0.0;
  NRyR_ready = 0;
  CaTOT_Cyto = 0;
  CaTOT_NSR = 0;
  CaTOT_TRPN = 0;

  CaJSR = FRU_states[index_frustates_CaJSR];

  CaTOT_JSR = FRU_states[index_frustates_CaJSR] * (1.0 + CSQNtot / (KmCSQN + FRU_states[index_frustates_CaJSR]) );

  CaTOT_SS = 0.0;

    CaTOT_SS += FRU_states[1]
                * (1.0 + BSRtot / (KBSR + FRU_states[1]) + BSLtot / (KBSL + FRU_states[1]));

  CaSS = 0.0;
  double N_RyR_Open_Cleft = Get_RyR_Open_Cleft();
  
  PRyR_Open = Get_RyR_Open();
  NRyR_ready = Get_RyR_Ready();
  PLType_Open = Get_LCC_Open();
  PNorm_Mode = Get_LCC_NormMode();
  PnotVinact = Get_LCC_Vinact();
  
  CaSS  += FRU_states[1];
  
  JRyR += JRyRmax * N_RyR_Open_Cleft * (FRU_states[index_frustates_CaJSR] - FRU_states[1]);
    
  PIto2_Open += Ito2_state;
    
  // The following values are sums, not averages
  num_stat[0] += CaSS;
  num_stat[1] += CaJSR;
  num_stat[2] += JRyR;
  num_stat[3] += PRyR_Open;
  num_stat[4] += NRyR_ready;
  num_stat[5] += CaTOT_SS;
  num_stat[6] += CaTOT_JSR;
  num_stat[7] += PNorm_Mode;
  num_stat[8] += PnotVinact;
  num_stat[9] += PLType_Open;
  num_stat[10] += PIto2_Open;
  num_stat[11] += CaTOT_Cyto;
  num_stat[12] += CaTOT_NSR;
  num_stat[13] += CaTOT_TRPN;
}


