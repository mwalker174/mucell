/*       ----------------------------------------------------

				 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
				 Copyright 2003, The Johns Hopkins University
				 School of Medicine. All rights reserved.

				 Name of Program: Local Control Model
				 Version: Documented Version, C
				 Date: November 2003

				 --------------------------------------------------

				 simfru_local.c - This subroutine runs the Monte Carlo simulation
				 from time0 to timef for all of the stochastic elements
				 within a single release unit (indexed by iFRU).
				 It runs in incremental time steps determined based
				 on the transition rates of currently occupied states.
				 On each step, the generation of random numbers are
				 used to determine into which state each stochastic
				 element has transitioned, if a transition has occured.
				 Local Ca2+ concentrations are integrated within these
				 time steps.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <cassert>

#include "StochModel.h"

double ReleaseUnit::Get_RyR_Open_Cleft()
{
  return (double)RyR_state;
}

double ReleaseUnit::Get_LCC_Open_Cleft()
{
  double N_open = 0;
    for (int j = 0; j < N_LCC_Active; j++) {
      if ((LCC_States[j] == O1_LType 
	|| LCC_States[j] == O2_LType) 
	  && (LCC_Vdep[j] == Oy_LType)) {
	N_open += 1.0;
      }
    }
    return N_open;
  //return N_LCC_Open;
}

void ReleaseUnit::simfru(const double time0,
                         const double timef,
                         const double FRUdep_states0[Nstates_FRUdep],
			 const double FRUdep_statesf[Nstates_FRUdep],
			 unsigned long mt[mtN+1], 
			 int &mti,
                         int beta_flag,
                         double Cao,
			 double ryr_kplus_scale)
{
  double FRU_states1[Nstates_FRU];
  double dFRU_states[Nstates_FRU];
  double k1[Nstates_FRU], y_1[Nstates_FRU];
  double Ito2_rates[Nstates_Ito2];
  double RyR_trans_rates[Nstates_RyR];
  double time_FRU, time_stepFRU;
  double FRUdep_states[Nstates_FRUdep];

  double LCC_rates[Max_LCCs_per_cleft][3];
  char LCC_index[Max_LCCs_per_cleft][3];
  double LCC_Vdep_rates[Max_LCCs_per_cleft];

  for (int i = 0; i < Nstates_FRUdep; i++) {
    FRUdep_states[i] = FRUdep_states0[i];
  }

  time_FRU = time0;
  // V       FRUdep_states[0]
  // CaSL     FRUdep_states[1]
  // CaNSR   FRUdep_states[2]
  // CaJSR   FRU_states[0]
  // CaSS1   FRU_states[1]
  // CaSS2   FRU_states[2]
  // CaSS3   FRU_states[3]
  // CaSS4   FRU_states[4)

  char bDone = 0;
  double Ri0, R_sum, alpha_frudep;

  while (!bDone) {

    //Sum of residuals at beginning of time step
    Ri0 = Ri;

    R_sum = fru_rates(FRUdep_states, LCC_rates, LCC_index, LCC_Vdep_rates, RyR_trans_rates, Ito2_rates, beta_flag, ryr_kplus_scale);
    
    //Time step
    time_stepFRU = timef - time_FRU;
    double tstep_max = (FRU_states[1] < 1.5e-3) ? 10e-3 : 5e-3;

    if (time_stepFRU > tstep_max) {
      time_stepFRU = tstep_max;
      bDone = 0;
    } else {
      bDone = 1;
    }

    //update residuals
    Ri += R_sum * time_stepFRU;

    //decrease time step if residual crosses zero
    if (Ri > 0) {
      //Determine time it took for reaction to fire
      time_stepFRU = time_stepFRU * Ri0 / (Ri0 - Ri);
    }

    //FRU states
    for(int i = 0; i < Nstates_FRU; i++) {
      FRU_states1[i] = FRU_states[i];
    }

    // local velocity field for Ca concentrations within the unit
    fcn_fru(FRU_states1, FRUdep_states, dFRU_states, beta_flag, Cao);

    for(int i = 0; i < Nstates_FRU; i++) { // First intermediate integration step (Trapezoidal Method)
      k1[i] = time_stepFRU * dFRU_states[i];
      y_1[i] = FRU_states1[i] + k1[i];
    }

    time_FRU = time_FRU + time_stepFRU;

    alpha_frudep = (time_FRU-time0)/(timef-time0);
    for (int i = 0; i < Nstates_FRUdep; i++) {
      FRUdep_states[i] = FRUdep_states0[i] + alpha_frudep*(FRUdep_statesf[i]-FRUdep_states0[i]);
    }

    // event occurs at delta t
    // local velocity field for Ca concentrations within the unit
    fcn_fru(y_1, FRUdep_states, dFRU_states, beta_flag, Cao);

    for(int i = 0; i < Nstates_FRU; i++) {	// Second intermediate integration step (Trapezoidal Method)
      FRU_states[i] = FRU_states[i] + (k1[i] + time_stepFRU * dFRU_states[i]) / 2.0;
    }

    if (FRU_states[1] < 1e-9) { 
      printf("Warning: CaSS < 1e-9, clamping!\n");
      FRU_states[1] = 1e-9;
    }

    /*
    time_FRU = time_FRU + time_stepFRU;

    alpha_frudep = (time_FRU-time0)/(timef-time0);
    for (int i = 0; i < Nstates_FRUdep; i++) {
      FRUdep_states[i] = FRUdep_states0[i] + alpha_frudep*(FRUdep_statesf[i]-FRUdep_states0[i]);
      }
    
    // event occurs at delta t
    // local velocity field for Ca concentrations within the unit
    fcn_fru(FRU_states, FRUdep_states, dFRU_states, beta_flag, Cao);

    for(int i = 0; i < Nstates_FRU; i++) {	// Second intermediate integration step (Trapezoidal Method)
      FRU_states[i] += time_stepFRU * dFRU_states[i];
      }*/

    //effect reactions if residual crosses zero
    if (Ri > 0) {
      //Not done with CaRU simulation
      bDone = 0;

      //Generate uniform random number on (0,1)
      double runif;
      MersenneTwister_fast(&mti, mt, 1, &runif);

      //fire reaction
      fru_fire_rxn(runif * R_sum, LCC_rates, LCC_index, LCC_Vdep_rates, RyR_trans_rates, Ito2_rates);

      //Generate new exponential random variable and assign to residual
      MersenneTwister_fast(&mti, mt, 1, &runif);
      Ri = log(runif);

    } //if fire rxn

  }
}
