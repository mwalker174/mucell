/*       ----------------------------------------------------

				 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
				 Copyright 2003, The Johns Hopkins University
				 School of Medicine. All rights reserved.

				 Name of Program: Local Control Model
				 Version: Documented Version, C
				 Date: November 2003

				 --------------------------------------------------

				 fru_rates_local.c - This subroutine returns the individual rates
				 for leaving the current state, the sum of those rates,
				 and an index array that specifies the destination
				 state that corresponds to each exit rate. All of
				 these quantities are returned for each RyR and
				 the LType channel and the V-dep inactivation gate
				 of the LType channel and the Ito2 channel
				 Also returned is the maximum of the sum of rates
				 leaving any current state.
*/

#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>

#include "StochModel.h"

double ReleaseUnit::fru_rates(const double FRUdep_states[Nstates_FRUdep],
			      double LCC_rates[Max_LCCs_per_cleft][3],
			      char LCC_index[Max_LCCs_per_cleft][3],
			      double LCC_Vdep_rates[Max_LCCs_per_cleft],
                              double RyR_rates[Nstates_RyR],
                              double Ito2_rates[Nstates_Ito2],
                              int beta_flag,
			      double ryr_kplus_scale)
{
  // real	rnum
  double CaSS;

  const double fL = 0.85; // transition	rate into open state (1/ms)
  const double gL = 2.0; //transition rate	out	of open	state (1/ms)
  const double fLprime = 0.005;	// transition rate into	Ca mode	open state (1/ms)
  const double gLprime = 7.0; // transition	rate out of	Ca mode	open state (1/ms)
  const double bL = 1.9356;		// mode	transition parameter
  const double bL2 = bL * bL;
  const double bL3 = bL * bL * bL;
  const double bL4 = bL * bL * bL * bL;
  const double kappa = 1;
  const double kappa2 = kappa * kappa;
  const double kappa3 = kappa * kappa * kappa;
  const double kappa4 = kappa * kappa * kappa * kappa;
  const double aL = 2.0; //	mode transition	parameter
  const double aL2 = aL * aL;
  const double aL3 = aL * aL * aL;
  const double aL4 = aL * aL * aL * aL;
   double omega = 0.83 * 2.0 * 1.3 * 0.01; // mode transition parameter	(1/ms)
  const double gPhosph = 0.049 * 1.75; //REFLECTS CHANGE IN GATING WHEN CHANNEL IS PHOSPHORYLATED

  const double alphacf = 4.0 * 1.2 * 0.416;
  const double betacf = 4.0 * 0.45 * 0.049;
  const double gammacf = 0.83 * 1.9 * 1.3 * 0.31 * 7.5 * 0.09233;	// (ms-1 mM-1) 

  double alpha, beta, alpha_prime, beta_prime, gamma_rate;
  double C0_to_C1, C1_to_C2, C2_to_C3, C3_to_C4;
  double C1_to_C0, C2_to_C1, C3_to_C2, C4_to_C3;
  double CCa0_to_CCa1, CCa1_to_CCa2;
  double CCa2_to_CCa3, CCa3_to_CCa4;
  double CCa1_to_CCa0, CCa2_to_CCa1;
  double CCa3_to_CCa2, CCa4_to_CCa3;
  double C0_to_CCa0, C1_to_CCa1, C2_to_CCa2;
  double C3_to_CCa3, C4_to_CCa4;

  const double CCa0_to_C0	= kappa4 * omega;		// = omega
  const double CCa1_to_C1	= kappa3 * omega / bL;	// = omega/bL
  const double CCa2_to_C2	= kappa2 * omega / bL2;	// = omega/bL^2
  const double CCa3_to_C3	= kappa * omega / bL3;	// = omega/bL^3
  const double CCa4_to_C4	= omega / bL4;	// = omega/bL^4


  double yCa_inf,	tau_yCa;
  //Moved to set_FRUdep_states:
  //const double yCa_frac = 0.4;	// asymptotic value	for	fraction of	LCCs that
  // voltage-inactivate at depolarized potentials
  //	int	Ito2_state[Nclefts_FRU][NFRU_sim_max];
  //	double Ito2_exitrate[Nclefts_FRU];
  const double KdIto2 = 0.1502;	 //	(mM)
  const double kbIto2 = 2.0; // (ms-1)
  const double kfIto2 = kbIto2 / KdIto2;	// (ms-1 mM-1)

  double rate_sum = 0.0;
  CaSS = FRU_states[1];
  //double CaSL = FRUdep_states[index_frudep_Cai];

  alpha =	alphacf	* FRUdep_states[index_frudep_exp_alpha]; //exp(0.012 * (V - 35.0));
  beta = betacf *	FRUdep_states[index_frudep_exp_beta]; //exp(-0.05 * (V - 35.0));
  alpha_prime	= aL * alpha;
  beta_prime = beta / (bL * kappa);

  //const double gammakd1 = 8e-3; //8e-3; //8e-3; //TODO (4e-3; works with 5:1 ratio 12/9/15)
  //gamma_rate =	  gammacf  * CaSS * (gammakd1/(CaSS+gammakd1)); // * gammakd1 / (CaSS + gammakd1);

  //gamma_rate = 0.08*CaSS; //0.2*CaSS;
  //gamma_rate = 40e-3*pow(CaSS,0.5); //0.1*pow(CaSS,0.9); //0.05*pow(CaSS,0.75);
  //gamma_rate = 3e-3 * CaSS * CaSS * CaSS / (pow(1e-3,3.0) + CaSS*CaSS*CaSS);
  gamma_rate = 2.2e-3 * CaSS / (3e-3 + CaSS); //3e-3 //15e-3

  for (int j = 0; j < N_LCC_Active; j++) {

    LCC_rates[j][1] = 0;
    LCC_rates[j][2] = 0;
    LCC_index[j][1] = 0;
    LCC_index[j][2] = 0;
    
    switch (LCC_States[j]) {
      
    case 1:	// if LType	is in state	1
      C0_to_C1 = 4.0 * alpha;
      C0_to_CCa0 = gamma_rate;
      LCC_rates[j][0] = C0_to_C1;		// rate	from 1 to 2
      LCC_rates[j][1] = C0_to_CCa0;		// rate	from 1 to 7
      LCC_index[j][0] = 2;			// index(2)	state into which rates(2) takes	you	= 2
      LCC_index[j][1] = 7;			// index(3)	state into which rates(3) takes	you	= 7
      break;
      
    case 2:
      C1_to_C2 = 3.0 * alpha;
      C1_to_C0 =		beta;
      C1_to_CCa1 = aL * gamma_rate;
      LCC_rates[j][0] = C1_to_C0;
      LCC_rates[j][1] = C1_to_C2;
      LCC_rates[j][2] = C1_to_CCa1;
      LCC_index[j][0] = 1;
      LCC_index[j][1] = 3;
      LCC_index[j][2] = 8;
      break;
      
    case 3:
      C2_to_C1 = 2.0 * beta;
      C2_to_C3 = 2.0 * alpha;
      C2_to_CCa2 = aL2 * gamma_rate;
      LCC_rates[j][0] = C2_to_C1;
      LCC_rates[j][1] = C2_to_C3;
      LCC_rates[j][2] = C2_to_CCa2;
      LCC_index[j][0] = 2;
      LCC_index[j][1] = 4;
      LCC_index[j][2] = 9;
      break;
      
    case 4:
      C3_to_C4 =		alpha;
      C3_to_C2 = 3.0 * beta;
      C3_to_CCa3 = aL3 * gamma_rate;
      LCC_rates[j][0] = C3_to_C2;
      LCC_rates[j][1] = C3_to_C4;
      LCC_rates[j][2] = C3_to_CCa3;
      LCC_index[j][0] = 3;
      LCC_index[j][1] = 5;
      LCC_index[j][2] = 10;
      break;
	
    case 5:
      C4_to_C3 = 4.0 * beta;
      C4_to_CCa4 = aL4 * gamma_rate;
      LCC_rates[j][0] = C4_to_C3;
      LCC_rates[j][1] = fL;
      LCC_rates[j][2] = C4_to_CCa4;
      LCC_index[j][0] = 4;
      LCC_index[j][1] = 6;
      LCC_index[j][2] = 11;
      break;
      
    case 6:
      if(LCC_Mode2[j] == 1) {
	LCC_rates[j][0] = gL * gPhosph;
      } else if(LCC_Mode2[j] == 0) {
	LCC_rates[j][0] = gL;
      } else {
	fprintf(stderr, "Unknown LCC_Mode2 state %d in LCC local\n",LCC_Mode2[j]);
      }
      LCC_index[j][0] = 5;
      break;
      
    case 7:
      //		CCa0_to_C0 = omega	// constant
      CCa0_to_CCa1 = 4.0 * alpha_prime;
      LCC_rates[j][0] = CCa0_to_C0;
      LCC_rates[j][1] = CCa0_to_CCa1;
      LCC_index[j][0] = 1;
      LCC_index[j][1] = 8;
      break;
      
    case 8:
      //		CCa1_to_C1 = omega/bL  // constant
      CCa1_to_CCa2 = 3.0 * alpha_prime;
      CCa1_to_CCa0 =		beta_prime;
      LCC_rates[j][0] = CCa1_to_CCa0;
      LCC_rates[j][1] = CCa1_to_C1;
      LCC_rates[j][2] = CCa1_to_CCa2;
      LCC_index[j][0] = 7;
      LCC_index[j][1] = 2;
      LCC_index[j][2] = 9;
      break;
      
    case 9:
      //		CCa2_to_C2 = omega/bL2	// constant
      CCa2_to_CCa3 = 2.0 * alpha_prime;
      CCa2_to_CCa1 = 2.0 * beta_prime;
      LCC_rates[j][0] = CCa2_to_CCa1;
      LCC_rates[j][1] = CCa2_to_C2;
      LCC_rates[j][2] = CCa2_to_CCa3;
      LCC_index[j][0] = 8;
      LCC_index[j][1] = 3;
      LCC_index[j][2] = 10;
      break;
      
    case 10:
      //		CCa3_to_C3 = omega/bL3	// constant
      CCa3_to_CCa4 =		alpha_prime;
      CCa3_to_CCa2 = 3.0 * beta_prime;
      LCC_rates[j][0] = CCa3_to_CCa2;
      LCC_rates[j][1] = CCa3_to_C3;
      LCC_rates[j][2] = CCa3_to_CCa4;
      LCC_index[j][0] = 9;
      LCC_index[j][1] = 4;
      LCC_index[j][2] = 11;
      break;
      
    case 11:
      //		CCa4_to_C4 = omega/bL4	// constant
      CCa4_to_CCa3 = 4.0 * beta_prime;
      LCC_rates[j][0] = CCa4_to_CCa3;
      LCC_rates[j][1] = CCa4_to_C4;
      LCC_rates[j][2] = fLprime;
      LCC_index[j][0] = 10;
      LCC_index[j][1] = 5;
      LCC_index[j][2] = 12;
      break;
      
    case 12:
      LCC_rates[j][0] = gLprime;
      LCC_index[j][0] = 11;
      break;
      
    default: //	Unknown	state
      fprintf(stderr, "Unknown state %d in LCC  local\n",LCC_States[j]);
      break;
    }
    for (int i = 0; i < 3; i++) {
      rate_sum += LCC_rates[j][i];
    }
  }// for
  
  
  // Voltage dependent inactivation gate
  //double nu = 1; //1.5; //moved to set_FRUdep_states
  yCa_inf	= FRUdep_states[index_frudep_exp_inf]; //nu * ( yCa_frac / (1.0 + exp((V + 12.5) / 5.0)) + (1.0 - yCa_frac));
  //ORIGINALLY, MULTIPLIED TAU BY NU
  tau_yCa	=  FRUdep_states[index_frudep_exp_tau]; //60.0 + 340.0 / (1.0	+ exp((V + 30.0) / 12.0));

  for (int j = 0; j < N_LCC_Active; j++) {
    switch (LCC_Vdep[j]) {
    case Oy_LType:
      LCC_Vdep_rates[j]	= (1.0 - yCa_inf) / tau_yCa;
      break;
      
    case Cy_LType:
      LCC_Vdep_rates[j]	= yCa_inf / tau_yCa;
      break;
      
    default: //	Unknown	state
      fprintf(stderr, "Unknown state %d in LCC Vinact local\n",LCC_Vdep[j]);
      break;
    }
    rate_sum += LCC_Vdep_rates[j];
  }


  //RyR (Super-resolution)
  
  //2-State RyR Model
  double RyR_Rate_Open;
  double kplus =  0.1107e-3 * 1e6;
  double kminus = 0.5;
  double phi = 0.8025 + pow(FRU_states[index_frustates_CaJSR]/1.5,4);
  //double phi = 0.8025; 
  //double phi = FRU_states[index_frustates_CaJSR]/0.5;
  //double phi = 0.8025 + pow(FRU_states[index_frustates_CaJSR]/1.0,4);
  //double CJSR = FRU_states[index_frustates_CaJSR];
  //double phi = (CJSR < 0.2) ? 0 : 1;

  //Allosteric coupling (mean-field)
  //  	double astar = 0;//0.02*(49.0/28.0);
  //  	const double ecc = -0.92;
  //  	const double eoo = -0.85;
  //  	int N_O = 0;
  //  	for (icleft = 0;icleft<Nclefts_FRU;icleft++) {
  //  		for(i=0;i<NRyRs_per_cleft;i++)	{
  //  			if (RyR_state[i]==2) {
  //  				N_O++;
  //  			}
  //  		}
  //  	}
  //  	int N_C = NRyRs_per_cleft*Nclefts_FRU-N_O;
  //double Xoc = 1; //exp(-astar*(N_C*ecc - (N_O-1)*eoo));
  //double Xco = 1; //exp(-astar*(N_O*eoo - (N_C-1)*ecc));

  //TODO
  //  if (beta_flag) {
  kplus *= ryr_kplus_scale;
  //}
  if (HF_flag) {
    kplus *= HF_RyR_Scale;
  }
  
  double ryro = kplus * phi;
  double RyR_Rate_Close = kminus;
  double eta = 2.1;

  RyR_Rate_Open = ryro * pow(CaSS, eta);
  RyR_rates[0] = (NRyRs_per_cleft-RyR_state)*RyR_Rate_Open;
  RyR_rates[1] = RyR_state*RyR_Rate_Close;
  rate_sum += RyR_rates[0] + RyR_rates[1];
    
  //RyR rate constants
  /*const double k_iCa=0.2967;
  const double k_im = 0.0055;
  const double k_oCa = 6500;
  const double k_om = 0.0683;
  const double MinSR = 1.0;
  const double MaxSR = 15.0;
  const double EC50_SR = 0.45;
  const double H_SR = 2.5;
  double k_CaSRden = 1 + pow(EC50_SR/(FRU_states[index_frustates_CaJSR]),H_SR);
  double k_CaSR = MaxSR - ((MaxSR - MinSR)/k_CaSRden);
  double k_oSRCa = k_oCa/k_CaSR;
  double k_iSRCa = k_iCa*k_CaSR;*/

  /*double k_RItoR = k_im;
  double k_OtoR  = k_om;
  double k_ItoO  = k_im;
  double k_ItoRI = k_om;
  double k_RtoO, k_RItoI, k_RtoRI, k_OtoI;*/
  /*
  double k_RItoR = 0.01;
  double k_OtoR  = 0.5;
  double k_ItoO  = 0;
  double k_ItoRI = 1;
  double k_RItoI = 1;
  double k_RtoRI = 0;
  double k_OtoI  = pow(0.65,23) / (pow(0.65,23) + pow(FRU_states[index_frustates_CaJSR],23));
  double k_RtoO;

  for(icleft = 0; icleft<Nclefts_FRU;icleft++) {

  	k_RtoO  = 0.1107e-3 * 1e6 * 10 * CaSS[icleft]*CaSS[icleft];
  	//k_RItoI = k_oSRCa*CaSS[icleft]*CaSS[icleft];
  	//k_RtoRI = k_iSRCa*CaSS[icleft];
  	//k_OtoI  = k_iSRCa*CaSS[icleft];

  	for(i=0;i<NRyRs_per_cleft;i++) {

  		switch (RyR_state[icleft][i]) {

  		case 1:
  		  rate_sum += k_RItoR + k_RItoI;
  	    RyR_rates[icleft][i][1] = k_RItoR;
  	    RyR_rates[icleft][i][2] = k_RItoI;

  	    RyR_length[icleft][i] = 3;
  	    //				RyR_index[icleft][i][0] = 1;
  	    RyR_index[icleft][i][1] = 2;
  	    RyR_index[icleft][i][2] = 4;
  	    break;

  		case 2:
  	    rate_sum += k_RtoRI + k_RtoO;
  	    RyR_rates[icleft][i][1] = k_RtoRI;
  	    RyR_rates[icleft][i][2] = k_RtoO;

  	    RyR_length[icleft][i] = 3;
  	    //				RyR_index[icleft][i][0] = 2;
  	    RyR_index[icleft][i][1] = 1;
  	    RyR_index[icleft][i][2] = 3;
  	    break;

  		case 3:
  	    rate_sum += k_OtoR + k_OtoI;
  	    RyR_rates[icleft][i][1] = k_OtoR;
  	    RyR_rates[icleft][i][2] = k_OtoI;

  			RyR_length[icleft][i] = 3;
  	    //				RyR_index[icleft][i][0] = 3;
  	    RyR_index[icleft][i][1] = 2;
  	    RyR_index[icleft][i][2] = 4;
  	    break;

  		case 4:
  	    rate_sum += k_ItoRI + k_ItoO;
  	    RyR_rates[icleft][i][1] = k_ItoRI;
  	    RyR_rates[icleft][i][2] = k_ItoO;

  	    RyR_length[icleft][i] = 3;
  	    //				RyR_index[icleft][i][0] = 4;
  	    RyR_index[icleft][i][1] = 1;
  	    RyR_index[icleft][i][2] = 3;
  	    break;

  		default: //	Unknown	state
  			printf("Unknown state %d in RyR %d cleft % d\n",RyR_state[icleft][i],i,icleft);
  			break;
      }//switch
    }//for
  }//for
  */
  // Ito2
  double ito2_rate_scale = 1.0;
  Ito2_rates[0] = ito2_rate_scale * kfIto2 * FRUdep_states[index_frudep_CaSL] * (NIto2_per_cleft - Ito2_state);
  Ito2_rates[1] = kbIto2 * Ito2_state;
  
  rate_sum += Ito2_rates[0] + Ito2_rates[1];
  
  return rate_sum;
}

