#include "StochModel.h"
#include <math.h>


//Determines which reaction fired based on rates and effects reaction
void ReleaseUnit::fru_fire_rxn(double rxn_rnd,
			       double LCC_rates[Max_LCCs_per_cleft][3],
			       char LCC_index[Max_LCCs_per_cleft][3],
			       double LCC_Vdep_rates[Max_LCCs_per_cleft],
                               double RyR_rates[Nstates_RyR],
                               double Ito2_rates[Nstates_Ito2])
{
  double a_sum = 0;
  int i;
  
  for (int j = 0; j < N_LCC_Active; j++) {
    
    //LCC states
    for (i = 0; i < 3; i++) {
      a_sum += LCC_rates[j][i];
      
      if (a_sum > rxn_rnd) {
	/*
	if ((LCC_index[j][i] == O1_LType 
	    || LCC_index[j][i] == O2_LType) 
	    && (LCC_Vdep[j] == Oy_LType)) {
	  N_LCC_Open++;
	} else if ((LCC_States[j] == O1_LType 
		   || LCC_States[j] == O2_LType) 
		   && (LCC_Vdep[j] == Oy_LType)) {
	  N_LCC_Open--;
	}
      */
	LCC_States[j] = LCC_index[j][i]; 
	return;
      }
    }
    
    //LCC v-dep inactivation
    a_sum += LCC_Vdep_rates[j];
    
    if (a_sum > rxn_rnd) {
      /*
      if ((LCC_States[j] == O1_LType 
	   || LCC_States[j] == O2_LType) 
	  && (LCC_Vdep[j] == Cy_LType)) {
	N_LCC_Open++;
      } else if ((LCC_States[j] == O1_LType 
		  || LCC_States[j] == O2_LType) 
		 && (LCC_Vdep[j] == Oy_LType)) {
	N_LCC_Open--;
      }
      */
      switch (LCC_Vdep[j]) {
      case Oy_LType:
        LCC_Vdep[j] = Cy_LType;
        return;
	
      case Cy_LType:
        LCC_Vdep[j] = Oy_LType;
        return;
	
      default: //	Unknown	state
        fprintf(stderr, "fru_fire_rxn(): Unknown state %d in LCC Vinact\n",LCC_States[j]);
        return;
      } //switch LType vinact
    } //if
    
  }
  
  //RyR
  for (i = 0; i < Nstates_RyR; i++) {
    a_sum += RyR_rates[i];
    
    if (a_sum > rxn_rnd) {
      if (i==0) {
	RyR_state++;
      } else {
	RyR_state--;
      }
      return;
    }
  }
  
  //Ito2
  for (i = 0; i < Nstates_Ito2; i++) {
    a_sum += Ito2_rates[i];
    
    if (a_sum > rxn_rnd) {
      if (i==0) {
	Ito2_state++;
      } else {
	Ito2_state--;
      }
      return;
    } //ito2 if
  } //for i
  
  
  //Warn if we never reach 1
  printf("fru_fire_rxn(): Warning: reaction probabilities did not sum to one. (sum = %f, rxn_rnd = %f)\n", a_sum, rxn_rnd);
}

