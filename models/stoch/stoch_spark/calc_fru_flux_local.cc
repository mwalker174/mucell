/*       ----------------------------------------------------

	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.

	 Name of Program: Local Control Model
	 Version: Documented Version, C
	 Date: November 2003

	 --------------------------------------------------

	 calc_fru_flux_local.c - This subroutine is called by rk4am preceding
	                         each call to fcn, and returns the values of
							 fluxes which cross release unit boundaries.
							 This routine sums over all the units to calculate
							 these fluxes (Jxfer,Jtr,ICa,Ito2)

	 A lot of work is moved to mpi_master
*/


#include <math.h>

#include "StochModel.h"

void ReleaseUnit::calc_fru_flux(double num[5])
{
  double ICa_numerator;
  double sum_CaSS_local;
  double sum_CaJSR_local;

  double OCa_numerator;
  int NIto2_Open;

  ICa_numerator = 0.0;
  OCa_numerator = 0;
  NIto2_Open = 0;

  double N_open_lcc = Get_LCC_Open_Cleft();

  //ICa_numerator = ICa_numerator + FRU_states[1]*exp_2VFRT-Cao*0.341;
  ICa_numerator = ICa_numerator + N_open_lcc*FRU_states[1] * LCC_GAMMA;
  OCa_numerator = OCa_numerator + N_open_lcc;

  NIto2_Open = NIto2_Open + Ito2_state;

  sum_CaSS_local = 0.0;
  sum_CaJSR_local = 0.0;

  sum_CaJSR_local += FRU_states[index_frustates_CaJSR];

  sum_CaSS_local += FRU_states[1];

  num[0] += sum_CaSS_local;
  num[1] += sum_CaJSR_local;
  num[2] += ICa_numerator;
  num[3] += OCa_numerator;
  num[4] += NIto2_Open;
}








