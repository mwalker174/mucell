/*       ----------------------------------------------------

	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.

	 Name of Program: Local Control Model
	 Version: Documented Version, C
	 Date: November 2003

	 --------------------------------------------------

	 parameters.h - constants, definitions and prototypes

*/
#ifndef _parameters_header
#define _parameters_header

// definitions of the model
//#define Nclefts_FRU 1 //4
#define Nstates_FRU 2 //(1+Nclefts_FRU)

#define NRyRs_per_cleft 48 
#define Max_LCCs_per_cleft 8 //10
#define NIto2_per_cleft 8
#define Nstates_RyR 2
#define Nstates_Ito2 2
//#define Nindepstates_LType 4

#define NRVseqs_per_cleft (2+NRyRs_per_cleft+1)
#define NRVseqs_per_FRU (NRVseqs_per_cleft*Nclefts_FRU)

#define NFRU_total 25000 //15000 //25000

//	//mt19937 parameters
#define mtN 624
#define OMP_THREAD_MAX 24
#define Nstat 15

#define V_CLAMP 0
#define LCC_GAMMA 1
#define ORPHAN_SCALE 100.0 //30.0
#define HF_flag 0
#define HF_NaCa 1.75
#define HF_SERCA 0.38
#define HF_Kv 0.34
#define HF_RyR_Scale 5.0

//#define clamp_nai 0
#define clamp_ki 1

/*------PHYSICAL constANTS----------------------------------------------------*/

#define Faraday 96.5 // Faraday's constant (C/mmol)
#define Temp 310.0 // absolute temperature (K)
#define Rgas 8.314 // ideal gas constant (J/[mol*K])
#define RT_over_F (Rgas*Temp/Faraday) // (Rgas*Temp/Faraday); // Rgas*Temp/Faraday (mV)

// Open and Closed states in various MC models

#define O1_LType 6
#define O2_LType 12
#define Oy_LType 2
#define Cy_LType 1
//#define O1_RyR 2
//#define C_RyR 1
//#define O_Ito2 2
//#define C_Ito2 1

//Random number generator functions
extern "C"
{
  double MersenneTwisterOne(int *mti, unsigned long[mtN + 1]);
  void MersenneTwister(int *, unsigned long[mtN + 1], const int, double *);
  void MersenneTwister_fast(int *, unsigned long [mtN + 1], const int , double *);
}

#endif
