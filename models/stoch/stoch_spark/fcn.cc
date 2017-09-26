/*       ----------------------------------------------------

	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.

	 Name of Program: Local Control Model
	 Version: Documented Version, C
	 Date: November 2003

	 --------------------------------------------------

	 fcn.c - This subroutine calculates the velocity field (state derivatives)
	         for states defined by the global model. Voltage clamp protocol is
                 implemented here as well. Currents and fluxes that cross release
		 unit boundaries are not calculated here, but have been passed in
		 as arguments (Jxfer,Jtr,ICa,Ito2).

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>

#include "StochModel.h"


void StochModel::fcn(const double time, std::vector<double> &state, std::vector<double> &Fstate, double current[Ncur],
                     const int keepc, const double Jxfer, const double Jtr, const double ICa, const double Ito2, const double Iext, const double JCa)
{
  //-------STATE VARIABLE AND STORAGE ARRAYS-------------------------------------

  //	double state[N];	// dynamic states
  //	double Fstate[N];	// state derivatives
  //	double current[Ncur];  // currents and fluxes
  //	int keepc; 			// flags when to store currents

  //-------LOCAL STATE VARIABLE ARRAYS-------------------------------------------

  double V;	// membane potential
  double mNa;	// INa activation
  double hNa;	// INa fast inactivation
  double jNa;	// INa slow inactivation

  double xKs;	// IKs activation
  double Nai;	// intracellular Na+ conc.
  double Ki;	// intracellular K+ conc.
  double Cai;	// intracellular Ca++ conc.
  double CaSL;	// subsarcolemmal Ca++ conc.
  double CaNSR;	// NSR Ca++ conc.
  double HTRPNCa; // fraction Ca++ bound high affinity troponin sites
  double LTRPNCa; // fraction Ca++ bound low affinity troponin sites

  double C0Kv43;   // Kv4.3 channel state C1
  double C1Kv43;   // Kv4.3 channel state C2
  double C2Kv43;   // Kv4.3 channel state C3
  double C3Kv43;   // Kv4.3 channel state C4
  double OKv43;    // Kv4.3 channel open state
  double CI0Kv43;  // Kv4.3 channel state I1
  double CI1Kv43;  // Kv4.3 channel state I2
  double CI2Kv43;  // Kv4.3 channel state I3
  double CI3Kv43;  // Kv4.3 channel state I4
  double OIKv43;   // Kv4.3 channel state I5
  double C0Kv14;   // Kv1.4 channel state C1
  double C1Kv14;   // Kv1.4 channel state C2
  double C2Kv14;   // Kv1.4 channel state C3
  double C3Kv14;   // Kv1.4 channel state C4
  double OKv14;    // Kv1.4 channel open state
  double CI0Kv14;  // Kv1.4 channel state I1
  double CI1Kv14;  // Kv1.4 channel state I2
  double CI2Kv14;  // Kv1.4 channel state I3
  double CI3Kv14;  // Kv1.4 channel state I4
  double OIKv14;   // Kv1.4 channel state I5
  //double CaTOT;	  // in fmoles
  double C1Herg;   //HERG channel state C1 (closed)
  double C2Herg;   //HERG channel state C2
  double C3Herg;   //HERG channel state C3
  double OHerg;    //HERG channel state O (Open)
  double IHerg;   //HERG channel state I1(Inactivated)

  //-------LOCAL STATE DERIVATIVE ARRAYS-----------------------------------------

  double dV;

  double dmNa;
  double dhNa;
  double djNa;
  double dxKs;
  double dNai;
  double dKi;
  double dCai;
  double dCaSL;
  double dCaNSR;
  double dHTRPNCa;
  double dLTRPNCa;
  double dC0Kv43;
  double dC1Kv43;
  double dC2Kv43;
  double dC3Kv43;
  double dOKv43;
  double dCI0Kv43;
  double dCI1Kv43;
  double dCI2Kv43;
  double dCI3Kv43;
  //double dOIKv43;
  double dC0Kv14;
  double dC1Kv14;
  double dC2Kv14;
  double dC3Kv14;
  double dOKv14;
  double dCI0Kv14;
  double dCI1Kv14;
  double dCI2Kv14;
  double dCI3Kv14;
  //double dOIKv14;
  double dCaTOT;
  double dC1Herg;
  double dC2Herg;
  double dC3Herg;
  double dOHerg;
  //double dIHerg;
  std::vector<double> dCa_spatial(NFRU,0);

  //-------LOCAL CURRENT AND FLUX ARRAYS-----------------------------------------

  double INa;	// Na+ current
  double IKr;    // rapid activating delayed rectifier K+ current
  double IKs;	// slow activating delayed rectifier K+ current
  double Ito1;	// transient outward K+ current
  double IK1;	// time independent K+ current
  double IKp;	// plateau K+ current
  std::vector<double> INaCa(NFRU,0);	// Na+/Ca++ exchanger current
  double INaK;	// Na+/K+ pump current
  double IpCa;	// sarcolemmal Ca++ pump current
  double ICab;	// Ca++ background current
  double INab;	// Na+ background current
  //	double ICa;	// L-type channel Ca++ current
  double Jup;	// SERCA pump flux
  double Jtrpn;	// troponin-Ca++ binding/unbinding flux
  //	double Jtr;	// Ca++ flux from NSR to JSR
  //	double Jxfer;	// Ca++ flux from SS to myoplasm
  double IKv43;     // Kv4.3 channel current
  double IKv14_K;   // Kv1.4 channel K+ current
  double IKv14_Na;  // Kv1.4 channel Na+ current
  double IKv14;     // Kv1.4 channel total current
  double Istim;	// Applied current stimulus
  //	double Ito2;	// transient outward Cl- current
  double Itot;   // Total membrane current

  //-------REVERSAL POTENTIALS---------------------------------------------------

  double ENa;	// reversal potential for INa,INab
  double EK;	// reversal potential for IKr,IKv43,IK1,IKp
  double EKs;	// reversal potential for IKs
  double ECa;	// reversal potential for ICab

  //-------MEMBRANE CURRENT PARAMETERS-------------------------------------------

  /*	peak IKr  conductance (mS/uF) */
  //const double GKr=0.024 // 0.023 in stable C++
  const double GKr = 0.029; // by Reza

  /*	peak IKs  conductance (mS/uF) */
  // const double GKs=0.0027134 // incorrect
  const double GKs = 0.035; //0.015702;  //TODO decreased by mark (Volders); // 0.035	G_Ks in Canine C++ code

  /*	peak IK1  conductance (mS/uF) */
  //UpK	const double GK1=2.8;
  const double GK1 = 3.0; // 2.73 by Reza

  /*  peak IKp  conductance (mS/uF) */
  const double GKp = 1.2 * 0.002216; // 0.001 in stable C++
  // const double GKp=0.001; // 0.001 in stable C++

  /*  peak INa  conductance (mS/uF) */
  const double GNa = 12.8;

  /*  scaling factor of Na+/Ca++ exchange (uA/uF) */
  const double kNaCa = 0.9 * 0.30;

  const double KmNa = 87.5; // Na+  half sat. constant for Na+/Ca++ exch. (mM)
  const double KmCa = 1.38; // Ca++ half sat. constant for Na+/Ca++ exch. (mM)
  const double KmK1 = 13.0; // Ca++ half sat. constant for IK1 (mM)
  const double ksat = 0.2; // Na+/Ca++ exch. sat. factor at negative potentials
  const double eta = 0.35; // controls voltage dependence of Na+/Ca++ exch.
  const double INaKmax = 1.3 * 0.693; // maximum Na+/K+ pump current (uA/uF)
  // const double INaKmax=1.0*0.693; // in Canine C++
  const double KmNai = 10.0; // Na+  half sat. constant for Na+/K+ pump (mM)
  const double KmKo = 1.5; // K+   half sat. constant for Na+/K+ pump (mM)
  const double IpCamax = 0.6 * 0.05; // maximum sarcolemmal Ca++ pump current (uA/uF)
  const double KmpCa = 0.0005; // half sat. constant for sarcolemmal Ca++ pump (mM)

  // max. background Ca++ current conductance (mS/uF)
  const double GCab = 3.3 * 7.684e-5; // 3e-5 in stable C++
  // const double GCab=20e-5; // 3e-5 in stable C++

  // max. background Na+  current conductance (mS/uF)
  //const double GNab=0.85*0.0031; // 0.000395 in stable C++
  const double GNab = 0.001; // 0.000395 in stable C++

  //-------SR PARAMETERS---------------------------------------------------------

  double Kfb = 0.26e-3; // foward half sat. constant for Ca++ ATPase (mM)
  const double Krb = 1.8; // backward half sat. constant for Ca++ ATPase (mM)
  const double KSR = 1.0; // scaling factor for Ca++ ATPase
  const double Nfb = 0.75; // foward cooperativity constant for Ca++ ATPase
  const double Nrb = 0.75; // backward cooperativity constant for Ca++ ATPase
  const double vmaxf = 1.53 * 137.0e-6; // Ca++ ATPase forward rate parameter (mM/ms)
  const double vmaxr = 1.53 * 137.0e-6; // Ca++ ATPase backward rate parameter (mM/ms)

  //-------Kv43 AND Kv14 CURRENT PARAMETERS FOR Ito1-----------------------------

  double KvScale = 1.13 * 1.03 * 1.55; // Scale factor for Kv4.3 and Kv1.4 currents

  if (HF_flag) {
    KvScale *= HF_Kv;
  }

  const double Kv43Frac = 0.77; // Fraction of Ito1 which is Kv4.3 current
  const double GKv43 = Kv43Frac * KvScale * 0.1; // Maximum conductance of Kv4.3 channel (mS/uF)
  const double PKv14 = (1.0 - Kv43Frac) * KvScale * 4.792933e-7; // Permeability of Kv1.4 channel (cm/s)

  //	NOTE: For step from -90mV to 40mV
  //	GKv43=0.1 and PKv14=4.792933e-7 ---> same peak current

  // alphaa0Kv43,aaKv43,betaa0Kv43,baKv43 // Voltage dependent
  // alphai0Kv43,aiKv43,betai0Kv43,biKv43 // rate parameters
  const double alphaa0Kv43 = 0.543708;
  const double aaKv43 = 0.028983;
  const double betaa0Kv43 = 0.080185;
  const double baKv43 = 0.0468437;
  const double alphai0Kv43 = 0.0498424;
  const double aiKv43 = 0.000373016;
  const double betai0Kv43 = 0.000819482;
  const double biKv43 = 0.00000005374;
  const double alphaa0Kv14 = 6.0 * 0.31551;
  const double aaKv14 = 0.00695;
  const double betaa0Kv14 = 6.0 * 0.001966;
  const double baKv14 = 0.08527;
  const double alphai0Kv14 = 6.0 * 0.0004938;
  const double betai0Kv14 = 6.0 * 0.0000176;

  // f1Kv43, f2Kv43, f3Kv43, f4Kv43  // Rate scaling factors
  const double f1Kv43 = 1.8936;
  const double f2Kv43 = 14.224647456;
  const double f3Kv43 = 158.574378389;
  const double f4Kv43 = 142.936645351;
  const double b1Kv43 = 6.77348;
  const double b2Kv43 = 15.6212705152;
  const double b3Kv43 = 28.7532603313;
  const double b4Kv43 = 524.576206679;
  const double f1Kv14 = 0.20005;
  const double f2Kv14 = 0.320280;
  const double f3Kv14 = 13.50909223;
  const double f4Kv14 = 1151.7651385;
  const double b1Kv14 = 2.2300;
  const double b2Kv14 = 12.000299;
  const double b3Kv14 = 5.3701338025;
  const double b4Kv14 = 5.2396395511;

  double alpha_act43, beta_act43;
  double alpha_inact43, beta_inact43;
  double C0Kv43_to_C1Kv43, C1Kv43_to_C2Kv43;
  double C2Kv43_to_C3Kv43, C3Kv43_to_OKv43;
  double CI0Kv43_to_CI1Kv43, CI1Kv43_to_CI2Kv43;
  double CI2Kv43_to_CI3Kv43, CI3Kv43_to_OIKv43;
  double C1Kv43_to_C0Kv43, C2Kv43_to_C1Kv43;
  double C3Kv43_to_C2Kv43, OKv43_to_C3Kv43;
  double CI1Kv43_to_CI0Kv43, CI2Kv43_to_CI1Kv43;
  double CI3Kv43_to_CI2Kv43, OIKv43_to_CI3Kv43;
  double C0Kv43_to_CI0Kv43, C1Kv43_to_CI1Kv43, C2Kv43_to_CI2Kv43;
  double C3Kv43_to_CI3Kv43, OKv43_to_OIKv43;
  double CI0Kv43_to_C0Kv43, CI1Kv43_to_C1Kv43, CI2Kv43_to_C2Kv43;
  double CI3Kv43_to_C3Kv43, OIKv43_to_OKv43;
  double alpha_act14, beta_act14;
  double alpha_inact14, beta_inact14;
  double C0Kv14_to_C1Kv14, C1Kv14_to_C2Kv14;
  double C2Kv14_to_C3Kv14, C3Kv14_to_OKv14;
  double CI0Kv14_to_CI1Kv14, CI1Kv14_to_CI2Kv14;
  double CI2Kv14_to_CI3Kv14, CI3Kv14_to_OIKv14;
  double C1Kv14_to_C0Kv14, C2Kv14_to_C1Kv14;
  double C3Kv14_to_C2Kv14, OKv14_to_C3Kv14;
  double CI1Kv14_to_CI0Kv14, CI2Kv14_to_CI1Kv14;
  double CI3Kv14_to_CI2Kv14, OIKv14_to_CI3Kv14;
  double C0Kv14_to_CI0Kv14, C1Kv14_to_CI1Kv14, C2Kv14_to_CI2Kv14;
  double C3Kv14_to_CI3Kv14, OKv14_to_OIKv14;
  double CI0Kv14_to_C0Kv14, CI1Kv14_to_C1Kv14, CI2Kv14_to_C2Kv14;
  double CI3Kv14_to_C3Kv14, OIKv14_to_OKv14;

  //-------RM HERG(+MiRP1) current parameters for IKr (10/00), H designates HERG

  double C1H_to_C2H, C2H_to_C1H;
  //	double C2H_to_C3H, C3H_to_C2H;  //Voltage independent rates
  double C3H_to_OH,  OH_to_C3H;
  double C3H_to_IH, IH_to_C3H;
  double OH_to_IH,  IH_to_OH;

  const double T_Const_HERG = 5.320000001; //Temp constant from 23 to 37C,with
  //Q10 of 3.3
  const double A0_HERG = 0.017147641733086;
  const double B0_HERG = 0.03304608038835;
  const double A1_HERG = 0.03969328381141;
  const double B1_HERG = -0.04306054163980;
  const double A2_HERG = 0.02057448605977;
  const double B2_HERG = 0.02617412715118;
  const double A3_HERG = 0.00134366604423;
  const double B3_HERG = -0.02691385498399;
  const double A4_HERG = 0.10666316491288;
  const double B4_HERG = 0.00568908859717;
  const double A5_HERG = 0.00646393910049;
  const double B5_HERG = -0.04536642959543;
  const double A6_HERG = 0.00008039374403;
  const double B6_HERG = 0.00000069808924;
  const double C2H_to_C3H = T_Const_HERG * 0.02608362043337;
  const double C3H_to_C2H = T_Const_HERG * 0.14832978132145;

  //-------MEMBRANE CURRENT AND PUMP ENTITIES------------------------------------

  double VF_over_RT, VFsq_over_RT;
  double fKo;
  double K1_inf, KpV;
  double fNaK, sigma;
  double alpha_m, beta_m, alpha_h, beta_h, alpha_j, beta_j;
  double xKs_inf, tau_xKs;
  double fb, rb, beta_i;
  double a1, a2, a3; //, a4, a5;
  double exp_VFRT, exp_etaVFRT;

  //printf("%g: Jxfer=%.20g ICa=%.20g Jtr=%.20g Ito2=%.20g\n",time,Jxfer,ICa,Jtr,Ito2);


  //-------MOVE STATES TO LOCAL VARIABLE ARRAYS----------------------------------

  //   STATE#   STATE NAME
  mNa = state[index_mNa];                 //       2      m
  hNa = state[index_hNa];                 //       3      h
  jNa = state[index_jNa];                 //       4      j

  Nai = state[index_Nai];			//	 5	Nai
  Ki = state[index_Ki];		//	 6	Ki
  Cai = state[index_Cai];			//	 7	Cai
  CaSL = state[index_CaSL];			//	 7	Cai
  CaNSR = state[index_CaNSR];		//	 8	CaNSR
  xKs = state[index_xKs];			//	10	xKs
  LTRPNCa = state[index_LTRPNCa];		// 	11	LTRPNCa
  HTRPNCa = state[index_HTRPNCa];		// 	12	HTRPNCa

  C0Kv43 = state[index_C0Kv43];		// 	13	Kv4.3 state C1
  C1Kv43 = state[index_C1Kv43];		// 	14	Kv4.3 state C2
  C2Kv43 = state[index_C2Kv43];		// 	15	Kv4.3 state C3
  C3Kv43 = state[index_C3Kv43];		// 	16	Kv4.3 state C4
  OKv43 = state[index_OKv43];		// 	17	Kv4.3 state O (open];
  CI0Kv43 = state[index_CI0Kv43];		// 	18	Kv4.3 state I1
  CI1Kv43 = state[index_CI1Kv43];		// 	19	Kv4.3 state I2
  CI2Kv43 = state[index_CI2Kv43];		// 	20	Kv4.3 state I3
  CI3Kv43 = state[index_CI3Kv43];		// 	21	Kv4.3 state I4
  //OIKv43=state[index_OIKv43];		// 	22	Kv4.3 state I5
  OIKv43 = 1.0 - CI0Kv43 - CI1Kv43 - CI2Kv43 - CI3Kv43 - OKv43 - C0Kv43 - C1Kv43 - C2Kv43 - C3Kv43;	// state[index_OIKv43)		// 	32	Kv1.4 state I5

  C0Kv14 = state[index_C0Kv14];		// 	23	Kv1.4 state C1
  C1Kv14 = state[index_C1Kv14];		// 	24	Kv1.4 state C2
  C2Kv14 = state[index_C2Kv14];		// 	25	Kv1.4 state C3
  C3Kv14 = state[index_C3Kv14];		// 	26	Kv1.4 state C4
  OKv14 = state[index_OKv14];		// 	27	Kv1.4 state O (open];
  CI0Kv14 = state[index_CI0Kv14];		// 	28	Kv1.4 state I1
  CI1Kv14 = state[index_CI1Kv14];		// 	29	Kv1.4 state I2
  CI2Kv14 = state[index_CI2Kv14];		// 	30	Kv1.4 state I3
  CI3Kv14 = state[index_CI3Kv14];		// 	31	Kv1.4 state I4
  //OIKv14=state[index_OIKv14];		// 	32	Kv1.4 state I5
  OIKv14 = 1.0 - CI0Kv14 - CI1Kv14 - CI2Kv14 - CI3Kv14 - OKv14 - C0Kv14 - C1Kv14 - C2Kv14 - C3Kv14;


  C1Herg = state[index_C1Herg];		//	34	HERG state C1
  C2Herg = state[index_C2Herg];		//	35	HERG state C2
  C3Herg = state[index_C3Herg];		//	36	HERG state C3
  OHerg = state[index_OHerg];		//	37	HERG state O1
  //IHerg=state[index_IHerg];		//	37	HERG state I1
  IHerg = 1.0 - OHerg - C1Herg - C2Herg - C3Herg;

  //------------------------- Membrane potential V -----------------

  V = state[index_V];	//	 1	V

  //Apply external current to the stimulus
  Istim = Iext / (Acap * 1e6);

  //------SET THE MEMBRANE POTENTIAL IF IN VOLTAGE CLAMP MODE--------------------

  //-------COMPUTE REVERSAL POTENTIALS-------------------------------------------

  ENa = RT_over_F * log(Nao / Nai);
  EK =  RT_over_F * log(Ko / Ki);
  a1 = Ko + 0.01833 * Nao;
  a2 = Ki + 0.01833 * Nai;
  EKs = RT_over_F * log(a1 / a2);

  //-------COMPUTE INa, IKr, IKs, Ito1, IK1, INab, IKp---------------------------

  //Saucerman IKs model
  /*double f_Ksavail = 0.8;
  double Xs05 = 2.95; //-35.6 for +ISO
  double g_Ks0 = 0.5; // mS/uF
  double p_Ca = -log(Cai)/log(10) + 3.0;
  double g_Ks = g_Ks0 * f_Ksavail * (0.057 + 0.19) / (1.0 + exp((-7.2+p_Ca)/0.6));
  a1 = 0.0000719*(V+30.0)/(1.0-exp(-0.148*(V+30)));
  a2 = 0.0001310*(V+30.0)/(exp(0.0687*(V+30))-1.0);
  double Xs_infinity = 1.0/(1.0+exp(-(V-Xs05)/16.7));
  tau_xKs = 1.0 / (a1 + a2);
  //printf("tau_xKs = %g\n",tau_xKs);
  dxKs = (Xs_infinity-xKs)/tau_xKs;*/

  //Hund-Rudy 2004 IKs
  /*double gksbar  = 0.1 * 0.0575*(0.433*(1+.6/(1+pow((0.000038/(Cai)),1.4))));  //TODO scaled
  double xsss    = 1/(1+exp(-(V-10.5)/24.7));
  double xs1tau  = 1/(0.0000761*(V+44.6)/(1-exp(-9.97*(V+44.6))) + 0.00036 *(V- 0.55)/(exp(0.128*(V - 0.55))-1));
  double xs2tau  = 2*xs1tau;
  dxKs = (xsss - xKs)/xs1tau;
  double xKs2 = state[index_xKs2];*/
  double dxKs2 = 0; //(xsss - xKs2)/xs2tau;

  // Make sure thatn mNa never goes negative
  if (mNa <= 0.0) mNa = 1.e-15;

  INa = GNa * (mNa * mNa * mNa) * hNa * jNa * (V - ENa);

  fKo = pow(Ko / 4.0, 0.5);
  IKr = GKr * fKo * OHerg * (V - EK); 
  IKs = GKs * (xKs * xKs) * (V - EKs);
  //IKs = g_Ks * (xKs * xKs) * (V - EKs);
  //IKs = gksbar * (xKs * xKs2) * (V - EKs);

  /*double GKur = 3.62 * 3.48e-3;
  double IKur = GKur*(1.0/(1.0+exp((36.0-V)/17.0)))*(V-EK);
  IKr += IKur;*/

  if (beta_flag) {
    IKs *= beta_iks_scale; //3.172248;
  }
  IKr *= ikr_scale;

  IKv43 = GKv43 * OKv43 * (V - EK);

  VF_over_RT = V / RT_over_F;
  VFsq_over_RT = (1000.0 * Faraday) * VF_over_RT;

  // Idea here is to reduce number of calls to exp
  exp_VFRT = exp(VF_over_RT);

  // a1 =  Ki*exp(VF_over_RT)-Ko
  //     a2 = exp(VF_over_RT)-1.0

  a1 =  Ki * exp_VFRT - Ko;
  a2 = exp_VFRT - 1.0;

  if (fabs(V) < 1.e-6) {
    IKv14_K = PKv14 * OKv14 * 1000.0 * Faraday * (Ki - Ko);
  } else {
    IKv14_K = PKv14 * OKv14 * VFsq_over_RT * (a1 / a2);
  }

  // a1 =  Nai*exp(VF_over_RT)-Nao
  a1 =  Nai * exp_VFRT - Nao;

  if (fabs(V) < 1.e-6) {
    IKv14_Na = 0.02 * PKv14 * OKv14 * 1000.0 * Faraday * (Nai - Nao);
  } else {
    IKv14_Na = 0.02 * PKv14 * OKv14 * VFsq_over_RT * (a1 / a2);
  }

  IKv14 = IKv14_K + IKv14_Na;

  Ito1 = ito1_scale * (IKv43 + IKv14);
  // Original IK1 implementation
  // K1_inf = 1.0/(2.0+exp(1.5/RT_over_F*(V-EK)))

  // New formulation endorsed by R. Mazhari. However, we are using GK1=3.0, not 2.73
  K1_inf = 1.0 / (0.94 + exp(1.76 / RT_over_F * (V - EK)));

  IK1 = ik1_scale * GK1 * Ko / (Ko + KmK1) * K1_inf * (V - EK);

  //Shannon IK1
  /*double gk1bar = 0.9*pow(Ko/5.4,0.5);
  double ak1 = 1.02/(1+exp(0.2385*(V-EK-59.215)));
  double bk1 = (0.49124*exp(0.08032*(V-EK+5.476)) + exp(0.06175*(V-EK-594.31)))/(1.0+exp(-0.5143*(V-EK+4.753)));
  double K1inf = ak1/(ak1+bk1);
  IK1 = gk1bar*K1inf*(V-EK);*/

  INab = GNab * (V - ENa); // orig

  KpV = 1.0 / (1.0 + exp((7.488 - V) / 5.98));
  IKp = GKp * KpV * (V - EK);

  //-------COMPUTE INaK, INaCa, ICab, IpCa --------------------------------------

  VF_over_RT = V / RT_over_F;

  sigma = (exp(Nao / 67.3) - 1.0) / 7.0;
  a1 = 1.0 + 0.1245 * exp(-0.1 * VF_over_RT);
  a2 = 0.0365 * sigma / exp_VFRT;
  // a2 = 0.0365*sigma*exp(-VF_over_RT)
  fNaK = 1.0 / (a1 + a2);
  a1 = Ko / (Ko + KmKo);
  a2 = 1.0 + pow(KmNai / Nai, 1.5);
  INaK = nak_scale * INaKmax * fNaK * (a1 / a2);

  exp_etaVFRT = exp(eta * VF_over_RT);

  //...........INa
  if(fabs(V + 47.13) <= 1.e-4) {
    alpha_m = 0.32 / (0.1 - 0.005 * (V + 47.13));
    //Taylor Approximation near 47.13mV
  } else {
    a1 = 0.32 * (V + 47.13);
    alpha_m = a1 / (1.0 - exp(-0.1 * (V + 47.13)));
  }

  beta_m = 0.08 * exp(-V / 11.0);

  if (V >= -90.0) {	// Tau_m < 0.0035ms when V < -90mV
    dmNa = alpha_m * (1.0 - mNa) - beta_m * mNa;
  } else {
    dmNa = 0.0;
    mNa = alpha_m / (alpha_m + beta_m); // steady state approx.
    state[index_mNa] = mNa;
  }

  // A slightly modified, continuous implementation
  alpha_h = 0.135 * exp((80 + V) / -6.8);

  if (V < -38.73809636838782) {
    // curves do not cross at -40mV, but close
    beta_h = 3.56 * exp(0.079 * V) + 310000 * exp(0.35 * V);
  } else {
    beta_h = 1 / (0.13 * (1 + exp((V + 10.66) / -11.1)));
  }

  if (V < -37.78) {
    a1 = -127140 * exp(0.2444 * V);
    a2 = 3.474E-5 * exp(-0.04391 * V);
    a3 = 1.0 + exp(0.311 * (V + 79.23));
    alpha_j = (a1 - a2) * (V + 37.78) / a3;
  } else {
    alpha_j = 0;
  }

  if (V < -39.82600037702883) {
    // curves do not cross at -40mV, but close
    beta_j = 0.1212 * exp(-0.01052 * V) / (1.0 + exp(-0.1378 * (V + 40.14)));
  } else {
    beta_j = 0.3 * exp(-2.535E-7 * V) / (1 + exp(-0.1 * (V + 32)));
  }

  dhNa = alpha_h * (1.0 - hNa) - beta_h * hNa;
  djNa = alpha_j * (1.0 - jNa) - beta_j * jNa;

  C1H_to_C2H = T_Const_HERG * A0_HERG * exp(B0_HERG * V);
  C2H_to_C1H = T_Const_HERG * A1_HERG * exp(B1_HERG * V);
  C3H_to_OH =  T_Const_HERG * A2_HERG * exp(B2_HERG * V);
  OH_to_C3H =  T_Const_HERG * A3_HERG * exp(B3_HERG * V);
  OH_to_IH =   T_Const_HERG * A4_HERG * exp(B4_HERG * V);
  IH_to_OH =   T_Const_HERG * A5_HERG * exp(B5_HERG * V);
  C3H_to_IH =  T_Const_HERG * A6_HERG * exp(B6_HERG * V);
  IH_to_C3H =  (OH_to_C3H * IH_to_OH * C3H_to_IH) / (C3H_to_OH * OH_to_IH);

  if (beta_flag) {
    OH_to_IH /= beta_ikr_rate_scale;
    C3H_to_IH /= beta_ikr_rate_scale;
    IH_to_OH *= beta_ikr_rate_scale;
    IH_to_C3H *= beta_ikr_rate_scale;
  }

  dC1Herg = C2H_to_C1H * C2Herg - C1H_to_C2H * C1Herg;
  a1 = C1H_to_C2H * C1Herg + C3H_to_C2H * C3Herg;
  a2 = (C2H_to_C1H + C2H_to_C3H) * C2Herg;
  dC2Herg = a1 - a2;
  a1 = C2H_to_C3H * C2Herg + OH_to_C3H * OHerg + IH_to_C3H * IHerg;
  a2 = (C3H_to_IH + C3H_to_OH + C3H_to_C2H) * C3Herg;
  dC3Herg = a1 - a2;
  a1 = C3H_to_OH * C3Herg + IH_to_OH * IHerg;
  a2 = (OH_to_C3H + OH_to_IH) * OHerg;
  dOHerg = a1 - a2;
  // IHerg "removed"
  // a1 = C3H_to_IH * C3Herg + OH_to_IH * OHerg
  // a2 = (IH_to_C3H + IH_to_OH) * IHerg
  // dIHerg = a1-a2


  //Mark's version (Based on Volders Circulation 2003)
  /*if (beta_flag) {
    xKs_inf = 1.0 / (1.0 + exp(-(V + 10.0) / 13.6));
  } else {
    xKs_inf = 1.0 / (1.0 + exp(-(V - 10.0) / 13.6));
  }
  
  if (beta_flag) {
    if (fabs(V) < 1.e-6) { // First order Taylor expansion
      a1 = (1.13398e-5) / 0.005;
      a2 = (1.31001e-5) / 0.0687;
    } else {
      a1 = 1.13398e-5 * (V) / (1.0 - exp(-0.005 * (V)));
      a2 = 1.31001e-5 * (V) / (exp(0.0687* (V)) - 1.0);
    }
  } else {
    if (fabs(V) < 1.e-6) { // First order Taylor expansion
      a1 = (3.60614e-5) / 0.148;
      a2 = (1.2291e-5) / 0.0687;
    } else {
      a1 = 3.60614e-5 * (V) / (1.0 - exp(-0.148 * (V)));
      a2 = 1.2291e-5 * (V) / (exp(0.0687 * (V)) - 1.0);
    }
    }*/
 
  double IKs_V_shift;
  if (beta_flag)  {
    IKs_V_shift = -35; //TODO
    IKs *= 1.4;
  } else {
    IKs_V_shift = 0;
  }
  //Original model
  xKs_inf = 1.0/(1.0+exp(-(V-(24.7+IKs_V_shift))/13.6));
  if (fabs(V-10)<1.e-6) { // First order Taylor expansion                                                                                                                                             
    a1 = (7.19e-5)/0.148;
    a2 = (1.31e-4)/0.0687;
  } else {
    a1 = 7.19e-5*(V-10.0-IKs_V_shift)/(1.0-exp(-0.148*(V-10.0-IKs_V_shift)));
    a2 = 1.31e-4*(V-10.0-IKs_V_shift)/(exp(0.0687*(V-10.0-IKs_V_shift))-1.0);
  }
    
  tau_xKs = 1.0 / (a1 + a2);
  dxKs = (xKs_inf - xKs) / tau_xKs;

  //-------COMPUTE CONCENTRATION AND VOLTAGE DERIVATIVES-------------------------


  if (Cai < 0) {
    Cai = 1e-6;
    //printf("Warning! Cai clamped to 0.\n");
  }
  
  //SERCA
  fb = pow(Cai / (serca_kf_scale * Kfb), Nfb);
  rb = pow(CaNSR / Krb, Nrb);
  Jup = KSR * (vmaxf * fb - vmaxr * rb) / (1.0 + fb + rb);
  
/*
  //Tran-Crampin SERCA Formulation
  const double K_D_SR_SERCA = 2240.0/1000.0; //mM
  double K_D_I_SERCA = 910.0/1000.0; //mM
  K_D_I_SERCA *= serca_kf_scale;
  const double A_P_SERCA = 150.0/1000.0; //mM
  double K_i = Cai / K_D_I_SERCA;
  K_i = K_i*K_i;
  double K_sr = CaNSR / K_D_SR_SERCA;
  K_sr = K_sr*K_sr;
  double D_cycle = 0.104217 + 17.923*K_sr + K_i*(1.75583e6 + K_sr*7.61673e6) + K_i*K_i*(6.08463e11 + K_sr*4.50544e11);
  double v_cycle = (K_i*K_i*3.24873e12 + K_i*(9.17846e6 - 11478.2*K_sr) - 0.329904*K_sr) / D_cycle;
  Jup = (2e-3) * v_cycle * A_P_SERCA; //mM / ms (Note it is given in uM/s in Williams et al 2011)
  */
  
  /*
  //Shannon SERCA Model
  double Q10SRCaP = 2.6;
  double Vmax_SRCaP = 2.86e-4;
  double Kmf2 = 0.246e-3;
  double Kmr2 = 1.7;
  double hillSRCaP = 1.787;
  fb = pow(Cai/(serca_kf_scale*Kmf2),hillSRCaP);
  rb = pow(CaNSR/Kmr2,hillSRCaP);
  Jup = Q10SRCaP*Vmax_SRCaP*(fb-rb)/(1.0+fb+rb);
  */
  
  Jup *= serca_vmax_scale;
  if (HF_flag) {
    Jup *= HF_SERCA;
  } 

  //Shannon NCX Model
  /*
  double Kda = 0.256e-3; //Shannon
  //double Kda = 0.11e-3; //Restrepo
  double Ka = 1.0/(1+pow(Kda/CaSL,3.0));
  double Vmax = 9.0; //A/F
  double QNCX = 1.57;
  double HNa = 3;
  double eta_ncx = 0.35;
  double ksat_ncx = 0.27;
  double KmCai = 3.59e-3;
  double KmCao = 1.3;
  double KmNai2 = 12.29;
  double KmNao = 87.5;
  INaCa = Ka*QNCX*Vmax*(exp(eta_ncx*VF_over_RT)*pow(Nai,HNa)*Cao - exp((eta_ncx-1)*VF_over_RT)*pow(Nao,HNa)*CaSL)
  	/((KmCai*pow(Nao,HNa)*(1+pow(Nai/KmNai2,HNa)) +
  		 pow(KmNao,HNa)*CaSL*(1+(CaSL/KmCai)) +
  		 KmCao*pow(Nai,HNa) +
  		 pow(Nai,HNa)*Cao +
  		 pow(Nao,HNa)*CaSL)*(1+ksat_ncx*exp((eta_ncx-1)*VF_over_RT)));
  */

  ECa = 0.5 * RT_over_F * log(Cao / Cai); // log=ln vai ei ??

  ICab = GCab * (V - ECa);
  IpCa = IpCamax * Cai / (KmpCa + Cai);

  if (beta_flag) {
    a1 = beta_ltrpn_koff_scale * kltrpn_minus * LTRPNCa;
  } else {
    a1 = kltrpn_minus * LTRPNCa;
  }

  dLTRPNCa = kltrpn_plus * Cai * (1.0 - LTRPNCa) - a1;

  a1 = khtrpn_minus * HTRPNCa;
  dHTRPNCa = khtrpn_plus * Cai * (1.0 - HTRPNCa) - a1;

  Jtrpn = LTRPNtot * dLTRPNCa + HTRPNtot * dHTRPNCa;

  double a3_ncx = 1.0 + ksat * exp_etaVFRT / exp_VFRT;
  double a4_ncx = KmCa + Cao;
  double a5_ncx = 5000.0 / (KmNa * KmNa * KmNa + Nao * Nao * Nao);
  double a1_ncx = exp_etaVFRT * (Nai * Nai * Nai) * Cao;

  double JSL = 0;
  double Jxfer_tot = 0;
  double INaCa_tot = 0;
  a1 = (Acap/NFRU_total) / (VSL * Faraday * 1000.0);
// #pragma omp parallel for default(shared) reduction(+:JSL,Jxfer_tot,INaCa_tot)
  for (int i = 0; i < NFRU; i++) {
    
    double ci = state[N_states+i];

    if (ci < 1e-9) { 
      //printf("Warning: CaSL < 1e-9, clamping!\n");
      state[N_states+i] = 1e-9;
      ci = 1e-9;
    }

    double a2_ncx = exp_etaVFRT / exp_VFRT * (Nao * Nao * Nao) * ci;
    INaCa[i] = ncx_scale * frac_ncx_sl * kNaCa * a5_ncx * (a1_ncx - a2_ncx) / (a4_ncx * a3_ncx);

    //double OScale = FRUs[i].isOrphan() ? ORPHAN_SCALE : 1.0;

    double jxfer_i = (FRUs[i].Get_CaSS_Avg() - ci) / tauxfer;

    //double jsl2cyto = ((Cai - ci) + (6-FRU_neighbs[i].size())*(Cai-ci)) / tauSL2cyto;
    double jsl2cyto = (Cai - ci) / tauSL2cyto;

    double Jdiff = 0;
    for (unsigned j = 0; j < FRU_neighbs[i].size(); j++) {
      double taugrid_i = taugrid*FRU_neighbs_distance[i][j];
      Jdiff += (state[N_states+FRU_neighbs[i][j]] - ci)/taugrid_i;
    }

    double a1_i = CMDNtot * KmCMDN / (pow(ci + KmCMDN, 2.0)); // + (165*1.1 / ((ci+1.1)*(ci+1.1)));
    //double a2_i = EGTAtot * KmEGTA / (pow(ci + KmEGTA, 2.0));
    double a3_i = BSLtot_SL * KBSL / (pow(ci + KBSL, 2.0));
    double beta_i_i = 1.0 / (1.0 + a1_i + a3_i);
    
    Fstate[N_states+i] = beta_i_i*(jxfer_i*(VSS/VSL) + jsl2cyto + Jdiff + INaCa[i]*a1); //Assign Fstate within omp loop for speed
    
    INaCa_tot += INaCa[i];
    Jxfer_tot += jxfer_i;
    JSL += jsl2cyto;
  }
  INaCa_tot /= NFRU;
  JSL *= VSL * NFRU_scale;
  Jxfer_tot *= NFRU_scale;

  a1 = exp_etaVFRT * (Nai * Nai * Nai) * Cao;
  a2 = exp_etaVFRT / exp_VFRT * (Nao * Nao * Nao) * Cai;
  a3 = 1.0 + ksat * exp_etaVFRT / exp_VFRT;
  double a4 = KmCa + Cao;
  double a5 = 5000.0 / (KmNa * KmNa * KmNa + Nao * Nao * Nao);
  double INaCa_cyto = ncx_scale * (1.0-frac_ncx_sl) * kNaCa * a5 * (a1 - a2) / (a4 * a3);

  INaCa_tot += INaCa_cyto;

  a1 = CMDNtot * KmCMDN / (pow(Cai + KmCMDN, 2.0));
  a2 = EGTAtot * KmEGTA / (pow(Cai + KmEGTA, 2.0));
  beta_i = 1.0 / (1.0 + a1 + a2);

  a1 = Acap / (Vmyo * Faraday * 1000.0);
  dCai = beta_i * ( -JSL / Vmyo - Jup - Jtrpn - (ICab + IpCa - JCa - 2.0*INaCa_cyto) * 0.5 * a1);

  a1 = CMDNtot * KmCMDN / (pow(CaSL + KmCMDN, 2.0));
  a2 = EGTAtot * KmEGTA / (pow(CaSL + KmEGTA, 2.0));
  a3 = BSLtot_SL * KBSL / (pow(CaSL + KBSL, 2.0));
  beta_i = 1.0 / (1.0 + a1 + a2 + a3);

  dCaSL = 0;

  if (!bEjectJup) { 
    dCaNSR = -Jtr * VJSR / VNSR + Jup * Vmyo / VNSR; 
  } else {
    dCaNSR = 0; //-Jtr * VJSR / VNSR; 
  }
  if (bClampNSR) {
    dCaNSR = 0;
  }

  a1 = Acap / (Faraday * 1000.0);
  dCaTOT =  -((ICab - 2.0 * INaCa_tot + IpCa) + ICa + JCa) * 0.5 * a1 * 1.0e6;

  //
  //Luo-Rudy 1994 non-selective Ca-dependent current (NOT IDENTIFIED IN CANINE)
  //

  /*
  const double R = 8314; // Universal Gas Constant (J/kmol*K) 
  const double frdy = 96485; // Faraday's Constant (C/mol) 
  const double temp = 310; // Temperature (K) 

  const double zna = 1; // Na valence 
  const double zk = 1; // K valence 

  const double pnsca = 0.000000175; // Permeability of channel to Na and K (cm/s) 
  const double kmnsca = 0.0012; // Half-saturation concentration of NSCa channel (mM) 
  const double gako = 0.75; // Activity coefficient of K 
  const double ganao = 0.75; // Activity coefficient of Na 
  const double gaki = 0.75; // Activity coefficient of K 
  const double ganai = 0.75; // Activity coefficient of Na 
  double ibarnsna = pnsca*zna*zna*((V*frdy*frdy)/(R*temp))
    *((ganai*Nai*exp((zna*V*frdy)/(R*temp))-ganao*Nao)/(exp((zna*V*frdy)/(R*temp))-1)); 
  double ibarnsk = pnsca*zk*zk*((V*frdy*frdy)/(R*temp))
  *((gaki*Ki*exp((zk*V*frdy)/(R*temp))-gako*Ko)/(exp((zk*V*frdy)/(R*temp))-1)); */

  double insna = 0; //ibarnsna/(1+pow(kmnsca/Cai,3));
  double insk = 0; //ibarnsk/(1+pow(kmnsca/Cai,3)); 

  //IK(ATP) as described in Ferrero 1996
  /*double gkatp = 3.8; //channels/um^2
  double ADPi = 15; //95; //15; //uM 
  double ATPi = 6.8; //4.3; //6.8;  //mM
  double Mgi = 0.5; //mM

  double Kmatp = 1e-3 * (35.8 + 17.9*pow(ADPi,0.256)); //mM
  double hatp = 1.3 + 0.74*exp(-0.09*ADPi);
  double fatp = 1.0/(1.0 + pow(ATPi/Kmatp,hatp));

  double gamma0 = 35.375*pow(Ko/5.4,0.24);
  double KhMg = (0.65/sqrt(Ko+5.0)) * exp(-2*0.32*VF_over_RT);
  double FM = 1.0/(1.0 + (Mgi/KhMg));
  double KhNa = 25.9*exp(-0.35*VF_over_RT);
  double FN = 1.0/(1.0 + pow(Nai/KhNa,2.0));
  double FT = 1.0; //No temperature effect
  double g0 = gamma0 * FM * FN * FT;
  //printf("%g %g %g; %g; %g %g; %g %g; %g\n", Kmatp, hatp, fatp, gamma0, KhMg, FM, KhNa, FN, g0);
  //printf("%g %g\n",KhNa,FN);

  double gkbaratp = 0.91 * gkatp*fatp*g0 * 1e-2; //pA/pF-mV 2 orders of magnitude too high?
  //printf("%g\n", gkbaratp);

  insk = gkbaratp*(V-EK); //TODO add new current for this
  */

  //IK(ATP) from Zhou 2009
  /*
  double ADPi = 15; //uM
  double ATPi = 6.8; //mM
  double sigma_atp = 0.6; //1/um^2
  double g0_atp = 30.95;
  double po_atp = 0.91;
  double HKATP = 0.0001; //Negative?
  double KMKATP = 0.56;
  double surf = 10e4; //um^2
  double KMATP = 35.8 + 17.9*pow(ADPi,KMKATP);
  double HATP = 1.3 + 0.74*exp(-HKATP*ADPi);
  double fATP = 1.0/(1.0+pow(ATPi/KMATP,HATP));

  insk = sigma_atp*g0_atp*po_atp*fATP*(V-EK);
  */
  

  /*
  //IK(Na) Na-dependent potassium channel (Not in canine, Saxena 1997)
  double gknabar = 0.12848e3; //uS/uF
  double PoV = 0.8 - (0.65/(1+exp((V+125)/15.0)));
  double PoNai = 0.85 / (1+pow(66.0/Nai,2.8));
  double insk = gknabar * PoV * PoNai * (Em - Ek);
  */

  if (clamp_nai == 0) {
    a1 = Acap / (Faraday * 1000.0 * Vmyo);
    dNai = -( INa + INab + 3.0 * (INaCa_tot + INaK) + IKv14_Na + insna) * a1;
  } else {
    dNai = 0;
  }

#if clamp_ki == 0
  a3 = IKr + IKs + IK1 + IKp + insk;
  dKi = -( a3 - 2.0 * INaK + (Ito1 - IKv14_Na) + Istim ) * a1;
#else
  dKi = 0;
#endif

  a1 = INa + ICa + IKr + IKs;
  a2 = IK1 + IKp + INaCa_tot + INaK + Ito1 + Ito2;
  a3 = IpCa + ICab + INab;

  Itot = a1 + a2 + a3 + Istim + JCa + insna +insk;

  if (VClamp_Flag) {
    dV = 0;
  } else {
    dV = -Itot;
  }

  //-------COMPUTE DERIVATIVES OF Kv4.3 CHANNEL STATES---------

  alpha_act43 = alphaa0Kv43 * exp(aaKv43 * V);
  beta_act43  = betaa0Kv43 * exp(-baKv43 * V);
  alpha_inact43 = alphai0Kv43 * exp(-aiKv43 * V);
  beta_inact43  = betai0Kv43 * exp(biKv43 * V);

  C0Kv43_to_C1Kv43 = 4.0 * alpha_act43;
  C1Kv43_to_C2Kv43 = 3.0 * alpha_act43;
  C2Kv43_to_C3Kv43 = 2.0 * alpha_act43;
  C3Kv43_to_OKv43  =      alpha_act43;

  CI0Kv43_to_CI1Kv43 = 4.0 * b1Kv43 * alpha_act43;
  CI1Kv43_to_CI2Kv43 = 3.0 * b2Kv43 * alpha_act43 / b1Kv43;
  CI2Kv43_to_CI3Kv43 = 2.0 * b3Kv43 * alpha_act43 / b2Kv43;
  CI3Kv43_to_OIKv43  =      b4Kv43 * alpha_act43 / b3Kv43;

  C1Kv43_to_C0Kv43 =      beta_act43;
  C2Kv43_to_C1Kv43 = 2.0 * beta_act43;
  C3Kv43_to_C2Kv43 = 3.0 * beta_act43;
  OKv43_to_C3Kv43  = 4.0 * beta_act43;

  CI1Kv43_to_CI0Kv43 =             beta_act43 / f1Kv43;
  CI2Kv43_to_CI1Kv43 = 2.0 * f1Kv43 * beta_act43 / f2Kv43;
  CI3Kv43_to_CI2Kv43 = 3.0 * f2Kv43 * beta_act43 / f3Kv43;
  OIKv43_to_CI3Kv43  = 4.0 * f3Kv43 * beta_act43 / f4Kv43;

  C0Kv43_to_CI0Kv43 = beta_inact43;
  C1Kv43_to_CI1Kv43 = f1Kv43 * beta_inact43;
  C2Kv43_to_CI2Kv43 = f2Kv43 * beta_inact43;
  C3Kv43_to_CI3Kv43 = f3Kv43 * beta_inact43;
  OKv43_to_OIKv43   = f4Kv43 * beta_inact43;

  CI0Kv43_to_C0Kv43 = alpha_inact43;
  CI1Kv43_to_C1Kv43 = alpha_inact43 / b1Kv43;
  CI2Kv43_to_C2Kv43 = alpha_inact43 / b2Kv43;
  CI3Kv43_to_C3Kv43 = alpha_inact43 / b3Kv43;
  OIKv43_to_OKv43   = alpha_inact43 / b4Kv43;

  a1 = (C0Kv43_to_C1Kv43 + C0Kv43_to_CI0Kv43) * C0Kv43;
  a2 = C1Kv43_to_C0Kv43 * C1Kv43 + CI0Kv43_to_C0Kv43 * CI0Kv43;
  dC0Kv43 = a2 - a1;

  a1 = (C1Kv43_to_C2Kv43 + C1Kv43_to_C0Kv43 + C1Kv43_to_CI1Kv43) * C1Kv43;
  a2 = C2Kv43_to_C1Kv43 * C2Kv43 + CI1Kv43_to_C1Kv43 * CI1Kv43 + C0Kv43_to_C1Kv43 * C0Kv43;
  dC1Kv43 = a2 - a1;

  a1 = (C2Kv43_to_C3Kv43 + C2Kv43_to_C1Kv43 + C2Kv43_to_CI2Kv43) * C2Kv43;
  a2 = C3Kv43_to_C2Kv43 * C3Kv43 + CI2Kv43_to_C2Kv43 * CI2Kv43 + C1Kv43_to_C2Kv43 * C1Kv43;
  dC2Kv43 = a2 - a1;

  a1 = (C3Kv43_to_OKv43 + C3Kv43_to_C2Kv43 + C3Kv43_to_CI3Kv43) * C3Kv43;
  a2 = OKv43_to_C3Kv43 * OKv43 + CI3Kv43_to_C3Kv43 * CI3Kv43 + C2Kv43_to_C3Kv43 * C2Kv43;
  dC3Kv43 = a2 - a1;

  a1 = (OKv43_to_C3Kv43 + OKv43_to_OIKv43) * OKv43;
  a2 = C3Kv43_to_OKv43 * C3Kv43 + OIKv43_to_OKv43 * OIKv43;
  dOKv43 = a2 - a1;

  a1 = (CI0Kv43_to_C0Kv43 + CI0Kv43_to_CI1Kv43) * CI0Kv43;
  a2 = C0Kv43_to_CI0Kv43 * C0Kv43 + CI1Kv43_to_CI0Kv43 * CI1Kv43;
  dCI0Kv43 = a2 - a1;

  a1 = (CI1Kv43_to_CI2Kv43 + CI1Kv43_to_C1Kv43 + CI1Kv43_to_CI0Kv43) * CI1Kv43;
  a2 = CI2Kv43_to_CI1Kv43 * CI2Kv43 + C1Kv43_to_CI1Kv43 * C1Kv43 + CI0Kv43_to_CI1Kv43 * CI0Kv43;
  dCI1Kv43 = a2 - a1;

  a1 = (CI2Kv43_to_CI3Kv43 + CI2Kv43_to_C2Kv43 + CI2Kv43_to_CI1Kv43) * CI2Kv43;
  a2 = CI3Kv43_to_CI2Kv43 * CI3Kv43 + C2Kv43_to_CI2Kv43 * C2Kv43 + CI1Kv43_to_CI2Kv43 * CI1Kv43;
  dCI2Kv43 = a2 - a1;

  a1 = (CI3Kv43_to_OIKv43 + CI3Kv43_to_C3Kv43 + CI3Kv43_to_CI2Kv43) * CI3Kv43;
  a2 = OIKv43_to_CI3Kv43 * OIKv43 + C3Kv43_to_CI3Kv43 * C3Kv43 + CI2Kv43_to_CI3Kv43 * CI2Kv43;
  dCI3Kv43 = a2 - a1;

  // OIKv43 is "removed"
  // a1 = (OIKv43_to_OKv43+OIKv43_to_CI3Kv43)*OIKv43;
  // a2 = OKv43_to_OIKv43*OKv43 + CI3Kv43_to_OIKv43*CI3Kv43;
  // dOIKv43 = a2 - a1;


  //-------COMPUTE DERIVATIVES OF Kv1.4 CHANNEL STATES---------

  alpha_act14 = alphaa0Kv14 * exp(aaKv14 * V);
  beta_act14  = betaa0Kv14 * exp(-baKv14 * V);
  alpha_inact14 = alphai0Kv14;
  beta_inact14  = betai0Kv14;

  C0Kv14_to_C1Kv14 = 4.0 * alpha_act14;
  C1Kv14_to_C2Kv14 = 3.0 * alpha_act14;
  C2Kv14_to_C3Kv14 = 2.0 * alpha_act14;
  C3Kv14_to_OKv14  =      alpha_act14;

  CI0Kv14_to_CI1Kv14 = 4.0 * b1Kv14 * alpha_act14;
  CI1Kv14_to_CI2Kv14 = 3.0 * b2Kv14 * alpha_act14 / b1Kv14;
  CI2Kv14_to_CI3Kv14 = 2.0 * b3Kv14 * alpha_act14 / b2Kv14;
  CI3Kv14_to_OIKv14  =      b4Kv14 * alpha_act14 / b3Kv14;

  C1Kv14_to_C0Kv14 =      beta_act14;
  C2Kv14_to_C1Kv14 = 2.0 * beta_act14;
  C3Kv14_to_C2Kv14 = 3.0 * beta_act14;
  OKv14_to_C3Kv14  = 4.0 * beta_act14;

  CI1Kv14_to_CI0Kv14 =             beta_act14 / f1Kv14;
  CI2Kv14_to_CI1Kv14 = 2.0 * f1Kv14 * beta_act14 / f2Kv14;
  CI3Kv14_to_CI2Kv14 = 3.0 * f2Kv14 * beta_act14 / f3Kv14;
  OIKv14_to_CI3Kv14  = 4.0 * f3Kv14 * beta_act14 / f4Kv14;

  C0Kv14_to_CI0Kv14 = beta_inact14;
  C1Kv14_to_CI1Kv14 = f1Kv14 * beta_inact14;
  C2Kv14_to_CI2Kv14 = f2Kv14 * beta_inact14;
  C3Kv14_to_CI3Kv14 = f3Kv14 * beta_inact14;
  OKv14_to_OIKv14   = f4Kv14 * beta_inact14;

  CI0Kv14_to_C0Kv14 = alpha_inact14;
  CI1Kv14_to_C1Kv14 = alpha_inact14 / b1Kv14;
  CI2Kv14_to_C2Kv14 = alpha_inact14 / b2Kv14;
  CI3Kv14_to_C3Kv14 = alpha_inact14 / b3Kv14;
  OIKv14_to_OKv14   = alpha_inact14 / b4Kv14;

  a1 = (C0Kv14_to_C1Kv14 + C0Kv14_to_CI0Kv14) * C0Kv14;
  a2 = C1Kv14_to_C0Kv14 * C1Kv14 + CI0Kv14_to_C0Kv14 * CI0Kv14;
  dC0Kv14 = a2 - a1;

  a1 = (C1Kv14_to_C2Kv14 + C1Kv14_to_C0Kv14 + C1Kv14_to_CI1Kv14) * C1Kv14;
  a2 = C2Kv14_to_C1Kv14 * C2Kv14 + CI1Kv14_to_C1Kv14 * CI1Kv14 + C0Kv14_to_C1Kv14 * C0Kv14;
  dC1Kv14 = a2 - a1;

  a1 = (C2Kv14_to_C3Kv14 + C2Kv14_to_C1Kv14 + C2Kv14_to_CI2Kv14) * C2Kv14;
  a2 = C3Kv14_to_C2Kv14 * C3Kv14 + CI2Kv14_to_C2Kv14 * CI2Kv14 + C1Kv14_to_C2Kv14 * C1Kv14;
  dC2Kv14 = a2 - a1;

  a1 = (C3Kv14_to_OKv14 + C3Kv14_to_C2Kv14 + C3Kv14_to_CI3Kv14) * C3Kv14;
  a2 = OKv14_to_C3Kv14 * OKv14 + CI3Kv14_to_C3Kv14 * CI3Kv14 + C2Kv14_to_C3Kv14 * C2Kv14;
  dC3Kv14 = a2 - a1;

  a1 = (OKv14_to_C3Kv14 + OKv14_to_OIKv14) * OKv14;
  a2 = C3Kv14_to_OKv14 * C3Kv14 + OIKv14_to_OKv14 * OIKv14;
  dOKv14 = a2 - a1;

  a1 = (CI0Kv14_to_C0Kv14 + CI0Kv14_to_CI1Kv14) * CI0Kv14;
  a2 = C0Kv14_to_CI0Kv14 * C0Kv14 + CI1Kv14_to_CI0Kv14 * CI1Kv14;
  dCI0Kv14 = a2 - a1;

  a1 = (CI1Kv14_to_CI2Kv14 + CI1Kv14_to_C1Kv14 + CI1Kv14_to_CI0Kv14) * CI1Kv14;
  a2 = CI2Kv14_to_CI1Kv14 * CI2Kv14 + C1Kv14_to_CI1Kv14 * C1Kv14 + CI0Kv14_to_CI1Kv14 * CI0Kv14;
  dCI1Kv14 = a2 - a1;

  a1 = (CI2Kv14_to_CI3Kv14 + CI2Kv14_to_C2Kv14 + CI2Kv14_to_CI1Kv14) * CI2Kv14;
  a2 = CI3Kv14_to_CI2Kv14 * CI3Kv14 + C2Kv14_to_CI2Kv14 * C2Kv14 + CI1Kv14_to_CI2Kv14 * CI1Kv14;
  dCI2Kv14 = a2 - a1;

  a1 = (CI3Kv14_to_OIKv14 + CI3Kv14_to_C3Kv14 + CI3Kv14_to_CI2Kv14) * CI3Kv14;
  a2 = OIKv14_to_CI3Kv14 * OIKv14 + C3Kv14_to_CI3Kv14 * C3Kv14 + CI2Kv14_to_CI3Kv14 * CI2Kv14;
  dCI3Kv14 = a2 - a1;

  // OIKv14 is "removed"
  // a1 = (OIKv14_to_OKv14+OIKv14_to_CI3Kv14)*OIKv14;
  // a2 = OKv14_to_OIKv14*OKv14 + CI3Kv14_to_OIKv14*CI3Kv14;
  // dOIKv14 = a2 - a1;


  //------MOVE LOCAL VARIABLES TO DERIVATIVE AND CURRENT ARRAYS------------------

  Fstate[index_V] = dV;
  Fstate[index_mNa] = dmNa;
  Fstate[index_hNa] = dhNa;
  Fstate[index_jNa] = djNa;
  Fstate[index_Nai] = dNai;
  Fstate[index_Ki] = dKi;
  Fstate[index_xKs] = dxKs;
  Fstate[index_C0Kv43] = dC0Kv43;
  Fstate[index_C1Kv43] = dC1Kv43;
  Fstate[index_C2Kv43] = dC2Kv43;
  Fstate[index_C3Kv43] = dC3Kv43;
  Fstate[index_OKv43] = dOKv43;
  Fstate[index_CI0Kv43] = dCI0Kv43;
  Fstate[index_CI1Kv43] = dCI1Kv43;
  Fstate[index_CI2Kv43] = dCI2Kv43;
  Fstate[index_CI3Kv43] = dCI3Kv43;
  Fstate[index_OIKv43] = -dCI0Kv43 - dCI1Kv43 - dCI2Kv43 - dCI3Kv43 - dOKv43 - dC0Kv43 - dC1Kv43 - dC2Kv43 - dC3Kv43;
  Fstate[index_C0Kv14] = dC0Kv14;
  Fstate[index_C1Kv14] = dC1Kv14;
  Fstate[index_C2Kv14] = dC2Kv14;
  Fstate[index_C3Kv14] = dC3Kv14;
  Fstate[index_OKv14] = dOKv14;
  Fstate[index_CI0Kv14] = dCI0Kv14;
  Fstate[index_CI1Kv14] = dCI1Kv14;
  Fstate[index_CI2Kv14] = dCI2Kv14;
  Fstate[index_CI3Kv14] = dCI3Kv14;
  Fstate[index_OIKv14] = -dCI0Kv14 - dCI1Kv14 - dCI2Kv14 - dCI3Kv14 - dOKv14 - dC0Kv14 - dC1Kv14 - dC2Kv14 - dC3Kv14;
  Fstate[index_CaTOT] = dCaTOT;
  Fstate[index_C1Herg] = dC1Herg;
  Fstate[index_C2Herg] = dC2Herg;
  Fstate[index_C3Herg] = dC3Herg;
  Fstate[index_OHerg] = dOHerg;
  Fstate[index_IHerg] = -dC1Herg - dC2Herg - dC3Herg - dOHerg;


  Fstate[index_Cai] = dCai;
  Fstate[index_CaNSR] = dCaNSR;
  Fstate[index_CaSL] = dCaSL;
  Fstate[index_LTRPNCa] = dLTRPNCa;
  Fstate[index_HTRPNCa] = dHTRPNCa;
  Fstate[index_xKs2] = dxKs2;


  if (keepc) {
    current[index_INa] = INa;
    current[index_IKr] = IKr;
    current[index_IKs] = IKs;
    current[index_Ito1] = Ito1;
    current[index_IK1] = IK1;
    current[index_IKp] = IKp;
    current[index_INaCa] = INaCa_tot;
    current[index_INaK] = INaK;
    current[index_IpCa] = IpCa;
    current[index_ICab] = ICab;
    current[index_INab] = INab;
    current[index_ICa] = ICa;
    current[index_JDHPR] = ICa * Acap / (-2.0 * Vmyo * Faraday * 1000.0);
    current[index_Jup] = Jup;
    current[index_Jtrpn] = Jtrpn ;
    current[index_Jtr] = Jtr * VJSR / Vmyo;
    current[index_Jxfer] = Jxfer_tot * VSS / Vmyo;
    current[index_IKv43] = IKv43;
    current[index_IKv14] = IKv14 ;
    current[index_IKv14_K] = IKv14_K ;
    current[index_IKv14_Na] = IKv14_Na;
    current[index_Ito2] = Ito2;
    current[index_Istim] = Istim;
    current[index_Itot] = Itot;
    current[index_Insca] = insna + insk;
  }
}


//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
/*
void StochModel::lastcall(double time, double Iext, double JCa)
{
  double F[N];
  const int dummytrue = 1;

  double FRUdep_states[Nstates_FRUdep];
  double Jxfer, Jtr, ICa, Ito2;

  set_FRUdep_states(FRUdep_states);

  distrib_simFRU(time, time, FRUdep_states, &Jxfer, &Jtr, &ICa, &Ito2);
  //send_calc_fru_flux(FRUdep_states, &Jxfer, &Jtr, &ICa, &Ito2);
  fcn(time, state, F, current, dummytrue, Jxfer, Jtr, ICa, Ito2, Iext, JCa);
}
*/
