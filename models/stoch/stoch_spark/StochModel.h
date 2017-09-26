#ifndef __STOCHMODEL_H
#define __STOCHMODEL_H

#include <stdio.h>
#include <string>
#include <iostream>
#include <map>
#include <list>

#include "ReleaseUnit.h"
#include "indices.h"

class StochModel
{

public:

  StochModel();
  StochModel(const StochModel& other);

  //Set parameters
  void Set_Parameters(std::map<std::string, double> params);
  void Read_Parameters(std::string filename);

  //Set random number generator seed (overwrites value from Set_Parameters)
  void Set_Seed(unsigned long s);

  void sgrnd(unsigned long, int *, unsigned long[mtN + 1]);
  double Get_Rand_Unif(unsigned long[mtN+1], int &);
  double generateGaussian(unsigned long[mtN+1], int &, double, double);

  //Initialize default states (run after reading parameters)
  void Initialize_Default_State();

  //Initializes values in Currents array
  void Initialize_Currents(const double Iext, const double JCa);

  //Initializes static parameters to ReleaseUnit class
  void Initialize_ReleaseUnit_Params();

  //Run to next time step and save states
  void Delta_V_Step(double dt, double V, double &dV);
  void Integrate_Iext(double dt, double Iext);
  void Integrate_Iext_Ca(double dt, double Iext, double JCa);

  //Run to next time step and do not save states
  void Delta_V_Peek(double dt, double I, double &dV);

  //Accessors
  double& Get_V();
  double& Get_Cai();
  double& Get_CaSL();
  void Set_CaSR(double c);
  double Get_Capacitance();
  double Get_SR_Load();
  double Get_LCC_Vinact();
  double Get_LCC_ModeCa();
  double Get_LCC_Mode2_Open();
  double Get_CaSS();
  double Get_CaJSR();
  double Get_mNa();
  double Get_hNa();
  double Get_jNa();
  double Get_Ca_Tot();
  double Get_RyR_Open();
  double Get_JRyR();
  double Get_INaCa();
  double Get_Ito2();
  double Get_Nai();

  //For single-cell simulations only (stochtest)
  double Get_T_Final();
  double Get_Period();
  double Get_T_End_Stim();
  double Get_I_Stim();
  double Get_VClamp_Flag();
  double Get_VClamp_ClampV();
  double Get_VClamp_HoldV();
  double Get_VClamp_Freq();
  double Get_VClamp_Duration();
  double Get_VClamp2_T1();
  double Get_VClamp2_T2();
  double Get_VClamp2_ClampV();

  //Returns default random seed
  int Get_Seed();

  //Read and write full model state (ASCII)
  void Write_State(std::ostream &os);
  void Read_State(std::istream &is);

  //Write states, currents, and "other states"
  void Write_Info_Header(std::ostream &os);
  void Write_Info(std::ostream &os, double t); //Writes time and all variables
  void Write_Info(std::ostream &os, std::list<std::string> vars); //Writes vars
  void Write_Grid(std::ostream &os, double t);

  //Write individual release unit states
  void Write_CRU_Header(std::ostream &os);
  void Write_CRU(std::ostream &os, double t, int iFRU);

  //Control functions
  int Toggle_Eject_Jup();
  int Toggle_Clamp_NSR();
  void Set_Cai(double c);

private:
  double rk54pd(double, double, double, double);
  void euler(double t, double tstep, double Iext);
  void euler_adaptive(double t, double tstep, double Iext);
  
  void fcn(const double, std::vector<double>&, std::vector<double>&, double[Ncur], const int, const double, const double, const double, const double, const double, const double);
  void lastcall(double, double, double);
  double read_next_double(FILE *);
  int read_next_int(FILE *);
  void initialize_mpi(int *, char **, int *, int *);
  void end_mpi(void);
  void parallel(void);
  void distrib_simFRU(double, double, double[Nstates_FRUdep], double[Nstates_FRUdep],
		      std::vector<double>&, std::vector<double>&,
		      double *, double *, double *, double *);
  void send_save_state(void);
  void parallel_save_state(void);
  void parallel_resume_state(void);
  void parallel_compute_simfru(void);
  void send_resume_state(void);
  void send_calc_fru_flux(double[Nstates_FRUdep], double *, double *, double *, double *);
  void send_calc_fru_avg();
  void parallel_calc_fru_avg(void);
  void parallel_calc_fru_flux(void);
  void set_FRUdep_states(double[Nstates_FRUdep]);
  void set_FRUdep_statesf(double[Nstates_FRUdep],std::vector<double>);

  // Release units
  std::vector<ReleaseUnit> FRUs;
  std::vector<ReleaseUnit> FRUs_hold;
  int Num_Orphans;

  //Release unit neighbors
  std::vector<std::vector<int> > FRU_neighbs;
  std::vector<std::vector<double> > FRU_neighbs_distance;

  //Solver error weights
  std::vector<double> errweight;

  //State variable arrays
  std::vector<double> state;
  double current[Ncur];
  double otherstates[Nother];

  //Time step size
  double oldstepsize;

  //Control variables
  int bEjectJup; //SERCA uptake ejected from cell instead
  int bClampNSR;

  //////////////////////////////
  //Parameters read in from file
  //////////////////////////////

  double Default_Cai;
  double Default_CaSR;
  double Default_Nai;
  double beta_flag;
  double mode2_frac;
  double ryr_kplus_scale;
  double serca_kf_scale;
  double ikr_scale;
  double serca_vmax_scale;
  double ncx_scale;
  double ik1_scale;
  double ito1_scale;
  double frac_orphan;
  double frac_ncx_sl;
  double Acap;
  double Cao; // extracellular Ca++ concentration (mM)
  double nak_scale;
  double NFRU_scale;
  int NFRU;
  int NFRU_X, NFRU_Y, NFRU_Z;
  double clamp_nai;

  unsigned long rand_seed;
  unsigned long mt[OMP_THREAD_MAX][mtN+1]; /* the array for the state vector */
  int mti[OMP_THREAD_MAX]; /* mti==mtN+1 means mt[mtN] is not initialized */
  unsigned long mt_hold[OMP_THREAD_MAX][mtN+1];
  int mti_hold[OMP_THREAD_MAX]; 

  //Single-cell simulations only (stochtest)
  double t_final;
  double period_stim;
  double t_end_stim;
  double I_stim;
  double VClamp_Flag;
  double VClamp_ClampV;
  double VClamp_HoldV;
  double VClamp_Freq;
  double VClamp_Duration;
  double VClamp2_T1, VClamp2_T2, VClamp2_ClampV;

  //////////////////////////////
  //Hard-coded static parameters
  //////////////////////////////
  static const double step_min;
  static const double step_max;
  static const double tolrk;

  //LCC active fraction
  static const double lcc_frac_active;

  //Ion concentrations
  static const double Nao; // extracellular Na+ concentration (mM)
  static const double Ko; // extracellular K+ concentration (mM)
  static const double Default_Ki;

  //Beta-Adrenergic conditions
  static const double beta_iks_scale;
  static const double beta_ikr_rate_scale;
  static const double beta_ltrpn_koff_scale;
  static const double beta_frac_active;

  // External physical constants

  // total troponin low affinity site conc. (mM)
  static const double LTRPNtot;
  // total troponin high affinity site conc. (mM)
  static const double HTRPNtot;
  // Ca++ on rate for troponin high affinity sites (1/[mM*ms])
  static const double khtrpn_plus;
  // Ca++ off rate for troponin high affinity sites (1/ms)
  static const double khtrpn_minus;
  // Ca++ on rate for troponin low affinity sites (1/[mM*ms])
  static const double kltrpn_plus;
  // Ca++ off rate for troponin low affinity sites (1/ms)
  static const double kltrpn_minus;
  // total myoplasmic calmodulin concentration (mM)
  static const double CMDNtot;
  // total myoplasmic EGTA concentration (mM)
  static const double EGTAtot;
  // Ca++ half sat. constant for calmodulin (mM)
  static const double KmCMDN;
  // Ca++ half sat. constant for EGTA (mM)
  static const double KmEGTA;

  /*------STANDARD IONIC CONCENTRATIONS-----------------------------------------*/

  static const double Clo; // extracellular Cl- concentration (mM)
  static const double Cli; // intracellular Cl- concentration (mM)

  /*------CELL GEOMETRY PARAMETERS----------------------------------------------*/

  /* NOTE: The assumption of 1uF/cm^2 holds for this model
  	 therefore Acap in cm^2 is equal to whole cell
  	 capacitance in uF.
  */

  //static const double Acap; // capacitive membrane area (cm^2)
  static const double Vmyo; // myoplasmic volume (uL)
  static const double VNSR; // junctional SR volume (uL)
  static const double PCl; //(cm/s) *uF

  //Subsarcolemmal compartment
  static const double VSL;
  static const double BSLtot_SL; // (mM) (SL Compartment)
  static const double tauSL2cyto; //(ms)

  //Static Const variables for ReleaseUnit class

  // JSR to subspace through a single RyR (1/ms)
  static const double JRyRmax;
  // NSR to JSR (ms)
  static const double tautr;
  // subspace to cytosol (ms)
  static const double tauxfer;
  // subspace to subspace (ms)
  static const double tauss2ss;
  // Grid diffusion
  static const double taugrid;

  static const double VJSR; // network SR volume (uL) //0.75 * 3e-12;
  static const double VSS; // subspace volume (uL)
  static const double PCa; //(cm/s) *uF

  static const double BSLtot; // (mM)
  static const double CSQNtot;
  static const double BSRtot;
  static const double KBSL;
  static const double KmCSQN;
  static const double KBSR;

};

#endif
