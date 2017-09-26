/*       ----------------------------------------------------

	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.

	 Name of Program: Local Control Model
	 Version: Documented Version, C
	 Date: November 2003

	 --------------------------------------------------

	 write_fru_props.c - Save properties of individual FRUs to files

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "StochModel.h"


// If the write_fru_props_flag is set, then this subroutine is called to write
// the data describing individual release unit variable to the output files

void StochModel::write_fru_props(FILE *filee[], FILE *file_avg[], FILE *file_avg_all,
                                 const double time,
                                 const double V,
                                 double Cai,
                                 double CaNSR,
                                 double FRU_states[NFRU][Nstates_FRU],
                                 int LType_state[NFRU][Nclefts_FRU][NLCCs_per_cleft][Nindepstates_LType],
                                 int RyR_state[NFRU][Nclefts_FRU][NRyRs_per_cleft],
                                 int Ito2_state[NFRU][Nclefts_FRU])
{
  double VF_over_RT, VFsq_over_RT;
  double Jtr, CaJSR;
  double JLType_tot1, JRyR_tot1, Jxfer_tot1, CaSS_avg1, PRyR_open1;
  double JLType_totALL, JRyR_totALL, Jtr_totALL, Jxfer_totALL;
  double CaSS_avgALL, CaJSR_avgALL, PRyR_openALL, num_avg;
  double CaSS[Nclefts_FRU], JLType[Nclefts_FRU], JRyR[Nclefts_FRU], Jxfer[Nclefts_FRU];
  int i, j, jcleft, NRyR_open[Nclefts_FRU];
  double Ito2[Nclefts_FRU], Ito2_tot1, Ito2_avgALL;
  double tprint;

  tprint = time;

  if (ts_sec) tprint = tprint * 1.e-3;

  VF_over_RT = V / RT_over_F;
  VFsq_over_RT = (1000.0 * Faraday) * VF_over_RT;
  JLType_totALL = 0.0;
  JRyR_totALL = 0.0;
  Jtr_totALL = 0.0;
  Jxfer_totALL = 0.0;
  CaSS_avgALL = 0.0;
  CaJSR_avgALL = 0.0;
  PRyR_openALL = 0.0;
  Ito2_tot1 = 0.0;
  Ito2_avgALL = 0.0;

  if (fru_end > NFRU) {
    fprintf(stderr, "Inappropriate fru_end (%d) in write_fru_props\n", fru_end);
    exit(-3);
  }

  if (fru_start < 0) {
    fprintf(stderr, "Inappropriate fru_start (%d) in write_fru_props\n", fru_start);
    exit(-4);
  }

  for(j = fru_start; j <= fru_end; j++) {
    JLType_tot1 = 0.0;
    JRyR_tot1 = 0.0;
    Jxfer_tot1 = 0.0;
    CaSS_avg1 = 0.0;
    PRyR_open1 = 0.0;
    Ito2_tot1 = 0.0;

    CaJSR = FRU_states[j][index_frustates_CaJSR];
    CaJSR_avgALL = CaJSR_avgALL + CaJSR;
    Jtr = (CaNSR - CaJSR) / tautr * VJSR / Vmyo;
    Jtr_totALL = Jtr_totALL + Jtr;

    for(jcleft = 0; jcleft < Nclefts_FRU; jcleft++) {
      CaSS[jcleft] = FRU_states[j][jcleft + 1];
      CaSS_avg1 += CaSS[jcleft];
      CaSS_avgALL += CaSS[jcleft];
      JLType[jcleft] = 0.0;
      int jlcc;

      for (jlcc = 0; jlcc < NLCCs_per_cleft; jlcc++) {

        if ( ((LType_state[j][jcleft][jlcc][index_LCC_states] == O1_LType) || (LType_state[j][jcleft][jlcc][index_LCC_states] == O2_LType))
             && LType_state[j][jcleft][jlcc][index_LCC_Vinact] == Oy_LType) {
          double jca = PCa * 4.0 * VFsq_over_RT * (CaSS[jcleft] * exp(2.0 * VF_over_RT) * lcc_gamma_ca - Cao * 0.341) / (exp(2.0 * VF_over_RT) - 1.0) / (-2.0 * Vmyo * Faraday * 1000.0);
          //     .				/(exp(2.0*VF_over_RT)-1.0)  //  /(-2.0*Vmyo*Faraday*1000.0)
          JLType[jcleft] += jca;
          JLType_tot1 += jca;
          JLType_totALL += jca;
        }
      }

      Ito2[jcleft] = 0.0;

      if (Ito2_state[j][jcleft] == O_Ito2) {
        Ito2[jcleft] = PCl * VFsq_over_RT * (Cli * exp(-VF_over_RT) - Clo) / (exp(-VF_over_RT) - 1.0);
        Ito2_tot1 += Ito2[jcleft];
        Ito2_avgALL += Ito2[jcleft];
      }

      NRyR_open[jcleft] = 0;

      for(i = 0; i < NRyRs_per_cleft; i++) {
        if ((RyR_state[j][jcleft][i] == O1_RyR) || (RyR_state[j][jcleft][i] == O2_RyR)
            || (RyR_state[j][jcleft][i] == O3_RyR)) {
          NRyR_open[jcleft] = NRyR_open[jcleft] + 1;
        }
      }

      PRyR_open1 = PRyR_open1 + (double)(NRyR_open[jcleft]);
      PRyR_openALL = PRyR_openALL + (double)(NRyR_open[jcleft]);
      JRyR[jcleft] = JRyRmax * (double)(NRyR_open[jcleft]) * (CaJSR - CaSS[jcleft]) * VSS / Vmyo;
      JRyR_tot1 = JRyR_tot1 + JRyR[jcleft];
      JRyR_totALL = JRyR_totALL + JRyR[jcleft];
      Jxfer[jcleft] = (CaSS[jcleft] - Cai) / tauxfer * VSS / Vmyo;
      Jxfer_tot1 = Jxfer_tot1 + Jxfer[jcleft];
      Jxfer_totALL = Jxfer_totALL + Jxfer[jcleft];

      fprintf(filee[j - fru_start], "%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %d\n",
              tprint, JLType[jcleft], JRyR[jcleft], Jtr, Jxfer[jcleft], CaSS[jcleft], CaJSR, Ito2[jcleft], NRyR_open[jcleft]);
      fflush(filee[j - fru_start]);
    }

    CaSS_avg1 = CaSS_avg1 / Nclefts_FRU;
    PRyR_open1 = PRyR_open1 / (double)(NRyRs_per_cleft * Nclefts_FRU);

    fprintf(file_avg[j - fru_start], "%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",
            tprint, JLType_tot1, JRyR_tot1, Jtr, Jxfer_tot1, CaSS_avg1, CaJSR, Ito2_tot1, PRyR_open1);
    fflush(filee[j - fru_start]);
  }

  num_avg = (double)(fru_end - fru_start + 1);
  CaSS_avgALL = CaSS_avgALL / num_avg / Nclefts_FRU;
  CaJSR_avgALL = CaJSR_avgALL / num_avg;
  PRyR_openALL = PRyR_openALL / (num_avg * (double)(NRyRs_per_cleft * Nclefts_FRU));
  fprintf(file_avg_all, "%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",
          tprint, JLType_totALL, JRyR_totALL, Jtr_totALL, Jxfer_totALL, CaSS_avgALL, CaJSR_avgALL, Ito2_avgALL, PRyR_openALL);
  fflush(file_avg_all);
}


void StochModel::open_fru_props(const int filenumber, FILE *filee[]) // Open output all files for local release unit data
{
  char filename[256];
  int j;

  for(j = fru_start; j <= fru_end; j++) {

    sprintf(filename, "%s%s_%d.%d.txt", output_dir, output_fruprops_file, j, filenumber);

    if ((filee[j - fru_start] = fopen(filename, "w+")) == NULL) {
      fprintf(stderr, "Cannot open file '%s'!\n", filename);
      exit(-5);
    }

    fprintf(filee[j - fru_start],
            "%s %s %s %s %s %s %s %s %s\n",
            "Time", "JLType", "JRyR", "Jtr", "Jxfer", "CaSS", "CaJSR", "Ito2", "NRyR_open");
  }
}

void StochModel::open_fru_props_avg(const int filenumber, FILE *filee[]) // Open output all files for local release unit data
{
  char filename[256];
  int j;

  for(j = fru_start; j <= fru_end; j++) {
    sprintf(filename, "%s%s_%d.%d.txt", output_dir, output_frupropsavg_file, j, filenumber);

    if ((filee[j - fru_start] = fopen(filename, "w+")) == NULL) {
      fprintf(stderr, "Cannot open file '%s'!\n", filename);
      exit(-5);
    }

    fprintf(filee[j - fru_start],
            "%s %s %s %s %s %s %s %s %s\n",
            "Time", "JLType_tot", "JRyR_tot", "Jtr", "Jxfer_tot", "CaSS_avg", "CaJSR", "Ito2", "PRyR_open");
  }
}

FILE *StochModel::open_fru_props_avgtot(const int filenumber) // Open output all files for local release unit data
{
  FILE *filee;
  char filename[256];

  sprintf(filename, "%s%s.%d.txt", output_dir, output_frupropsavgtot_file, filenumber);

  if ((filee = fopen(filename, "w+")) == NULL) {
    fprintf(stderr, "Cannot open file '%s'!\n", filename);
    exit(-5);
  }

  fprintf(filee,
          "%s %s %s %s %s %s %s %s %s\n",
          "Time", "JLType_totAll", "JRyR_totAll", "Jtr_totAll", "Jxfer_totAll", "CaSS_avgAll", "CaJSR_avgAll", "Ito2", "PRyR_openAll");

  return filee;
}

