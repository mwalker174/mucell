#include "ReleaseUnit.h"
#include "indices.h"
#include <cstring>

double ReleaseUnit::JRyRmax = 0;
double ReleaseUnit::tautr = 0;
double ReleaseUnit::tauxfer = 0;
double ReleaseUnit::tauss2ss = 0;
double ReleaseUnit::VJSR = 0;
double ReleaseUnit::VSS = 0;
double ReleaseUnit::PCa = 0;
double ReleaseUnit::BSLtot = 0;
double ReleaseUnit::CSQNtot = 0;
double ReleaseUnit::BSRtot = 0;
double ReleaseUnit::KBSL = 0;
double ReleaseUnit::KmCSQN = 0;
double ReleaseUnit::KBSR = 0;



ReleaseUnit::ReleaseUnit()
{
  bOrphan = 0;
  N_LCC_Active = 0;
}

ReleaseUnit::ReleaseUnit(const ReleaseUnit& other)
{

  std::memcpy(FRU_states, other.FRU_states, sizeof(double)*Nstates_FRU);
  std::memcpy(LCC_States, other.LCC_States, sizeof(char)*Max_LCCs_per_cleft);
  std::memcpy(LCC_Vdep, other.LCC_Vdep, sizeof(char)*Max_LCCs_per_cleft);
  std::memcpy(LCC_Mode2, other.LCC_Mode2, sizeof(char)*Max_LCCs_per_cleft);
  RyR_state = other.RyR_state;
  Ito2_state = other.Ito2_state;
  Ri = other.Ri;
  bOrphan = other.bOrphan;
  N_LCC_Active = other.N_LCC_Active;

}


double ReleaseUnit::Get_RyR_Open()
{
  return (double)RyR_state;
}

double ReleaseUnit::Get_RyR_Ready()
{

  return (double)(NRyRs_per_cleft-RyR_state);
}

double ReleaseUnit::Get_LCC_Open()
{
  return Get_LCC_Open_Cleft();
}

double ReleaseUnit::Get_LCC_Vinact()
{
  double sum = 0;

  for (int j = 0; j < N_LCC_Active; j++) {
    if (LCC_Vdep[j] == Oy_LType)
      sum += 1.0;
  }

  return sum;
}

double ReleaseUnit::Get_LCC_NormMode()
{
  double sum = 0;
    for (int j = 0; j < N_LCC_Active; j++) {
      if (LCC_States[j] <= 6)
	sum += 1.0;
    }

  return sum;
}

double ReleaseUnit::Get_Ito2_Open()
{
  return (double)Ito2_state;
}

double ReleaseUnit::Get_CaSS_Avg()
{
  return FRU_states[1];
}

double ReleaseUnit::Get_CaTOT_SS()
{
  double d = bOrphan ? ORPHAN_SCALE : 1; //Orphan
  return  d * FRU_states[1] * (1.0 + (BSRtot/d) / (KBSR + FRU_states[1]) + (BSLtot/d) / (KBSL + FRU_states[1]));
}

double ReleaseUnit::Get_CaTOT_JSR()
{
  return FRU_states[index_frustates_CaJSR] * (1.0 + CSQNtot / (KmCSQN + FRU_states[index_frustates_CaJSR]) );
}

double ReleaseUnit::Get_JRyR_Cleft()
{
    return JRyRmax * Get_RyR_Open() * (FRU_states[index_frustates_CaJSR] - FRU_states[1]);
}

double ReleaseUnit::Get_JRyR()
{
  return Get_JRyR_Cleft();
}


void ReleaseUnit::Set_Orphan(char c) {
  bOrphan = c;
}

char ReleaseUnit::isOrphan() {
  return bOrphan;
}


