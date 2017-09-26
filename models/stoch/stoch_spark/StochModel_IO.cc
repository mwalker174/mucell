#include "StochModel.h"
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <limits>
#include <omp.h>

void StochModel::Write_State(std::ostream &os)
{

  //int PREC = std::numeric_limits<double>::digits10+2; //Floating point decimal precision

  os << NFRU_X << std::endl;
  os << NFRU_Y << std::endl;
  os << NFRU_Z << std::endl;
  os << oldstepsize << std::endl;
  //os << rand_seed << std::endl;

  for (unsigned i = 0; i < state.size(); i++) { 
    os << state[i] << std::endl;
  }

  /*
  for (int j = 0; j < OMP_THREAD_MAX; j++) {
    for (int i = 0; i < mtN + 1; i++) {
      os << mt[j][i] << std::endl;
    }
  }

  for (int j = 0; j < OMP_THREAD_MAX; j++) {
    os << mti[j] << std::endl;
  }
  */

  for (int i = 0; i < NFRU; i++) {
    os << FRU_neighbs[i].size() << std::endl;
    for (unsigned j = 0; j < FRU_neighbs[i].size(); j++) {
      os << FRU_neighbs[i][j] << std::endl;
      os << FRU_neighbs_distance[i][j] << std::endl;
    }
  }

  for (int i = 0; i < NFRU; i++) {
    for (int j = 0; j < Nstates_FRU; j++) {
      os << FRUs[i].FRU_states[j] << std::endl;
    }
    
    os << FRUs[i].N_LCC_Active << std::endl;
    for (int m = 0; m < FRUs[i].N_LCC_Active; m++) {
	os << (int)FRUs[i].LCC_States[m] << std::endl;
	os << (int)FRUs[i].LCC_Vdep[m] << std::endl;
	os << (int)FRUs[i].LCC_Mode2[m] << std::endl;
    }
    
    os << (int)FRUs[i].RyR_state << std::endl;
    os << (int)FRUs[i].Ito2_state << std::endl;
    
    os << FRUs[i].Ri << std::endl;
    os << (int)FRUs[i].bOrphan << std::endl;
  }

}

void StochModel::Read_State(std::istream &is)
{

  std::string val;

  std::getline(is, val);
  NFRU_X = std::strtol(val.c_str(), NULL, 0);
  std::getline(is, val);
  NFRU_Y = std::strtol(val.c_str(), NULL, 0);
  std::getline(is, val);
  NFRU_Z = std::strtol(val.c_str(), NULL, 0);
  NFRU = NFRU_X*NFRU_Y*NFRU_Z;

  std::getline(is, val);
  oldstepsize = std::strtod(val.c_str(), NULL);

  //std::getline(is, val);
  //rand_seed = std::strtol(val.c_str(), NULL, 0);

  state = std::vector<double>(N_states+NFRU);
  for (unsigned i = 0; i < state.size(); i++) { 
    std::getline(is, val);
    state[i] = std::strtod(val.c_str(), NULL);
  }

  /*for (int j = 0; j < OMP_THREAD_MAX; j++) {
    for (int i = 0; i < mtN + 1; i++) {
      std::getline(is, val);
      mt[j][i] = (unsigned long)std::strtol(val.c_str(), NULL, 0);
    }
  }

  for (int j = 0; j < OMP_THREAD_MAX; j++) {
    std::getline(is, val);
    mti[j] = std::strtol(val.c_str(), NULL, 0);
    }*/

  FRU_neighbs = std::vector<std::vector<int> >(NFRU);
  FRU_neighbs_distance = std::vector<std::vector<double> >(NFRU);

  for (int i = 0; i < NFRU; i++) {
    std::getline(is, val);
    int n_neighbs = (int)std::strtol(val.c_str(), NULL, 0);
    FRU_neighbs[i] = std::vector<int>(n_neighbs);
    FRU_neighbs_distance[i] = std::vector<double>(n_neighbs);
    for (unsigned j = 0; j < FRU_neighbs[i].size(); j++) {
      std::getline(is, val);
      FRU_neighbs[i][j]= (int)std::strtol(val.c_str(), NULL, 0);
      std::getline(is, val);
      FRU_neighbs_distance[i][j]= std::strtod(val.c_str(), NULL);
    }
  }

  for (int i = 0; i < NFRU; i++) {
    for (int j = 0; j < Nstates_FRU; j++) {
      std::getline(is, val);
      FRUs[i].FRU_states[j] = std::strtod(val.c_str(), NULL);
    }
    
    std::getline(is, val);
    FRUs[i].N_LCC_Active = (int)std::strtol(val.c_str(), NULL, 0);
    for (int m = 0; m < FRUs[i].N_LCC_Active; m++) {
      std::getline(is, val);
      FRUs[i].LCC_States[m] = (char)std::strtol(val.c_str(), NULL, 0);
      std::getline(is, val);
      FRUs[i].LCC_Vdep[m] = (char)std::strtol(val.c_str(), NULL, 0);
      std::getline(is, val);
      FRUs[i].LCC_Mode2[m] = (char)std::strtol(val.c_str(), NULL, 0);
    }
      
    std::getline(is, val);
    FRUs[i].RyR_state =(int)std::strtol(val.c_str(), NULL, 0);
      
    std::getline(is, val);
    FRUs[i].Ito2_state = (int)std::strtol(val.c_str(), NULL, 0);

    std::getline(is, val);
    FRUs[i].Ri = std::strtod(val.c_str(), NULL);
      
    std::getline(is, val);
    FRUs[i].bOrphan = (int)std::strtol(val.c_str(), NULL, 0);
  }

  //Re-initialize current values
  Initialize_Currents(0, 0);

}

void StochModel::Write_Info_Header(std::ostream &os)
{

  os << "Time ";

  os << "V mNa hNa jNa Nai Ki xKs C0Kv43 C1Kv43 C2Kv43 C3Kv43 OKv43 CI0Kv43 CI1Kv43 CI2Kv43 CI3Kv43 OIKv43 C0Kv14 C1Kv14 C2Kv14 C3Kv14 OKv14 CI0Kv14 CI1Kv14 CI2Kv14 CI3Kv14 OIKv14 CaToT C1Herg C2Herg C3Herg OHerg IHerg Cai CaNSR CaSL LTRPNCa HTRPNCa xKs2 ";

  os << "INa IKr IKs Ito1 IK1 IKp INaCa INaK IpCa ICab INab ICa JDHPR Jup Jtrpn Jtr Jxfer IKv43 IKv14 IKv14_K IKv14_Na Ito2 Istim Itot Insca ";

  os << "CaSSavg CaJSRavg JRyRtot PRyR_open PRyR_ready PNorm_mode PnotVinact PLType_open CaToT2 PIto2_open CaJSRtot CaSStot CaSRtot CaSLtot";

  os << std::endl;

}

void StochModel::Write_Info(std::ostream &os, double t)
{

  send_calc_fru_avg();

  os << t * 1e-3 << " ";

  int i;

  for (i = 0; i < N_states; i++) {
    os << state[i] << " ";
  }

  for (i = 0; i < Ncur; i++) {
    os << current[i] << " ";
  }

  for (i = 0; i < Nother; i++) {
    os << otherstates[i] << " ";
  }

  os << std::endl;

}

void StochModel::Write_Grid(std::ostream &os, double t)
{
	double tprint;
	tprint = t * 1.e-3;

	//fwrite(&state[index_state_Cai],sizeof(double),NFRU,file);

	os << "# vtk DataFile Version 2.0\n";
	os << "Cardiac myocyte t=" << tprint << " sec\n";
	os << "ASCII\n";
	os << "DATASET STRUCTURED_POINTS\n";
	os << "DIMENSIONS " << NFRU_Z << " " << NFRU_Y << " " << NFRU_X << "\n";
	os << "ORIGIN 0 0 0\n";
	os << "SPACING 2 1 1\n";

	// Write states
	os << "POINT_DATA " << NFRU << "\n"; 
	os << "SCALARS CaSL float 1\n";
	os << "LOOKUP_TABLE default\n";

	int idx = 0;
	for (int i = 0; i < NFRU_X; i++) {
	  for (int j = 0; j < NFRU_Y; j++) {
	    for (int k = 0; k < NFRU_Z; k++) {
	      if (idx < NFRU) {
		os << state[N_states+idx]*1e3 << "\n";
		//os << FRUs[idx].Get_CaSS_Avg()*1e3 << "\n";
	      } else {
		os << "0\n";
	      }
	      idx++;
	    }
	  }
	}
	os << "\n";
	/*
	os << "SCALARS CaSS float 1\n";
	os << "LOOKUP_TABLE default\n";
	idx = 0;
	for (int i = 0; i < NFRU_X; i++) {
	  for (int j = 0; j < NFRU_Y; j++) {
	    for (int k = 0; k < NFRU_Z; k++) {
	      if (idx < NFRU) {
		//os << state[N_states+idx]*1e3 << "\n";
		os << FRUs[idx].Get_CaSS_Avg()*1e3 << "\n";
	      } else {
		os << "0\n";
	      }
	      idx++;
	    }
	  }
	}
	os << "\n";

	os << "SCALARS CaJSR float 1\n";
	os << "LOOKUP_TABLE default\n";
	idx = 0;
	for (int i = 0; i < NFRU_X; i++) {
	  for (int j = 0; j < NFRU_Y; j++) {
	    for (int k = 0; k < NFRU_Z; k++) {
	      if (idx < NFRU) {
		os << FRUs[idx].Get_CaTOT_JSR() << "\n";
	      } else {
		os << "0\n";
	      }
	      idx++;
	    }
	  }
	}
	os << "\n";
	*/
	/*
	os << "SCALARS dist_min float 1\n";
	os << "LOOKUP_TABLE default\n";
	idx = 0;
	for (int i = 0; i < NFRU_X; i++) {
	  for (int j = 0; j < NFRU_Y; j++) {
	    for (int k = 0; k < NFRU_Z; k++) {
	      if (idx < NFRU) {
		double min = 10;
		for (int m = 0; m < FRU_neighbs_distance[idx].size(); m++) {
		  if (min > FRU_neighbs_distance[idx][m]) {
		    min = FRU_neighbs_distance[idx][m];
		  }
		}
		os << min << "\n";
	      } else {
		os << "0\n";
	      }
	      idx++;
	    }
	  }
	}
	os << "\n";
	*/
}


void StochModel::Write_Info(std::ostream &os, std::list<std::string> vars)
{

  for (std::list<std::string>::iterator it = vars.begin(); it != vars.end(); ++it) {
    if (!(*it).compare("V")) {
      os << state[index_V] << " ";
    } else if (!(*it).compare("mNa")) {
      os << state[index_mNa] << " ";
    } else if (!(*it).compare("hNa")) {
      os << state[index_hNa] << " ";
    } else if (!(*it).compare("jNa")) {
      os << state[index_jNa] << " ";
    } else if (!(*it).compare("Nai")) {
      os << state[index_Nai] << " ";
    } else if (!(*it).compare("Ki")) {
      os << state[index_Ki] << " ";
    } else if (!(*it).compare("xKs")) {
      os << state[index_xKs] << " ";
    } else if (!(*it).compare("C0Kv43")) {
      os << state[index_C0Kv43] << " ";
    } else if (!(*it).compare("C1Kv43")) {
      os << state[index_C1Kv43] << " ";
    } else if (!(*it).compare("C2Kv43")) {
      os << state[index_C2Kv43] << " ";
    } else if (!(*it).compare("C3Kv43")) {
      os << state[index_C3Kv43] << " ";
    } else if (!(*it).compare("OKv43")) {
      os << state[index_OKv43] << " ";
    } else if (!(*it).compare("CI0Kv43")) {
      os << state[index_CI0Kv43] << " ";
    } else if (!(*it).compare("CI1Kv43")) {
      os << state[index_CI1Kv43] << " ";
    } else if (!(*it).compare("CI2Kv43")) {
      os << state[index_CI2Kv43] << " ";
    } else if (!(*it).compare("CI3Kv43")) {
      os << state[index_CI3Kv43] << " ";
    } else if (!(*it).compare("OIKv43")) {
      os << state[index_OIKv43] << " ";
    } else if (!(*it).compare("C0Kv14")) {
      os << state[index_C0Kv14] << " ";
    } else if (!(*it).compare("C1Kv14")) {
      os << state[index_C1Kv14] << " ";
    } else if (!(*it).compare("C2Kv14")) {
      os << state[index_C2Kv14] << " ";
    } else if (!(*it).compare("C3Kv14")) {
      os << state[index_C3Kv14] << " ";
    } else if (!(*it).compare("OKv14")) {
      os << state[index_OKv14] << " ";
    } else if (!(*it).compare("CI0Kv14")) {
      os << state[index_CI0Kv14] << " ";
    } else if (!(*it).compare("CI1Kv14")) {
      os << state[index_CI1Kv14] << " ";
    } else if (!(*it).compare("CI2Kv14")) {
      os << state[index_CI2Kv14] << " ";
    } else if (!(*it).compare("CI3Kv14")) {
      os << state[index_CI3Kv14] << " ";
    } else if (!(*it).compare("OIKv14")) {
      os << state[index_OIKv14] << " ";
    } else if (!(*it).compare("CaTOT")) {
      os << state[index_CaTOT] << " ";
    } else if (!(*it).compare("C1Herg")) {
      os << state[index_C1Herg] << " ";
    } else if (!(*it).compare("C2Herg")) {
      os << state[index_C2Herg] << " ";
    } else if (!(*it).compare("C3Herg")) {
      os << state[index_C3Herg] << " ";
    } else if (!(*it).compare("OHerg")) {
      os << state[index_OHerg] << " ";
    } else if (!(*it).compare("IHerg")) {
      os << state[index_IHerg] << " ";
    } else if (!(*it).compare("Cai")) {
      os << state[index_Cai] << " ";
    } else if (!(*it).compare("CaNSR")) {
      os << state[index_CaNSR] << " ";
    } else if (!(*it).compare("CaSL")) {
      os << state[index_CaSL] << " ";
    } else if (!(*it).compare("LTRPNCa")) {
      os << state[index_LTRPNCa] << " ";
    } else if (!(*it).compare("HTRPNCa")) {
      os << state[index_HTRPNCa] << " ";

    } else if (!(*it).compare("CaSSavg")) {
      os << otherstates[index_CaSSavg] << " ";
    } else if (!(*it).compare("CaJSRavg")) {
      os << otherstates[index_CaJSRavg] << " ";
    } else if (!(*it).compare("JRyRtot")) {
      os << otherstates[index_JRyRtot] << " ";
    } else if (!(*it).compare("PRyR_Open")) {
      os << otherstates[index_PRyR_Open] << " ";
    } else if (!(*it).compare("PRyR_ready")) {
      os << otherstates[index_PRyR_ready] << " ";
    } else if (!(*it).compare("PNorm_Mode")) {
      os << otherstates[index_PNorm_Mode] << " ";
    } else if (!(*it).compare("PnotVinact")) {
      os << otherstates[index_PnotVinact] << " ";
    } else if (!(*it).compare("PLType_Open")) {
      os << otherstates[index_PLType_Open] << " ";
    } else if (!(*it).compare("CaTOT2")) {
      os << otherstates[index_CaTOT2] << " ";
    } else if (!(*it).compare("PIto2_Open")) {
      os << otherstates[index_PIto2_Open] << " ";
    } else if (!(*it).compare("CaJSRtot")) {
      os << otherstates[index_CaJSRtot] << " ";
    } else if (!(*it).compare("CaSStot")) {
      os << otherstates[index_CaSStot] << " ";
    } else if (!(*it).compare("CaSRtot")) {
      os << otherstates[index_CaSRtot] << " ";
    } else if (!(*it).compare("CaSLtot")) {
      os << otherstates[index_CaSLtot] << " ";

    } else if (!(*it).compare("INa")) {
      os << current[index_INa] << " ";
    } else if (!(*it).compare("IKr")) {
      os << current[index_IKr] << " ";
    } else if (!(*it).compare("IKs")) {
      os << current[index_IKs] << " ";
    } else if (!(*it).compare("Ito1")) {
      os << current[index_Ito1] << " ";
    } else if (!(*it).compare("IK1")) {
      os << current[index_IK1] << " ";
    } else if (!(*it).compare("IKp")) {
      os << current[index_IKp] << " ";
    } else if (!(*it).compare("INaCa")) {
      os << current[index_INaCa] << " ";
    } else if (!(*it).compare("INaK")) {
      os << current[index_INaK] << " ";
    } else if (!(*it).compare("IpCa")) {
      os << current[index_IpCa] << " ";
    } else if (!(*it).compare("ICab")) {
      os << current[index_ICab] << " ";
    } else if (!(*it).compare("INab")) {
      os << current[index_INab] << " ";
    } else if (!(*it).compare("ICa")) {
      os << current[index_ICa] << " ";
    } else if (!(*it).compare("JDHPR")) {
      os << current[index_JDHPR] << " ";
    } else if (!(*it).compare("Jup")) {
      os << current[index_Jup] << " ";
    } else if (!(*it).compare("Jtrpn")) {
      os << current[index_Jtrpn] << " ";
    } else if (!(*it).compare("Jtr")) {
      os << current[index_Jtr] << " ";
    } else if (!(*it).compare("Jxfer")) {
      os << current[index_Jxfer] << " ";
    } else if (!(*it).compare("IKv43")) {
      os << current[index_IKv43] << " ";
    } else if (!(*it).compare("IKv14")) {
      os << current[index_IKv14] << " ";
    } else if (!(*it).compare("IKv14_K")) {
      os << current[index_IKv14_K] << " ";
    } else if (!(*it).compare("IKv14_Na")) {
      os << current[index_IKv14_Na] << " ";
    } else if (!(*it).compare("Ito2")) {
      os << current[index_Ito2] << " ";
    } else if (!(*it).compare("Istim")) {
      os << current[index_Istim] << " ";
    } else if (!(*it).compare("Itot")) {
      os << current[index_Itot] << " ";
    } else {
      printf("Warning: unrecognized output variable: %s\n", (*it).c_str());
      os << "NA ";
    }
  }

}

void StochModel::Write_CRU_Header(std::ostream &os)
{

  os << "Time NLCC_open NRyR_open CaJSR CaSS" << std::endl;

}

void StochModel::Write_CRU(std::ostream &os, double t, int iFRU)
{

  os << t << " ";
  os << FRUs[iFRU].Get_LCC_Open() << " ";
  os << FRUs[iFRU].Get_RyR_Open() << " ";
  os << FRUs[iFRU].FRU_states[index_frustates_CaJSR] << " ";
  os << FRUs[iFRU].Get_CaSS_Avg() << " ";
  os << std::endl;

}
