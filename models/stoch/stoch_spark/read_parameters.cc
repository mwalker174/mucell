#include "StochModel.h"
#include "spirit_wrapper.h"

void StochModel::Read_Parameters(std::string filename)
{

  try {
    spirit_wrapper::SpiritReader reader;
    reader.read_file(filename);
    const spirit_wrapper::SpiritReaderObject& object = reader.get_object("stoch");

    object.get_prop("Default_Cai", Default_Cai);
    object.get_prop("Default_CaSR", Default_CaSR);
    object.get_prop("Default_Nai", Default_Nai);
    object.get_prop("beta_flag", beta_flag);
    object.get_prop("mode2_frac", mode2_frac);
    object.get_prop("ryr_kplus_scale", ryr_kplus_scale);
    object.get_prop("serca_kf_scale", serca_kf_scale);
    object.get_prop("ikr_scale", ikr_scale);
    object.get_prop("serca_vmax_scale", serca_vmax_scale);
    object.get_prop("ik1_scale", ik1_scale);
    object.get_prop("ncx_scale", ncx_scale);
    object.get_prop("ito1_scale", ito1_scale);
    object.get_prop("frac_orphan", frac_orphan);
    object.get_prop("frac_ncx_sl", frac_ncx_sl);
    object.get_prop("Acap", Acap);
    object.get_prop("Cao", Cao);
    object.get_prop("nak_scale", nak_scale);
    object.get_prop("NFRU_X", NFRU_X);
    object.get_prop("NFRU_Y", NFRU_Y);
    object.get_prop("NFRU_Z", NFRU_Z);
    object.get_prop("clamp_nai", clamp_nai);
  
    object.get_prop("t_final", t_final);
    object.get_prop("period_stim", period_stim);
    object.get_prop("t_end_stim", t_end_stim);
    object.get_prop("I_stim", I_stim);

    object.get_prop("VClamp_Flag", VClamp_Flag);
    object.get_prop("VClamp_ClampV", VClamp_ClampV);
    object.get_prop("VClamp_HoldV", VClamp_HoldV);
    object.get_prop("VClamp_Freq", VClamp_Freq);
    object.get_prop("VClamp_Duration", VClamp_Duration);

    object.get_prop("VClamp2_T1", VClamp2_T1);
    object.get_prop("VClamp2_T2", VClamp2_T2);
    object.get_prop("VClamp2_ClampV", VClamp2_ClampV);
    
    NFRU = NFRU_X*NFRU_Y*NFRU_Z;

    NFRU_scale = ((double)NFRU_total)/((double)NFRU);

    double temp;
    object.get_prop("rand_seed", temp);
    rand_seed = (unsigned long)temp;

  } catch(exceptions::ExceptionBase e) {
    std::cout << e.get_message() << std::endl;
  }

  Initialize_ReleaseUnit_Params();

}

