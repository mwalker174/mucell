#include "stoch.h"
#include "ExceptionBase.h"

namespace ryr
{

void Stoch::set_output(const std::vector<std::string>& save)
{

  ModelOutput<Stoch>* container[] = {new CaOutput(), new SRLoadOutput(), new CaSLOutput(), new LCCVinactOutput(), new LCCModeCaOutput(), new LCCMode2OpenOutput(), new CaSSOutput(), new CaJSROutput(), new mNaOutput(), new jNaOutput(), new hNaOutput(), new CaTotOutput(), new RyROpenOutput(), new JRyROutput(), new INaCaOutput(), new Ito2Output(), new NaiOutput()}; //memory leak so far
  std::map<std::string, ModelOutput<Stoch>*> m;
  m["Ca_i"] = container[0];
  m["SR_Load"] = container[1];
  m["CaSL"] = container[2];
  m["LCC_Vinact"] = container[3];
  m["LCC_ModeCa"] = container[4];
  m["LCC_Mode2_Open"] = container[5];
  m["CaSS"] = container[6];
  m["CaJSR"] = container[7];
  m["mNa"] = container[8];
  m["jNa"] = container[9];
  m["hNa"] = container[10];
  m["Ca_Tot"] = container[11];
  m["RyR_Open"] = container[12];
  m["JRyR"] = container[13];
  m["INaCa"] = container[14];
  m["Ito2"] = container[15];
  m["Nai"] = container[16];

  unsigned n_save = save.size();

  for(unsigned ii = 0; ii < n_save; ii++) {
    std::map<std::string, ModelOutput<Stoch>*>::iterator it = m.find(save[ii]);

    if(it == m.end()) {
      throw exceptions::WrongOption() << "Saving option " << save[ii]
                                      << " is not defined for Stoch model!";
    }

    output_vector_.push_back(it->second);
  }

}

void Stoch::output_header(std::ostream& file)
{

  unsigned n_output_vector = output_vector_.size();

  for(unsigned ii = 0; ii < n_output_vector; ii++) {
    output_vector_[ii]->output_header(file);
    file << " ";
  }

}

void Stoch::output_data(std::ostream& file)
{

  unsigned n_output_vector = output_vector_.size();

  for(unsigned ii = 0; ii < n_output_vector; ii++) {
    output_vector_[ii]->output_data(file, *this);
    file << " ";
  }
}

}
