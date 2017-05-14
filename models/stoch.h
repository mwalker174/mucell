#ifndef STOCH_ION_HH
#define STOCH_ION_HH

#include "IonicModel.h"
#include "stoch/stoch_spark/StochModel.h"

#include <string>
#include <vector>

namespace ryr
{

class Stoch : public IonicModel
{
public:

  virtual void read_parameters(std::map<std::string, double>& parameters) {

    model_.Set_Parameters(parameters);
    model_.Set_Seed(model_.Get_Seed());
    model_.Initialize_Default_State();
    model_.Initialize_Currents(0, 0);

  }

  virtual void set_id(unsigned cell_id) {
    unsigned seed = model_.Get_Seed() + cell_id * cell_id * 7 + cell_id * 92748;
    model_.Set_Seed(seed);
  }

  /// returns dv/dt
  virtual double integrate(double time, double dt, double V) {
    double dV;
    model_.Delta_V_Step(dt, V, dV);
    return dV / dt;
  }

  virtual double get_default_voltage() {
    return model_.Get_V();
  }

  virtual IonicModel* clone() {

    Stoch* new_model = new Stoch();
    *new_model = *this;
    return new_model;

  }

  void write_restart(std::ostream& o) {
    model_.Write_State(o);
  }
  void read_restart(std::istream& i) {
    model_.Read_State(i);
  }

  bool share_calcium() const {
    return false; //true; //Enable intercellular Ca2+ diffusion
  }

  double& get_calcium() {
    return model_.Get_Cai();
  }

  double& get_cai() {
    return model_.Get_Cai();
  }

  virtual double& get_v() {
    return model_.Get_V();
  }

  double get_capacitance() {
    return model_.Get_Capacitance();
  }

  double get_sr_load() {
    return model_.Get_SR_Load();
  }

  double get_casl() {
    return model_.Get_CaSL();
  }

  double get_lcc_vinact() {
    return model_.Get_LCC_Vinact();
  }

  double get_lcc_modeca() {
    return model_.Get_LCC_ModeCa();
  }

  double get_lcc_mode2_open() {
    return model_.Get_LCC_Mode2_Open();
  }

  double get_cass() {
    return model_.Get_CaSS();
  }

  double get_cajsr() {
    return model_.Get_CaJSR();
  }

  double get_mna() {
    return model_.Get_mNa();
  }

  double get_jna() {
    return model_.Get_jNa();
  }

  double get_hna() {
    return model_.Get_hNa();
  }

  double get_ca_tot() {
    return model_.Get_Ca_Tot();
  }

  double get_ryr_open() {
    return model_.Get_RyR_Open();
  }

  double get_jryr() {
    return model_.Get_JRyR();
  }

  double get_inaca() {
    return model_.Get_INaCa();
  }

  double get_ito2() {
    return model_.Get_Ito2();
  }

  double get_nai() {
    return model_.Get_Nai();
  }

  //output
  virtual void set_output(const std::vector<std::string>& save);
  virtual void output_header(std::ostream& file);
  virtual void output_data(std::ostream& file);


private:
  StochModel model_;
  std::vector<ModelOutput<Stoch>*> output_vector_;
};



class CaOutput : public ModelOutput<Stoch>
{

public:

  virtual void output_header(std::ostream& file) {
    file << "Ca_i";
  }

  virtual void output_data(std::ostream& file, Stoch& model) {
    file << model.get_cai();
  }

};

class SRLoadOutput : public ModelOutput<Stoch>
{

public:

  virtual void output_header(std::ostream& file) {
    file << "SR_Load";
  }

  virtual void output_data(std::ostream& file, Stoch& model) {
    file << model.get_sr_load();
  }

};

class CaSLOutput : public ModelOutput<Stoch>
{

public:

  virtual void output_header(std::ostream& file) {
    file << "CaSL";
  }

  virtual void output_data(std::ostream& file, Stoch& model) {
    file << model.get_casl();
  }

};

class LCCVinactOutput : public ModelOutput<Stoch>
{

public:

  virtual void output_header(std::ostream& file) {
    file << "LCC_Vinact";
  }

  virtual void output_data(std::ostream& file, Stoch& model) {
    file << model.get_lcc_vinact();
  }

};

class LCCModeCaOutput : public ModelOutput<Stoch>
{

public:

  virtual void output_header(std::ostream& file) {
    file << "LCC_ModeCa";
  }

  virtual void output_data(std::ostream& file, Stoch& model) {
    file << model.get_lcc_modeca();
  }

};

class LCCMode2OpenOutput : public ModelOutput<Stoch>
{

public:

  virtual void output_header(std::ostream& file) {
    file << "LCC_Mode2_Open";
  }

  virtual void output_data(std::ostream& file, Stoch& model) {
    file << model.get_lcc_mode2_open();
  }

};

class CaSSOutput : public ModelOutput<Stoch>
{

public:

  virtual void output_header(std::ostream& file) {
    file << "CaSS";
  }

  virtual void output_data(std::ostream& file, Stoch& model) {
    file << model.get_cass();
  }

};

class CaJSROutput : public ModelOutput<Stoch>
{

public:

  virtual void output_header(std::ostream& file) {
    file << "CaJSR";
  }

  virtual void output_data(std::ostream& file, Stoch& model) {
    file << model.get_cajsr();
  }

};

class mNaOutput : public ModelOutput<Stoch>
{

public:

  virtual void output_header(std::ostream& file) {
    file << "mNa";
  }

  virtual void output_data(std::ostream& file, Stoch& model) {
    file << model.get_mna();
  }

};

class jNaOutput : public ModelOutput<Stoch>
{

public:

  virtual void output_header(std::ostream& file) {
    file << "jNa";
  }

  virtual void output_data(std::ostream& file, Stoch& model) {
    file << model.get_jna();
  }

};

class hNaOutput : public ModelOutput<Stoch>
{

public:

  virtual void output_header(std::ostream& file) {
    file << "hNa";
  }

  virtual void output_data(std::ostream& file, Stoch& model) {
    file << model.get_hna();
  }

};

class CaTotOutput : public ModelOutput<Stoch>
{

public:

  virtual void output_header(std::ostream& file) {
    file << "Ca_Tot";
  }

  virtual void output_data(std::ostream& file, Stoch& model) {
    file << model.get_ca_tot();
  }

};

class JRyROutput : public ModelOutput<Stoch>
{

public:

  virtual void output_header(std::ostream& file) {
    file << "JRyR";
  }

  virtual void output_data(std::ostream& file, Stoch& model) {
    file << model.get_jryr();
  }

};

class INaCaOutput : public ModelOutput<Stoch>
{

public:

  virtual void output_header(std::ostream& file) {
    file << "INaCa";
  }

  virtual void output_data(std::ostream& file, Stoch& model) {
    file << model.get_inaca();
  }

};

class Ito2Output : public ModelOutput<Stoch>
{

public:

  virtual void output_header(std::ostream& file) {
    file << "Ito2";
  }

  virtual void output_data(std::ostream& file, Stoch& model) {
    file << model.get_ito2();
  }

};

class RyROpenOutput : public ModelOutput<Stoch>
{

public:

  virtual void output_header(std::ostream& file) {
    file << "RyR_Open";
  }

  virtual void output_data(std::ostream& file, Stoch& model) {
    file << model.get_ryr_open();
  }

};

class NaiOutput : public ModelOutput<Stoch>
{

public:

  virtual void output_header(std::ostream& file) {
    file << "Nai";
  }

  virtual void output_data(std::ostream& file, Stoch& model) {
    file << model.get_nai();
  }

};

class StochFactory : public IonicFactory
{
public:
  virtual IonicModel* create_model() {
    return new Stoch();
  }
};

}


#endif
