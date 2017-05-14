#ifndef TT06_RRG_HH
#define TT06_RRG_HH

#include "IonicModel.h"

#include <string>
#include <vector>
#include <iostream>

class TT06_RRG
{
public:

  TT06_RRG() {
    states_.resize(20);
  }
  TT06_RRG(int cellType);
  void init(int cellType) {
    initConsts(cellType);
    initStates(cellType);
    defaultVoltage_ = states_[0];
  }
  double defaultVoltage();
  double computeRates(double dt, double V, double iStim, bool integrate = false);

  void write_restart(std::ostream& o) {
    for(unsigned ii = 0; ii < 20; ii++) {
      o << states_[ii] << " ";
    }
  }


  void read_restart(std::istream& i) {
    for(unsigned ii = 0; ii < 20; ii++) {
      i >> states_[ii];
    }
  }


  double& get_calcium() {
    return states_[3];
  }

  double& get_v() {
    return states_[0];
  }

  double get_capacitance() {
    return 3.2808e-5 * 2.0; //uF
  }

  double get_casr() {
    return states_[17];
  }


private:

  void initConsts(int cellType);
  void initStates(int cellType);



  static double constants_[53];
  int s_switch_;
  double defaultVoltage_;
  double g_Ks_;  // formerly CONSTANTS[15]
  double g_Kr_;  // formerly CONSTANTS[14]
  double g_to_;  // formerly CONSTANTS[20]
  double P_NaK_; // formerly CONSTANTS[21]
  double g_NaL_;
  std::vector<double> states_;
};


namespace ryr
{
class TT06 : public IonicModel
{
public:

  virtual void read_parameters(std::map<std::string, double>& parameters) {
    unsigned id = parameters["id"];
    model_.init(id);
  }

  /// returns dv/dt
  virtual double integrate(double time, double dt, double V) {
    return model_.computeRates(dt, V, 0.0, true);
  }

  virtual double get_default_voltage() {
    return model_.defaultVoltage();
  }

  virtual IonicModel* clone() {

    TT06* new_model = new TT06();
    *new_model = *this;
    return new_model;

  }

  void write_restart(std::ostream& o) {
    model_.write_restart(o);
  }


  void read_restart(std::istream& i) {
    model_.read_restart(i);
  }

  bool share_calcium() const {
    return true;
  }

  double& get_calcium() {
    return model_.get_calcium();
  }

  double& get_v() {
    return model_.get_v();
  }

  double get_capacitance() {
    return model_.get_capacitance();
  }

  double get_casr() {
    return model_.get_casr();
  }

  //output
  virtual void set_output(const std::vector<std::string>& save);
  virtual void output_header(std::ostream& file);
  virtual void output_data(std::ostream& file);


private:
  TT06_RRG model_;
  std::vector<ModelOutput<TT06>*> output_vector_;
};

class TT06Factory : public IonicFactory
{
public:
  virtual IonicModel* create_model() {
    return new TT06();
  }
};


}


#endif
