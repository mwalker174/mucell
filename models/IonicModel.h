#ifndef __IONIC_MODEL_RYR__
#define __IONIC_MODEL_RYR__

#include <map>
#include <string>
#include <vector>

namespace ryr
{


template<class T>
class ModelOutput
{

public:

  virtual void output_header(std::ostream& file) = 0;
  virtual void output_data(std::ostream& file, T& model) = 0;

};

class IonicModel
{

public:
  virtual void read_parameters(std::map<std::string, double>& parameters) = 0;

  virtual void set_id(unsigned cell_id) {
  }

  virtual double integrate(double time, double dt, double V) = 0;
  virtual double get_default_voltage() = 0;
  virtual IonicModel* clone() = 0;
  virtual void write_restart(std::ostream& o) = 0;
  virtual void read_restart(std::istream& i) = 0;
  virtual bool share_calcium() const = 0;
  virtual double& get_calcium() = 0;
  virtual double& get_v() = 0;
  virtual double get_capacitance() = 0;

  //output
  virtual void set_output(const std::vector<std::string>& save) = 0;
  virtual void output_header(std::ostream& file) = 0;
  virtual void output_data(std::ostream& file) = 0;

};


class IonicFactory
{
public:
  virtual IonicModel* create_model() = 0;
};


}

#endif
