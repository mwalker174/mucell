#ifndef _____CELL______RYR
#define _____CELL______RYR

#include <list>
#include <vector>
#include <set>
#include <iostream>

#include "models/TT06_RRG.h"


namespace ryr
{

namespace cell
{
template <class C>
class Cell : public C
{

public:
  void init_cell(double dt) {
    this->init_cell_base(dt);
  }

  double& dv() {
    return this->dv_base();
  }

};

class CellBase
{
public:
  unsigned& matrix_index() {
    return matrix_index_;
  }

  unsigned& original_index() {
    return original_index_;
  }

  double& v() {
    return v_;
  }

  double& Ca() {
    return Ca_;
  }

  double& dv_base() {
    return dv_;
  }

  bool& share_calcium() {
    return share_calcium_;
  }

protected:
  unsigned param_id_;
  double v_, dv_, Ca_;
  bool share_calcium_;

private:
  unsigned original_index_, matrix_index_;
};


class MainCellBase : public CellBase
{

protected:

  double calcalate_dv(double dt) {
    return 0;
  }

public:

  //output
  void set_output(const std::vector<std::string>& save) {
    model_->set_output(save);
  }

  void output_header(std::ostream& file) {
    model_->output_header(file);
  }

  void output_data(std::ostream& file) {
    model_->output_data(file);
  }

  unsigned& get_param_id() {
    return param_id_;
  }

  MainCellBase() {
  }

  void clone_model(IonicModel* model) {

    model_ = model->clone();
    v_ = model_->get_default_voltage();
    Ca_ = model_->get_calcium();
    share_calcium_ = model_->share_calcium();
    model_->set_id(this->original_index());
  }

  void add_cell(CellBase* cell, double res) {
    adj_cell_container_.push_back(std::make_pair(cell, res));
  }

  void integrate(double time, double dt) {
    v_ += model_->integrate(time, dt, v_) * dt;
  }

  void get_n_matrix_entry(std::vector<unsigned>& e);
  void calc_matrix_entry_operator_spliting(double dt, unsigned offset, std::vector<double>& mat_value);
  void calc_rhs_operator_spliting(double dt, double& value);
  void integrate() {
  }

  void write_restart(std::ostream& o) {
    o << v_ << " ";
    model_->write_restart(o);
  }

  void read_restart(std::istream& i) {
    i >> v_;
    model_->read_restart(i);
  }

  void update_calcium();
  void update_calcium_from_model();
  void calc_calcium_explicit_method(double resistance_factor, double& dCa, double& dv);

  void set_id(unsigned cell_id) {
    model_->set_id(cell_id);
  }

  void set_resistance() {
    double cap = model_->get_capacitance();

    for(std::list<std::pair<CellBase*, double> >::iterator it = adj_cell_container_.begin();
        it != adj_cell_container_.end(); it++) {
      it->second *= cap;
    }
  }

private:
  std::list<std::pair<CellBase*, double> > adj_cell_container_;
  //  std::set<unsigned> matrix_index_;
  //std::vector<unsigned> cell_entry_index_;
  IonicModel *model_;
};

class HaloCellBase : public CellBase
{
protected:

private:

};

typedef Cell<MainCellBase> MainCell;
typedef Cell<HaloCellBase> HaloCell;

}

}


#endif
