#include "Cell.h"
//#include <cstdlib>
#include <ostream>
#include <iostream>
#include <cmath>

namespace ryr
{

namespace cell
{

void MainCellBase::
get_n_matrix_entry(std::vector<unsigned>& index)
{

  index.clear();
  unsigned n_adj = adj_cell_container_.size();
  index.resize(n_adj + 1);

  index[0] = this->matrix_index();
  unsigned count = 1;

  for(std::list<std::pair<CellBase*, double> >::const_iterator it = adj_cell_container_.begin();
      it != adj_cell_container_.end(); it++) {
    index[count++] = it->first->matrix_index();
  }

  //	exit(0);
}


void MainCellBase::
update_calcium()
{
  model_->get_calcium() = Ca_;
}


void MainCellBase::
update_calcium_from_model()
{
  Ca_ = model_->get_calcium();
}


void MainCellBase::
calc_calcium_explicit_method(double resistance_factor, double& dCa, double& dv)
{

  if(!model_->share_calcium()) {
    dCa = 0;
    dv = 0;
    return;
  }

  double result = 0;
  //double result2 = 0;

  for(std::list<std::pair<CellBase*, double> >::const_iterator it = adj_cell_container_.begin();
      it != adj_cell_container_.end(); it++) {
    if(it->first->share_calcium()) {
      result += ((it->first->v() - v_) - 13.3541 * std::log(Ca_ / it->first->Ca())) / resistance_factor / it->second;
      //result2 += ( - 13.3541 * std::log(Ca_ / it->first->Ca())) / resistance_factor / it->second;
    }
  }

  dCa = result / (2 * 25.84 * 96.485 * 1000.0); //  / (z * Vcell * F)
  dv = result * 1e-9 / 153.4; //result2;

}


void MainCellBase::
calc_matrix_entry_operator_spliting(double dt, unsigned offset, std::vector<double>& mat_value)
{

  unsigned index = 1;
  mat_value[offset] = 1.0 / dt;

  for(std::list<std::pair<CellBase*, double> >::const_iterator it = adj_cell_container_.begin();
      it != adj_cell_container_.end(); ++it) {
    mat_value[offset + index] = -0.5 / it->second;
    mat_value[offset] += 0.5 / it->second;
    index++;
  }

}

void MainCellBase::
calc_rhs_operator_spliting(double dt, double& value)
{

  value = 1.0 / dt;

  for(std::list<std::pair<CellBase*, double> >::const_iterator it = adj_cell_container_.begin();
      it != adj_cell_container_.end(); ++it) {
    value -= 0.5 / it->second;
  }

  value *= this->v();

  for(std::list<std::pair<CellBase*, double> >::const_iterator it = adj_cell_container_.begin();
      it != adj_cell_container_.end(); ++it) {
    value += 0.5 / it->second * it->first->v();
  }

}


}

}
