#ifndef __PROBLEM__RYR
#define __PROBLEM__RYR

#include "Mesh.h"
#include "helpers/Helpers.h"
#include "Stimul.h"
#include "Matrix.h"
#include "models/ModelFactory.h"
#include "pio/PIO.h"

#include <fstream>

namespace ryr
{

class Problem
{

public:

  Problem(): matrix_(NULL), time_(0), linear_solver_(NULL) {}
  virtual ~Problem() {

    delete matrix_;
    delete linear_solver_;
  }

  void init(mpi::MPIComm* comm);
  void solve_operator_splitting();

protected:

  void create_matrix_map();

  void create_matrix_operator_splitting();
  void iteration_operator_splitting();
  void calc_rhs_operator_splitting(std::vector<double>& rhs);

  void process_stimuli();

  void integrate_ionic_models();
  //  void init_ionic_models();

  void add_ionic_models();

  void save_to_disc();

  void mpi_read();
  void mpi_write(double time);
  void read_restart(std::istream& file);
  void write_restart(std::ostream& file);
  void update_calcium();

private:

  void close_probe_files();

  factory::ModelFactory factory_;

  std::string probe_prefix_;
  std::list<cell::MainCell*> probe_cells_;
  std::list<std::ofstream*> probe_files_;
  double save_dt_;



  mesh::Mesh mesh_;
  std::list<Stimul> stimul_;
  double dt_, duration_, time_;


  sparse_matrix::HypreMatrix *matrix_;
  mpi::MPIComm* comm_;

  std::vector<int> n_cols_per_row_;
  std::vector<int> matrix_cols_;
  std::vector<double> matrix_vals_;

  linear_solver::LinearSolver* linear_solver_;


  //restart
  bool save_restart_, read_restart_;
  double restart_interval_, restart_start_, restart_end_;
  std::string restart_file_name_, restart_prefix_;
  double calcium_resistance_factor_;


	bool single_output_;
	pio::PIO pio;

};


}
#endif
