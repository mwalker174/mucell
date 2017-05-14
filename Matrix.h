#ifndef __MATRIX_LINEAR_SOLVER__
#define __MATRIX_LINEAR_SOLVER__

#include "helpers/Helpers.h"

#include "_hypre_utilities.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_krylov.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_mv.h"
#include "_hypre_parcsr_ls.h"
#include "_hypre_parcsr_mv.h"


#include <vector>

namespace ryr
{
namespace sparse_matrix
{
class HypreMatrix
{
public:
  HypreMatrix(): matrix_ij_(NULL),
    lower_row_(0),
    upper_row_(0) {
  }

  virtual ~HypreMatrix() {
    if(matrix_ij_)
      HYPRE_IJMatrixDestroy(matrix_ij_);
  }

  void clean() {
    if(matrix_ij_)
      HYPRE_IJMatrixDestroy(matrix_ij_);
  }

  void create_matrix(unsigned lower_row, unsigned upper_row,
                     unsigned lower_column, unsigned upper_column,
                     int* n_cols_per_row,
                     const int* matrix_cols, const double* matrix_vals,
                     MPI_Comm comm);


  hypre_ParCSRMatrix* matrix() {
    return (hypre_ParCSRMatrix*)matrix_;
  }

  const unsigned& get_first_row() const {
    return lower_row_;
  }

  const unsigned& get_last_row() const {
    return upper_row_;
  }

  hypre_ParVector* get_par_vector() {
    return hypre_ParVectorInRangeOf((hypre_ParCSRMatrix *)matrix_);
  }

  void
  print() {
    HYPRE_IJMatrixPrint(matrix_ij_, "matrix.txt");
  }


private:
  HYPRE_ParCSRMatrix matrix_;
  HYPRE_IJMatrix matrix_ij_;
  unsigned lower_row_;
  unsigned upper_row_;
};
}



namespace linear_solver
{
class LinearSolver
{
public:
  LinearSolver(mpi::MPIComm* comm): cg_solver_(NULL), paras_(NULL) {
    comm_ = comm;
  }
  virtual ~LinearSolver() {
    clean_memory();
  }

  void  set_matrix(sparse_matrix::HypreMatrix& matrix);
  void solve(sparse_matrix::HypreMatrix& matrix,
             const std::vector<double>& b, std::vector<double>& x);
  void precond_solve(HYPRE_Matrix matrix, HYPRE_Vector b, HYPRE_Vector x);
protected:

  void clean_memory();
  void solve_block(HYPRE_Matrix matrix, HYPRE_Vector b, HYPRE_Vector x);
  void setup_SPAI(sparse_matrix::HypreMatrix& matrix, double threshold, int n_levels);
  void setup_CG(sparse_matrix::HypreMatrix& matrix);

private:
  mpi::MPIComm* comm_;
  HYPRE_Solver cg_solver_;
  HYPRE_Solver paras_;
  hypre_ParVector *hypre_vector_, *b_vector_;
};


}


}


#endif
