#include "Matrix.h"
#include "ExceptionBase.h"

#include <sstream>
#include <cstring>
#include <cassert>

namespace ryr
{
namespace sparse_matrix
{
void
HypreMatrix::
create_matrix(unsigned lower_row, unsigned upper_row,
              unsigned lower_column, unsigned upper_column,
              int* n_cols_per_row,
              const int* matrix_cols, const double* matrix_vals,
              MPI_Comm comm)
{
  lower_row_ = lower_row;
  upper_row_ = upper_row;

  HYPRE_IJMatrixCreate(comm,
                       lower_row,    upper_row,
                       lower_column, upper_column,
                       &(matrix_ij_));

  // Choose a parallel csr format storage
  HYPRE_IJMatrixSetObjectType(matrix_ij_, HYPRE_PARCSR);
  HYPRE_IJMatrixInitialize(matrix_ij_);

  const unsigned hypre_nrow_local = upper_row - lower_row + 1;
  const unsigned hypre_first_row = lower_row;

  std::vector<int> row_map(hypre_nrow_local);

  for (unsigned ii = 0; ii < hypre_nrow_local; ii++) {
    row_map[ii] = hypre_first_row + ii;
  }

  HYPRE_IJMatrixSetValues(matrix_ij_,
                          hypre_nrow_local,
                          n_cols_per_row,
                          &row_map[0],
                          matrix_cols,
                          matrix_vals);

  HYPRE_IJMatrixAssemble(matrix_ij_);
  HYPRE_IJMatrixGetObject(matrix_ij_, (void **) & (matrix_));

  std::ostringstream message;
  int err = HYPRE_GetError();

  if (err) {
    throw exceptions::ExceptionBase() << "Hypre matrix setup failed";
  }
}

}

static	void print_vector(hypre_ParVector *s);

static
void print_vector(const std::vector<double>& s)
{
  unsigned n = s.size();

  for(unsigned ii = 0; ii < n; ii++) {
    std::cout << s[ii] << std::endl;
  }

}


static
void print_vector(hypre_ParVector *s)
{

  hypre_Vector* b_local_vector = hypre_ParVectorLocalVector((hypre_ParVector *)s);
  double* b_vector_data = hypre_VectorData(b_local_vector);

  unsigned n = b_local_vector->size;

  for(unsigned ii = 0; ii < n; ii++) {
    std::cout << b_vector_data[ii] << std::endl;
  }
}



namespace linear_solver
{
extern "C"
{

  HYPRE_Int
  block_setup
  (HYPRE_Solver solver, HYPRE_Matrix matrix, HYPRE_Vector b, HYPRE_Vector x)
  {
    return 0;
  }


  HYPRE_Int
  block_solve
  (HYPRE_Solver solver, HYPRE_Matrix matrix, HYPRE_Vector b, HYPRE_Vector x)
  {
    LinearSolver* bls = reinterpret_cast<LinearSolver*>(solver);
    bls->precond_solve(matrix, b, x);
    return 0;
  }
}



void LinearSolver::
precond_solve(HYPRE_Matrix matrix, HYPRE_Vector b, HYPRE_Vector x)
{

  //hypre_Vector* b_local_vector = hypre_ParVectorLocalVector((hypre_ParVector *)b);
  //double* b_vector_data = hypre_VectorData(b_local_vector);
  //hypre_Vector* x_local_vector = hypre_ParVectorLocalVector((hypre_ParVector *)x);
  //double* x_vector_data = hypre_VectorData(x_local_vector);

  //  std::memcpy(x_vector_data, b_vector_data, sizeof(double) * b_local_vector->size);
  //	HYPRE_PCGSolve(cg_solver_, (HYPRE_Matrix) matrix->matrix(), (HYPRE_Vector)b_vector_, (HYPRE_Vector)hypre_vector_);
  //print_vector((hypre_ParVector *)b);
  HYPRE_ParaSailsSolve(paras_, (HYPRE_ParCSRMatrix)matrix, (HYPRE_ParVector)b, (HYPRE_ParVector)x);
}


void LinearSolver::
setup_SPAI(sparse_matrix::HypreMatrix& matrix, double threshold, int n_levels)
{

  MPI_Comm comm = comm_->comm();
  HYPRE_ParaSailsCreate(comm, &paras_);
  HYPRE_ParaSailsSetParams(paras_, threshold, n_levels);

  int sym = 1;
  HYPRE_ParaSailsSetSym(paras_, sym);
  HYPRE_ParaSailsSetup(paras_,
                       (HYPRE_ParCSRMatrix) matrix.matrix(),
                       (HYPRE_ParVector) hypre_vector_,
                       (HYPRE_ParVector) hypre_vector_);

  //HYPRE_IJMatrix inv_matrix;
  //HYPRE_ParaSailsBuildIJMatrix(paras, &inv_matrix);
  //HYPRE_IJMatrixGetObject(inv_matrix, (void **) &A_inv_);

}

void LinearSolver::
setup_CG(sparse_matrix::HypreMatrix& matrix)
{

  const double tolerance = 0.00001;
  const double absolute_tolerance = 0.000001;
  const double krylov_print_level = 0;
  const unsigned max_iter = 1000;

  HYPRE_ParCSRPCGCreate(comm_->comm(), &cg_solver_);
  //HYPRE_PCGSetTol(cg_solver_, tolerance);
  HYPRE_PCGSetAbsoluteTol(cg_solver_, absolute_tolerance);
  HYPRE_PCGSetLogging(cg_solver_, 0);
  HYPRE_PCGSetPrintLevel(cg_solver_, krylov_print_level);
  HYPRE_PCGSetMaxIter(cg_solver_, max_iter);
  HYPRE_PCGSetPrecond(cg_solver_,
                      (HYPRE_PtrToSolverFcn) block_solve,
                      (HYPRE_PtrToSolverFcn) block_setup,
                      reinterpret_cast<HYPRE_Solver>(this));

  HYPRE_PCGSetup(cg_solver_,
                 (HYPRE_Matrix) matrix.matrix(),
                 (HYPRE_Vector) hypre_vector_,
                 (HYPRE_Vector) hypre_vector_);
}


void LinearSolver::
clean_memory()
{
  if(cg_solver_) {
    hypre_ParVectorDestroy(b_vector_);
    hypre_ParVectorDestroy(hypre_vector_);
    HYPRE_ParCSRPCGDestroy(cg_solver_);
    HYPRE_ParaSailsDestroy(paras_);
  }

}

void LinearSolver::
set_matrix(sparse_matrix::HypreMatrix& matrix)
{

  //matrix.print();
  clean_memory();
  hypre_vector_ = matrix.get_par_vector();
  b_vector_ = matrix.get_par_vector();
  this->setup_CG(matrix);
  double threshold = 0.5;
  int n_levels = 0;
  this->setup_SPAI(matrix, threshold, n_levels);

}


void LinearSolver::
solve(sparse_matrix::HypreMatrix& matrix,
      const std::vector<double>& b, std::vector<double>& x)
{



  hypre_Vector* b_local_vector = hypre_ParVectorLocalVector((hypre_ParVector *)b_vector_);
  double* b_vector_data = hypre_VectorData(b_local_vector);
  std::memcpy(b_vector_data, &b[0], sizeof(double) * b_local_vector->size);

  HYPRE_PCGSolve(cg_solver_, (HYPRE_Matrix) matrix.matrix(), (HYPRE_Vector)b_vector_, (HYPRE_Vector)hypre_vector_);

  hypre_Vector* x_local_vector = hypre_ParVectorLocalVector((hypre_ParVector *)hypre_vector_);
  double* x_vector_data = hypre_VectorData(x_local_vector);
  std::memcpy(&x[0], x_vector_data, sizeof(double) * x_local_vector->size);
}


}
}
