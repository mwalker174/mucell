#include "parser/spirit_wrapper.h"
#include "Problem.h"
#include "helpers/Helpers.h"

int main(int argc, char **argv)
{

  ryr::mpi::MPIComm comm;
  comm.init(&argc, &argv);
  ryr::info::info.init(comm);

  try {

    ryr::Problem p;

    p.init(&comm);
    p.solve_operator_splitting();

  } catch(exceptions::ExceptionBase e) {
    std::cout << e.get_message() << std::endl;
  }

  comm.finalize();
}

