#include "parser/spirit_wrapper.h"
#include "PIO.h"
#include "helpers/Helpers.h"

int main(int argc, char **argv)
{

  ryr::mpi::MPIComm comm;
  comm.init(&argc, &argv);
  ryr::info::info.init(comm);

  try {

		pio::PIO p; p.set_file_name("test.txt"); p.init(comm.comm());
		for(unsigned ii = 0; ii < 10; ii++){
			p.get_stream() << ii + 10 * comm.rank() << "T";
			p.next(ii + 10 * comm.rank());
		} 
		p.commit();
  
  } catch(exceptions::ExceptionBase e) {
    std::cout << e.get_message() << std::endl;
  }

  comm.finalize();
}

