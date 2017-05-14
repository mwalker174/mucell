#include "Cell.h"
#include "helpers/Helpers.h"
#include "Stimul.h"
#include "models/ModelFactory.h"
#include "ExceptionBase.h"

#include <cassert>
#include <vector>
#include <map>


namespace ryr
{

namespace mesh
{
class Mesh
{
public:
  typedef std::vector<cell::HaloCell*> HaloVector;
  typedef std::vector<cell::MainCell*> HaloedVector;
  typedef std::map<unsigned, HaloVector > HaloContainer;
  typedef std::map<unsigned, HaloedVector > HaloedContainer;

  void init(mpi::MPIComm* comm, double res) {
    comm_ = comm;
    res_ = res;
  }

  void create_matrix_index();
  void exchange_data();
  void exchange_calcium_info_data();
  void read_from_file(std::string file_name);

  void fill_stimuli_cell_pointers(std::list<Stimul>& stim);

  unsigned get_first_dof() {
    return bump_;
  }

  unsigned get_n_cell() {
    return cells_.size();
  }

  std::vector<cell::MainCell*>& get_cells() {
    return cells_;
  }

  void init_cell_params(const std::vector<IonicModel*>& clones) {

    unsigned n_cells = cells_.size();

    for(unsigned ii = 0; ii < n_cells; ii++) {
      unsigned id = cells_[ii]->get_param_id();

      if(id >= clones.size()) {
        throw exceptions::WrongModelId() << "Id " << id << " is outside model types";
      }

      unsigned cell_id = ii + comm_->rank() * 1e4; //todo
      cells_[ii]->clone_model(clones[id]);
      cells_[ii]->set_resistance();
    }

  }

  void write_restart(std::ostream& o) {

    unsigned n_cells = cells_.size();

    for(unsigned ii = 0; ii < n_cells; ii++) {
      o << cells_[ii]->original_index() << std::endl;
      std::ostringstream s;
      cells_[ii]->write_restart(s);
      long long pos = s.tellp();
      o << pos << " ";
      o.write((char*)s.str().c_str(), pos);
    }

  }

  void read_restart(std::istream& i);


  void find_probe(const std::vector<unsigned>& cell_ids, std::list<cell::MainCell*>& out);

  void update_calcium_from_model() {

    unsigned n_cell = this->get_n_cell();

    for(unsigned ii = 0; ii < n_cell; ii++) {
      if(cells_[ii]->share_calcium()) {
        cells_[ii]->update_calcium_from_model();
      }
    }

  }

private:
  std::vector<cell::MainCell*> cells_;
  std::vector<cell::HaloCell*> halo_cells_;

  HaloedContainer haloed_container_;
  HaloContainer halo_container_;

  mpi::MPIComm* comm_;

  unsigned bump_;
  double res_;

};
}

}
