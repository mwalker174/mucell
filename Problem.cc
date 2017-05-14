#include "Problem.h"
#include "parser/spirit_wrapper.h"
#include "mpi_io_stream.h"
#include "models/stoch.h"

#include <iomanip>
#include <fstream>

using spirit_wrapper::SpiritReaderObject;

namespace ryr
{


void
Problem::
read_restart(std::istream& file)
{
  file >> time_;
  mesh_.read_restart(file);
}


void
Problem::
write_restart(std::ostream& file)
{
  if(comm_->rank() == 0) {
    file << time_ << " ";
  }

  mesh_.write_restart(file);
}



void
Problem::
mpi_read()
{

  MPI_Comm comm = comm_->comm();
  //input_output::mpi_istream some_file(comm, restart_file_name_);
  std::ifstream some_file(restart_file_name_.c_str());
  this->read_restart(some_file);
  some_file.close();

}


void Problem::
mpi_write(double time)
{

  std::stringstream s;
  s << std::setprecision(9);
  this->write_restart(s);
  long pos = s.tellp();
  char name[1000];
  sprintf(name, "%g", time);
  MPI_Comm comm = comm_->comm();
  unsigned long offset = input_output::mpi_io_streambuffer::get_writing_offset(comm, pos);
  input_output::mpi_ostream some_file(comm, restart_prefix_ + name);
  some_file.disable_delimiter();
  some_file.set_write_position(offset);
  some_file << s.rdbuf();
  some_file.close();

}



void Problem::
create_matrix_map()
{
  unsigned n_cell = mesh_.get_n_cell();
  n_cols_per_row_.resize(n_cell + 1, 0);
  std::vector<cell::MainCell*>& cells = mesh_.get_cells();

  for(unsigned ii = 0; ii < n_cell; ii++) {
    std::vector<unsigned> e;
    cells[ii]->get_n_matrix_entry(e);
    n_cols_per_row_[ii] = e.size();
    matrix_cols_.insert(matrix_cols_.end(), e.begin(), e.end());
  }

  matrix_vals_.resize(matrix_cols_.size());
}


void Problem::
update_calcium()
{

  unsigned n_cell = mesh_.get_n_cell();
  std::vector<cell::MainCell*>& cells = mesh_.get_cells();
  std::vector<double> dCa(n_cell);
  std::vector<double> dv(n_cell);

  for(unsigned ii = 0; ii < n_cell; ii++) {
    cells[ii]->calc_calcium_explicit_method(calcium_resistance_factor_, dCa[ii], dv[ii]);
  }

  for(unsigned ii = 0; ii < n_cell; ii++) {
    if(cells[ii]->share_calcium()) {
      cells[ii]->Ca() += dt_ * dCa[ii];
      cells[ii]->update_calcium();
      cells[ii]->v() += dt_ * dv[ii];
    }
  }

}


void Problem::
create_matrix_operator_splitting()
{

  delete matrix_;

  std::vector<cell::MainCell*>& cells = mesh_.get_cells();
  unsigned n_cell = cells.size();
  unsigned offset = 0;

  for(unsigned ii = 0; ii < n_cell; ii++) {
    cells[ii]->calc_matrix_entry_operator_spliting(dt_, offset, matrix_vals_);
    offset += n_cols_per_row_[ii];
  }

  matrix_ = new sparse_matrix::HypreMatrix();
  unsigned low = mesh_.get_first_dof();
  unsigned high = mesh_.get_n_cell() + low - 1;
  matrix_->create_matrix(low, high,
                         low, high,
                         &n_cols_per_row_[0],
                         &matrix_cols_[0], &matrix_vals_[0],
                         comm_->comm());

}


void Problem::
calc_rhs_operator_splitting(std::vector<double>& rhs)
{

  unsigned n_cell = mesh_.get_n_cell();
  std::vector<cell::MainCell*>& cells = mesh_.get_cells();

  for(unsigned ii = 0; ii < n_cell; ii++) {
    cells[ii]->calc_rhs_operator_spliting(dt_, rhs[ii]);
  }

}


void Problem::
solve_operator_splitting()
{

  this->create_matrix_operator_splitting();
  linear_solver_->set_matrix(*matrix_);
  mesh_.exchange_data();

  if(read_restart_) {
    this->mpi_read();
    mesh_.update_calcium_from_model();
  }

  mesh_.exchange_calcium_info_data();

  EventTimer restart_timer;
	unsigned index = 0;
  while(time_ < duration_) {
    this->iteration_operator_splitting();


    if(save_restart_ && time_ >= restart_start_ && time_ <= restart_end_ &&
       restart_timer.check_event(dt_, restart_interval_, time_ - restart_start_)) {
      this->mpi_write(time_);
    }

    //save to file
    double save_time = time_;
    save_time -= int(time_ / save_dt_) * save_dt_;

    if(save_time < dt_) {
			if(single_output_){
				std::string name = probe_prefix_ + conversion::to_string(index++, 5);
				pio.set_file_name(name);
			}			
      this->save_to_disc();
      info::info << "time " << time_ << std::endl;
    }
  }

	if(!single_output_){
		this->close_probe_files();
	}

}


void Problem::
save_to_disc()
{

	info::info << "saving... " << std::endl;
  std::list<cell::MainCell*>::const_iterator probe_cell_it = probe_cells_.begin();
	
	if(single_output_){
		std::string name = probe_prefix_ + conversion::to_string((*probe_cell_it)->original_index());
		std::vector<cell::MainCell*>& cells = mesh_.get_cells();
		unsigned n_cell = cells.size();
		for(unsigned ii = 0; ii < n_cell; ii++) {
			pio.get_stream() << cells[ii]->v() << " ";
			cells[ii]->output_data(pio.get_stream());
			pio.next(cells[ii]->original_index());
		}
		pio.commit();
	} else {
		for(std::list<std::ofstream*>::iterator it = probe_files_.begin();
				it != probe_files_.end(); ++it) {
			**it << (*probe_cell_it)->v() << " ";
			(*probe_cell_it)->output_data(**it);
			**it << std::endl;
			++probe_cell_it;
		}
	}

	
	info::info << "saving done " << std::endl;
}

void Problem::
close_probe_files()
{

  for(std::list<std::ofstream*>::iterator it = probe_files_.begin();
      it != probe_files_.end(); ++it) {
    (*it)->close();
    delete (*it);
  }

}



void Problem::
iteration_operator_splitting()
{

  this->process_stimuli();

  unsigned n_cell = mesh_.get_n_cell();
  std::vector<cell::MainCell*>& cells = mesh_.get_cells();
  std::vector<double> rhs(n_cell), x(n_cell);
  this->calc_rhs_operator_splitting(rhs);
  linear_solver_->solve(*matrix_, rhs, x);

  for(unsigned ii = 0; ii < n_cell; ii++) {
    cells[ii]->v() = x[ii];
  }

  this->integrate_ionic_models();
  mesh_.update_calcium_from_model();

  this->update_calcium();
  mesh_.exchange_data();

  time_ += dt_;
}


void Problem::
process_stimuli()
{

  for(std::list<Stimul>::iterator it = stimul_.begin();
      it != stimul_.end(); ++it) {
    if(time_ > it->end_) {
      it = stimul_.erase(it);
    }
  }

  for(std::list<Stimul>::iterator it = stimul_.begin();
      it != stimul_.end(); ++it) {
    it->update_voltage(time_, dt_);
  }

}

/*
void Problem::
init_ionic_models()
{

std::vector<cell::MainCell*>& cells = mesh_.get_cells();
unsigned n_cell = cells.size();

for(unsigned ii = 0; ii < n_cell; ii++) {
  cells[ii]->init_cell(dt_);
}

}
*/

void Problem::
integrate_ionic_models()
{

  std::vector<cell::MainCell*>& cells = mesh_.get_cells();
  unsigned n_cell = cells.size();
#pragma omp parallel for
  for(unsigned ii = 0; ii < n_cell; ii++) {
    cells[ii]->integrate(time_, dt_);
  }

}


void Problem::
add_ionic_models()
{

  std::map<std::string, IonicFactory*> container;
  container["TT06"] = new TT06Factory();
  //Place new models here
  container["Stoch"] = new StochFactory();
  factory_.set_model_factories(container);

}


void Problem::
init(mpi::MPIComm* comm)
{

  this->add_ionic_models();

  // reading options
  comm_ = comm;

  spirit_wrapper::SpiritReader reader;
  reader.read_file("object.data");

  const SpiritReaderObject& object = reader.get_object("main");

  double r_factor;
  object.get_prop("R_factor", r_factor);
  object.get_prop("Ca_R_factor", calcium_resistance_factor_);
  mesh_.init(comm, r_factor);

  object.get_prop("dt", dt_);
  object.get_prop("duration", duration_);
  std::string mesh;
  object.get_prop("mesh", mesh);
  {
    const SpiritReaderObject& object = reader.get_object(mesh);
    std::string file_name;
    object.get_prop("fileName", file_name);
    //input_output::mpi_istream stream(comm->comm(), file_name);
    mesh_.read_from_file(file_name);
  }

  //stimulation info
  std::string stimulation;
  object.get_prop("stimulation", stimulation);
  {

    const SpiritReaderObject& object = reader.get_object(stimulation);
    std::vector<std::string> objects;
    object.get_prop("objects", objects);
    unsigned n_objects = objects.size();

    for(unsigned ii = 0; ii < n_objects; ii++) {
      Stimul s;
      const SpiritReaderObject& object = reader.get_object(objects[ii]);
      object.get_prop("cellIds", s.cell_id_);
      object.get_prop("start", s.start_);
      object.get_prop("end", s.end_);
      object.get_prop("duration", s.duration_);
      object.get_prop("interval", s.interval_);
      object.get_prop("I", s.current_);
      stimul_.push_back(s);
    }

    mesh_.fill_stimuli_cell_pointers(stimul_);
  }

  //restart info
  std::string restart;
  object.get_prop("restart", restart);
  {

    const SpiritReaderObject& object = reader.get_object(restart);
    object.get_prop("saveRestart", save_restart_);

    if(save_restart_) {

      object.get_prop("start", restart_start_);
      object.get_prop("end", restart_end_);
      object.get_prop("interval", restart_interval_);
      object.get_prop("prefix", restart_prefix_);

    }

    object.get_prop("readRestart", read_restart_);
    object.get_prop("restartFileName", restart_file_name_);

  }



  std::vector<std::string> cells;
  object.get_prop("cell_types", cells);

  std::vector<std::map<std::string, double> > param_map;
  std::vector<unsigned> model_type;
  {

    unsigned n_cells = cells.size();
    std::vector<IonicModel*> im(n_cells);
    param_map.resize(n_cells);

    for(unsigned ii = 0; ii < n_cells; ii++) {

      const SpiritReaderObject& object = reader.get_object(cells[ii]);
      std::vector<std::string> params;
      object.get_prop("params", params);

      unsigned n_params = params.size();
      std::map<std::string, double> ref;

      for(unsigned jj = 0; jj < n_params; jj++) {
        double value;
        object.get_prop(params[jj], value);
        ref[params[jj]] = value;
      }

      std::string type;
      object.get_prop("type", type);
      im[ii] = factory_.create_model(type);
      im[ii]->read_parameters(ref);

      //creating output header
      std::string output;
      object.get_prop("output", output);
      {
        const SpiritReaderObject& object = reader.get_object(output);
        std::vector<std::string> vars;
        object.get_prop("variables", vars);
        im[ii]->set_output(vars);
				MPI_Comm mpi_comm = comm_->comm();
				pio.init(mpi_comm);
      }

    }

    //creating one instance of every ionic model with different param sets
    mesh_.init_cell_params(im);

    for(unsigned ii = 0; ii < n_cells; ii++) {
      delete im[ii];
    }

  }

  std::string output;
  object.get_prop("output", output);
  {

    const SpiritReaderObject& object = reader.get_object(output);
    std::vector<unsigned> cell_ids;
    object.get_prop("cellIds", cell_ids);
    object.get_prop("prefix", probe_prefix_);
    object.get_prop("dt", save_dt_);
		object.get_prop("single_file", single_output_);
		
    //finding pointer to cells in the mesh and opening files to write probe
    mesh_.find_probe(cell_ids, probe_cells_);
    unsigned n_probe = probe_cells_.size();

    probe_files_.resize(n_probe);
    std::list<cell::MainCell*>::const_iterator probe_cell_it = probe_cells_.begin();

    unsigned cell_index = 0;

		if(single_output_){
			/*for(;;++probe_cell_it) {
				std::string name = probe_prefix_ + conversion::to_string((*probe_cell_it)->original_index());
				*it = new std::ofstream(name.c_str());
				**it << "# Vm ";
				(*probe_cell_it)->output_header(**it);
				**it << "\n";				
				}*/
		} else {
			for(std::list<std::ofstream*>::iterator it = probe_files_.begin();
					it != probe_files_.end(); ++it) {
				std::string name = probe_prefix_ + conversion::to_string((*probe_cell_it)->original_index());
				*it = new std::ofstream(name.c_str());
				**it << "# Vm ";
				(*probe_cell_it)->output_header(**it);
				**it << "\n";
				++probe_cell_it;
			}
		}

  }


  this->create_matrix_map();

  linear_solver_ = new linear_solver::LinearSolver(comm);

  //  this->init_ionic_models();

}



}
