#include "Mesh.h"
#include "ExceptionBase.h"
#include "mpi_io_stream.h"

#include <stdio.h>
#include <mpi.h>
#include <ptscotch.h>

#include <iostream>
#include <algorithm>
#include <set>
#include <fstream>
#include <cassert>

namespace ryr
{

namespace mesh
{


static
unsigned
get_n_local(unsigned size, unsigned rank, unsigned n_vertex)
{
  unsigned local = n_vertex / size;
  unsigned reminder = n_vertex % size;

  unsigned n_local = local;

  if(rank < reminder)
    n_local += 1;

  return n_local;
}

static
void
get_local_indices(unsigned size, unsigned rank, unsigned n_vertex, unsigned& start, unsigned& n_local)
{
  unsigned local = n_vertex / size;
  unsigned reminder = n_vertex % size;

  n_local = local;

  if(rank < reminder)
    n_local += 1;

  start = rank * local;

  if(rank < reminder)
    start += rank;
  else
    start += reminder;
}


static void
partition_with_scotch(mpi::MPIComm* mpi_comm, unsigned n_vertex,  const std::vector<std::vector<unsigned> > &vertex_adj,
                      std::vector<SCOTCH_Num> &vertex_rank_buffer);


void Mesh::
read_from_file(std::string file_name)
{

	
	std::map<unsigned, std::list<std::pair<unsigned, cell::MainCell*> > > haloed;
	std::map<unsigned, std::list<std::pair<unsigned, cell::HaloCell*> > > halo;
			

	input_output::mpi_istream fstream(comm_->comm(), file_name.c_str());
//	std::ifstream fstream(file_name.c_str());

  const unsigned max_string_size = 10000;
  char s[max_string_size];

  unsigned n_cell;
  fstream >> n_cell;
  std::vector<SCOTCH_Num> vertex_rank_buffer;
  std::vector<unsigned> local_vertex;
  unsigned rank = comm_->rank();


  unsigned vertex_counter = 0;
  {

    unsigned n_local_cell, start;
    get_local_indices(comm_->n_proc(), rank, n_cell, start, n_local_cell);


    for(unsigned ii = 0; ii < start + 1; ii++) {
      fstream.getline(s, max_string_size);
    }

    std::vector<std::vector<unsigned> > vertex_adj(n_local_cell);

    for(unsigned ii = 0; ii < n_local_cell; ii++) {

      unsigned n_adj, param;
      fstream >> param >> n_adj;

      std::vector<unsigned>& a = vertex_adj[ii];
      a.resize(n_adj);
      double blah;

      for(unsigned jj = 0; jj < n_adj; jj++) {
        fstream >> a[jj] >> blah;
      }
    }

    partition_with_scotch(comm_, n_cell,
                          vertex_adj,
                          vertex_rank_buffer);

    //count number of local cells


    for(unsigned ii = 0; ii < n_cell; ii++) {
      if(unsigned(vertex_rank_buffer[ii]) == rank)
        vertex_counter++;
    }

    //creating array of local cells
    local_vertex.resize(vertex_counter);
    unsigned index = 0;

    for(unsigned ii = 0; ii < n_cell; ii++) {
      if(static_cast<unsigned>(vertex_rank_buffer[ii]) == rank)
        local_vertex[index++] = ii;
    }

    std::sort(local_vertex.begin(), local_vertex.end());
  }

  fstream.close();
	input_output::mpi_istream stream(comm_->comm(), file_name.c_str());
  stream >> n_cell;
  {

    unsigned cell_index = 0;
    std::vector<cell::CellBase*> visited(n_cell, NULL);
    cells_.resize(vertex_counter, NULL);

    for(unsigned ii = 0; ii < vertex_counter; ii++) {
      unsigned index = local_vertex[ii];
      visited[index] = cells_[ii] = new cell::MainCell();
      cells_[ii]->original_index() = index;
    }

    unsigned n_size = local_vertex.size();
    std::sort(local_vertex.begin(), local_vertex.end());

    unsigned line_index = 0;

    for(unsigned ii = 0; ii < n_size; ii++) {
      unsigned n_skip = local_vertex[ii] - line_index;

      for(unsigned ii = 0; ii < n_skip + 1; ii++) {
        stream.getline(s, max_string_size);
        line_index++;
      }

      //creaing the cell
      cell::MainCell* main_cell = cells_[cell_index];

      unsigned n_adj, param;
      stream >> param >> n_adj;
      main_cell->get_param_id() = param;
      std::vector<std::pair<unsigned, double> > sort_vector(n_adj);

      for(unsigned jj = 0; jj < n_adj; jj++) {
        unsigned adj_cell;
        double res;
        stream >> adj_cell >> res;
        res *= res_;
        std::pair<unsigned, double>& sr = sort_vector[jj];
        sr.first = adj_cell;
        sr.second = res;
      }

      std::sort(sort_vector.begin(), sort_vector.end());


      for(unsigned jj = 0; jj < n_adj; jj++) {
        unsigned adj_cell;
        double res;
        std::pair<unsigned, double>& sr = sort_vector[jj];
        adj_cell = sr.first;
        res = sr.second;
        /// \todo sort vectors

        unsigned process = static_cast<unsigned>(vertex_rank_buffer[adj_cell]);

        //checking if this cell belongs to current process
        if(process == rank) {
          main_cell->add_cell(visited[adj_cell], res);
        } else {

          cell::CellBase*& vis = visited[adj_cell];

          if(vis == NULL) { //new halo cell
            cell::HaloCell* h = new cell::HaloCell();
            vis = h;
            vis->original_index() = adj_cell;
            halo_cells_.push_back(h);
						//  halo_container_[process].push_back(h);
            //haloed_container_[process].push_back(main_cell);
						halo[process].push_back(std::make_pair(adj_cell, h));						
          }
					haloed[process].push_back(std::make_pair(main_cell->original_index(), main_cell));
          main_cell->add_cell(vis, res);
        }
      }

      cell_index++;
			// line_index++;
    }
  }
  stream.close();

	// fill haloed and halo lists
	//first we sort and extract unique values than set to final halo and haloed arrays
	{

		std::map<unsigned, std::list<std::pair<unsigned, cell::HaloCell*> > >::iterator halo_it = halo.begin();
		for(; halo_it!= halo.end(); ++halo_it){
			halo_it->second.sort();
			std::list<std::pair<unsigned, cell::HaloCell*> >::iterator end_it = 
				std::unique(halo_it->second.begin(), halo_it->second.end());      
			halo_it->second.resize(std::distance(halo_it->second.begin(), end_it));
			for(std::list<std::pair<unsigned, cell::HaloCell*> >::iterator it = halo_it->second.begin(); it != halo_it->second.end();
					++it){
				halo_container_[halo_it->first].push_back(it->second);
			}
		}


		std::map<unsigned, std::list<std::pair<unsigned, cell::MainCell*> > >::iterator haloed_it = haloed.begin();		
		for(; haloed_it!= haloed.end(); ++haloed_it){
			haloed_it->second.sort();
			std::list<std::pair<unsigned, cell::MainCell*> >::iterator end_it = 
				std::unique(haloed_it->second.begin(), haloed_it->second.end());      
			haloed_it->second.resize(std::distance(haloed_it->second.begin(), end_it));
			for(std::list<std::pair<unsigned, cell::MainCell*> >::iterator it = haloed_it->second.begin(); it != haloed_it->second.end();
					++it){
				haloed_container_[haloed_it->first].push_back(it->second);
			}
		}
	}


  info::info << "Matrix partitioned" << std::endl;
  this->create_matrix_index();
  info::info << "Matrix index created" << std::endl;
}


struct ExData {
  double v, dv, Ca;
};


void Mesh::
exchange_data()
{


  MPI_Comm comm = comm_->comm();
  unsigned my_rank = comm_->rank();
  unsigned nproc = comm_->n_proc();

  unsigned n_haloed = haloed_container_.size();
  unsigned n_halo = halo_container_.size();

  std::vector<MPI_Request> send_req(n_haloed + n_halo);

  unsigned index = 0;
  std::vector<std::vector<ExData> > data_container(haloed_container_.size());

  for(HaloedContainer::const_iterator it = haloed_container_.begin();
      it != haloed_container_.end(); ++it) {
    const HaloedVector& halo_vector = it->second;
    unsigned n_halo_vector = halo_vector.size();
    std::vector<ExData>& ex = data_container[index];
    ex.resize(n_halo_vector);

    for(unsigned ii = 0; ii < n_halo_vector; ii++) {
      ex[ii].v = halo_vector[ii]->v();
      ex[ii].dv = halo_vector[ii]->dv();
      ex[ii].Ca = halo_vector[ii]->Ca();
    }

    MPI_Isend(&ex[0], n_halo_vector * sizeof(ExData), MPI_BYTE, it->first, 0,
              comm, &send_req[index++]);
  }

  unsigned vector_index = 0;
  std::vector<std::vector<ExData> > data(halo_container_.size());

  for(HaloContainer::const_iterator it = halo_container_.begin();
      it != halo_container_.end(); ++it) {
    const HaloVector& halo_vector = it->second;
    unsigned n_halo_vector = halo_vector.size();
    std::vector<ExData>& d =  data[vector_index++];
    d.resize(n_halo_vector);
    MPI_Irecv(&d[0], n_halo_vector * sizeof(ExData), MPI_BYTE, it->first, 0,
              comm, &send_req[index++]);
  }

  MPI_Waitall(index, &send_req[0], MPI_STATUSES_IGNORE);

  vector_index = 0;

  for(HaloContainer::const_iterator it = halo_container_.begin();
      it != halo_container_.end(); ++it) {
    std::vector<ExData>& d =  data[vector_index];
    const HaloVector& halo_vector = it->second;
    unsigned n_halo_vector = halo_vector.size();

    for(unsigned ii = 0; ii < n_halo_vector; ii++) {
      cell::HaloCell* cell = halo_vector[ii];
      cell->v() = d[ii].v;
      cell->dv() = d[ii].dv;
      cell->Ca() = d[ii].Ca;
    }

    vector_index++;
  }


}


void Mesh::
create_matrix_index()
{

  unsigned n_cell = cells_.size();

  MPI_Comm comm = comm_->comm();
  unsigned my_rank = comm_->rank();
  unsigned n_proc = comm_->n_proc();
  {
    unsigned n_send = n_proc - my_rank - 1;
    std::vector<MPI_Request> send_req;

    if(n_send > 0) {
      send_req.resize(n_send);

      for (unsigned p = my_rank + 1; p < n_proc; p++) {
        MPI_Isend(&n_cell, 1, MPI_UNSIGNED, p, 0,
                  comm, &send_req[p - my_rank - 1]);
      }
    }

    std::vector<unsigned> n_eqn_on_proc;
    std::vector<MPI_Request> recv_req;

    if(my_rank != 0) {
      recv_req.resize(my_rank);
      n_eqn_on_proc.resize(my_rank);

      for (unsigned p = 0; p < my_rank; p++) {
        MPI_Irecv(&n_eqn_on_proc[p], 1, MPI_UNSIGNED, p, 0,
                  comm, &recv_req[p]);
      }

      MPI_Waitall(my_rank, &recv_req[0], MPI_STATUSES_IGNORE);
    }

    if(n_send > 0)
      MPI_Waitall(n_send, &send_req[0], MPI_STATUSES_IGNORE);

    bump_ = 0;

    for (unsigned p = 0; p < my_rank; p++) {
      bump_ += n_eqn_on_proc[p];
    }

    for(unsigned ii = 0; ii < n_cell; ii++) {
      cells_[ii]->matrix_index() = ii + bump_;
    }
  }

  //exchanging indexes
  {

    unsigned n_haloed = haloed_container_.size();
    unsigned n_halo = halo_container_.size();

    std::vector<std::vector<unsigned> > halo_index(n_halo);
    std::vector<std::vector<unsigned> > haloed_index(n_haloed);

    //filling sending vector with matrix indexes for haloed nodes
    {
      unsigned index = 0;

      for(HaloedContainer::const_iterator it = haloed_container_.begin();
          it != haloed_container_.end(); ++it) {
        std::vector<unsigned>& v = haloed_index[index++];
        const HaloedVector& hv = it->second;
        unsigned h_v = hv.size();
        v.resize(h_v);

        for(unsigned ii = 0; ii < h_v; ii++) {
          v[ii] = hv[ii]->matrix_index();
        }
      }
    }

    std::vector<MPI_Request> send_req(n_haloed + n_halo);

    unsigned index = 0;

    for(HaloedContainer::const_iterator it = haloed_container_.begin();
        it != haloed_container_.end(); ++it) {

      std::vector<unsigned>& v = haloed_index[index];
      MPI_Isend(&v[0], v.size(), MPI_UNSIGNED, it->first, 0,
                comm, &send_req[index++]);
    }

    unsigned vector_index = 0;

    for(HaloContainer::iterator it = halo_container_.begin();
        it != halo_container_.end(); ++it) {
      std::vector<unsigned>& v = halo_index[vector_index++];
      v.resize(it->second.size());
      MPI_Irecv(&v[0], v.size(), MPI_UNSIGNED, it->first, 0,
                comm, &send_req[index++]);
    }

    MPI_Waitall(index, &send_req[0], MPI_STATUSES_IGNORE);


    {
      unsigned index = 0;

      for(HaloContainer::const_iterator it = halo_container_.begin();
          it != halo_container_.end(); ++it) {
        std::vector<unsigned>& v = halo_index[index++];
        const HaloVector& hv = it->second;
        unsigned n_v = v.size();
        assert(hv.size() == n_v);

        for(unsigned ii = 0; ii < n_v; ii++) {
          hv[ii]->matrix_index() = v[ii];
        }
      }
    }

  }
}


static void
partition_with_scotch(mpi::MPIComm* mpi_comm,
                      unsigned n_vertex,
                      const std::vector<std::vector<unsigned> > &vertex_adj,
                      std::vector<SCOTCH_Num> &vertex_rank_buffer)
{

  MPI_Comm comm = mpi_comm->comm();
  unsigned my_rank = mpi_comm->rank();
  unsigned n_proc = mpi_comm->n_proc();

  // calculating number of edges for scotch routine
  unsigned edge_total = 0;

  unsigned n_local = vertex_adj.size();

  for(unsigned ii = 0; ii < n_local; ii++) {
    const std::vector<unsigned>& v = vertex_adj[ii];
    edge_total += v.size();
  }


  // sum number of edges on communicator
  unsigned global_edge_total = 0;
  MPI_Allreduce(&edge_total, &global_edge_total, 1, MPI_UNSIGNED, MPI_SUM, comm);

  SCOTCH_randomReset();

  //initializing scotch graph
  SCOTCH_Strat strat;
  SCOTCH_Dgraph grafdat;

  if (SCOTCH_dgraphInit (&grafdat, comm) != 0) {
    throw exceptions::ExceptionBase() << "graph init failed";
  }

  info::info << "setup scotch input" << std::endl;
  //fill array to prepare input to build graph
  std::vector<SCOTCH_Num> *vertloctab =  new std::vector<SCOTCH_Num>(n_local + 1, 0);
  std::vector<SCOTCH_Num> *edgeloctab = new std::vector<SCOTCH_Num>(edge_total);

  unsigned index = 0;

  for(unsigned ii = 0; ii < n_local; ii++) {
    const std::vector<unsigned>& v = vertex_adj[ii];
    unsigned degree = v.size();
    (*vertloctab)[ii + 1] = (*vertloctab)[ii] + degree;

    for(std::vector<unsigned>::const_iterator it = v.begin(); it != v.end(); it++) {
      (*edgeloctab)[index] = (*it);
      index++;
    }
  }


  //finally build graph
  SCOTCH_Num baseval = 0, vertlocnbr = n_local, vertlocmax = n_local, edgelocnbr = edge_total;
  SCOTCH_Num edgelocsiz = edgelocnbr;
  SCOTCH_Num *vendloctab = &(*vertloctab)[1], *veloloctab = NULL,
              *vlblocltab = NULL, *edloloctab = NULL, *edgegsttab = NULL;

  info::info << "building scotch graph ..." << std::endl;

  if(SCOTCH_dgraphBuild (&grafdat,
                         baseval, vertlocnbr, vertlocmax, &(*vertloctab)[0], vendloctab, veloloctab,
                         vlblocltab, edgelocnbr, edgelocsiz, &(*edgeloctab)[0], edgegsttab, edloloctab)) {
    throw exceptions::ExceptionBase() << "graph building failed";
  }


  info::info << "checking graph ..." << std::endl;

  if(SCOTCH_dgraphCheck(&grafdat)) {
    throw exceptions::ExceptionBase() << "graph inconsistent";
  }


  std::vector<SCOTCH_Num> num(n_local);

  if(SCOTCH_stratInit(&strat)) {
    throw exceptions::ExceptionBase() << "loading of default strategy failed";
  }

  //SCOTCH_randomReset();

  info::info << "partitioning ...." << std::endl;

  if(SCOTCH_dgraphPart(&grafdat,
                       n_proc,
                       &strat,
                       &num[0])) {
    throw exceptions::ExceptionBase() << "main procedure of partitioning failed";
  }

	info::info << "partitioning done" << std::endl;

  delete vertloctab;
  delete edgeloctab;

  SCOTCH_dgraphExit(&grafdat);
  SCOTCH_stratExit(&strat);


  std::vector<int> displs(n_proc);
  std::vector<int> rcounts(n_proc);

  unsigned displacement = 0;

  for (unsigned ii = 0; ii < n_proc; ii++) {
    rcounts[ii] = get_n_local(n_proc, ii, n_vertex);
    displs[ii] = displacement;
    displacement += rcounts[ii];
  }

  vertex_rank_buffer.resize(displacement);
  MPI_Allgatherv(&num[0], n_local, MPI_UNSIGNED, &vertex_rank_buffer[0], &rcounts[0], &displs[0], MPI_UNSIGNED, comm);

}


void Mesh::
fill_stimuli_cell_pointers(std::list<Stimul>& stim)
{

  // the goal is to make the join of two containers below using cell_id (first argument in pair) and
  // set pointers in std::vector<cell::MainCell*> of Stimul

  // cell_id, index in std::vector<cell::MainCell*>, Stimul*
  std::list<std::pair<unsigned, std::pair<unsigned, Stimul*> > > stimul;
  // cell_id, cell::MainCell*
  std::vector<std::pair<unsigned, cell::MainCell*> > local_cells;


  // filling local_cells container and sort
  {
    unsigned n_cells = cells_.size();
    local_cells.resize(cells_.size());

    for(unsigned ii = 0; ii < n_cells; ii++) {
      local_cells[ii].first = cells_[ii]->original_index();
      local_cells[ii].second = cells_[ii];
    }

    std::sort(local_cells.begin(), local_cells.end());
  }


  //filling second container and sort
  {
    unsigned n_stim = stim.size();
    unsigned ii = 0;

    for(std::list<Stimul>::iterator  it = stim.begin();
        it != stim.end();	++it) {


      //copy all cells from current stimulus to temp vector
      std::vector<unsigned>& s  = it->cell_id_;
      unsigned n_s = s.size();
      it->cells_.resize(n_s, NULL);
      std::vector<std::pair<unsigned, std::pair<unsigned, Stimul*> > > temp_container(n_s);

      for(unsigned jj = 0; jj < n_s; jj++) {
        temp_container[jj].first = s[jj];
        temp_container[jj].second.first = jj;
        temp_container[jj].second.second = &(*it);
      }

      stimul.insert(stimul.end(), temp_container.begin(), temp_container.end());
      ii++;
    }

    stimul.sort();
  }

  //join two arrays using modified algorithm from stl
  {
    std::list<std::pair<unsigned, std::pair<unsigned, Stimul*> > >::iterator first1 = stimul.begin(), last1 = stimul.end();
    std::vector<std::pair<unsigned, cell::MainCell*> >::iterator first2 = local_cells.begin(), last2 = local_cells.end();

    while (first1 != last1 && first2 != last2) {
      if (first1->first < first2->first) ++first1;
      else if (first2->first < first1->first) ++first2;
      else {
        //assigning cell::MainCell* pointer to element of Stimul std::vector<cell::MainCell*>
        first1->second.second->cells_[first1->second.first] = first2->second;
        ++first1;
        ++first2;
      }
    }
  }
}


void Mesh::
find_probe(const std::vector<unsigned>& cell_ids, std::list<cell::MainCell*>& out)
{

  std::vector<std::pair<unsigned, cell::MainCell*> > cell;
  unsigned n_cell = cells_.size();
  cell.resize(n_cell);

  for(unsigned ii = 0; ii < n_cell; ii++) {
    cell[ii].first = cells_[ii]->original_index();
    cell[ii].second = cells_[ii];
  }

  std::sort(cell.begin(), cell.end());
  std::vector<unsigned> ids(cell_ids);
  std::sort(ids.begin(), ids.end());

  {
    std::vector<std::pair<unsigned, cell::MainCell*> >::iterator first1 = cell.begin(), last1 = cell.end();
    std::vector<unsigned>::iterator first2 = ids.begin(), last2 = ids.end();

    while (first1 != last1 && first2 != last2) {
      if (first1->first < *first2) ++first1;
      else if (*first2 < first1->first) ++first2;
      else {
        //assigning cell::MainCell* pointer to element of Stimul std::vector<cell::MainCell*>
        out.push_back(first1->second);
        ++first1;
        ++first2;
      }
    }
  }

}


void Mesh::
read_restart(std::istream& i)
{

  unsigned n_cells = cells_.size();

  for(unsigned ii = 0; ii < n_cells; ii++) {
    unsigned index = cells_[ii]->original_index();
    long long length;
    unsigned read_index;

    do {
      double v;
      i >> read_index;
      i >> length;

      char c;
      i.ignore(1, ' ');

      if(read_index != index) {
        std::vector<char> buffer(length);
        i.read(&buffer[0], length);
      }

    } while(i.good() && (read_index != index));

    if(read_index != index)
      throw exceptions::BadFile() << "Restart file is corrupted";

    std::vector<char> buffer(length);
    i.read(&buffer[0], length);

    std::string buffer_string(&buffer[0], length);
    std::istringstream s(buffer_string);
    cells_[ii]->read_restart(s);
  }

  this->exchange_data();

}




void Mesh::
exchange_calcium_info_data()
{

  MPI_Comm comm = comm_->comm();
  unsigned my_rank = comm_->rank();
  unsigned nproc = comm_->n_proc();

  unsigned n_haloed = haloed_container_.size();
  unsigned n_halo = halo_container_.size();

  std::vector<MPI_Request> send_req(n_haloed + n_halo);

  unsigned index = 0;
  std::vector<std::vector<char> > data_container(haloed_container_.size());

  for(HaloedContainer::const_iterator it = haloed_container_.begin();
      it != haloed_container_.end(); ++it) {
    const HaloedVector& halo_vector = it->second;
    unsigned n_halo_vector = halo_vector.size();
    std::vector<char>& ex = data_container[index];
    ex.resize(n_halo_vector);

    for(unsigned ii = 0; ii < n_halo_vector; ii++) {
      ex[ii] = halo_vector[ii]->share_calcium();
    }

    MPI_Isend(&ex[0], n_halo_vector * sizeof(char), MPI_BYTE, it->first, 0,
              comm, &send_req[index++]);
  }

  unsigned vector_index = 0;
  std::vector<std::vector<char> > data(halo_container_.size());

  for(HaloContainer::const_iterator it = halo_container_.begin();
      it != halo_container_.end(); ++it) {
    const HaloVector& halo_vector = it->second;
    unsigned n_halo_vector = halo_vector.size();
    std::vector<char>& d =  data[vector_index++];
    d.resize(n_halo_vector);
    MPI_Irecv(&d[0], n_halo_vector * sizeof(char), MPI_BYTE, it->first, 0,
              comm, &send_req[index++]);
  }

  MPI_Waitall(index, &send_req[0], MPI_STATUSES_IGNORE);

  vector_index = 0;

  for(HaloContainer::const_iterator it = halo_container_.begin();
      it != halo_container_.end(); ++it) {
    std::vector<char>& d =  data[vector_index];
    const HaloVector& halo_vector = it->second;
    unsigned n_halo_vector = halo_vector.size();

    for(unsigned ii = 0; ii < n_halo_vector; ii++) {
      cell::HaloCell* cell = halo_vector[ii];
      cell->share_calcium() = static_cast<bool>(d[ii]);
    }

    vector_index++;
  }


}




}
}
