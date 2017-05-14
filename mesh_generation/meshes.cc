#include "meshes.h"
#include <iostream>

void
generate_mesh(std::string file_name, unsigned n_cables, unsigned n_per_cable, unsigned diameter_dad, unsigned dimensions)
{

  const unsigned type_norm = 1;
  const unsigned type_dad = 0;

  unsigned n_cell = n_per_cable * n_cables;
  unsigned i_center = n_per_cable/2;
  unsigned j_center = n_cables/2;

  std::ofstream file(file_name.c_str());
  file << n_cell << std::endl;

  if (dimensions == 2) {
    unsigned r_dad = diameter_dad/2;
    for(unsigned ii = 0; ii < n_per_cable; ii++) {
      for (unsigned jj = 0; jj < n_cables; jj++) {

	unsigned cell_id = get_cell_id(ii,jj,n_per_cable,n_cables);
	unsigned cell_type;
	if ((ii-i_center)*(ii-i_center) + (jj-j_center)*(jj-j_center) <= r_dad*r_dad) {
	    cell_type = type_dad;
	} else {
	  cell_type = type_norm;
        }
	
	unsigned n_neighb = 0;
	if (ii > 0)
	  n_neighb++;
	if (ii < n_per_cable-1)
	  n_neighb++;
	if (jj > 0)
	  n_neighb++;
	if (jj < n_cables-1)
	  n_neighb++;
	
	file << cell_type << " " << n_neighb;
	add_neighbor(file, 1.0, ii, jj-1, n_per_cable, n_cables);
	add_neighbor(file, 1.0, ii, jj+1, n_per_cable, n_cables);
	add_neighbor(file, 1.0, ii-1, jj, n_per_cable, n_cables);
	add_neighbor(file, 1.0, ii+1, jj, n_per_cable, n_cables);
	file << std::endl;
	
      }
    }
  } else if (dimensions == 1) {

    for(unsigned ii = 0; ii < n_per_cable; ii++) {
      for (unsigned jj = 0; jj < n_cables; jj++) {

	unsigned cell_id = get_cell_id(ii,jj,n_per_cable,n_cables);
	unsigned cell_type;
	unsigned pad0 = (n_per_cable-diameter_dad)/2;
	unsigned pad1 = n_per_cable-diameter_dad-pad0;
	if (ii >= pad0 && ii < n_per_cable-pad1) {
	  cell_type = type_dad;
	} else {
	  cell_type = type_norm;
	}
	
	unsigned n_neighb = 0;
	if (ii > 0)
	  n_neighb++;
	if (ii < n_per_cable-1)
	  n_neighb++;
	
	file << cell_type << " " << n_neighb;
	add_neighbor(file, 1.0, ii-1, jj, n_per_cable, n_cables);
	add_neighbor(file, 1.0, ii+1, jj, n_per_cable, n_cables);
	file << std::endl;
	
      }
    }
  } else {
    std::cerr << "Error: invalid dimensions value, must be 1 or 2." << std::endl;
  }

  file.close();

}

unsigned get_cell_id(unsigned ix, unsigned iy, unsigned n_per_cable, unsigned n_cables) {

  if (ix < 0 || ix >= n_per_cable || iy < 0 || iy >= n_cables)
    return -1;

  return ix*n_cables + iy;

}

void add_neighbor(std::ofstream &f, double resist, unsigned ix, unsigned iy, unsigned n_per_cable, unsigned n_cables) {

  unsigned neighb_id = get_cell_id(ix, iy, n_per_cable, n_cables);
  if (neighb_id != -1) {
    f << " " << neighb_id << " " << resist;
  }

}
