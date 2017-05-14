#ifndef __MESHES__
#define __MESHES__

#include <string>
#include <fstream>

void generate_mesh(std::string file_name, unsigned n_cables, unsigned n_per_cable, unsigned diameter_dad, unsigned dimensions);
unsigned get_cell_id(unsigned ix, unsigned iy, unsigned n_x, unsigned n_y);
void add_neighbor(std::ofstream &f, double resist, unsigned ix, unsigned iy, unsigned n_x, unsigned n_y);

#endif
