#include "parser/spirit_wrapper.h"
#include "meshes.h"

int main(int argc, char **argv)
{

  try {

    spirit_wrapper::SpiritReader reader;
    reader.read_file("object.data");
    const spirit_wrapper::SpiritReaderObject& object = reader.get_object("main");

    std::string file_name;
    object.get_prop("fileName", file_name);
    unsigned n_cables, n_per_cable, diameter_dad, dimensions;
    object.get_prop("n_cables", n_cables);
    object.get_prop("n_per_cable", n_per_cable);
    object.get_prop("diameter_dad", diameter_dad);
    object.get_prop("dimensions", dimensions);
    generate_mesh(file_name, n_cables, n_per_cable, diameter_dad, dimensions);

  } catch(exceptions::ExceptionBase e) {
    std::cout << e.get_message() << std::endl;
  }

}

