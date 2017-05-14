#include "ModelFactory.h"
#include "ExceptionBase.h"

namespace ryr
{

namespace factory
{


IonicModel*
ModelFactory::
create_model(std::string model)
{

  std::map<std::string, IonicFactory*>::iterator it = container_.find(model);

  if(it == container_.end()) {
    throw exceptions::ModelNotFound() << "Ionic model " << model << " is not supported.";
  }

  return it->second->create_model();

}



}

}
