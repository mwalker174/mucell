#ifndef __MODEL_FACTORY__
#define __MODEL_FACTORY__

#include "IonicModel.h"

#include <map>

namespace ryr
{

namespace factory
{

class ModelFactory
{
public:
  IonicModel* create_model(std::string type);

  void set_model_factories(const std::map<std::string, IonicFactory*>& container) {
    container_ = container;
  }

private:
  std::map<std::string, IonicFactory*> container_;
};

}

}
#endif
