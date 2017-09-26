#ifndef __EXCEPTION_BASE__
#define __EXCEPTION_BASE__

#include <string>
#include <sstream>

namespace exceptions
{

class ExceptionBase
{

public:
  ExceptionBase() {}
  ExceptionBase(std::string message) {
    message_ << message;
  }


  template <typename Type>
  ExceptionBase & operator << (const Type & value) {
    message_ << value;
    return *this;
  }

  ExceptionBase(ExceptionBase& a) {
    message_ << a.message_.str();
  }

  ExceptionBase &operator = (ExceptionBase& a) {
    return a;
  }

  virtual std::string get_message() {
    return std::string("Error: ") + message_.str();
  }

private:

  std::stringstream message_;

};

}

#endif
