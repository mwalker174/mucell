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

  //ExceptionBase &operator = (ExceptionBase& a) {
  //  return a;
  //}

  virtual std::string get_message() {
    return std::string("Error: ") + message_.str();
  }

protected:

  std::stringstream message_;

};

template<class T>
class MainException : public exceptions::ExceptionBase
{
public:
  MainException() {}
  MainException(MainException<T>& a) {
    message_ << a.message_.str();
  }
  template <typename Type>
  MainException<T> & operator << (const Type & value) {
    message_ << value;
    return *this;
  }
};

class BadFileBase
{
};
typedef exceptions::MainException<BadFileBase> BadFile;


class FileNotExistsBase
{
};
typedef exceptions::MainException<FileNotExistsBase> FileNotExists;


class ModelNotFoundBase
{
};
typedef exceptions::MainException<ModelNotFoundBase> ModelNotFound;


class WrongModelIdBase
{
};
typedef exceptions::MainException<WrongModelIdBase> WrongModelId;


class WrongOptionBase
{
};
typedef exceptions::MainException<WrongOptionBase> WrongOption;


}

#endif
