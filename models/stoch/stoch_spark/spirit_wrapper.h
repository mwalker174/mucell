#ifndef __SPIRIT_WRAPPER__
#define __SPIRIT_WRAPPER__

#include <string>
#include <map>
#include <sstream>
#include <iostream>
#include <vector>


#include "ExceptionBase.h"


namespace spirit_wrapper
{

class ParserException : public exceptions::ExceptionBase
{

};


class PathString : public std::string
{
public:
  PathString(std::string s) : std::string(s) {
  }
  PathString() : std::string() {
  }
};

/// \short Single object (text in {}) from text IO
///
/// This class is returned by SpiritReader::get_object
class SpiritReaderObject
{

private:

  friend class SpiritReader;
  static void parse_string_vector(std::string& stream, std::vector<std::string>& vec);
  template <typename T>
  static void parse_string_brackets(std::string& stream, std::vector<std::vector<T> >& vec);
  template <typename T>
  static void get_unsigned_list(std::string vec, std::vector<T> &value);

  template <typename T>
  inline void get_prop_base(const std::string& prop_name, T &value) const {

    std::map<std::string, std::string>::const_iterator it = props_.find(prop_name);

    if(it == props_.end()) {
      throw ParserException() << "Property " << prop_name << " does not exists in object " << name_;
    }

    std::istringstream stream(it->second);
    stream >> value;

    if(stream.fail()) {
      throw ParserException() << "Type conversion for " << prop_name << " with value " << it->second << " failed";
    }
  }

public:

  typedef PathString PathType;

  template <typename T>
  inline void get_prop(const std::string& prop_name, T &value) const {

    this->get_prop_base(prop_name, value);

  }

  template <typename T>
  inline void set_prop(const std::string& prop_name, const T &value) {

    props_[prop_name] = value;

  }

  std::map<std::string, std::string>& get_props() {
    return props_;
  }

  const std::map<std::string, std::string>& get_props() const {
    return props_;
  }

  std::string get_name() const {
    return name_;
  }

private:

  enum TYPES {NUMBER, STRING, VECTOR, VECTOR_STRING};

  std::map<std::string, std::string> types_;
  std::map<std::string, std::string> props_;
  std::string name_;

};


template <>
inline void SpiritReaderObject::
get_prop_base(const std::string& prop_name, std::string &value) const
{

  std::map<std::string, std::string>::const_iterator it = props_.find(prop_name);

  if(it == props_.end()) {
    throw ParserException() << "Property " << prop_name << " does not exists in object " << name_;
  }

  value = it->second;

}


template <>
inline void SpiritReaderObject::get_prop(const std::string& prop_name,
    std::string &value) const
{
  this->get_prop_base(prop_name, value);
}


template <typename T>
void put_vector_to_string(std::ostream& os, const std::vector<std::vector<T> >& value);


template <>
inline void SpiritReaderObject::
set_prop(const std::string& prop_name, const std::string& value)
{

  types_[prop_name] = STRING;
  props_[prop_name] = std::string("\"") + value + std::string("\"");

}


template <>
inline void SpiritReaderObject::
set_prop(const std::string& prop_name, const std::vector<std::string> &value)
{

  std::ostringstream stream;
  unsigned n_value = value.size();

  for(unsigned ii = 0; ii < n_value; ii++) {
    stream << value[ii];

    if(ii != n_value - 1)
      stream << ",";
  }

  props_[prop_name] = stream.str();
  types_[prop_name] = VECTOR_STRING;

}


template <>
inline void SpiritReaderObject::
set_prop(const std::string& prop_name, const unsigned &value)
{

  std::ostringstream stream;
  stream << value;
  props_[prop_name] = stream.str();
  types_[prop_name] = NUMBER;

}

template <>
inline void SpiritReaderObject::
set_prop(const std::string& prop_name, const std::vector<std::vector<double> > &value)
{

  std::ostringstream stream;
  put_vector_to_string(stream, value);
  props_[prop_name] = stream.str();
  types_[prop_name] = VECTOR;
}


template <>
inline void SpiritReaderObject::
set_prop(const std::string& prop_name, const std::vector<std::vector<unsigned> > &value)
{

  std::ostringstream stream;
  put_vector_to_string(stream, value);
  props_[prop_name] = stream.str();
  types_[prop_name] = VECTOR;
}



/*!
	Specialization of get_prop that is a placeholder, in case we need to do folder conversions for some OS.
	This specialization should be used to read all paths even it is identical to original function now
*/
template <>
inline void SpiritReaderObject::get_prop(const std::string& prop_name,
    SpiritReaderObject::PathType &value) const
{
  std::string string_value;
  this->get_prop_base(prop_name, string_value);
  value = string_value;
}


template <>
inline void SpiritReaderObject::get_prop(const std::string& prop_name, bool &value) const
{
  std::string boolean;
  this->get_prop(prop_name, boolean);
  value = (boolean == "true");

  if(!value && boolean != "false") {
    throw ParserException() << "Boolean type conversion for " << prop_name << " with value " << boolean << " failed";
  }
}



template <>
inline void SpiritReaderObject::get_prop(const std::string& prop_name,
    std::vector<std::vector<unsigned> > &value) const
{
  std::string vec;
  this->get_prop(prop_name, vec);
  parse_string_brackets(vec, value);
}


template <>
inline void SpiritReaderObject::get_prop(const std::string& prop_name,
    std::vector<std::vector<double> > &value) const
{
  std::string vec;
  this->get_prop(prop_name, vec);
  parse_string_brackets(vec, value);
}


template <>
inline void SpiritReaderObject::get_prop(const std::string& prop_name,
    std::vector<unsigned> &value) const
{

  std::string vec;
  this->get_prop(prop_name, vec);
  get_unsigned_list(vec, value);

}

template <>
inline void SpiritReaderObject::get_prop(const std::string& prop_name,
    std::vector<double> &value) const
{

  std::string vec;
  this->get_prop(prop_name, vec);
  get_unsigned_list(vec, value);

}



template <>
inline void SpiritReaderObject::get_prop(const std::string& prop_name,
    std::vector<std::string> &value) const
{

  std::string vec;
  this->get_prop(prop_name, vec);
  parse_string_vector(vec, value);

}

/// \brief Simple analog of PropertyTree from the Boost library. This class will be a wrapper around PropertyTree later. Currently, Spirit is used to parse the input files in "unwrapped" JSON format.
/*!
  ~~~~~~~~~~~~~~~{.py}
  main{
       problemType = heart;
       problemObject = heartProblem; //std::string
  }


  heartProblem{
       typeInfo = problemObject; // object type, optional

       meshObject = finiteElementMesh; // child object, see this object below; use std::string
  }

  finiteElementMesh{
       typeInfo = meshObject;

       node_file = "canine.node" // string can be used within "" or without them; use std::string type
       double_var = 1.1; // double
       unsigned_var = 1; // unsigned
       vector_of_unsigned = 1,3,4,2; // std::vector<unsigned>
       vector_of_double = 1.1,3.3,4.1,2.9; // std::vector<double>
       vector_of_string = aa,ccc,mmm; // no spaces
       vector_of_vectors_double = [1.1,2.2],[3.3,9.9,5.5]; // std::vector<std::vector<double> >, no spaces
       vector_of_vectors_unsigned = [1,2],[3,9,5],[4,5]; // std::vector<std::vector<double> >, no spaces
       path_var = "folder/folder"; // SpiritReaderObject::PathType

  }

  ~~~~~~~~~~~~~~~
 */

class SpiritReader
{
public:
  /// Reading input file, should be called before getting SpiritReaderObject (JSON node objects)
  /*!
    \param file_name a file name of the input file
  */
  virtual void read_file(std::string file_name);
  void read_stream(std::istream& stream);
  void read_string(std::string& stream);

  void write_string(std::string& stream);

  /*!
    \param object_name a file name of the input file
    \return SpiritReaderObject that is used to access
  */
  const SpiritReaderObject& get_object(std::string object_name) const {

    std::map<std::string, SpiritReaderObject>::const_iterator it = container_.find(object_name);

    if(it == container_.end())
      throw ParserException() << "Object " << object_name << " not found";

    return it->second;

  }

  SpiritReaderObject& set_object(std::string object_name) {

    return container_[object_name];

  }


  std::map<std::string, SpiritReaderObject>& get_container() {
    return container_;
  }

  void print_objects();

private:

  std::map<std::string, SpiritReaderObject> container_;

};


}

#endif
