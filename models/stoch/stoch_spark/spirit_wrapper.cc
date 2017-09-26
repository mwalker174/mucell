#include "spirit_wrapper.h"


#include <fstream>
#include <iterator>
#include <typeinfo>
#include <vector>

// spirit headers
//#define BOOST_SPIRIT_DEBUG
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>


namespace spirit_wrapper
{

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

using qi::int_;
using qi::char_;
using qi::_1;
using ascii::space;
using phoenix::ref;
using boost::spirit::eol;
using boost::spirit::lit;
using boost::spirit::standard::string;
using boost::phoenix::push_back;
using boost::spirit::lexeme;

template <typename Iterator>
struct object_parser
: qi::grammar<Iterator, std::map<std::string, std::string>()> {
  object_parser()
    : object_parser::base_type(query) {

  query = space >>
          *(comment >> object_name >> space  >> qi::lit('{') >> object_content >> qi::lit('}') >> space);

  space = *boost::spirit::ascii::space;
  object_name = *(qi::char_ - '{' - (+boost::spirit::ascii::space));
  object_content = *(qi::char_ - '}');
  comment = *(space >> char_('/') >> char_('/') >> *(qi::char_ - eol) >> eol >> space);

  //BOOST_SPIRIT_DEBUG_NODE(object_name);
  //		BOOST_SPIRIT_DEBUG_NODE(object_content);
  //BOOST_SPIRIT_DEBUG_NODE(value);
  //BOOST_SPIRIT_DEBUG_NODE(query);
  //BOOST_SPIRIT_DEBUG_NODE(first);
}
qi::rule<Iterator, std::map<std::string, std::string>()> query;
qi::rule<Iterator, std::string()> object_name, object_content;
qi::rule<Iterator> space, comment;
};


template <typename Iterator>
struct comma_string_parser
: qi::grammar<Iterator, std::vector<std::string>()> {
  comma_string_parser()
    : comma_string_parser::base_type(query) {

  query = object % ',';

  object = space >> *(qi::char_ - ',') >> space;
  space = *boost::spirit::ascii::space;
  //BOOST_SPIRIT_DEBUG_NODE(object);
  //BOOST_SPIRIT_DEBUG_NODE(query);
}
qi::rule<Iterator, std::vector<std::string>()> query;
qi::rule<Iterator, std::string()> object;
qi::rule<Iterator> space;
};



template <typename Iterator>
struct comma_bracket_string_parser
: qi::grammar<Iterator, std::vector<std::string>()> {
  comma_bracket_string_parser()
    : comma_bracket_string_parser::base_type(query) {

  query = space >> ('[' >> object >> *(qi::lit(']') >> qi::lit(',') >> space >> comment >> qi::lit('[') >> object) >> qi::lit(']')) >> space;
  space = *boost::spirit::ascii::space;
  object = space >> *((qi::char_ - ']' - '[') >> comment) >> space;
  comment = *(space >> char_('/') >> char_('/') >> *(qi::char_ - eol) >> eol >> space);
  //BOOST_SPIRIT_DEBUG_NODE(query);
  //BOOST_SPIRIT_DEBUG_NODE(object);
}
qi::rule<Iterator, std::vector<std::string>()> query;
qi::rule<Iterator, std::string()> object;
qi::rule<Iterator> space, comment;
//qi::rule<Iterator> brackets;
};


template <typename Iterator>
struct keys_and_values
: qi::grammar<Iterator, std::map<std::string, std::string>()> {
  keys_and_values()
    : keys_and_values::base_type(query) {

  query = space >> *( comment >> pair >> space >> *(qi::lit(';') >> space));
  space = *boost::spirit::ascii::space;
  pair  = key >> -('=' >> value);
  //  key   = space >> +(qi::char_ - '=' - '}' - (+boost::spirit::ascii::space)) >> space;
  key   = space >> +(qi::char_ - '=' - '}' - (+boost::spirit::ascii::space))  >> space;
  //*qi::lit('"') >>
  value = space >>
          (qi::lit('"') >> +((qi::char_) - '"') >> qi::lit('"')) |
          +(space >> (qi::char_ - ';') >> comment) >> space;
  comment = *(space >> char_('/') >> char_('/') >> *(qi::char_ - eol) >> eol >> space);


  //BOOST_SPIRIT_DEBUG_NODE(comment);
  BOOST_SPIRIT_DEBUG_NODE(value);
  //BOOST_SPIRIT_DEBUG_NODE(query);
  //BOOST_SPIRIT_DEBUG_NODE(first);
}
qi::rule<Iterator, std::map<std::string, std::string>()> query;
qi::rule<Iterator, std::pair<std::string, std::string>()> pair;
qi::rule<Iterator, std::string()> key, value;
qi::rule<Iterator> space, comment;
};


void
SpiritReader::
read_stream(std::istream& stream)
{

  std::string s;
  {

    std::ostringstream sout;
    std::copy(std::istreambuf_iterator<char>(stream),
              std::istreambuf_iterator<char>(),
              std::ostreambuf_iterator<char>(sout));
    s = sout.str();
  }

  this->read_string(s);

}


void
SpiritReaderObject::
parse_string_vector(std::string& stream, std::vector<std::string>& vec)
{

  comma_string_parser<std::string::iterator> kv;
  std::string::iterator beg = stream.begin();
  std::string::iterator en = stream.end();
  boost::spirit::qi::parse(beg, en, kv, vec);

}


template <typename T>
void
SpiritReaderObject::
get_unsigned_list(std::string vec, std::vector<T> &value)
{

  comma_string_parser<std::string::iterator> kv;
  std::string::iterator beg = vec.begin();
  std::string::iterator en = vec.end();
  std::vector<std::string> vals;
  boost::spirit::qi::parse(beg, en, kv, vals);

  unsigned n_vec = vals.size();
  T val;

  for(unsigned ii = 0; ii < n_vec; ii++) {
    std::stringstream ss(vals[ii]);
    ss >> val;
    value.push_back(val);
  }

}

template
void
SpiritReaderObject::
parse_string_brackets(std::string& stream, std::vector<std::vector<double> >& vec);

template
void
SpiritReaderObject::
parse_string_brackets(std::string& stream, std::vector<std::vector<unsigned> >& vec);


template <typename T>
void
SpiritReaderObject::
parse_string_brackets(std::string& stream, std::vector<std::vector<T> >& vec)
{

  comma_bracket_string_parser<std::string::iterator> kv;
  std::string::iterator beg = stream.begin();
  std::string::iterator en = stream.end();
  std::vector<std::string> brackets;
  boost::spirit::qi::parse(beg, en, kv, brackets);
  unsigned n_brackets = brackets.size();
  vec.resize(n_brackets);

  for(unsigned ii = 0; ii < n_brackets; ii++) {

    std::vector<T>& v = vec[ii];
    get_unsigned_list(brackets[ii], v);

  }

}


void
SpiritReader::
read_string(std::string& stream)
{

  std::map<std::string, std::string> object_content;
  object_parser<std::string::iterator> kv;

  std::string::iterator beg = stream.begin();
  std::string::iterator en = stream.end();
  boost::spirit::qi::parse(beg, en, kv, object_content);

  for(std::map<std::string, std::string>::iterator it = object_content.begin(); it != object_content.end();
      it++) {

    SpiritReaderObject spr;
    spr.name_ = it->first;
    std::string::iterator beg = it->second.begin();
    std::string::iterator en = it->second.end();
    keys_and_values<std::string::iterator> kv;
    std::map<std::string, std::string>& props = spr.get_props();
    boost::spirit::qi::parse(beg, en, kv, props);
    container_[it->first] = spr;

  }
}


void
SpiritReader::
read_file(std::string file_name)
{

  std::ifstream file(file_name.c_str());
  this->read_stream(file);
  file.close();

}


void
SpiritReader::
print_objects()
{

  for(std::map<std::string, SpiritReaderObject>::iterator it = container_.begin(); it != container_.end();
      it++) {

    std::cout << "Object: " << it->first << std::endl;
    std::map<std::string, std::string>& props = it->second.get_props();

    for(std::map<std::string, std::string>::iterator it2 = props.begin(); it2 != props.end();
        it2++) {
      if(std::string::npos != it2->second.find(',')) {
        if(std::string::npos != it2->second.find('[')) {
          std::cout << "  prop: " << it2->first << " ";
          std::vector<std::vector<unsigned> > v;
          SpiritReaderObject::parse_string_brackets(it2->second, v);
          unsigned n_v = v.size();

          for(unsigned ii = 0; ii < n_v; ii++) {
            std::vector<unsigned>& h = v[ii];
            unsigned n_h = h.size();

            for(unsigned jj = 0; jj < n_h; jj++) {
              std::cout << h[jj] << "LL";
            }

            std::cout << "GG";
          }

          std::cout << std::endl;
        } else {
          std::cout << "  prop: " << it2->first << " ";
          std::vector<std::string> v;
          SpiritReaderObject::parse_string_vector(it2->second, v);
          unsigned n_v = v.size();

          for(unsigned jj = 0; jj < n_v; jj++) {
            std::cout << v[jj] << "LL";
          }

          std::cout << std::endl;
        }
      }

      std::cout << "  prop: " << it2->first << " value: |" << it2->second << "|" << std::endl;
    }
  }

}


template <typename T>
void put_vector_to_string(std::ostream& os, const std::vector<std::vector<T> >& value)
{

  unsigned n_value = value.size();

  for(unsigned ii = 0; ii < n_value; ii++) {
    if(ii) {
      os << ",";
    }

    const std::vector<T>& v = value[ii];
    unsigned n_v = v.size();

    if(n_v)
      os << "[" << v[0];

    for(unsigned jj = 1; jj < n_v; jj++) {
      os << "," << v[jj];
    }

    os << "]";
  }
}


template
void put_vector_to_string(std::ostream& os, const std::vector<std::vector<unsigned> >& value);
template
void put_vector_to_string(std::ostream& os, const std::vector<std::vector<double> >& value);

void
SpiritReader::
write_string(std::string& s)
{
  std::ostringstream stream;

  for(std::map<std::string, SpiritReaderObject>::const_iterator it = container_.begin();
      container_.end() != it; ++it) {
    stream << it->first << "{" << std::endl << std::endl;
    const std::map<std::string, std::string>& s = it->second.get_props();

    for(std::map<std::string, std::string>::const_iterator sit = s.begin();
        sit != s.end(); ++sit) {
      stream << "  " << sit->first << " = " << sit->second << ";" << std::endl;
    }

    stream << std::endl << "}" << std::endl << std::endl;
  }

  s = stream.str();
}


}
