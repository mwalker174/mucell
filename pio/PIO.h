#ifndef __PIO_MODULE__
#define __PIO_MODULE__

#include <string>
#include <mpi.h>
#include <list>
#include <vector>
#include <sstream>

namespace pio
{

class StringData
{

public:

	void receive_data(MPI_Comm comm, std::string file_name);
	void send_data(MPI_Comm comm);
	void add_data(unsigned index, std::string s);
	void clear();

private:

	static const unsigned limit_ = 10000;
	typedef std::list<std::pair<unsigned, std::string> > list_type;
	list_type data_;
	MPI_Comm comm_;

	unsigned create_buffer(std::vector<char>& buffer, list_type::iterator begin_iterator, list_type::iterator end_iterator);
	unsigned get_n_data();	
	void extract_list_from_buffer(const std::vector<char>& buffer, list_type::iterator& it);

};



class PIO 
{

public:

	void init(MPI_Comm comm);

	void set_file_name(std::string file_name);
	void reset();

	PIO& operator<<(unsigned& ii);
	PIO& operator<<(double& ii);

	void next(unsigned index);
	void commit();
	void close();

	std::ostringstream& get_stream(){
		return stream_;
	}

private: 

	std::ostringstream stream_;

	StringData sd_;
	std::string file_name_;
	MPI_Comm comm_;

};



}

#endif
