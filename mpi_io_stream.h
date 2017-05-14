#ifndef CARDIAC_MPI_IO_STREAM
#define CARDIAC_MPI_IO_STREAM

#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <set>
#include <mpi.h>
#include <stdio.h>
#include <limits.h>
#include <map>

#define CARDIAC_FILE_IO
//#define CARDIAC_MPI_IO
#define CARDIAC_MMAP_IO

namespace ryr
{

namespace input_output
{

#ifdef CARDIAC_MMAP_IO
class pio
{
private:
  MPI_Comm Comm, group_comm;
  int rank, size, n_group, group_rank;
  static const int n_writer = 10000;//64 * 8 * 8;
  MPI_Group group, all_communicator;
  int fd;
  long file_offset, file_size;
  bool write_at_position;
  bool reading;
  std::map<int, bool> processors_done_writing;
  bool append;
  long initial_position_;
protected:

  virtual void get_child_range(std::pair<int, int>& range);

  virtual bool is_writer() {
    return (rank % n_writer == 0);
  }

  virtual int my_writer() {
    return (int(rank / n_writer)) * n_writer;
  }

  void synchronize_file_size();

public:

  pio(): file_offset(0), write_at_position(false), reading(0), append(0), initial_position_(0)  {}
  void general_open(MPI_Comm comm, std::string file_name,
                    std::ios_base::openmode mode);
  void general_close();

  void set_view(long offset) {
    file_offset = offset;
  }

  void initialize_child_processors() {
    if(is_writer()) {
      std::pair<int, int> range;
      get_child_range(range);
      int n_range = range.second - range.first;

      if(reading)
        for(int ii = 0; ii < n_range; ii++)
          processors_done_writing[ii] = false;
      else
        for(int ii = 0; ii < n_range + 1; ii++)
          processors_done_writing[ii] = false;
    }
  }

  void set_write_at_position(unsigned long writing_position) {
    file_offset = writing_position;
    write_at_position = true;
    initialize_child_processors();
  }

  long mpi_read_all(char* buffer, long count);
  void mpi_write_shared(char* buffer, long count);
  bool mpi_write(char* buffer, long count);
  long get_file_size();
  void flash();
};
#endif

#ifdef PELOTON
const int CACHE_FILE_SIZE = 1024000;
#else
const int CACHE_FILE_SIZE = 10240000;
#endif

class mpi_io_streambuffer :  public std::basic_streambuf<char, std::char_traits<char> >
{
  std::vector<char> _inputbuffer;
  std::vector<char> _outputbuffer;
#ifdef CARDIAC_MPI_IO
  std::vector<char> _assync_buffer;
  bool first_assync_writing;
#endif
  MPI_Request request;
  /// \short _PBeg, _PEnd, _GBeg, _GCur, _GEnd are pointer to the current position for output,
  /// end of the output buffer, beginning of the input buffer, current position in the output
  /// buffer, end of the output buffer, respectively.
  char *_PBeg, *_PEnd, *_GBeg, *_GCur, *_GEnd;
  /// \short Pointer to the delimiter in the output buffer. Atomic MPI output operation is from
  /// beginnig of the buffer to the delimiter.
  char *_PDel;

  typedef std::char_traits<char> _Tr;
  MPI_Comm comm;
#ifdef CARDIAC_MPI_IO
  MPI_File fh;
#endif
#ifdef CARDIAC_MMAP_IO
  pio Pio;
#endif
  bool _eof;
  bool is_read;
  bool is_write;

  /// \short If true, no real output. After the file is closed total_size
  /// variable is equal to ouput size
  bool simulate_output;
  /// \short If simulate_output is true, no real output. After the file is closed total_size
  /// variable is equal to ouput size
  unsigned long total_size;

  /// \short If true, MPI_File_write_at is used starting from position which is set by
  /// mpi_io_streambuffer::set_write_position()
  bool need_write_at_position;
  /// \short If need_write_at_position is true, MPI_File_write_at is used
  /// starting from position which is set by mpi_io_streambuffer::set_write_position()
  unsigned long write_position;
  /// \short This variable is used to disable delimiters. See mpi_io_streambuffer::set_delimiter()
  bool use_delimiter;

private:
  /// \short Flush all data to the file. This function is used when the file is closed.
  void flush_output_buffer();
#ifdef CARDIAC_MPI_IO
  long file_size, file_offset;
  long get_file_size() {
    MPI_Offset total_len;
    MPI_File_get_size(fh, &total_len);
    return total_len;
  }
#endif

protected:
  /// \short See documentation to the standard basic_streambuffer class.
  virtual int fill_buff();
  /// \short See documentation to the standard basic_streambuffer class.
  virtual int uflow();
  /// \short See documentation to the standard basic_streambuffer class.
  virtual int underflow();
  /// \short See documentation to the standard basic_streambuffer class.
  virtual int overflow(int = _Tr::eof());
  /// \short See documentation to the standard basic_streambuffer class.
  virtual pos_type seekoff(off_type _Off,
                           std::ios_base::seekdir _Way,
                           std::ios_base::openmode _Which = std::ios_base::in | std::ios_base::out
                          );
  /// \short See documentation to the standard streambuffer class.
  virtual pos_type seekpos(pos_type _Sp,
                           std::ios_base::openmode _Which = std::ios_base::in | std::ios_base::out);

  /// \short Update Get and Put pointers usign pointers of our class.
  void update_pointers() {
    if(is_read)
      this->setg(_GBeg, _GCur, _GEnd);

    if(is_write)
      this->setp(_PBeg, _PEnd);
  }

public:

  void read_all();

  /// \short Constructor.
  mpi_io_streambuffer(MPI_Comm comm, std::string file_name, std::ios_base::openmode mode = std::ios_base::in, bool simulate = false);
  void open(MPI_Comm comm, std::string file_name, std::ios_base::openmode mode = std::ios_base::in,
            bool simulate = false);
  void close() {
    if(is_write) {
      flush_output_buffer();

      if(!simulate_output)
        MPI_Barrier(comm);
    }

#ifdef CARDIAC_MPI_IO

    if(!simulate_output)
      MPI_File_close(&fh);

#endif

#ifdef CARDIAC_MMAP_IO

    if(!simulate_output)
      Pio.general_close();

#endif
  }

  void set_delimiter() {
    _PDel = this->pptr();
  }

  /// \short Calculate offset of the process for writing based on
  /// the size of writing data for each processor
  static unsigned long get_writing_offset(MPI_Comm mpi_comm, unsigned long my_size) {
    int rank, size;
    MPI_Comm_size( mpi_comm, &size );
    MPI_Comm_rank( mpi_comm, &rank );

    std::vector<unsigned long > buf(size);

    MPI_Allgather(
      &my_size,
      1,
      MPI_UNSIGNED_LONG,
      &buf[0],
      1,
      MPI_UNSIGNED_LONG,
      mpi_comm
    );

    unsigned long offset = 0;

    for(int ii = 0; ii < rank; ii++)
      offset += buf[ii];

    return offset;
  }

  /// \short Calculate offset of the process for writing based on
  /// the size of writing data for each processor. This size is calculated in the simulated output
  unsigned long get_offset() {
    return get_writing_offset(comm, total_size);
  }

  /// \short Set offset for writing for the current processor
  void set_write_position(unsigned long my_write_position) {
#ifdef CARDIAC_MPI_IO
    char native[] = "native";
    MPI_File_set_view(fh, my_write_position, MPI_CHAR, MPI_CHAR, native, MPI_INFO_NULL);
    _assync_buffer.resize(CACHE_FILE_SIZE);
#endif

#ifdef CARDIAC_MMAP_IO
    Pio.set_write_at_position(my_write_position);
#endif
    this->write_position = 0;
    need_write_at_position = true;
  }

  /// \short Disable delimiter. This is used when writing with offset rather then MPI_File_write_shared
  void disable_delimiter() {
    use_delimiter = false;
  }

};


/// \short Wrapper for mpi_io_streambuffer to provide stl basic_istream functionality
class mpi_istream : public std::basic_istream<char, std::char_traits<char> >
{
  mpi_io_streambuffer _streambuffer;
public:
  mpi_istream(MPI_Comm comm, std::string file_name, std::ios_base::openmode mode = std::ios_base::in) :
    std::basic_istream<char, std::char_traits<char> >(&_streambuffer), _streambuffer(comm, file_name, mode) {
    clear();
  }

  void open(MPI_Comm comm, std::string file_name, std::ios_base::openmode mode = std::ios_base::in) {
    _streambuffer.open(comm, file_name, mode);
  }

  void close() {
    _streambuffer.close();
  }

  void read_all() {
    _streambuffer.read_all();
  }

  bool is_open() {
    return true;
  }

};

/// \short Wrapper for mpi_io_streambuffer to provide stl basic_ostream functionality
class mpi_ostream : public std::basic_ostream<char, std::char_traits<char> >
{
  mpi_io_streambuffer _streambuffer;
public:
  mpi_ostream(MPI_Comm comm, std::string file_name, bool simulate = false, std::ios_base::openmode mode = std::ios_base::out) :
    std::basic_ostream<char, std::char_traits<char> >(&_streambuffer),
    _streambuffer(comm, file_name, mode, simulate) {
    clear();
  }
  void open(MPI_Comm comm, std::string file_name, std::ios_base::openmode mode = std::ios_base::out) {
    _streambuffer.open(comm, file_name, mode);
  }
  void close() {
    _streambuffer.close();
  }
  void set_delimiter() {
    _streambuffer.set_delimiter();
  }
  unsigned long get_offset() {
    return _streambuffer.get_offset();
  }
  void set_write_position(unsigned long my_write_position) {
    _streambuffer.set_write_position(my_write_position);
  }
  void disable_delimiter() {
    _streambuffer.disable_delimiter();
  }
};


}

}
#endif
