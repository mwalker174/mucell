#include "mpi_io_stream.h"
#include "ExceptionBase.h"
#include "helpers/Helpers.h"

#ifdef CARDIAC_MMAP_IO
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <numeric>
#endif


namespace ryr
{

namespace input_output
{

mpi_io_streambuffer::pos_type mpi_io_streambuffer::seekpos(
  mpi_io_streambuffer::pos_type _Sp,
  std::ios_base::openmode _Which)
{
  return this->seekoff(_Sp, std::ios_base::beg, _Which);
}

mpi_io_streambuffer::pos_type mpi_io_streambuffer::seekoff(mpi_io_streambuffer::off_type _Off,
    std::ios_base::seekdir _Way,
    std::ios_base::openmode _Which
                                                          )
{
  /*
  int whence = MPI_SEEK_SET;

  switch(_Way) {
  case std::ios_base::beg:
    whence = MPI_SEEK_SET;
    break;
  case std::ios_base::cur:
    whence = MPI_SEEK_CUR;
    break;
  case std::ios_base::end:
    whence = MPI_SEEK_END;
    break;
  default:
    ;
  };

  #ifdef CARDIAC_MPI_IO
  int ierr;

  ierr = MPI_File_seek(fh, _Off, whence);

  if(ierr != MPI_SUCCESS)
    throw  std::ios_base::failure("");

  #endif
  */
  if(is_write) {
    _PBeg = &_outputbuffer[0];
    unsigned size = _outputbuffer.size();
    _PEnd = &_outputbuffer[size - 1];
  }

  if(is_read) {
    _GBeg = _GCur = &_inputbuffer[0], _GEnd = &_inputbuffer[0];
  }

  update_pointers();
#ifdef CARDIAC_MPI_IO
  MPI_Offset offset;
  MPI_File_get_position(fh, &offset);
#endif

#ifdef CARDIAC_MMAP_IO
  mpi_io_streambuffer::pos_type offset = 0;
#endif
  return offset;
}

void
mpi_io_streambuffer::open(MPI_Comm comm,
                          std::string file_name,
                          std::ios_base::openmode mode,
                          bool simulate)
{



  this->use_delimiter = false;

  this->need_write_at_position = false;
  this->write_position = 0;

  this->total_size = 0;
  this->simulate_output = simulate;

  this->_eof = false;
  this->comm = comm;

  /*
    MPI_MODE_RDONLY --- read only,
    MPI_MODE_RDWR --- reading and writing,
    MPI_MODE_WRONLY --- write only,
    MPI_MODE_CREATE --- create the file if it does not exist,
    MPI_MODE_EXCL --- error if creating file that already exists,
    MPI_MODE_DELETE_ON_CLOSE --- delete file on close,
    MPI_MODE_UNIQUE_OPEN --- file will not be concurrently opened elsewhere,
    MPI_MODE_SEQUENTIAL --- file will only be accessed sequentially,
    MPI_MODE_APPEND --- set initial position of all file pointers to end of file.
  */

  int amode = MPI_MODE_RDWR;

  if(
    !((mode & std::ios_base::in) ||
      (mode & std::ios_base::out))
  ) {
    throw exceptions::ExceptionBase() << "wrong file reading options, should be ::in or ::out";
  }


  if(
    ((mode & std::ios_base::in) &&
     (mode & std::ios_base::out))
  )
    throw exceptions::ExceptionBase() << "simulteneous reading and writing is not supported";



  is_read = is_write = true;

  if(!(mode & std::ios_base::in)) {
    is_read = false;
    amode = MPI_MODE_WRONLY;
    amode = amode | MPI_MODE_CREATE;

    if(mode & std::ios_base::app)
      amode = amode | MPI_MODE_APPEND;
  }

  if(!(mode & std::ios_base::out)) {
    is_write = false;
    amode = MPI_MODE_RDONLY;
  }

  if(mode & std::ios_base::app)
    amode = amode | MPI_MODE_APPEND;

  if(!simulate) {
#ifdef CARDIAC_MPI_IO
    this->first_assync_writing = true;
    this->file_offset = 0;
    int ierr = MPI_File_open(comm, const_cast<char*>(file_name.c_str()), amode, MPI_INFO_NULL, &fh);


    if(ierr != MPI_SUCCESS) {
      throw exceptions::ExceptionBase() << "Unable to open file";
    }

    file_size = get_file_size();
#endif

#ifdef CARDIAC_MMAP_IO
    Pio.general_open(comm, file_name, mode);
#endif
  }



  if(is_read) {
    _inputbuffer.resize(CACHE_FILE_SIZE);
    _GBeg = _GCur = _GEnd = &_inputbuffer[0];
  }

  if(is_write) {
    if(_outputbuffer.size() == 0) {
      _outputbuffer.resize(CACHE_FILE_SIZE);
    }

    unsigned size = _outputbuffer.size();
    _PBeg = &_outputbuffer[0];
    _PEnd = &_outputbuffer[size - 1];
    _PDel = NULL;
  }

  update_pointers();
}

mpi_io_streambuffer::mpi_io_streambuffer(MPI_Comm comm, std::string file_name, std::ios_base::openmode mode, bool simulate)
{
  this->open(comm, file_name, mode, simulate);
}

int mpi_io_streambuffer::uflow()
{

  if(gptr() && gptr() < egptr()) {
    int res = _Tr::to_int_type(*gptr());
    _GCur++;
    update_pointers();
    return std::char_traits<char>::not_eof(res);
  }

  int fill = fill_buff();
  _GCur++;
  update_pointers();
  return std::char_traits<char>::not_eof(fill);
}

int mpi_io_streambuffer::underflow()
{

  if(gptr() && gptr() < egptr()) {
    return _Tr::to_int_type(*gptr());
  }

  return std::char_traits<char>::not_eof(fill_buff());
}

void mpi_io_streambuffer::read_all()
{
  while(!_eof)
    fill_buff();
}

int mpi_io_streambuffer::fill_buff()
{

  if(_eof)
    return _Tr::eof();

  char c;
  int count;

#ifdef CARDIAC_MPI_IO
  count = CACHE_FILE_SIZE;

  //if(file_offset + count > file_size)
  //  count = file_size - file_offset;

  MPI_Status status;
  int ierr;
  std::memset( &status, 0xff, sizeof(MPI_Status));
  ierr = MPI_File_read_all(fh, &_inputbuffer[0], count,
                           MPI_CHAR, &status);

  file_offset += count;

  if(ierr != MPI_SUCCESS)
    throw  std::ios_base::failure("");

  if(MPI_SUCCESS != MPI_Get_count(&status, MPI_CHAR, &count))
    throw  std::ios_base::failure("");

#endif

#ifdef CARDIAC_MMAP_IO
  count = Pio.mpi_read_all(&_inputbuffer[0], CACHE_FILE_SIZE);
#endif


  if(count == 0) {
    _eof = true;
    return _Tr::eof();
  }

  if(count < CACHE_FILE_SIZE)
    _eof = true;

  _GEnd = &_inputbuffer[count - 1] + 1;
  _GCur = _GBeg = &_inputbuffer[0];
  update_pointers();

  c = *gptr();
  return _Tr::to_int_type(c);
}


void mpi_io_streambuffer::flush_output_buffer()
{

  char* pointer_zero_element = &_outputbuffer[0];
  int count = this->pptr() - pointer_zero_element;

  if(simulate_output)
    total_size += count;
  else {
    if(need_write_at_position) {

#ifdef CARDIAC_MPI_IO
      MPI_Status status;

      if(!first_assync_writing)
        MPI_Wait(&request, &status);

      if(count)
        MPI_File_write(fh, &_outputbuffer[0],
                       count, MPI_CHAR,
                       &status);

#endif

#ifdef CARDIAC_MMAP_IO
      Pio.mpi_write(&_outputbuffer[0], count);
      Pio.flash();
#endif
      write_position += unsigned(count);
    } else
#ifdef CARDIAC_MPI_IO
    {
      MPI_Status status;
      MPI_File_write_shared(fh, &_outputbuffer[0], count, MPI_CHAR, &status);
    }

#endif

#ifdef CARDIAC_MMAP_IO
    Pio.mpi_write_shared(&_outputbuffer[0], count);
#endif

  }
}

int mpi_io_streambuffer::overflow(int c)
{

  unsigned size = _outputbuffer.size();

  if(simulate_output) { // if we simulate output, no writing is needed. We only count size
    char* pointer_zero_element = &_outputbuffer[0];

    if(!use_delimiter)
      _PDel = this->pptr() - 1;

    int count = _PDel - pointer_zero_element + 1;
    total_size += count;
    _PBeg = &_outputbuffer[0];
    _PEnd = &_outputbuffer[size - 1];
    update_pointers();
    return _Tr::not_eof(c);
  }


  //We always set the last element of the vector as the end of the output buffer
  // to keep space for single character passed to overflow()
  _outputbuffer[size - 1] = c;

  //if we don't use delimiters then ouput the whole buffer
  //by setting _PDel to the end of the buffer
  if(!use_delimiter)
    _PDel = this->pptr() - 1;

  //We are usign delimiters. Therefore we output buffer from beginning to the last set delimiter.
  //Then we copy the rest of the buffer to the beginning of the buffer and update pointers.
  if(_PDel) {
    char* pointer_zero_element = &_outputbuffer[0];
    int count = _PDel - pointer_zero_element + 2;
    int pptr_offset = this->pptr() - pointer_zero_element + 1;




    if(need_write_at_position) {
#ifdef CARDIAC_MPI_IO
      MPI_Status status;

      if(!first_assync_writing)
        MPI_Wait(&request, &status);

      first_assync_writing = false;
      std::copy(&_outputbuffer[0], &_outputbuffer[0] + count, &_assync_buffer[0]);

      MPI_File_iwrite(fh, &_assync_buffer[0],
                      count, MPI_CHAR,
                      &request);
#endif

#ifdef CARDIAC_MMAP_IO
      Pio.mpi_write(&_outputbuffer[0], count);
#endif

      write_position += unsigned(count);
    } else

#ifdef CARDIAC_MPI_IO
    {
      MPI_Status status;
      MPI_File_write_shared(fh, &_outputbuffer[0], count, MPI_CHAR, &status);
    }

#endif

#ifdef CARDIAC_MMAP_IO
    Pio.mpi_write_shared(&_outputbuffer[0], count);
#endif

    std::copy(&_outputbuffer[0] + count, &_outputbuffer[0] + size, &_outputbuffer[0]);
    pptr_offset -= count;

    _PBeg = &_outputbuffer[pptr_offset];
    _PDel = NULL;
  } else {
    //If delimiter is not set then we have to keep buffering.
    //Therefore, we increase the buffer size and update all put pointers
    char* pointer_zero_element = &_outputbuffer[0];
    int pptr_offset = this->pptr() - pointer_zero_element;
    _outputbuffer.resize(size + CACHE_FILE_SIZE);
#ifdef CARDIAC_MPI_IO
    MPI_Status status;

    if(need_write_at_position) {
      if(!first_assync_writing)
        MPI_Wait(&request, &status);

      _assync_buffer.resize(size + CACHE_FILE_SIZE);
    }

#endif
    _PBeg = &_outputbuffer[pptr_offset];
    size = _outputbuffer.size();
  }

  _PEnd = &_outputbuffer[size - 1];
  update_pointers();
  //writing
  return _Tr::not_eof(c);

}


#ifdef CARDIAC_MMAP_IO

void pio::flash()
{

  bool writer = is_writer();

  // writing post actions
  if(writer) {
    while(mpi_write(NULL, 0));
  } else
    mpi_write(NULL, 0);

}

bool pio::mpi_write(char* buffer, long count)
{
  bool writer = is_writer();
  bool childs_have_output = false;

  if(writer) {

    std::pair<int, int> range;
    this->get_child_range(range);
    int n_range = range.second - range.first;
    std::vector<long> count_and_offset((n_range + 1) * 2);

    unsigned non_empty_processor_counter = 0;

    if(n_range != 0) {


      std::vector<MPI_Request> req(n_range);
      std::vector<MPI_Status> stat(n_range);

      //reading size of the data and offset for each child process
      non_empty_processor_counter = 0;

      for(std::map<int, bool>::const_iterator it = processors_done_writing.begin();
          it != processors_done_writing.end(); it++) {

        if(it->first != n_range) {

          if(!it->second) {
            int ii = it->first;
            MPI_Irecv(
              &count_and_offset[non_empty_processor_counter * 2],
              2,
              MPI_LONG,
              ii + range.first,
              0,
              Comm,
              &req[non_empty_processor_counter]);
            non_empty_processor_counter++;
          }
        }

        if(non_empty_processor_counter != 0)
          MPI_Waitall(
            non_empty_processor_counter,
            &req[0],
            &stat[0]);
      }
    }

    //adding writer to the array to process within the loop for childs
    count_and_offset[non_empty_processor_counter * 2] = count;
    count_and_offset[non_empty_processor_counter * 2 + 1] = file_offset;

    // setting the value to add current process to the map
    processors_done_writing[n_range] = (buffer == NULL);

    //      processors_done_writing_copy.insert(processors_done_writing.size(), std::pair<int, bool>(rank, count == 0));
    non_empty_processor_counter = 0;

    for(std::map<int, bool>::iterator it = processors_done_writing.begin();
        it != processors_done_writing.end(); it++) {

      if(!it->second) {
        int  ii = it->first;
        long halo_count = count_and_offset[non_empty_processor_counter * 2];

        if((halo_count == 0) && (ii != n_range)) {
          it->second = true;
          non_empty_processor_counter++;
          continue;
        }

        long halo_offset = count_and_offset[non_empty_processor_counter * 2 + 1];
        non_empty_processor_counter++;

#ifdef  CARDIAC_FILE_IO
        std::vector<char> vector_data(halo_count);
        char* data = &vector_data[0];
#else
        long end = halo_count + halo_offset - 1;

        if(lseek(fd, end + initial_position_, SEEK_SET) == -1) {
          throw exceptions::ExceptionBase() << "lseek failed";
        }

        if(write(fd, "", 1) != 1) {
          throw exceptions::ExceptionBase() << "write failed";
        }


        long page_size = sysconf(_SC_PAGE_SIZE);
        long page_offset = (long((halo_offset) / page_size)) * page_size;
        long page_reminder = (halo_offset) % page_size;
        char* data_map = (char*)mmap(NULL, halo_count + page_reminder, PROT_READ | PROT_WRITE, MAP_SHARED, fd, page_offset);

        if (data_map == MAP_FAILED) {
          throw exceptions::ExceptionBase() << "map failed";
        }

        char* data = data_map + page_reminder;
#endif

        if(ii != n_range) {
          MPI_Request request;
          MPI_Status status;
          MPI_Irecv(
            data,
            halo_count,
            MPI_CHAR,
            ii + range.first,
            0,
            Comm,
            &request);
          childs_have_output = true;
          MPI_Wait(&request, &status);
        } else if(count != 0) {
          std::memcpy(data, buffer, count);
          file_offset += count;
        }

#ifdef  CARDIAC_FILE_IO

        if(halo_count) {
          if(lseek(fd, halo_offset + initial_position_, SEEK_SET) == -1) {
            throw exceptions::ExceptionBase() << "lseek failed";
          }

          if(write(fd, data, halo_count) != halo_count) {
            throw exceptions::ExceptionBase() << "write failed";
          }
        }

#else
        munmap(data_map, halo_count + page_reminder);
#endif
      }
    }
  } else {
    int my = my_writer();
    MPI_Request request_size, request_buffer;
    MPI_Status status;
    long count_and_offset[2] = {count, file_offset};
    MPI_Isend(
      &count_and_offset,
      2,
      MPI_LONG,
      my,
      0,
      Comm,
      &request_size);

    if(count != 0) {
      MPI_Isend(
        buffer,
        count,
        MPI_CHAR,
        my,
        0,
        Comm,
        &request_buffer);
      file_offset += count;
      MPI_Wait(&request_size, &status);
      MPI_Wait(&request_buffer, &status);
    } else
      MPI_Wait(&request_size, &status);
  }

  return childs_have_output;

}

////////////////
void pio::mpi_write_shared(char* buffer, long count)
{
  bool writer = is_writer();

  if(writer) {
    std::pair<int, int> range;
    this->get_child_range(range);
    int n_range = range.second - range.first;

    if(n_range != 0) {
      std::vector<long> count_receive(n_range);
      std::vector<MPI_Request> req(n_range);
      std::vector<MPI_Status> stat(n_range);

      for(int ii = 0; ii < n_range; ii++) {

        MPI_Irecv(
          &count_receive[ii],
          1,
          MPI_LONG,
          ii + range.first,
          0,
          Comm,
          &req[ii]);
      }

      MPI_Waitall(
        n_range,
        &req[0],
        &stat[0]);
      long zero = 0;
      long sum = std::accumulate(count_receive.begin(), count_receive.end(), zero);
      std::vector<long> rbuf(n_group);
      MPI_Allgather(&sum, 1, MPI_LONG, &rbuf[0], n_group, MPI_LONG, group_comm);
      int offset = std::accumulate(&rbuf[0], &rbuf[group_rank], zero);
      int end = offset + rbuf[group_rank] - 1;
      int total = std::accumulate(rbuf.begin(), rbuf.end(), zero);

      if(lseek(fd, end + file_offset + initial_position_, SEEK_SET) == -1) {
        throw exceptions::ExceptionBase() << "lseek failed";
      }

      if(write(fd, "", 1) != 1) {
        throw exceptions::ExceptionBase() << "write failed";
      }


      long page_size = sysconf(_SC_PAGE_SIZE);
      long page_offset = ((long)((offset + file_offset) / page_size)) * page_size;
      long page_reminder = (offset + file_offset) % page_size;
      char* data_map = (char*)mmap(NULL, sum + page_reminder, PROT_READ | PROT_WRITE, MAP_SHARED, fd, page_offset);

      if (data_map == MAP_FAILED) {
        throw exceptions::ExceptionBase() << "map failed";
      }


      char* data = data_map + page_reminder;

      for(int ii = 0; ii < n_range; ii++) {
        MPI_Irecv(
          data,
          count_receive[ii],
          MPI_CHAR,
          ii + range.first,
          0,
          Comm,
          &req[ii]);
      }

      MPI_Waitall(
        n_range,
        &req[0],
        &stat[0]);
      munmap(data_map, sum + page_reminder);
      file_offset += total;
    }

  } else {
    int my = my_writer();
    MPI_Request request_size, request_buffer;
    MPI_Status status;
    MPI_Isend(
      &count,
      1,
      MPI_LONG,
      my,
      0,
      Comm,
      &request_size);
    MPI_Isend(
      buffer,
      count,
      MPI_CHAR,
      my,
      0,
      Comm,
      &request_buffer);
    MPI_Wait(&request_size, &status);
    MPI_Wait(&request_buffer, &status);
  }
}


long pio::mpi_read_all(char* buffer, long count)
{

  //calculating actual number of bytes to read
  if(file_offset + count >= file_size)
    count = file_size - file_offset;

  if(count == 0) {
    return 0;
  }

  bool writer = is_writer();

  if(writer) {

    std::pair<int, int> range;
    this->get_child_range(range);
    int n_range = range.second - range.first;

    std::vector<MPI_Request> req;
    std::vector<MPI_Status> stat;
    unsigned processor_counter = 0;

    if(n_range != 0) {

      req.resize(n_range);
      stat.resize(n_range);
      std::vector<int> done(n_range);


      //checking if any processor called close() and finish reading
      processor_counter = 0;

      for(std::map<int, bool>::const_iterator it = processors_done_writing.begin();
          it != processors_done_writing.end(); it++) {

        if(!it->second) {
          int ii = it->first;
          MPI_Irecv(
            &done[processor_counter],
            1,
            MPI_INT,
            ii + range.first,
            0,
            Comm,
            &req[processor_counter]);
          processor_counter++;
        }
      }

      MPI_Waitall(
        processor_counter,
        &req[0],
        &stat[0]);


      // setting flag for the processors that finished reading
      processor_counter = 0;
      unsigned processor_reading = 0;

      for(std::map<int, bool>::iterator it = processors_done_writing.begin();
          it != processors_done_writing.end(); it++) {
        if(!it->second) {
          if(done[processor_counter]) {
            it->second = true;
            processor_reading++;
          }

          processor_counter++;
        }
      }

      //if all processors are done including the writer, return count = 0
      if((!buffer) && (processor_counter == processor_reading))
        return 0;
    } else if(!buffer) return 0;

#ifdef  CARDIAC_FILE_IO
    std::vector<char> vector_data(count);
    char* data = &vector_data[0];

    char* data_tmp = data;
    int read_count = count, actual_read_count = 0;

    do {
      data_tmp += actual_read_count;
      actual_read_count = read(fd, data_tmp, read_count);
      read_count -= actual_read_count;
    } while ((actual_read_count != -1) && (read_count > 0));

    if(actual_read_count == -1) {
      throw exceptions::ExceptionBase() << "read failed";
    }

    data = &vector_data[0];

#else
    //assigning pointer to the file
    long page_size = sysconf(_SC_PAGE_SIZE);
    long page_offset = long(file_offset / page_size) * page_size;
    long page_reminder = file_offset % page_size;

    char* data_map = (char*)mmap(NULL, count + page_reminder, PROT_READ, MAP_SHARED, fd, page_offset);

    if (data_map == MAP_FAILED) {
      throw exceptions::ExceptionBase() << "map failed";
    }


    char* data = data_map + page_reminder;
#endif
    file_offset += count;

    //sending only to the processors which are waiting for reading
    if(n_range != 0) {
      processor_counter = 0;

      for(std::map<int, bool>::const_iterator it = processors_done_writing.begin();
          it != processors_done_writing.end(); it++) {
        if(!it->second) {
          int ii = it->first;
          MPI_Isend(
            data,
            count,
            MPI_CHAR,
            ii + range.first,
            0,
            Comm,
            &req[processor_counter]);
          processor_counter++;
        }
      }
    }

    if(buffer)
      std::memcpy(buffer, data, count * sizeof(char));

    if(n_range != 0)
      MPI_Waitall(
        processor_counter,
        &req[0],
        &stat[0]);

#ifndef CARDIAC_FILE_IO
    munmap(data_map, count + page_reminder);
#endif

  } else {

    MPI_Request request_buffer;
    MPI_Request request_done;
    MPI_Status status;
    int my = my_writer();

    int done = 0;
    MPI_Isend(
      &done,
      1,
      MPI_INT,
      my,
      0,
      Comm,
      &request_done);

    MPI_Irecv(
      buffer,
      count,
      MPI_CHAR,
      my,
      0,
      Comm,
      &request_buffer);

    file_offset += count;
    MPI_Wait(&request_done, &status);
    MPI_Wait(&request_buffer, &status);
  }

  return count;
}

void pio::synchronize_file_size()
{
  bool writer = this->is_writer();

  if(writer) {

    std::pair<int, int> range;
    this->get_child_range(range);
    int n_range = range.second - range.first;

    if(n_range != 0) {

      MPI_Request req[n_range];

      for(int ii = 0; ii < n_range; ii++) {

        MPI_Isend(
          &file_size,
          1,
          MPI_LONG,
          ii + range.first,
          0,
          Comm,
          &req[ii]);
      }

      MPI_Waitall(n_range, req, MPI_STATUSES_IGNORE);
    }

  } else {

    MPI_Request request;
    MPI_Status status;
    int my = my_writer();

    MPI_Irecv(
      &file_size,
      1,
      MPI_LONG,
      my,
      0,
      Comm,
      &request);
    MPI_Wait(&request, &status);

  }
}

////////////////////////
void pio::general_open(MPI_Comm comm, std::string file_name, std::ios_base::openmode mode)
{
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  Comm = comm;

  /*
  MPI_Comm_group(comm, &all_communicator);

  n_group = (size-1)/n_writer+1;
  std::vector<int> ranks(n_group);
  for(unsigned ii=0; ii< n_group; ii++)
    ranks[ii] = ii*n_writer;

  MPI_Group_incl(all_communicator, n_group, &ranks[0], &group);
  MPI_Comm_create(comm, group, &group_comm);
  MPI_Comm_rank(group_comm, &group_rank);
  */

  reading = mode & std::ios_base::in;

  if(is_writer()) {
    if(!reading) {
      int oflags = O_RDWR | O_CREAT;

      if(mode & std::ios_base::app)
        oflags = oflags | O_APPEND;
      else
        oflags = oflags | O_TRUNC;

      if ((fd = open (file_name.c_str(), oflags, S_IREAD | S_IWRITE)) < 0) {
        throw exceptions::ExceptionBase() << "can't open file for writing";
      }

      if((mode & std::ios_base::app) && (initial_position_ = (lseek(fd, 0, SEEK_CUR)) == -1)) {
        throw exceptions::ExceptionBase() << "lseek failed";
      }

    } else {
      info::info << "file opened " << file_name << std::endl;

      if ((fd = open (file_name.c_str(), O_RDONLY)) < 0) {
        throw exceptions::ExceptionBase() << "can't open file for reading";
      }

      initialize_child_processors();
    }

    if(reading)
      file_size = get_file_size();

    initialize_child_processors();
  }

  if(reading)
    synchronize_file_size();
}

/////////////////////////
void pio::general_close()
{
  //  MPI_Group_free(&group);
  // MPI_Group_free(&all_communicator);

  bool writer = this->is_writer();

  // reading post actions
  if( reading ) {
    if(!writer) { //child sends to the writer that it is done
      if(file_size  != file_offset) {
        MPI_Request request_buffer;
        MPI_Status status_buffer;
        int my = my_writer();

        int done = 1;
        MPI_Isend(
          &done,
          1,
          MPI_INT,
          my,
          0,
          Comm,
          &request_buffer);
        MPI_Wait(&request_buffer, &status_buffer);
      }

    } else {
      while(mpi_read_all(NULL, CACHE_FILE_SIZE) != 0);
    }
  }

  MPI_Barrier(Comm);

  //  std::cout << file_size <<  " " << file_offset << " " << reading<< std::endl;
  if(is_writer())
    close(fd);
}

void pio::get_child_range(std::pair<int, int>& range)
{
  range.first = my_writer() + 1;
  range.second = (int(rank / n_writer) + 1) * n_writer;
  range.first = std::min(size, range.first);
  range.second = std::min(size, range.second);

}

long pio::get_file_size()
{
  struct stat fileStat;

  if(fstat(fd, &fileStat) < 0) {
    throw exceptions::ExceptionBase() << "can't get file size";
  }


  return fileStat.st_size;

}

#endif
}

}
