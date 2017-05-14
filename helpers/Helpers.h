#ifndef __HELPERS__RYR__
#define __HELPERS__RYR__

#include <mpi.h>
#include <iostream>
#include <cstdio>
#include <sstream>
#include <cmath>

namespace ryr
{

namespace conversion
{

inline	std::string to_string(unsigned value)
{
  std::ostringstream os;
  os << value;
  return os.str();
}

inline	std::string to_string(unsigned value, unsigned n_zeros)
{
  char s[1000];
  std::sprintf(s, "%0*d", n_zeros, value);
  return s;
}

inline	std::string to_string(double value, unsigned n_zeros)
{
  char s[1000];
  std::sprintf(s, "%0*d", n_zeros, value);
  return s;
}


}


namespace mpi
{

class MPIComm
{
public:
  void init(int *argc, char ***argv ) {
    MPI_Init(argc, argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &comm_);
    int rank, size;
    MPI_Comm_size(comm_, &size);
    MPI_Comm_rank(comm_, &rank);

    size_ = size;
    rank_ = rank;
  }
  MPI_Comm comm() const {
    return comm_;
  }
  unsigned n_proc() const {
    return size_;
  }
  unsigned rank() const {
    return rank_;
  }

  void finalize() {
    MPI_Finalize();
  }
private:
  MPI_Comm comm_;
  unsigned rank_, size_;
};

}


namespace info
{

class Nullstream : public std::ostream
{
public:
  Nullstream(): std::ios(0), std::ostream(0) {}
};




class Info
{

private:

  std::ostream *stream_;
  Nullstream nullstream_;
public:

  Info() {}
  void init(const mpi::MPIComm& comm) {

    if (comm.rank() != 0) {
      stream_ = &nullstream_;
    } else {
      stream_ = &std::cout;
    }

  }

  template<class _Tp>
  std::ostream &operator<<(_Tp argument) {
    *stream_ << argument;
    return *stream_;
  }

  std::ostream* &stream() {
    return stream_;
  }

  std::ostream &operator<<(std::ostream & (*f)(std::ostream &)) {
    return f(*stream_);
  }

};


extern Info info;

}

class EventTimer
{

public:

  EventTimer() {
    step_last = -1;
    dt_err = 1e-3;
  }

  bool check_event(double dt_solver, double dt_event, double t) {
    int step1 = (int)std::floor((t - dt_solver * dt_err) / dt_event);
    int step2 = (int)std::floor((t + dt_solver * dt_err) / dt_event);

    if (step1 > step_last) {
      step_last = step1;
      return true;
    } else if (step2 > step_last) {
      step_last = step2;
      return true;
    }

    return false;
  }

private:

  int step_last;
  double dt_err;

};

}

#endif
