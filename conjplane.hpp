#ifndef CONJPLANE_HPP
#define CONJPLANE_HPP

#include <complex>
#include <vector>

#include "fftw_mpi_3darray.hpp"
#include "smatrix.hpp"
#include "timestep.hpp"

class ConjPlane
{

  const int mpi_size;
  const int id;
  const MPI_Comm comm;

  const int cNz;
  
  bool equal_flag;
  int divider;
  const int nreq;
  MPI_Request *block_requests;
  MPI_Request *line_requests;
  int loopstart;         // 1 or 0 depending on if id == 0 or not

  
  // arrays for transferring via mpi
  sMatrix<std::complex<double>> block0,block1,block2,block3;
  sMatrix<std::complex<double>> line0,line1,line2,line3;

  void lefthalf_equal(TimeStep &,int);
  void lefthalf_unequal(TimeStep &,int);
  void righthalf_equal(TimeStep &,int);
  void righthalf_unequal(TimeStep &,int);
  void middle_odd_equal(TimeStep &,int);
  void middle_odd_unequal(TimeStep &,int);
  void middle_even_equal(TimeStep &,int);
  void middle_even_unequal(TimeStep &,int);


  void line_lefthalf_equal(TimeStep &,int);
  void line_lefthalf_unequal(TimeStep &,int);
  void line_righthalf_equal(TimeStep &,int);
  void line_righthalf_unequal(TimeStep &,int);
  void line_middle_odd_equal(TimeStep &,int);
  void line_middle_odd_unequal(TimeStep &,int);
  void line_middle_even_equal(TimeStep &,int);
  void line_middle_even_unequal(TimeStep &,int);

  
  
public:

  ConjPlane(int,int,MPI_Comm,std::vector<int>,int);
  ~ConjPlane();

  bool is_equal() {
    return equal_flag;
  }

  int get_divider() {
    return divider;
  }
  
  void single(TimeStep &,int);
  
  void first(TimeStep &,int);
  void lefthalf(TimeStep &,int);
  void righthalf(TimeStep &,int);
  void middle_odd(TimeStep &,int);
  void middle_even(TimeStep &,int);
  void last(TimeStep &,int);

  void line_single(TimeStep &,int);
  
  void line_first(TimeStep &,int);
  void line_lefthalf(TimeStep &,int);
  void line_righthalf(TimeStep &,int);
  void line_middle_odd(TimeStep &,int);
  void line_middle_even(TimeStep &,int);
  void line_last(TimeStep &,int);


};

#endif
