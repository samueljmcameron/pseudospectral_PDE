#ifndef PSPDE_CONJPLANE_HPP
#define PSPDE_CONJPLANE_HPP

#include <complex>
#include <vector>

#include "fftw_mpi_3darray.hpp"
#include "smatrix.hpp"
#include "integrator.hpp"

namespace psPDE{
class ConjPlane
{

  const int mpi_size;
  const int id;
  const MPI_Comm comm;

  const int cNy;
  
  bool equal_flag;
  int divider;
  const int nreq;
  MPI_Request *block_requests;
  MPI_Request *line_requests;
  int loopstart;         // 1 or 0 depending on if id == 0 or not

  
  // arrays for transferring via mpi
  std::vector<sMatrix<std::complex<double>>> blocks;
  std::vector<sMatrix<std::complex<double>>> lines;

  void lefthalf_equal(Integrator &,int);
  void lefthalf_unequal(Integrator &,int);
  void righthalf_equal(Integrator &,int);
  void righthalf_unequal(Integrator &,int);
  void middle_odd_equal(Integrator &,int);
  void middle_odd_unequal(Integrator &,int);
  void middle_even_equal(Integrator &,int);
  void middle_even_unequal(Integrator &,int);


  void line_lefthalf_equal(Integrator &,int);
  void line_lefthalf_unequal(Integrator &,int);
  void line_righthalf_equal(Integrator &,int);
  void line_righthalf_unequal(Integrator &,int);
  void line_middle_odd_equal(Integrator &,int);
  void line_middle_odd_unequal(Integrator &,int);
  void line_middle_even_equal(Integrator &,int);
  void line_middle_even_unequal(Integrator &,int);

  
  
public:

  ConjPlane(int,int,MPI_Comm,std::vector<int>,int);
  ~ConjPlane();

  bool is_equal() {
    return equal_flag;
  }

  int get_divider() {
    return divider;
  }
  
  void single(Integrator &,int);
  
  void first(Integrator &,int);
  void lefthalf(Integrator &,int);
  void righthalf(Integrator &,int);
  void middle_odd(Integrator &,int);
  void middle_even(Integrator &,int);
  void last(Integrator &,int);

  void line_single(Integrator &,int);
  
  void line_first(Integrator &,int);
  void line_lefthalf(Integrator &,int);
  void line_righthalf(Integrator &,int);
  void line_middle_odd(Integrator &,int);
  void line_middle_even(Integrator &,int);
  void line_last(Integrator &,int);


};

};
#endif
