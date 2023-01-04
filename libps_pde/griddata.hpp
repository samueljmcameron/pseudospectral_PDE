#ifndef PSPDE_GRIDDATA_HPP
#define PSPDE_GRIDDATA_HPP

#include "fftw_mpi_3darray.hpp"
#include <fftw3-mpi.h>

namespace psPDE {
class gridData {
public:
  
  gridData(const MPI_Comm &, ptrdiff_t, ptrdiff_t,
	   ptrdiff_t, Domain &);

  ~gridData();
  void reverseFlat(int, int &, int &, int &);
  void reverseFlatFourier(int, int &, int &, int &);
  
  
  psPDE::fftw_MPI_3Darray<double> phi;
  psPDE::fftw_MPI_3Darray<double> nonlinear;
  
  Domain &domain;
  Domain ft_domain;
  Domain ft_mod_domain;
  
  psPDE::fftw_MPI_3Darray<std::complex<double>> ft_phi;
  psPDE::fftw_MPI_3Darray<std::complex<double>> ft_nonlinear;
  psPDE::fftw_MPI_3Darray<std::complex<double>> ft_noise;
  
  psPDE::fftw_MPI_3Darray<double> modulus;
  
  fftw_plan forward_phi, backward_phi;
  fftw_plan forward_nonlinear, backward_nonlinear;

};

}

#endif
