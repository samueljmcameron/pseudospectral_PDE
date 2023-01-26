#ifndef PSPDE_GRID_HPP
#define PSPDE_GRID_HPP

#include <string>
#include <memory>
#include <array>

#include "fftw_mpi_3darray.hpp"
#include <fftw3-mpi.h>

namespace psPDE {
class Grid {
public:

  
  Grid(const std::vector<std::string> &,const MPI_Comm &);

  ~Grid();


  std::array<ptrdiff_t,3> boxgrid; // global number of grid points {Nx,Ny,Nz}
  std::array<ptrdiff_t,3> ft_boxgrid; // global number of grid points {Nx,Ny,Nz}
  
  const MPI_Comm &comm;

  std::unique_ptr<psPDE::fftw_MPI_3Darray<double>> phi;
  std::unique_ptr<psPDE::fftw_MPI_3Darray<double>> nonlinear;
  
  std::unique_ptr<psPDE::fftw_MPI_3Darray<std::complex<double>>> ft_phi;
  std::unique_ptr<psPDE::fftw_MPI_3Darray<std::complex<double>>> ft_nonlinear;
  std::unique_ptr<psPDE::fftw_MPI_3Darray<std::complex<double>>> ft_noise;

  std::unique_ptr<psPDE::fftw_MPI_3Darray<double>> ft_phi_modulus;
  std::unique_ptr<psPDE::fftw_MPI_3Darray<double>> ft_nonlinear_modulus;
  std::unique_ptr<psPDE::fftw_MPI_3Darray<double>> ft_noise_modulus;
  
  
  fftw_plan forward_phi, backward_phi;
  fftw_plan forward_nonlinear, backward_nonlinear;
  
  
  
  void populate(const std::vector<std::string> &);


  


private:

    
  void create_concentration(int , int , int);
  void create_noise(int, int, int);
  void constant_phi(double, double, int);
  void constant_noise(double,double,int);
  
};

}

#endif
