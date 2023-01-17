#ifndef PSPDE_GRID_HPP
#define PSPDE_GRID_HPP

#include <string>
#include <memory>

#include "domain.hpp"
#include "fftw_mpi_3darray.hpp"
#include <fftw3-mpi.h>

namespace psPDE {
class Grid {
public:

  
  Grid(const std::vector<std::string> &,const std::vector<std::string> &,
       int, int,MPI_Comm );

  ~Grid();


  Domain domain;
  std::array<ptrdiff_t,3> boxgrid; // global number of grid points {Nx,Ny,Nz}
  std::array<ptrdiff_t,3> ft_boxgrid; // global number of grid points {Nx,Ny,Nz}
  
  const MPI_Comm comm;

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
  
  // grid spacings
  
  double dx() { return domain.period[0]/boxgrid[0];};
  double dy() { return domain.period[1]/boxgrid[1];};
  double dz() { return domain.period[2]/boxgrid[2];};

  // fourier grid spacings
  
  double dqx() { return 2*3.141592/domain.period[0];};
  double dqy() { return 2*3.141592/domain.period[2];};  
  double dqz() { return 2*3.141592/domain.period[1];};

  
  
  void populate(const std::vector<std::string> &);
  


private:

    
  void create_concentration(int , int , int);
  void create_noise(int, int, int);
  void constant_phi(double, double, int);
  void constant_noise(double,double,int);

  void set_real_subdomain();
};

}

#endif
