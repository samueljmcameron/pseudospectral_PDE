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

  /* concentration (model B) arrays */
  std::unique_ptr<psPDE::fftw_MPI_3Darray<double>> phi;
  std::unique_ptr<psPDE::fftw_MPI_3Darray<double>> chempot;
  std::unique_ptr<psPDE::fftw_MPI_3Darray<double>> gradphi_x,gradphi_y,gradphi_z;
  
  std::unique_ptr<psPDE::fftw_MPI_3Darray<std::complex<double>>> ft_phi;
  std::unique_ptr<psPDE::fftw_MPI_3Darray<std::complex<double>>> ft_chempot;
  std::unique_ptr<psPDE::fftw_MPI_3Darray<std::complex<double>>> ft_noise;
  std::unique_ptr<psPDE::fftw_MPI_3Darray<std::complex<double>>> ft_gradphi_x,ft_gradphi_y,ft_gradphi_z;

  fftw_plan forward_phi, backward_phi;
  fftw_plan forward_chempot, backward_chempot;
  
  fftw_plan backward_gradphi_x;
  fftw_plan backward_gradphi_y;
  fftw_plan backward_gradphi_z;
  
  /* velocity arrays */


  
  std::unique_ptr<psPDE::fftw_MPI_3Darray<std::complex<double>>> ft_Znoise_x,ft_Znoise_y,ft_Znoise_z;
  std::unique_ptr<psPDE::fftw_MPI_3Darray<std::complex<double>>> ft_vtherm_x,ft_vtherm_y,ft_vtherm_z;


  std::unique_ptr<psPDE::fftw_MPI_3Darray<double>> vtherm_x,vtherm_y,vtherm_z;
  std::unique_ptr<psPDE::fftw_MPI_3Darray<double>> v_x, v_y, v_z,vtherm_dot_gradphi;

  std::unique_ptr<psPDE::fftw_MPI_3Darray<std::complex<double>>> ft_v_x,ft_v_y,ft_v_z;
  std::unique_ptr<psPDE::fftw_MPI_3Darray<std::complex<double>>> ft_vtherm_dot_gradphi;
  std::unique_ptr<psPDE::fftw_MPI_3Darray<std::complex<double>>> ft_gradphi_tilde_x,ft_gradphi_tilde_y,ft_gradphi_tilde_z;
  


  
  
  void populate(const std::vector<std::string> &);


  


private:

    
  void create_concentration(int , int , int);
  void create_noise(int, int, int);
  void constant_phi(double, double, int);
  void constant_noise(double,double,int);
  void create_modelH(int, int, int);

  
};

}

#endif
