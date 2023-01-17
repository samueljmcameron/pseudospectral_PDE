#include <random>
#include <iostream>

#include "grid.hpp"

#include "iovtk.cpp"

using namespace psPDE;


Grid::Grid(const std::vector<std::string> &v_line,MPI_Comm comm)
  : comm(comm),phi(nullptr),nonlinear(nullptr),
    ft_phi(nullptr),ft_nonlinear(nullptr),ft_noise(nullptr)
    
{


  if (v_line.size() < 4) 
    throw std::invalid_argument("Bad arguments for grid_setup.");

  int Nx = std::stoi(v_line[0]);
  int Ny = std::stoi(v_line[1]);
  int Nz = std::stoi(v_line[2]);

  boxgrid = {Nx,Ny,Nz};
  // since using transposed fftw methods, need ft_boxgrid to swap Nz and Ny
  ft_boxgrid = {Nx,Nz,Ny};

  int iarg = 3;
  
  while (iarg < v_line.size()) {
  
    if (v_line[iarg] == "hybrid") {
      iarg += 1;

      if (iarg == v_line.size())
	throw std::invalid_argument("Need to specify hybrid styles in grid_style.");    


    } else if (v_line[iarg] == "concentration") {
      create_concentration(Nx,Ny,Nz);
      iarg += 1;
    } else if (v_line[iarg] == "noise") {
      create_noise(Nx,Ny,Nz);      
      iarg += 1;
    } else {
      throw std::invalid_argument("Incompatible style type in grid_style.");
    }
  }
  


}


void Grid::populate(const std::vector<std::string> &v_line)
/*
  Read in parameter information. v_line should either have the format

  "constant grid_style average variance seed"

  or

  "read grid_style filename"

  ... more formats can be added later.
*/
{

  
  if (v_line[0] == "constant") {
    
    if (v_line.size() != 5)
      throw std::invalid_argument("Bad arguments (creating constant grid) for grid_populate.");
    
    
    double average = std::stod(v_line[2]);
    double variance = std::stod(v_line[3]);
    int seed = std::stoi(v_line[4]);
    
    if (v_line[1] == "concentration") {
      constant_phi(average,variance, seed);
    } else if (v_line[1] == "noise") {
      constant_noise(average,variance, seed);
    } else {
      throw std::invalid_argument("Bad arguments for grid_setup.");
    }
    
    
  } else if (v_line[0] == "read") {
    
    
    if (v_line.size() != 3)
      throw std::invalid_argument("Bad arguments (reading in data) for grid_populate.");

    if (v_line[1] == "concentration" && phi) {
    
      try {
	ioVTK::readVTKImageData({phi.get()},v_line[2]);
      } catch (std::runtime_error &err) {
	throw std::invalid_argument("invalid grid_setup populateialization file.");
      }
    } else {
      throw std::invalid_argument("Bad arguments for grid_setup.");
    }
  }



}

void Grid::create_concentration(int Nx, int Ny, int Nz)
{


  if (!phi) 
    phi = std::make_unique<fftw_MPI_3Darray<double>
			   >(comm,"concentration",Nx,Ny,Nz);

  if (!nonlinear) 
    nonlinear = std::make_unique<fftw_MPI_3Darray<double>
				 >(comm,"chemPot",Nx,Ny,Nz);
  if (!ft_phi)
    ft_phi = std::make_unique<fftw_MPI_3Darray<std::complex<double>>
			      >(comm,"ft_concentration",Nx,Nz,Ny);

  if (!ft_nonlinear)
    ft_nonlinear = std::make_unique<fftw_MPI_3Darray<std::complex<double>>
				    >(comm,"ft_chemPot",Nx,Nz,Ny);

  forward_phi = fftw_mpi_plan_dft_r2c_3d(Nz,Ny,Nx,phi->data(),
					 reinterpret_cast<fftw_complex*>
					 (ft_phi->data()),
					 comm, FFTW_MPI_TRANSPOSED_OUT);
  
  backward_phi = fftw_mpi_plan_dft_c2r_3d(Nz,Ny,Nx,reinterpret_cast<fftw_complex*>
					  (ft_phi->data()), phi->data(),
					  comm,FFTW_MPI_TRANSPOSED_IN);
  
  forward_nonlinear = fftw_mpi_plan_dft_r2c_3d(Nz,Ny,Nx,
					       nonlinear->data(),
					       reinterpret_cast<fftw_complex*>
					       (ft_nonlinear->data()),
					       comm, FFTW_MPI_TRANSPOSED_OUT);
  
  backward_nonlinear = fftw_mpi_plan_dft_c2r_3d(Nz,Ny,Nx,
						reinterpret_cast<fftw_complex*>
						(ft_nonlinear->data()),
						nonlinear->data(),comm,
						FFTW_MPI_TRANSPOSED_IN);
  return;
  
}

void Grid::create_noise(int Nx, int Ny, int Nz)
{

  if (!ft_noise)
    ft_noise = std::make_unique<fftw_MPI_3Darray<std::complex<double>>
				>(comm,"ft_noise",Nx,Nz,Ny);
  return;
}

void Grid::constant_phi(double average, double variance, int seed)
{

  if (!phi)
    throw std::runtime_error("Cannot create constant phi.");
  std::uniform_real_distribution<double> real_dist(-0.5,0.5);
    
  std::mt19937 gen;
    
  gen.seed(seed);
	  
  for (int i = 0; i < phi->Nz(); i++) 
    for (int j = 0; j < phi->Ny(); j++) 
      for (int k = 0; k < phi->Nx(); k++) 
	(*phi)(i,j,k) = average + variance*real_dist(gen);
  
  return;
}

void Grid::constant_noise(double average, double variance, int seed)
{
  if (!ft_noise)
    throw std::runtime_error("Cannot create constant phi.");
  
  std::uniform_real_distribution<double> real_dist(-0.5,0.5);
    
  std::mt19937 gen;
    
  gen.seed(seed);
	  
  for (int i = 0; i < ft_noise->Nz(); i++) 
    for (int j = 0; j < ft_noise->Ny(); j++) 
      for (int k = 0; k < ft_noise->Nx(); k++) {
	(*ft_noise)(i,j,k).real(average + variance*real_dist(gen));
	(*ft_noise)(i,j,k).imag(variance*real_dist(gen));
      }
  
  return;
}



Grid::~Grid() {

  
  fftw_destroy_plan(forward_phi);
  fftw_destroy_plan(backward_phi);
  fftw_destroy_plan(forward_nonlinear);
  fftw_destroy_plan(backward_nonlinear);

}


