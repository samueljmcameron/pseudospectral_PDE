#include <iostream>

#include "griddata.hpp"

using namespace psPDE;


gridData::gridData(const MPI_Comm &comm, ptrdiff_t Nx, ptrdiff_t Ny,
		   ptrdiff_t Nz, Domain &domain_o)
  : phi(comm,"concentration",Nx,Ny,Nz,domain_o),nonlinear(phi,"chempotential"),
    domain(domain_o),ft_domain(domain_o.fourier_transpose(Nx,Ny,Nz)),
    ft_mod_domain(domain_o.fourier_mod_transpose(Nx,Ny,Nz)),
    ft_phi(phi.make_fourier_transpose(ft_domain)),
    ft_nonlinear(nonlinear.make_fourier_transpose(ft_domain)),
    ft_noise(ft_phi,"ft_noise"),
    modulus(ft_phi.make_fourier_mod(ft_mod_domain))
{

  forward_phi = fftw_mpi_plan_dft_r2c_3d(Nz,Ny,Nx,phi.data(),
					 reinterpret_cast<fftw_complex*>
					 (ft_phi.data()),
					 comm, FFTW_MPI_TRANSPOSED_OUT);
  
  backward_phi = fftw_mpi_plan_dft_c2r_3d(Nz,Ny,Nx,reinterpret_cast<fftw_complex*>
					  (ft_phi.data()), phi.data(),
					  comm,FFTW_MPI_TRANSPOSED_IN);
  
  forward_nonlinear = fftw_mpi_plan_dft_r2c_3d(Nz,Ny,Nx,
					       nonlinear.data(),
					       reinterpret_cast<fftw_complex*>
					       (ft_nonlinear.data()),
					       comm, FFTW_MPI_TRANSPOSED_OUT);
  
  backward_nonlinear = fftw_mpi_plan_dft_c2r_3d(Nz,Ny,Nx,
						reinterpret_cast<fftw_complex*>
						(ft_nonlinear.data()),
						nonlinear.data(),comm,
						FFTW_MPI_TRANSPOSED_IN);
  

}




gridData::~gridData() {

  
  fftw_destroy_plan(forward_phi);
  fftw_destroy_plan(backward_phi);
  fftw_destroy_plan(forward_nonlinear);
  fftw_destroy_plan(backward_nonlinear);

  

}


void gridData::reverseFlat(int gridindex, int &i, int &j, int &k)
{
  k = gridindex % phi.Nx();
  j = (gridindex / phi.Nx()) % phi.Ny();
  i = (gridindex / phi.Nx())/phi.Ny();

}


void gridData::reverseFlatFourier(int gridindex, int &i, int &j, int &k)
{
  k = gridindex % ft_phi.Nx();
  j = (gridindex / ft_phi.Nx()) % ft_phi.Ny();
  i = (gridindex / ft_phi.Nx())/ft_phi.Ny();

}
