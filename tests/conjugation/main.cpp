
#include <fstream>
#include <iostream>

#include <fftw3-mpi.h>

#include "conjugate_trig.hpp"
#include "fftw_mpi_3darray.hpp"
#include "domain.hpp"

int main()
{

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;
  
  int ierr = MPI_Init(NULL,NULL);


  if (ierr != MPI_SUCCESS) {
    std::cout << "Error starting MPI program. Terminating." << std::endl;
    MPI_Abort(MPI_COMM_WORLD,ierr);
  }
  int mpi_size, id;

  ierr = MPI_Comm_size(comm,&mpi_size);
  ierr = MPI_Comm_rank(comm,&id);

  fftw_mpi_init();

  const int globalN = 10;
  const double L = 100.0;

  psPDE::Domain domain(L,L,L,-L/2,-L/2,-L/2,id,mpi_size);

  psPDE::Domain ft_domain(domain.fourier_transpose(globalN,globalN,globalN));
  
  psPDE::fftw_MPI_3Darray<std::complex<double>>
    ft_noise(comm,"ft_noise",globalN,globalN,globalN,ft_domain);

  
  psPDE::ConjugateTrig nUp(comm,mpi_size,id,&ft_noise);

  // conjugate things appropriately
  nUp.update();

  
  std::ofstream realfile,imagfile;

  realfile.open("real_" + std::to_string(id) + std::string(".txt"));
  imagfile.open("imag_" + std::to_string(id) + std::string(".txt"));
  
  for (int i = 0; i < ft_noise.Nz(); i ++) 
    for (int j = 0; j < ft_noise.Ny(); j++) 
      for (int k = 0; k < ft_noise.Nx(); k++) {
	realfile << ft_noise(i,j,k).real() << " ";
	imagfile << ft_noise(i,j,k).imag() << " ";
      }


  fftw_mpi_cleanup();

  ierr = MPI_Finalize();

  
  return 0;
}
