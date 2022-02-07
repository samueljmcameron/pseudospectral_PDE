#include <iostream>
#include <cmath>
#include <complex>


#include "fftw_mpi_3darray.hpp"
#include "integrator.hpp"
#include "conjplane.hpp"
#include "randompll.hpp"
#include "griddata.hpp"
#include "iovtk.hpp"
#include "solutionparams.hpp"
#include "timestep.hpp"

int main()
{

  GridData realspace(64,100.0);
  GridData fourier = realspace.fft_grid();

  fourier.transpose_yz();
  

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
  
  fftw_MPI_3Darray<double> phi(comm,"concentration",realspace);
  fftw_MPI_3Darray<double> nonlinear(comm,"chempotential",realspace);

  int baseseed = 129480;
  RandomPll rpll(comm,id,baseseed,mpi_size);
  
  double dt = 1e-4;


  std::vector<std::string> splitvec;

  splitvec.push_back("mobility");
  splitvec.push_back("1.0");
  splitvec.push_back("gamma");
  splitvec.push_back("100.0");
  splitvec.push_back("temp");
  splitvec.push_back("1.0");
  splitvec.push_back("chi");
  splitvec.push_back("2.5");
  splitvec.push_back("volFH");
  splitvec.push_back("0.01");

  SolutionParams solparams(splitvec);

  if (id == 0) {
    solparams.printall();
  }
  

  
  Integrator integrator(comm,fourier,rpll.get_processor_seed(),solparams,dt);



  double volfrac = 0.3;
  double variance = 0.0;





  fftw_plan forward_phi, backward_phi;
  fftw_plan forward_nonlinear, backward_nonlinear;

  
  forward_phi = fftw_mpi_plan_dft_r2c_3d(realspace.get_Nz(),realspace.get_Ny(),
					 realspace.get_Nx(),
					 phi.data(),
					 reinterpret_cast<fftw_complex*>
					 (integrator.ft_phi.data()),
					 comm, FFTW_MPI_TRANSPOSED_OUT);
  
  backward_phi = fftw_mpi_plan_dft_c2r_3d(realspace.get_Nz(),realspace.get_Ny(),
					  realspace.get_Nx(),
					  reinterpret_cast<fftw_complex*>
					  (integrator.ft_phi.data()),
					  phi.data(),comm,FFTW_MPI_TRANSPOSED_IN);

  forward_nonlinear = fftw_mpi_plan_dft_r2c_3d(realspace.get_Nz(),realspace.get_Ny(),
					       realspace.get_Nx(),
					       nonlinear.data(),
					       reinterpret_cast<fftw_complex*>
					       (integrator.ft_nonlinear.data()),
					       comm, FFTW_MPI_TRANSPOSED_OUT);
  
  backward_nonlinear = fftw_mpi_plan_dft_c2r_3d(realspace.get_Nz(),realspace.get_Ny(),
						realspace.get_Nx(),
						reinterpret_cast<fftw_complex*>
						(integrator.ft_nonlinear.data()),
						nonlinear.data(),comm,
						FFTW_MPI_TRANSPOSED_IN);
  

  integrator.initialize(phi,volfrac,variance);

  int numsteps = 20;
  int startstep = 0;

  double t = 0;

  std::string fname = std::string("data/test") + std::string("_p") + std::to_string(id) ;
  std::string fname_p = fname + std::string("_") +  std::to_string(startstep) +  std::string(".vti");

  bool restart_flag = false;
  
  
  if (restart_flag) {
    ioVTK::readVTKImageData({&phi},fname_p,realspace);
  } else {
    ioVTK::writeVTKImageData(fname_p,{&phi},realspace);
    
  }

  TimeStep timestep(comm,mpi_size,id,fourier.get_Nz(),fourier.get_Ny());
  
  for (int it = 1+startstep; it <= numsteps+startstep; it ++) {
    t += integrator.get_dt();
    std::cout << "id " << id << " is on t = " << t << std::endl;
    integrator.nonlinear(nonlinear,phi); // compute nl(t) given phi(t)
    
    fftw_execute(forward_phi);
    fftw_execute(forward_nonlinear);

    timestep.update(t,integrator);    

    fftw_execute(backward_phi); // get phi(t+dt)

    if (it % 1 == 0) {

      fname_p = fname + std::string("_") +  std::to_string(it) +  std::string(".vti");
      
      ioVTK::writeVTKImageData(fname_p,{&phi},realspace);
    }
  }

  fftw_destroy_plan(forward_phi);
  fftw_destroy_plan(backward_phi);
  fftw_destroy_plan(forward_nonlinear);
  fftw_destroy_plan(backward_nonlinear);
    
  fftw_mpi_cleanup();
  ierr = MPI_Finalize();
  
  return 0;
}

void check_results(Integrator &integrator, const int id,const int divider,const int nx)
{
  using namespace std::complex_literals;
  
  // check the resulting blocks to see if everything is done correctly
  int cNy = integrator.ft_phi.axis_size(0);
  int Nz = integrator.ft_phi.axis_size(1);
  int Ny = integrator.ft_phi.get_globalNz();
  int complex_local = integrator.ft_phi.get_local0start();

  int loopstart;
  if (id == 0 ) {
    loopstart = 1;
  } else {
    loopstart = 0;
  }

  


  std::complex diff = 0.0;


  if (id != divider) {

    if (id == 0) {
      int ny = 0;
      for (int j = 1; j < Nz/2; j++) {
	diff += std::abs(1.0*(ny+complex_local)+1i*(1.0*j) - integrator.ft_phi(ny,j,nx));
      }
      for (int j = Nz/2+1; j < Nz; j++) {
	diff += std::abs(1.0*(ny+complex_local)-1i*(1.0*(Nz-j)) - integrator.ft_phi(ny,j,nx));
      }
    }
	
    for (int i = loopstart; i < cNy; i++) {
      for (int j = 1; j < Nz/2; j++) {
      	diff += std::abs(1.0*(i+complex_local)+1i*(1.0*j) - integrator.ft_phi(i,j,nx));
      }

      for (int j = Nz/2+1; j < Nz; j++) {
	diff += std::abs(1.0*(Ny-i-complex_local)-1i*(1.0*(Nz-j)) - integrator.ft_phi(i,j,nx));
	
      }
      
    }
    if (id < divider) {
      
      for (int i = loopstart; i < cNy; i++) {
	diff += std::abs(1.0*(i+complex_local)+1i*(1.0*0) - integrator.ft_phi(i,0,nx));
	diff += std::abs(1.0*(Ny-i-complex_local)-1i*(1.0*Nz/2) - integrator.ft_phi(i,Nz/2,nx));
      }
    } else {
      for (int i = loopstart; i < cNy; i++) {
	diff += std::abs(1.0*(i+complex_local)+1i*(1.0*Nz/2) - integrator.ft_phi(i,Nz/2,nx));
	diff += std::abs(1.0*(Ny-i-complex_local)-1i*(1.0*0) - integrator.ft_phi(i,0,nx));
      }      
    }

    std::cout << "id = " << id << " diff = " << diff << std::endl;
  } else {
    
    for (int i = loopstart; i < cNy; i++) {
      for (int j = 1; j < Nz/2; j++) {
	diff += std::abs(1.0*(i+complex_local)+1i*(1.0*j) - integrator.ft_phi(i,j,nx));
      }
      for (int j = Nz/2+1; j < Nz; j++) {
	diff += std::abs(1.0*(Ny-i-complex_local)-1i*(1.0*(Nz-j)) - integrator.ft_phi(i,j,nx));
      }
    }

    for (int i = loopstart; i < Ny/2-complex_local; i++) {
      diff += std::abs(1.0*(i+complex_local)+1i*(1.0*0) - integrator.ft_phi(i,0,nx));
      diff += std::abs(1.0*(Ny-i-complex_local)-1i*(1.0*Nz/2) - integrator.ft_phi(i,Nz/2,nx));
    }
    for (int i = Ny/2-complex_local+1; i < cNy; i++) {
      diff += std::abs(1.0*(i+complex_local)+1i*(1.0*Nz/2) - integrator.ft_phi(i,Nz/2,nx));
      diff += std::abs(1.0*(Ny-i-complex_local)-1i*(1.0*0) - integrator.ft_phi(i,0,nx));
    }

    std::cout << "id = " << id << " diff = " << diff << std::endl;
  }
  return;
}
