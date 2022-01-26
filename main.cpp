#include <iostream>
#include <cmath>
#include <complex>


#include "fftw_mpi_3darray.hpp"
#include "timestep.hpp"
#include "conjplane.hpp"
#include "randompll.hpp"
#include "griddata.hpp"
#include "writevtk.hpp"

int main()
{


  const int fullNx = 64;
  const int fullNy = 64;
  const int fullNz = 64;

  double Lx = 100.0;
  double Ly = 100.0;
  double Lz = 100.0;


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
  
  fftw_MPI_3Darray<double> phi(comm,"concentration",fullNz,fullNy,fullNx);
  fftw_MPI_3Darray<double> nonlinear(comm,"chempotential",fullNz,fullNy,fullNx);

  int baseseed = 129480;
  RandomPll rpll(comm,id,baseseed,mpi_size);
  
  double mobility = 1.0;
  double gamma = 1.0;
  double temp = 1.0;
  double chi = 2.5;
  double volFH = 0.01;
  double dt = 1e-4;
  
  TimeStep timestep(comm,fullNy,fullNz,fullNx,Ly,Lz,Lx,rpll.get_processor_seed(),mobility,
		    gamma,temp,chi,volFH,dt);



  double volfrac = 0.3;
  double variance = 0.1;



  
  GridData gridstart;
  gridstart.dx = Lx/fullNx;
  gridstart.dy = Ly/fullNy;
  gridstart.dz = Lz/fullNz;
  gridstart.Ox = -1*(fullNx-1)/2*gridstart.dx;
  gridstart.Oy = -1*(fullNy-1)/2*gridstart.dy;
  gridstart.Oz = -1*(fullNz-1)/2*gridstart.dz;



  fftw_plan forward_phi, backward_phi;
  fftw_plan forward_nonlinear, backward_nonlinear;

  
  forward_phi = fftw_mpi_plan_dft_r2c_3d(fullNz,fullNy,fullNx,
					 phi.data(),
					 reinterpret_cast<fftw_complex*>
					 (timestep.ft_phi.data()),
					 comm, FFTW_MPI_TRANSPOSED_OUT);
  
  backward_phi = fftw_mpi_plan_dft_c2r_3d(fullNz,fullNy,fullNx,
					  reinterpret_cast<fftw_complex*>
					  (timestep.ft_phi.data()),
					  phi.data(),comm,FFTW_MPI_TRANSPOSED_IN);

  forward_nonlinear = fftw_mpi_plan_dft_r2c_3d(fullNz,fullNy,fullNx,
					       nonlinear.data(),
					       reinterpret_cast<fftw_complex*>
					       (timestep.ft_nonlinear.data()),
					       comm, FFTW_MPI_TRANSPOSED_OUT);
  
  backward_nonlinear = fftw_mpi_plan_dft_c2r_3d(fullNz,fullNy,fullNx,
						reinterpret_cast<fftw_complex*>
						(timestep.ft_nonlinear.data()),
						nonlinear.data(),comm,
						FFTW_MPI_TRANSPOSED_IN);
  

  timestep.initialize(phi,volfrac,variance);
  
  std::vector<int> all_cNys(mpi_size);

  for (int i = 0; i < mpi_size; i++) {

    if (i == id) {
      all_cNys[i] = timestep.ft_phi.axis_size(0);

    }
    
    MPI_Bcast(&all_cNys[i],1,MPI_INT,i,comm);

    

  }

  int cNy = timestep.ft_phi.axis_size(0);
  int cNz = timestep.ft_phi.axis_size(1);
  int cNx = timestep.ft_phi.axis_size(2);
  int local0start = timestep.ft_phi.get_local0start();

  int numsteps = 1000;

  double t = 0;
  {
    std::vector<fftw_MPI_3Darray<double>*> scalars;
  
    scalars.push_back(&phi);
    
    std::string fname = std::string("data/test") + std::string("_p") + std::to_string(id) ;
    std::string fname_p = fname + std::string("_") +  std::to_string(0) +  std::string(".vti");

    writeVTK::writeVTKImageData(fname_p,scalars,gridstart);
  }


  
  for (int it = 1; it <= numsteps; it ++) {
    t = it*timestep.get_dt();
    std::cout << "id " << id << " is on t = " << t << std::endl;
    timestep.nonlinear(nonlinear,phi); // compute nl(t) given phi(t)
    
    fftw_execute(forward_phi);
    fftw_execute(forward_nonlinear);
    


    // update the bulk of the system which has no constraints
    
    for (int ny = 0; ny < cNy; ny ++) {
      
      for (int nz = 0; nz < cNz; nz ++) {
	
	for (int nx = 1; nx < cNx-1; nx++) {
	  timestep.update(ny,nz,nx);
	}
      }
    }
    
    

    // update the two planes at nx = 0 and nx = Nx-1
    // which are constrained to ensure that phi remains real
    
    ConjPlane cpcl(mpi_size,id,comm,all_cNys,fullNz);
    
    for (int nx = 0; nx < cNx; nx += cNx-1) {
      
      if (mpi_size == 1) {
	
	cpcl.single(timestep,nx);
	cpcl.line_single(timestep,nx);
	
      } else {
	
	if (cpcl.is_equal() && id == 0) {
	  cpcl.first(timestep,nx);
	  cpcl.line_first(timestep,nx);
	} else if (id < cpcl.get_divider()) {
	  cpcl.lefthalf(timestep,nx);
	  cpcl.line_lefthalf(timestep,nx);
	} else if (id == cpcl.get_divider()) {
	  if (mpi_size % 2 == 0) {
	    cpcl.middle_even(timestep,nx);
	    cpcl.line_middle_even(timestep,nx);
	  } else {
	    cpcl.middle_odd(timestep,nx);
	    cpcl.line_middle_odd(timestep,nx);
	  }
	} else if (!cpcl.is_equal() && id == mpi_size-1) {
	  cpcl.last(timestep,nx);
	  cpcl.line_last(timestep,nx);
	} else {
	  cpcl.righthalf(timestep,nx);
	  cpcl.line_righthalf(timestep,nx);
	}
      }
    }
    
    
    // update the eight points where ft_phi must be real
    if (id == 0) {
      int nx = cNx-1;
      timestep.update_real(0,cNz/2,0);
      timestep.update_real(0,cNz/2,nx);
      timestep.update_real(0,0,nx);
      timestep.ft_phi(0,0,0) = timestep.ft_phi(0,0,0)/(1.0*fullNx*fullNy*fullNz);
    }
    if (timestep.ft_phi.get_local0start() <= fullNy/2
	&& timestep.ft_phi.get_local0start() + cNy -1 >= fullNy/2) {
      int ny = fullNy/2 - local0start;
      int nx = cNx-1;
      timestep.update_real(ny,0,0);
      timestep.update_real(ny,cNz/2,0);
      timestep.update_real(ny,cNz/2,nx);
      timestep.update_real(ny,0,nx);
    }

    fftw_execute(backward_phi); // get phi(t+dt)

    std::vector<fftw_MPI_3Darray<double>*> scalars;
  
    scalars.push_back(&phi);
    
    std::string fname = std::string("data/test") + std::string("_p") + std::to_string(id) ;
    std::string fname_p = fname + std::string("_") +  std::to_string(it) +  std::string(".vti");

    writeVTK::writeVTKImageData(fname_p,scalars,gridstart);
  }

  fftw_destroy_plan(forward_phi);
  fftw_destroy_plan(backward_phi);
  fftw_destroy_plan(forward_nonlinear);
  fftw_destroy_plan(backward_nonlinear);
    
  fftw_mpi_cleanup();
  ierr = MPI_Finalize();
  
  return 0;
}

void check_results(TimeStep &timestep, const int id,const int divider,const int nx)
{
  using namespace std::complex_literals;
  
  // check the resulting blocks to see if everything is done correctly
  int cNy = timestep.ft_phi.axis_size(0);
  int Nz = timestep.ft_phi.axis_size(1);
  int Ny = timestep.ft_phi.get_globalNz();
  int complex_local = timestep.ft_phi.get_local0start();

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
	diff += std::abs(1.0*(ny+complex_local)+1i*(1.0*j) - timestep.ft_phi(ny,j,nx));
      }
      for (int j = Nz/2+1; j < Nz; j++) {
	diff += std::abs(1.0*(ny+complex_local)-1i*(1.0*(Nz-j)) - timestep.ft_phi(ny,j,nx));
      }
    }
	
    for (int i = loopstart; i < cNy; i++) {
      for (int j = 1; j < Nz/2; j++) {
      	diff += std::abs(1.0*(i+complex_local)+1i*(1.0*j) - timestep.ft_phi(i,j,nx));
      }

      for (int j = Nz/2+1; j < Nz; j++) {
	diff += std::abs(1.0*(Ny-i-complex_local)-1i*(1.0*(Nz-j)) - timestep.ft_phi(i,j,nx));
	
      }
      
    }
    if (id < divider) {
      
      for (int i = loopstart; i < cNy; i++) {
	diff += std::abs(1.0*(i+complex_local)+1i*(1.0*0) - timestep.ft_phi(i,0,nx));
	diff += std::abs(1.0*(Ny-i-complex_local)-1i*(1.0*Nz/2) - timestep.ft_phi(i,Nz/2,nx));
      }
    } else {
      for (int i = loopstart; i < cNy; i++) {
	diff += std::abs(1.0*(i+complex_local)+1i*(1.0*Nz/2) - timestep.ft_phi(i,Nz/2,nx));
	diff += std::abs(1.0*(Ny-i-complex_local)-1i*(1.0*0) - timestep.ft_phi(i,0,nx));
      }      
    }

    std::cout << "id = " << id << " diff = " << diff << std::endl;
  } else {
    
    for (int i = loopstart; i < cNy; i++) {
      for (int j = 1; j < Nz/2; j++) {
	diff += std::abs(1.0*(i+complex_local)+1i*(1.0*j) - timestep.ft_phi(i,j,nx));
      }
      for (int j = Nz/2+1; j < Nz; j++) {
	diff += std::abs(1.0*(Ny-i-complex_local)-1i*(1.0*(Nz-j)) - timestep.ft_phi(i,j,nx));
      }
    }

    for (int i = loopstart; i < Ny/2-complex_local; i++) {
      diff += std::abs(1.0*(i+complex_local)+1i*(1.0*0) - timestep.ft_phi(i,0,nx));
      diff += std::abs(1.0*(Ny-i-complex_local)-1i*(1.0*Nz/2) - timestep.ft_phi(i,Nz/2,nx));
    }
    for (int i = Ny/2-complex_local+1; i < cNy; i++) {
      diff += std::abs(1.0*(i+complex_local)+1i*(1.0*Nz/2) - timestep.ft_phi(i,Nz/2,nx));
      diff += std::abs(1.0*(Ny-i-complex_local)-1i*(1.0*0) - timestep.ft_phi(i,0,nx));
    }

    std::cout << "id = " << id << " diff = " << diff << std::endl;
  }
  return;
}
