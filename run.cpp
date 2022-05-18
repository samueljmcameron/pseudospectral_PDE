
#include "fftw_mpi_3darray.hpp"
#include "integrator.hpp"
//#include "linker.hpp"
#include "conjplane.hpp"
#include "randompll.hpp"
#include "griddata.hpp"
#include "iovtk.hpp"
#include "timestep.hpp"

#include "run.hpp"

void run(GlobalParams gp, SolutionParams solparams) {

  
  fftw_MPI_3Darray<double> phi(gp.comm,"concentration",gp.realspace);
  fftw_MPI_3Darray<double> nonlinear(gp.comm,"chempotential",gp.realspace);


  RandomPll rpll(gp.comm,gp.id,gp.seed,gp.mpi_size);
  
  Integrator integrator(gp.comm,gp.fourier,rpll.get_processor_seed(),solparams,gp.dt);

  fftw_plan forward_phi, backward_phi;
  fftw_plan forward_nonlinear, backward_nonlinear;

  
  forward_phi = fftw_mpi_plan_dft_r2c_3d(gp.realspace.get_Nz(),gp.realspace.get_Ny(),
					 gp.realspace.get_Nx(),
					 phi.data(),
					 reinterpret_cast<fftw_complex*>
					 (integrator.ft_phi.data()),
					 gp.comm, FFTW_MPI_TRANSPOSED_OUT);
  
  backward_phi = fftw_mpi_plan_dft_c2r_3d(gp.realspace.get_Nz(),gp.realspace.get_Ny(),
					  gp.realspace.get_Nx(),
					  reinterpret_cast<fftw_complex*>
					  (integrator.ft_phi.data()),
					  phi.data(),gp.comm,FFTW_MPI_TRANSPOSED_IN);

  forward_nonlinear = fftw_mpi_plan_dft_r2c_3d(gp.realspace.get_Nz(),gp.realspace.get_Ny(),
					       gp.realspace.get_Nx(),
					       nonlinear.data(),
					       reinterpret_cast<fftw_complex*>
					       (integrator.ft_nonlinear.data()),
					       gp.comm, FFTW_MPI_TRANSPOSED_OUT);
  
  backward_nonlinear = fftw_mpi_plan_dft_c2r_3d(gp.realspace.get_Nz(),gp.realspace.get_Ny(),
						gp.realspace.get_Nx(),
						reinterpret_cast<fftw_complex*>
						(integrator.ft_nonlinear.data()),
						nonlinear.data(),gp.comm,
						FFTW_MPI_TRANSPOSED_IN);
  

  integrator.initialize(phi,gp.volFrac,gp.variance);

  std::vector<double> X_i = {0.0,0.0,5.0};
  double t = gp.starttime;

  std::string prefix = gp.dump_file + std::string("_p") + std::to_string(gp.id) ;
  std::string fname_p = prefix + std::string("_") +  std::to_string(gp.startstep) +  std::string(".vti");

  std::string collection_name = prefix + std::string(".pvd");
  
  std::string complexprefix = gp.dump_file + std::string("_complex")
    + std::string("_p") + std::to_string(gp.id) ;

  std::string complexfname_p = complexprefix + std::string("_") +  std::to_string(gp.startstep) +  std::string(".vti");

  std::string complexcollection_name = complexprefix + std::string(".pvd");
  
  fftw_MPI_3Darray<double> modulus(gp.comm,integrator.ft_phi.get_name()+std::string("_mod"),
				   gp.fourier.get_positiveNx_grid());

  int running_average_count = 0;
  if (gp.restart_flag) {
    
    ioVTK::restartVTKcollection(collection_name);
    ioVTK::restartVTKcollection(complexcollection_name);
    
    ioVTK::readVTKImageData({&phi},fname_p,gp.realspace);
    
    for (int i = 0; i < integrator.ft_phi.axis_size(0); i++) {
      for (int j = 0; j < integrator.ft_phi.axis_size(1); j++) {
	for (int k = 0; k < integrator.ft_phi.axis_size(2); k++) {
	  modulus(i,j,k) = 0.0;
	}
      }
    }

    std::cout << "restarting!" << std::endl;

  } else {

    ioVTK::writeVTKcollectionHeader(collection_name);
    ioVTK::writeVTKcollectionHeader(complexcollection_name);

    fftw_execute(forward_phi);
    integrator.ft_phi.mod(modulus);

    double norm = 1.0/(integrator.ft_phi.grid.get_Nx()*integrator.ft_phi.grid.get_Ny()
		       *integrator.ft_phi.grid.get_Nz());
    
    for (int i = 0; i < integrator.ft_phi.axis_size(0); i++) {
      for (int j = 0; j < integrator.ft_phi.axis_size(1); j++) {
	for (int k = 0; k < integrator.ft_phi.axis_size(2); k++) {
	  integrator.ft_phi(i,j,k) = integrator.ft_phi(i,j,k)*norm;
	  modulus(i,j,k) = modulus(i,j,k)*norm;
	}
      }
    }
    if (gp.id == 0) {
      modulus(0,0,0) = 0.0;
      std::cout << "modulus at zero is = " << modulus(0,0,0) << std::endl;
    }

    fftw_execute(backward_phi);
    
    ioVTK::writeVTKImageData(fname_p,{&phi},gp.realspace);
    ioVTK::writeVTKImageData(complexfname_p,{&modulus},modulus.grid);
    
    ioVTK::writeVTKcollectionMiddle(collection_name,fname_p,t);
    ioVTK::writeVTKcollectionMiddle(complexcollection_name,complexfname_p,t);
  }


  std::vector<double> free_energy_derivative;
  
  TimeStep timestep(gp.comm,gp.mpi_size,gp.id,integrator.ft_phi.axis_size(0),
		    integrator.ft_phi.axis_size(1));
  
  for (int it = 1+gp.startstep; it <= gp.steps+gp.startstep; it ++) {
    t += integrator.get_dt();

    integrator.nonlinear(nonlinear,phi,X_i); // compute nl(t) given phi(t)
    
    fftw_execute(forward_phi);
    fftw_execute(forward_nonlinear);

    
    timestep.update(t,integrator);

    integrator.ft_phi.running_mod(modulus);
    running_average_count += 1;

    if (it % gp.dump_every == 0) {
      std::cout << "id " << gp.id << " saving at t = " << t << std::endl;
      fname_p = prefix + std::string("_") +  std::to_string(it) +  std::string(".vti");
      complexfname_p = complexprefix + std::string("_") +  std::to_string(it) +  std::string(".vti");
      modulus /= running_average_count;
      if (gp.id == 0) {
	modulus(0,0,0) = 0.0;
	std::cout << "modulus at zero is = " << modulus(0,0,0) << std::endl;
      }
      
      running_average_count = 0;

      ioVTK::writeVTKImageData(complexfname_p,{&modulus},modulus.grid);
      ioVTK::writeVTKcollectionMiddle(complexcollection_name,complexfname_p,t);
      fftw_execute(backward_phi); // get phi(t+dt)

      integrator.initialize(modulus,0,0);

      ioVTK::writeVTKImageData(fname_p,{&phi},phi.grid);
      ioVTK::writeVTKcollectionMiddle(collection_name,fname_p,t);
    } else {
      fftw_execute(backward_phi); // get phi(t+dt)
    }

    //    link.wrap_X_i(X_i,phi.grid.get_Lx(),phi.grid.get_Ly(),phi.grid.get_Lz());
    //    free_energy_derivative
    //      = link.free_energy_derivative(X_i,integrator.linker_phi,
    //				    integrator.linker_derivative,phi);

    
  }

  ioVTK::writeVTKcollectionFooter(collection_name);
  ioVTK::writeVTKcollectionFooter(complexcollection_name);
  
  fftw_destroy_plan(forward_phi);
  fftw_destroy_plan(backward_phi);
  fftw_destroy_plan(forward_nonlinear);
  fftw_destroy_plan(backward_nonlinear);

  return;
}
