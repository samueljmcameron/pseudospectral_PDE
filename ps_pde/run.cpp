#include <fstream>
#include <chrono>
#include <random>

#include "fftw_mpi_3darray.hpp"
#include "integrator.hpp"
#include "conjplane.hpp"
#include "randompll.hpp"
#include "griddata.hpp"
#include "iovtk.hpp"
#include "timestep.hpp"
#include "input.hpp"

#include "run.hpp"

void X_i_of_t(std::vector<double> & X_i,
	      const std::vector<double> & dFdX_i,double dt,
	      double radius, double viscosity,double temp,
	      double randx, double randy, double randz,
	      int noiseon)
{
  double fric = 6*M_PI*viscosity*radius;
  X_i[0] += -dt*dFdX_i[0]/fric + noiseon*sqrt(24*temp/fric*dt)*randx;
  X_i[1] += -dt*dFdX_i[1]/fric + noiseon*sqrt(24*temp/fric*dt)*randy;
  X_i[2] += -dt*dFdX_i[2]/fric + noiseon*sqrt(24*temp/fric*dt)*randz;
  return;
}


void run(psPDE::GlobalParams gp, psPDE::SolutionParams solparams,
	 std::vector<std::vector<double>> &X_is,
	 std::vector<double> & radii, std::vector<double> & viscosities) {

  
  psPDE::fftw_MPI_3Darray<double> phi(gp.comm,"concentration",gp.realspace);
  psPDE::fftw_MPI_3Darray<double> nonlinear(gp.comm,"chempotential",gp.realspace);


  psPDE::RandomPll rpll(gp.comm,gp.id,gp.seed,gp.mpi_size);
  
  psPDE::Integrator integrator(gp.comm,gp.fourier,rpll.get_processor_seed(),solparams,gp.dt);

  fftw_plan forward_phi, backward_phi;
  fftw_plan forward_nonlinear, backward_nonlinear;

  // can use the global seed, as it is just used to generate seeds but not actually
  //  used as a seed itself.
  std::mt19937 dropgen(gp.seed);
  std::uniform_real_distribution<double> real_dist(-0.5,0.5);
  
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
  

  
  std::vector<std::vector<double>> free_energy_derivs;

  
  double t = gp.starttime;

  std::string prefix = gp.dump_file + std::string("_p") + std::to_string(gp.id) ;
  std::string fname_p = prefix + std::string("_") +  std::to_string(gp.startstep) +  std::string(".vti");

  std::string collection_name = prefix + std::string(".pvd");
  
  std::string complexprefix = gp.dump_file + std::string("_complex")
    + std::string("_p") + std::to_string(gp.id) ;

  std::string complexfname_p = complexprefix + std::string("_") +  std::to_string(gp.startstep) +  std::string(".vti");

  std::string complexcollection_name = complexprefix + std::string(".pvd");
  
  psPDE::fftw_MPI_3Darray<double> modulus(gp.comm,integrator.ft_phi.get_name()+std::string("_mod"),
					  gp.fourier.get_positiveNx_grid());


  
  int running_average_count = 0;

  std::ofstream myfile;

  
  if (gp.restart_flag) {


    if (X_is.size() > 0) 
      throw std::runtime_error("Cannot specify additional nucleation sites when "
			       + std::string("restarting a file. Use read instead."));


    psPDE::input::read_in_nuclei_properties(radii,viscosities,gp.nucs_to_keep,
					    true,gp.thermo_file,
					    gp.comm,gp.id);
    
    psPDE::input::put_in_vectors(X_is,gp.nucs_to_keep,gp.thermo_file,gp.comm,
				 gp.id,gp.starttime);

    
    psPDE::ioVTK::restartVTKcollection(collection_name,gp.comm);
    psPDE::ioVTK::restartVTKcollection(complexcollection_name,gp.comm);
    
    psPDE::ioVTK::readVTKImageData({&phi},fname_p);
    
    for (int i = 0; i < integrator.ft_phi.axis_size(0); i++) {
      for (int j = 0; j < integrator.ft_phi.axis_size(1); j++) {
	for (int k = 0; k < integrator.ft_phi.axis_size(2); k++) {
	  modulus(i,j,k) = 0.0;
	}
      }
    }

    if (gp.id == 0) {
      myfile.open(gp.thermo_file,std::ios::app);
      
    }

  } else {

    psPDE::ioVTK::writeVTKcollectionHeader(collection_name);
    psPDE::ioVTK::writeVTKcollectionHeader(complexcollection_name);


    if (gp.read_flag) {
      psPDE::ioVTK::readVTKImageData({&phi},gp.read_dump_file);
      
      psPDE::input::read_in_nuclei_properties(radii,viscosities,gp.nucs_to_keep,
					      gp.all_nucs_flag,gp.read_thermo_file,
					      gp.comm,gp.id);
      
      psPDE::input::put_in_vectors(X_is,gp.nucs_to_keep,gp.read_thermo_file,gp.comm,
				   gp.id,gp.starttime);

	
    } else { 
      integrator.initialize(phi,gp.volFrac,gp.variance);
    }

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
    }

    fftw_execute(backward_phi);
    
    psPDE::ioVTK::writeVTKImageData(fname_p,{&phi},gp.realspace);
    psPDE::ioVTK::writeVTKImageData(complexfname_p,{&modulus},modulus.grid);
    
    psPDE::ioVTK::writeVTKcollectionMiddle(collection_name,fname_p,t);
    psPDE::ioVTK::writeVTKcollectionMiddle(complexcollection_name,complexfname_p,t);


    
    if (gp.id == 0) {
      myfile.open(gp.thermo_file);

      myfile << "# nucnum \t radius \t viscosity " << std::endl;
      for (int index =  0; index < viscosities.size(); index++ ) {
	myfile << "# " << index << " \t " << radii.at(index) << " \t "
	       << viscosities.at(index) << std::endl;
      }
      
      
      myfile << "# t ";
      for (unsigned index = 0; index < X_is.size() ; index ++) {
	myfile << "\t (X_" << index << ")_x " << "(X_" << index << ")_y"
	       << "(X_" << index << ")_z";
      }
      myfile << "\t F(X) ";
      for (unsigned index = 0; index < X_is.size() ; index ++) {
	myfile << "\t (dF/dX_" << index << ")_x " << "(dF/dX_" << index << ")_y"
	       << "(dF/dX_" << index << ")_z";
      }
      myfile << std::endl;

  }

    
  }



  using Clock = std::chrono::steady_clock;
  using namespace std::literals;
  auto constexpr chronoitvl = 1.0s/60.0;
  using duration = std::chrono::duration<double>;
  using time_point = std::chrono::time_point<Clock,duration>;

  duration tt_integral = 0s;

  duration tt_freeenergy = 0s;

  
  psPDE::TimeStep timestep(gp.comm,gp.mpi_size,gp.id,integrator.ft_phi.axis_size(0),
		    integrator.ft_phi.axis_size(1));

  double free_energy;
  for (int it = 1+gp.startstep; it <= gp.steps+gp.startstep; it ++) {


    auto current_time = Clock::now();

    free_energy_derivs = integrator.nonlinear(nonlinear,phi,X_is,free_energy); // compute nl(t) given phi(t)

    tt_freeenergy += Clock::now() - current_time;

    if (gp.id == 0 && it % gp.thermo_every == 0) {
      myfile << t;

      for (unsigned index = 0; index < X_is.size() ; index ++) {

	myfile << "\t " << X_is.at(index).at(0) << "\t "
	       << X_is.at(index).at(1) << "\t "
	       << X_is.at(index).at(2);
      }

      myfile << "\t " << free_energy;
      for (unsigned index = 0; index < free_energy_derivs.size() ; index ++) {

	myfile << "\t " << free_energy_derivs.at(index).at(0) << "\t "
	       << free_energy_derivs.at(index).at(1) << "\t "
	       << free_energy_derivs.at(index).at(2);
      }
      

      myfile << std::endl;
      
    }

    t += integrator.get_dt();

    for (unsigned index = 0; index < X_is.size(); index ++ ) {
      double randx = real_dist(dropgen);
      double randy = real_dist(dropgen);
      double randz = real_dist(dropgen);
      X_i_of_t(X_is[index],free_energy_derivs[index],integrator.get_dt(),radii[index],
	       viscosities[index],solparams.temp,randx,randy,randz,gp.X_i_noise);
    }
    
    
    fftw_execute(forward_phi);
    fftw_execute(forward_nonlinear);

    current_time = Clock::now();
    timestep.update(t,integrator);
    tt_integral = Clock::now() - current_time;

    
    integrator.ft_phi.running_mod(modulus);
    running_average_count += 1;


    if (std::isnan(phi(0,0,0))) {

      if (gp.id == 0) {
	myfile.close();
      }
      
      psPDE::ioVTK::writeVTKcollectionFooter(collection_name);
      psPDE::ioVTK::writeVTKcollectionFooter(complexcollection_name);
      throw std::runtime_error("Solution concentration diverged at t = " +
			       std::to_string(t));
      
    }


    if (it % gp.dump_every == 0) {
      std::cout << "id " << gp.id << " saving at t = " << t << std::endl;
      fname_p = prefix + std::string("_") +  std::to_string(it) +  std::string(".vti");
      complexfname_p = complexprefix + std::string("_") +  std::to_string(it) +  std::string(".vti");
      modulus /= running_average_count;
      if (gp.id == 0) {
	modulus(0,0,0) = 0.0;
      }
      
      running_average_count = 0;

      psPDE::ioVTK::writeVTKImageData(complexfname_p,{&modulus},modulus.grid);
      psPDE::ioVTK::writeVTKcollectionMiddle(complexcollection_name,complexfname_p,t);
      fftw_execute(backward_phi); // get phi(t+dt)

      integrator.initialize(modulus,0,0);

      psPDE::ioVTK::writeVTKImageData(fname_p,{&phi},phi.grid);
      psPDE::ioVTK::writeVTKcollectionMiddle(collection_name,fname_p,t);
    } else {
      fftw_execute(backward_phi); // get phi(t+dt)
    }

    //    link.wrap_X_i(X_i,phi.grid.get_Lx(),phi.grid.get_Ly(),phi.grid.get_Lz());
    //    free_energy_derivative
    //      = link.free_energy_derivative(X_i,integrator.linker_phi,
    //				    integrator.linker_derivative,phi);


  }

  if (gp.id == 0) {
    myfile.close();
  }
  
  psPDE::ioVTK::writeVTKcollectionFooter(collection_name);
  psPDE::ioVTK::writeVTKcollectionFooter(complexcollection_name);
  
  fftw_destroy_plan(forward_phi);
  fftw_destroy_plan(backward_phi);
  fftw_destroy_plan(forward_nonlinear);
  fftw_destroy_plan(backward_nonlinear);


  const double tot_integral = tt_integral/chronoitvl;
  const double tot_freeenergy = tt_freeenergy/chronoitvl;

  std::cout << "total integration time on process " << gp.id << " is " << tot_integral
	    << std::endl;

  std::cout << "total free energy time on process " << gp.id << " is " << tot_freeenergy
	    << std::endl;

  

  return;
}
