#include <fstream>


#include "fftw_mpi_3darray.hpp"
#include "integrator.hpp"
//#include "linker.hpp"
#include "conjplane.hpp"
#include "randompll.hpp"
#include "griddata.hpp"
#include "iovtk.hpp"
#include "timestep.hpp"

#include "run.hpp"

void X_1_of_t(std::vector<double> & X_1,
	      const std::vector<double> & dFdX_1,double dt)
{
  double fric = 100.0;
  X_1[0] += -dt*dFdX_1[0]*fric;
  X_1[1] += -dt*dFdX_1[1]*fric;
  X_1[2] += -dt*dFdX_1[2]*fric;
  return;
}


void X_2_of_t(std::vector<double> & X_2,
	      const std::vector<double> & dFdX_2,double dt)
{
  double fric = 100.0;
  X_2[0] += -dt*dFdX_2[0]*fric;
  X_2[1] += -dt*dFdX_2[1]*fric;
  X_2[2] += -dt*dFdX_2[2]*fric;
  return;
  
}

std::string getLastLine(std::string filename);

void put_in_vectors(std::vector<std::vector<double>> & X_is,
		    std::string filename, MPI_Comm comm,int mpi_id);


void run(GlobalParams gp, SolutionParams solparams,
	 std::vector<std::vector<double>> &X_is) {

  
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

  
  std::vector<std::vector<double>> free_energy_derivs;

  //  std::vector<std::vector<double>> free_energy_derivative;
  
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

  std::ofstream myfile;

  
  if (gp.restart_flag) {


    put_in_vectors(X_is,gp.thermo_file,gp.comm,gp.id);
    
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

    // insert file read here!



    

    if (gp.id == 0) {
      myfile.open(gp.thermo_file,std::ios::app);
      
    }

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
    }

    fftw_execute(backward_phi);
    
    ioVTK::writeVTKImageData(fname_p,{&phi},gp.realspace);
    ioVTK::writeVTKImageData(complexfname_p,{&modulus},modulus.grid);
    
    ioVTK::writeVTKcollectionMiddle(collection_name,fname_p,t);
    ioVTK::writeVTKcollectionMiddle(complexcollection_name,complexfname_p,t);


    
    if (gp.id == 0) {
      myfile.open(gp.thermo_file);
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




  
  TimeStep timestep(gp.comm,gp.mpi_size,gp.id,integrator.ft_phi.axis_size(0),
		    integrator.ft_phi.axis_size(1));

  double free_energy;
  for (int it = 1+gp.startstep; it <= gp.steps+gp.startstep; it ++) {


    free_energy_derivs = integrator.nonlinear(nonlinear,phi,X_is,free_energy); // compute nl(t) given phi(t)


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

    X_1_of_t(X_is[0],free_energy_derivs[0],integrator.get_dt());
    X_2_of_t(X_is[1],free_energy_derivs[1],integrator.get_dt());
    
    
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

  if (gp.id == 0) {
    myfile.close();
  }
  
  ioVTK::writeVTKcollectionFooter(collection_name);
  ioVTK::writeVTKcollectionFooter(complexcollection_name);
  
  fftw_destroy_plan(forward_phi);
  fftw_destroy_plan(backward_phi);
  fftw_destroy_plan(forward_nonlinear);
  fftw_destroy_plan(backward_nonlinear);


  return;
}



std::string getLastLine(std::string filename)
{

  std::ifstream myfile (filename);
  std::string lastline;

  if (myfile) {

    // assumes that last line is just a newline char, hence the -2 below (vs -1)
    myfile.seekg(-2,myfile.end);

    char ch;

    myfile.get(ch);

    while (ch != '\n' && myfile.tellg() >1) {
    
      myfile.seekg(-2,myfile.cur);
      myfile.get(ch);
    }




    std::getline(myfile,lastline);

  } else {
    throw std::runtime_error("File " + std::string(filename) + " does not exist");
  }

  myfile.close();
  return lastline;
}


void put_in_vectors(std::vector<std::vector<double>> & X_is,
		    std::string filename, MPI_Comm comm,int mpi_id)
{

  std::string finalline;
  int signal = 0;
  if (mpi_id == 0) {

    try {
      finalline = getLastLine(filename);
    }
    catch (const std::runtime_error & error) {
      signal = -1;
    }

  }

  MPI_Bcast(&signal,1,MPI_INT,0,comm);

  if (signal == -1) {
    throw std::runtime_error("File " + std::string(filename) + " does not exist");    
  }

  std::istringstream iss(finalline);
  
  std::string subs;
  if (mpi_id == 0) {
    iss >> subs;
  }

  for (auto &cmp : X_is) {
    for (int i = 0; i < 3; i++) {
      if (mpi_id == 0) {
	iss >> cmp.at(i);
      }

      MPI_Bcast(&cmp.at(i),1,MPI_DOUBLE,0,comm);
    }
  }

  return;
  
}
