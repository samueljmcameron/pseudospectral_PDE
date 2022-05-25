#include <iostream>
#include <map>
#include <complex>
#include <chrono>
#include <cstring>


//#include "fftw_mpi_3darray.hpp"
//#include "integrator.hpp"
//#include "conjplane.hpp"
//#include "randompll.hpp"
//#include "griddata.hpp"
//#include "iovtk.hpp"
#include <fftw3-mpi.h>
#include "solutionparams.hpp"
//#include "timestep.hpp"
#include "globalparams.hpp"
#include "run.hpp"

int main(int argc, char **argv)
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

  std::ifstream infile,nucfile;

  std::map<std::string,std::string> variables;

  std::string simulation_type = "run";

  int iarg = 1;  
  while(iarg < argc) {
    if (strcmp(argv[iarg],"-in") == 0) {
      if (iarg+1 == argc) {
	std::cerr << "Error: input flag '-in' specified, but no file given."
		  << std::endl;
	return EXIT_FAILURE;
      }
      infile.open(argv[iarg+1]);
      iarg += 2;
      
    } else if (strcmp(argv[iarg],"-nuc") == 0) {
      if (iarg+1 == argc) {
	std::cerr << "Error: input flag '-nuc' specified, but no file given."
		  << std::endl;
	return EXIT_FAILURE;
      }
      nucfile.open(argv[iarg+1]);
      iarg += 2;
      
    } else if (strcmp(argv[iarg],"-var") == 0) {
      
      if (iarg + 2 >= argc) {
	std::cerr << "Error: invalid command line variable specification."
		  << std::endl;
	return EXIT_FAILURE;
      }
      variables[argv[iarg+1]] = argv[iarg+2];
      iarg += 3;
    } else {
      std::cerr << "Error: invalid command line variable specification."
		<< std::endl;
      return EXIT_FAILURE;
    }
  }

  
  if (not infile.is_open()) {
    std::cerr << "Error: need to specify input file." << std::endl;
    return EXIT_FAILURE;
  }

  if (not nucfile.is_open()) {
    std::cerr << "Error: need to specify input file." << std::endl;
    return EXIT_FAILURE;
  }

  std::string line;
  std::vector<std::string> splitvec;
    
  GlobalParams gp(comm,id,mpi_size,infile,variables,line);



  line = line.substr(0,line.find_first_of("#"));
  splitvec = input::split_line(line);
  

  while (std::getline(infile,line) &&
	 splitvec[0] != "build_solution") {

    line = line.substr(0,line.find_first_of("#"));
    splitvec = input::split_line(line);
    
  }

  splitvec.erase(splitvec.begin());

  for (auto &c : splitvec)
    input::convertVariable(c,variables);
  
  
  SolutionParams solparams(splitvec);


  if (id == 0) {
    gp.printall();
    solparams.printall();
  }

  std::vector<std::vector<double>> X_is;
  

  double X_x,X_y,X_z;
  while(nucfile >> X_x >> X_y >> X_z) {
    X_is.push_back({X_x,X_y,X_z});
  }

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();  

  if (simulation_type == "run") {
    std::cout << "Running simulation of solution." << std::endl;
    run(gp,solparams,X_is);
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Run time = "
	    << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1e6
	    << "seconds." << std::endl;  



  fftw_mpi_cleanup();

  ierr = MPI_Finalize();
  
  return 0;
}


