#include <iostream>
#include <map>
#include <complex>
#include <chrono>
#include <cstring>
#include <fstream>

#include <fftw3-mpi.h>

#include "input.hpp"
#include "run.hpp"

#include "grid.hpp"
#include "domain.hpp"
#include "conjugate_volfrac.hpp"
#include "fixgrid_floryhuggins.hpp"

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

  std::ifstream infile;

  std::string nucfilename;

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


  std::string line;
  std::vector<std::string> v_line;

  while (std::getline(infile,line)) {
    
    if (line == "" || line[0] == '#') continue;

    input::convertVariables(line,variables);

    v_line = input::split_line(line);


    if (v_line[0] != "domain") {
      std::cerr << "error, first line of input file needs to start with domain."
		<< std::endl;
      MPI_Abort(comm,1);
      return EXIT_FAILURE;
    }

    v_line.erase(v_line.begin());

    break;
    
  }

  psPDE::Domain domain(id,mpi_size,v_line);
  
  while (std::getline(infile,line)) {

    if (line == "" || line[0] == '#') continue;

    input::convertVariables(line,variables);
    
    v_line = input::split_line(line);
    
    if (v_line[0] != "grid_style") {
      std::cerr << "error, second line of input file needs to start with grid_style."
		<< std::endl;
      MPI_Abort(comm,1);
      return EXIT_FAILURE;
    }

    v_line.erase(v_line.begin());

    break;
    
  }


  // need to wrap grid in braces to get 
  {
    psPDE::Grid grid(v_line,comm);
  
    while (std::getline(infile,line)) {
      
      if (line == "" || line[0] == '#') continue;
      
      input::convertVariables(line,variables);

      input::replacePercentages(line,id);
    
      v_line = input::split_line(line);
      
      
      if (v_line[0] != "grid_populate") {
	std::cerr << "error, third line of input file needs to start with grid_populate."
		  << std::endl;
	MPI_Abort(comm,1);
	return EXIT_FAILURE;
      }
      
      v_line.erase(v_line.begin());
      
      break;
      
    }

    grid.populate(v_line);

    psPDE::ConjugateVolFrac conjvfrac(domain,grid);

    while (std::getline(infile,line)) {
      
      if (line == "" || line[0] == '#') continue;

      input::convertVariables(line,variables);

      v_line = input::split_line(line);

      if (v_line[0] != "conjugate") {
	std::cerr << "error, fourth line of input file needs to start with conjugate."
		  << std::endl;
	MPI_Abort(comm,1);
	return EXIT_FAILURE;
      }
      v_line.erase(v_line.begin());

      if (v_line[0] != "volfrac") {
	std::cerr << "error, fourth line of input file needs to start with conjugate volfrac."
		  << std::endl;
	MPI_Abort(comm,1);
	return EXIT_FAILURE;
      }
      v_line.erase(v_line.begin());

      break;
    }

    conjvfrac.readCoeffs(v_line);

    psPDE::FixGridFloryHuggins fxgridFH;

    while (std::getline(infile,line)) {
      
      if (line == "" || line[0] == '#') continue;

      input::convertVariables(line,variables);

      v_line = input::split_line(line);

      if (v_line[0] != "fixgrid") {
	std::cerr << "error, fifth line of input file needs to start with fixgrid."
		  << std::endl;
	MPI_Abort(comm,1);
	return EXIT_FAILURE;
      }

      if (v_line[0] != "fixgrid") {
	std::cerr << "error, fifth line of input file needs to start with fixgrid."
		  << std::endl;
	MPI_Abort(comm,1);
	return EXIT_FAILURE;
      }
      v_line.erase(v_line.begin());

      if (v_line[0] != "floryhuggins") {
	std::cerr << "error, fifth line of input file needs to start with fixgrid floryhuggins."
		  << std::endl;
	MPI_Abort(comm,1);
	return EXIT_FAILURE;
      }
      v_line.erase(v_line.begin());

      break;
      
    }

    fxgridFH.readCoeffs(v_line);


    double dt;
    while (std::getline(infile,line)) {
      
      if (line == "" || line[0] == '#') continue;

      input::convertVariables(line,variables);

      v_line = input::split_line(line);

      if (v_line[0] != "dt") {

	std::cerr << "error, sixth line of input file needs to start with dt."
		  << std::endl;
	MPI_Abort(comm,1);
	return EXIT_FAILURE;
      }
      v_line.erase(v_line.begin());

      dt = std::stod(v_line[0]);
      
      break;
    }


    conjvfrac.reset_dt(dt);

    int Nsteps;

    while (std::getline(infile,line)) {
      
      if (line == "" || line[0] == '#') continue;

      input::convertVariables(line,variables);

      v_line = input::split_line(line);

      if (v_line[0] != "run") {    

	std::cerr << "error, sixth line of input file needs to start with dt."
		  << std::endl;
	MPI_Abort(comm,1);
	return EXIT_FAILURE;

      }

      v_line.erase(v_line.begin());

      Nsteps = std::stoi(v_line[0]);

      break;

    }


    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();  

    std::cout << "Running simulation of solution." << std::endl;
    
    run(domain,grid,conjvfrac,fxgridFH,dt,Nsteps);
    
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Run time = "
	      << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1e6
	      << "seconds." << std::endl;  

    

    
  }


  fftw_mpi_cleanup();

  ierr = MPI_Finalize();
  
  return 0;
}
