#include <iostream>
#include <cmath>

#include "globalparams.hpp"

using namespace psPDE;


GlobalParams::GlobalParams(const MPI_Comm comm, const int id, const int mpi_size,
			   std::ifstream& input,
			   std::map<std::string,std::string> const& varMap,
			   std::string& lastline) :

  restart_flag(false),startstep(0),starttime(0.0),
  steps(100),seed(129480),dump_every(100),thermo_every(100),
  dt(1e-4),volFrac(0.3),variance(0.0),dump_file("dump"),
  thermo_file("thermo"),comm(comm), id(id),mpi_size(mpi_size),
  realspace(64,64,64,100.0,100.0,100.0),fourier(realspace)
{


  double Lx = realspace.get_Lx();
  double Ly = realspace.get_Ly();
  double Lz = realspace.get_Lz();

  int Nx = realspace.get_Nx();
  int Ny = realspace.get_Ny();
  int Nz = realspace.get_Nz();

  std::string line;  
  if (input) {

    
    std::vector<std::string> splitvec;

    
    
    while(std::getline(input,line)) {
      
      std::size_t found = line.find_first_of("#");
      line = line.substr(0,found);
      
      if (line == "")  continue;
      
      
      splitvec = input::split_line(line);
      
      // if build_solution, then global variable definitions are done
      if (splitvec[0] == "build_solution") {
	break;
      }
      
      if (splitvec.size() == 2) {
	input::convertVariable(splitvec[1],varMap);
	
	if (splitvec[0] == "steps") {
	  input::isInt(splitvec[1],steps,splitvec[0]);
	} else if (splitvec[0] == "dt") {
	  input::isDouble(splitvec[1],dt,splitvec[0]);
	} else if (splitvec[0] == "dump_every") {
	  input::isInt(splitvec[1],dump_every,splitvec[0]);
	} else if (splitvec[0] == "thermo_every") {
	  input::isInt(splitvec[1],thermo_every,splitvec[0]);	  
	} else if (splitvec[0] == "dump_file") {
	  dump_file = splitvec[1];
	} else if (splitvec[0] == "thermo_file") {
	  thermo_file = splitvec[1];
	} else if (splitvec[0] == "seed") {
	  input::isInt(splitvec[1],seed,splitvec[0]);
	} else {
	  throw std::runtime_error("Error: invalid input file.");
	}


      } else if (splitvec.size() == 3) {
	input::convertVariable(splitvec[1],varMap);
	input::convertVariable(splitvec[2],varMap);	
	if (splitvec[0] == "volFrac") {
	  input::isDouble(splitvec[1],volFrac,splitvec[0]);
	  input::isDouble(splitvec[2],variance,splitvec[0]);
	} else	if (splitvec[0] == "restart") {
	  input::isInt(splitvec[1],startstep,splitvec[0]);
	  input::isDouble(splitvec[2],starttime,splitvec[0]);
	  restart_flag = true;
	} else {
	  throw std::runtime_error("Error: invalid input file.");
	}
	
	
      } else if (splitvec.size() == 4) {
	input::convertVariable(splitvec[1],varMap);
	input::convertVariable(splitvec[2],varMap);
	input::convertVariable(splitvec[3],varMap);
	
	if (splitvec[0] == "boxgrid") {
	  input::isInt(splitvec[1],Nx,splitvec[0]);
	  input::isInt(splitvec[2],Ny,splitvec[0]);
	  input::isInt(splitvec[3],Nz,splitvec[0]);
	} else if (splitvec[0] == "boxdims") {

	  input::isDouble(splitvec[1],Lx,splitvec[0]);
	  input::isDouble(splitvec[2],Ly,splitvec[0]);
	  input::isDouble(splitvec[3],Lz,splitvec[0]);
	} else {
	  throw std::runtime_error("Error: invalid input file.");
	}
      } else {
	throw std::runtime_error("Error: invalid input file.");
      }
      
    }
    
  }
  
  realspace = GridData(Nz,Ny,Nx,Lz,Ly,Lx);
  fourier = realspace.fft_grid();


  lastline = line;

  // later might make it an option to have the fft not output the transposed input,
  // but for now just set this flag to always being true
  fft_is_transposed = true;

  if (fft_is_transposed) {
    fourier.transpose_yz();
  }
  
  return; 
}


void GlobalParams::printall()
{
  std::cout << std::endl << std::endl;
  std::cout << "Global parameters in this simulation: " << std::endl << std::endl;

  std::cout << "steps: " << steps << "."
	    << std::endl;
  std::cout << "volFrac: " << volFrac << "."
	    << std::endl;
  std::cout << "variance: " << variance << "."
	    << std::endl;
  std::cout << "dt: " << dt << "."
	    << std::endl;
  std::cout << "dump_every: " << dump_every << "."
	    << std::endl;
  std::cout << "dump_file: " << dump_file << "."
	    << std::endl;

  std::cout << "thermo_every: " << thermo_every << "."
	    << std::endl;
  std::cout << "thermo_file: " << thermo_file << "."
	    << std::endl;

  std::cout << "seed: " << seed << "."
	    << std::endl;
  
  std::cout << "boxgrid: (" << realspace.get_Nx() << "," << realspace.get_Ny() << ","
	    << realspace.get_Nz() << ")."
	    << std::endl;
  std::cout << "boxdims: (" << realspace.get_Lx() << ","
	    << realspace.get_Ly() << ","  << realspace.get_Lz()
	    << ")."
	    << std::endl;

  std::cout << "boxorigin: (" << realspace.get_Ox() << ","
	    << realspace.get_Oy() << ","  << realspace.get_Oz()
	    << "). "  << std::endl;

  std::cout << "fouriergrid (just to confirm): (" << fourier.get_Nx() << "," << fourier.get_Ny()
	    << "," << fourier.get_Nz() << ")."
	    << std::endl;
  std::cout << "fourierdims (just to confirm): (" << fourier.get_Lx() << ","
	    << fourier.get_Ly() << ","  << fourier.get_Lz()
	    << ")."
	    << std::endl;

  std::cout << "fourierorigin: (" << fourier.get_Ox() << ","
	    << fourier.get_Oy() << ","  << fourier.get_Oz()
	    << "). "  << std::endl;


  std::cout << "Restart: " << (restart_flag ? "yes" : "no" ) << std::endl;
  std::cout << "startstep: " << startstep << std::endl << std::endl;
}
