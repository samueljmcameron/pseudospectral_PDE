#include <iostream>
#include <cmath>

#include "globalparams.hpp"

GlobalParams::GlobalParams(const MPI_Comm comm, const int id, const int mpi_size,
			   std::ifstream& input,
			   std::map<std::string,std::string> const& varMap,
			   std::string& lastline) :
  steps(default_steps),dt(default_dt),volFrac(default_volFrac),
  dump_every(default_dump_every), dump_file(default_dump_file),
  id(id),comm(comm),mpi_size(mpi_size),variance(default_variance),
  restart_flag(false),starttime(0.0),startstep(0),seed(default_seed),
  realspace(default_grid,default_length),
  fourier(GridData(default_grid,default_length).fft_grid())
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
      
      // check if keyword is global parameter, assume if not then
      // no longer defining global params
      if (pset.find(splitvec[0]) == pset.end()) {
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
	} else if (splitvec[0] == "dump_file") {
	  dump_file = splitvec[1];
	} else if (splitvec[0] == "seed") {
	  input::isInt(splitvec[1],seed,splitvec[0]);
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


	  
	}
      } else {
	throw std::runtime_error("Error: invalid input file.");
      }
      
    }
    
  }
  
  realspace.reset(Nz,Ny,Nx,Lz,Ly,Lx);
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

  std::cout << "steps: " << steps << ". (Default is " << default_steps << ".)"
	    << std::endl;
  std::cout << "volFrac: " << volFrac << ". (Default is " << default_volFrac << ".)"
	    << std::endl;
  std::cout << "variance: " << variance << ". (Default is " << default_variance << ".)"
	    << std::endl;
  std::cout << "dt: " << dt << ". (Default is " << default_dt << ".)"
	    << std::endl;
  std::cout << "dump_every: " << dump_every << ". (Default is " << default_dump_every
	    << ".)"   << std::endl;
  std::cout << "dump_file: " << dump_file << ". (Default is " << default_dump_file << ".)"
	    << std::endl;
  std::cout << "seed: " << seed << ". (Default is " << default_seed << ".)"
	    << std::endl;
  
  
  std::cout << "boxgrid: (" << realspace.get_Nx() << "," << realspace.get_Ny() << ","
	    << realspace.get_Nz() << "). (Default is " << default_grid << ".)"
	    << std::endl;
  std::cout << "boxdims: (" << realspace.get_Lx() << ","
	    << realspace.get_Ny() << ","  << realspace.get_Nz()
	    << "). (Default is " << default_length << ".)"
	    << std::endl;

  std::cout << "fouriergrid (just to confirm): (" << fourier.get_Nx() << "," << fourier.get_Ny()
	    << "," << fourier.get_Nz() << "). (Default is " << default_grid << ".)"
	    << std::endl;
  std::cout << "fourierdims (just to confirm): (" << fourier.get_Lx() << ","
	    << fourier.get_Ny() << ","  << fourier.get_Nz()
	    << "). (Default is " << default_length << ".)"
	    << std::endl;

  
  std::cout << "Restart: " << restart_flag << std::endl;
  std::cout << "startstep: " << startstep << std::endl << std::endl;
}
