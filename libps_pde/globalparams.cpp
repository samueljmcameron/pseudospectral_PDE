#include <iostream>
#include <cmath>

#include "globalparams.hpp"

using namespace psPDE;


GlobalParams::GlobalParams(const MPI_Comm comm, const int id, const int mpi_size,
			   std::ifstream& input,
			   std::map<std::string,std::string> const& varMap,
			   std::string& lastline,
			   std::vector<std::string> EO_globals) :

  restart_flag(false),read_flag(false),all_nucs_flag(false),startstep(0),
  starttime(0.0),X_i_noise(0),steps(100),seed(129480),dump_every(100),
  thermo_every(100),dt(1e-4),volFrac(-1),variance(-1),dump_file(""),
  thermo_file(""),comm(comm), id(id),mpi_size(mpi_size),
  domain(100.0,100.0,100.0,-50.0,-50.0,-50.0,id,mpi_size),
  Nx(64),Ny(64),Nz(64)
{
  double Lx = domain.period[0];
  double Ly = domain.period[1];
  double Lz = domain.period[2];

  std::string line;  
  if (input) {

    
    std::vector<std::string> splitvec;

    
    
    while(std::getline(input,line)) {
      
      std::size_t found = line.find_first_of("#");
      line = line.substr(0,found);
      
      if (line == "")  continue;
      
      
      splitvec = input::split_line(line);
      
      int break_global = 0;
      // if build_solution, then global variable definitions are done
      for (auto &eo_glob : EO_globals)
	if (splitvec[0] == eo_glob) {
	  break_global = 1;
	  break;
	}

      if (break_global) break;
      
      if (splitvec[0] == "read") {
	input::convertVariable(splitvec.at(1),varMap);
	input::convertVariable(splitvec.at(2),varMap);
	input::convertVariable(splitvec.at(3),varMap);
	input::convertVariable(splitvec.at(4),varMap);

	read_flag = true;
	
	input::isDouble(splitvec.at(1),starttime,splitvec[0]);

	
	read_dump_file = splitvec.at(2);
	
	// replace any % with the processor number
	auto tmppos = read_dump_file.find("%");
	if (tmppos != read_dump_file.npos)
	  read_dump_file.replace(tmppos,1,std::to_string(id));
	
	// no need to replace % with processor number as thermo files are
	//  always on proc 0
	read_thermo_file = splitvec.at(3);

	// read in number of nucleation sites that are being kept and which ones
	int number_of_nucs_to_keep = 0;
	try {
	  input::isInt(splitvec.at(4),number_of_nucs_to_keep,splitvec[0]);
	} catch (std::runtime_error &err) {
	  if (splitvec.at(4) == "all")
	    all_nucs_flag = true;
	  else
	    throw std::runtime_error("Invalid input file.");

	}
	
	for (int index = 0; index < number_of_nucs_to_keep; index ++) {
	  int nucnum;
	  input::convertVariable(splitvec.at(5+index),varMap);
	  
	  input::isInt(splitvec.at(5+index),nucnum,splitvec[0]);
	  nucs_to_keep.push_back(nucnum);
	}


      } else if (splitvec.size() == 2) {
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
	} else if (splitvec[0] == "nucleation_noise") {
	  if (splitvec[1] == "on") 
	    X_i_noise = 1;	   
	  else if (splitvec[1] == "off") 
	    X_i_noise = 0;
	  else
	    throw std::runtime_error("Error: invalid input file.");
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

  domain = Domain(Lx,Ly,Lz,-Lx/2,-Ly/2,-Lz/2,id,mpi_size);



  lastline = line;


  if (read_flag && restart_flag) 
    throw std::runtime_error("Cannot restart and read simultaneously (choose one).");

  if (dump_file == "")
    throw std::runtime_error("Must specify dump_file.");

  if (thermo_file == "")
    throw std::runtime_error("Must specify thermo_file.");

  if (volFrac <= 0 && !read_flag) 
    throw std::runtime_error("Must specify volume fraction.");

  if (volFrac > 0 && read_flag)
    std::cerr << "Warning: volFrac specification is ignored when reading in data."
	      << std::endl;
  
  return; 
}


void GlobalParams::printall()
{
  std::cout << std::endl << std::endl;
  std::cout << "Global parameters in this simulation: " << std::endl << std::endl;

  std::cout << "steps: " << steps << "."
	    << std::endl;
  if (!read_flag) {
    std::cout << "volFrac: " << volFrac << "."
	      << std::endl;
    std::cout << "variance: " << variance << "."
	      << std::endl;
  }
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
  
  std::cout << "boxdims: (" << domain.period[0] << ","
	    << domain.period[1] << ","  << domain.period[2]
	    << ")."
	    << std::endl;

  std::cout << "boxorigin: (" << domain.boxlo[0] << ","
	    << domain.boxlo[1] << ","  << domain.boxlo[2]
	    << "). "  << std::endl;

  if (restart_flag) 
    std::cout << "Restarting from step " << startstep << " at time " << starttime
	      << "." << std::endl;
  if (read_flag) {
    std::cout << "Reading from files " << read_dump_file << " and " << read_thermo_file
	      << " to start the simulation at time " << starttime << "." << std::endl;
  }

}
