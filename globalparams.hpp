#ifndef GLOBALPARAMS_HPP
#define GLOBALPARAMS_HPP

#include <fstream>
#include <set>
#include <map>
#include <mpi.h>

#include "input.hpp"
#include "griddata.hpp"

class GlobalParams {
private:

  std::set<std::string> pset {"steps","volFrac",
      "dt", "dump_every", "dump_file", "thermo_every",
      "thermo_file","boxgrid", "boxdims",
      "seed","restart"};
  
  const int default_steps = 100;
  const int default_grid = 64;
  const double default_length = 100.0;
  const double default_dt = 1e-4;
  const int default_dump_every = 100;
  const int default_thermo_every = 100;
  const std::string default_dump_file{"dump"};
  const std::string default_thermo_file{"thermo"};
  const double default_volFrac = 0.3;
  const double default_variance = 0.0;
  const int default_seed = 129480;


  bool fft_is_transposed;
public:
  GlobalParams(const MPI_Comm, const int, const int,std::ifstream&,
	       std::map<std::string,std::string> const &,
	       std::string&);

  bool restart_flag;
  int startstep;
  double starttime;
  
  int steps,seed,dump_every,thermo_every;
  double dt, volFrac,variance;
  std::string dump_file,thermo_file;
  
  const MPI_Comm comm;
  const int id, mpi_size;

  GridData realspace;
  GridData fourier;
  


  void printall();

  
  
};

#endif
