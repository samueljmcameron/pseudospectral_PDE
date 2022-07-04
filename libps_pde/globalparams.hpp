#ifndef PSPDE_GLOBALPARAMS_HPP
#define PSPDE_GLOBALPARAMS_HPP

#include <fstream>
#include <set>
#include <map>
#include <mpi.h>

#include "input.hpp"
#include "griddata.hpp"
namespace psPDE {
class GlobalParams {
private:

  std::set<std::string> pset {"steps","volFrac",
      "dt", "dump_every", "dump_file", "thermo_every",
      "thermo_file","boxgrid", "boxdims",
      "seed","restart"};
  
  int default_steps;
  int default_grid;
  double default_length;
  double default_dt;
  int default_dump_every;
  int default_thermo_every;
  std::string default_dump_file;
  std::string default_thermo_file;
  double default_volFrac;
  double default_variance;
  int default_seed;


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


};
#endif
