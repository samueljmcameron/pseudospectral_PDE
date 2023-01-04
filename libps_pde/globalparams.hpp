#ifndef PSPDE_GLOBALPARAMS_HPP
#define PSPDE_GLOBALPARAMS_HPP

#include <fstream>
#include <set>
#include <map>
#include <mpi.h>
#include <vector>

#include "input.hpp"
#include "domain.hpp"

namespace psPDE {
class GlobalParams {
private:

  bool fft_is_transposed;
public:
  GlobalParams(const MPI_Comm, const int, const int,std::ifstream&,
	       std::map<std::string,std::string> const &,
	       std::string&,
	       std::vector<std::string> EO_globals = {"build_solution"});

  bool restart_flag,read_flag,all_nucs_flag;
  int startstep;
  double starttime;
  int X_i_noise;
  
  int steps,seed,dump_every,thermo_every;
  double dt, volFrac,variance;
  std::string dump_file,thermo_file;
  std::string read_dump_file, read_thermo_file;
  std::vector<int> nucs_to_keep;
  
  const MPI_Comm comm;
  const int id, mpi_size;

  int Nx,Ny,Nz;

  Domain domain,ft_domain;

  void printall();

  
  
};


};
#endif
