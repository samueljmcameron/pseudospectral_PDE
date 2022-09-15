#ifndef PSPDE_INPUT_HPP
#define PSPDE_INPUT_HPP

#include <string>
#include <vector>
#include <map>
#include <fftw3-mpi.h>

namespace psPDE {
namespace input {
  std::vector<std::string> split_line(std::string&);
  bool isInt(std::string&,int&,std::string);
  bool isDouble(std::string&,double&,std::string);
  void convertVariable(std::string &,
		       std::map<std::string, std::string> const&);
  void put_in_vectors(std::vector<std::vector<double>> & X_is,
		      const std::vector<int> & nucs_to_keep,
		      std::string filename, MPI_Comm comm,int mpi_id,
		      double starttime);

  void put_in_vectors(std::vector<std::vector<double>> & X_is,
		      std::string filename, MPI_Comm comm,int mpi_id,
		      double starttime);

  
  std::string getLastLine(std::string filename);
  std::string getMatchingLine(std::string filename,double starttime);
}
};
#endif
