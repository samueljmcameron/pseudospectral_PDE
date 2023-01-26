#ifndef INPUT_HPP
#define INPUT_HPP

#include <string>
#include <vector>
#include <map>
#include <mpi.h>

namespace input {
  std::vector<std::string> split_line(std::string&);
  void replacePercentages(std::string &, int);

  void convertVariables(std::string &,
			std::map<std::string, std::string> const&);

  void replace_with_new_seed(std::vector<std::string> &,
			     const std::string &,int,int,
			     MPI_Comm )  ;
};
#endif
