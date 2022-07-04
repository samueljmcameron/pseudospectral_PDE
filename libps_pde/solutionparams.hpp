#ifndef PSPDE_SOLUTIONPARAMS_HPP
#define PSPDE_SOLUTIONPARAMS_HPP

#include <vector>

#include "input.hpp"

namespace psPDE {
class SolutionParams {
  double default_mobility;
  double default_volFH;
  double default_gamma;
  double default_chi;
  double default_temp;
  double default_chi_LP;
  double default_chi_LL;
  double default_nucmax;
  double default_nucwidth;
  
public:
  SolutionParams();

  SolutionParams(std::vector<std::string> splitvec);


  double mobility, volFH, gamma,temp, chi;
  double chi_LP,chi_LL, nucmax, nucwidth;
  
  void printall();  
};
};
#endif
