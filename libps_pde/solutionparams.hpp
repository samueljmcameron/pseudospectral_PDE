#ifndef PSPDE_SOLUTIONPARAMS_HPP
#define PSPDE_SOLUTIONPARAMS_HPP

#include <vector>

#include "input.hpp"

namespace psPDE {
class SolutionParams {
  
public:
  SolutionParams();

  SolutionParams(std::vector<std::string> splitvec);


  double mobility, volFH, gamma,temp, chi;
  double chi_LP,chi_LL, nucwidth;
  double viscosity;
  
  void printall();  
};
};
#endif
