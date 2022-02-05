#ifndef SOLUTIONPARAMS_HPP
#define SOLUTIONPARAMS_HPP

#include <vector>

#include "input.hpp"

class SolutionParams {
  const double default_mobility = 1.0;
  const double default_volFH = 0.01;
  const double default_gamma = 100.0;
  const double default_chi = 2.5;
  const double default_temp = 1.0;
  
public:
  SolutionParams(std::vector<std::string> splitvec);


  double mobility, volFH, gamma,temp, chi;
  
  void printall();  
};

#endif
