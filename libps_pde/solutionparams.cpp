#include <iostream>

#include "solutionparams.hpp"

using namespace psPDE;


SolutionParams::SolutionParams()
  : mobility(0.01),volFH(0.01),gamma(10.0),
    temp(4.114),chi(2.5),chi_LP(2.5),chi_LL(2.5),nucmax(0.9),
    nucwidth(1.0) {}



SolutionParams::SolutionParams(std::vector<std::string> splitvec)
  : mobility(0.01),volFH(0.01),gamma(10.0),
    temp(4.114),chi(2.5),chi_LP(2.5),chi_LL(2.5),nucmax(0.9),
    nucwidth(1.0)
{
  int nargs = splitvec.size();
  int iarg = 0;

  
  while (iarg < nargs) {
    if (splitvec[iarg] == "mobility") {
      input::isDouble(splitvec[iarg+1],mobility,splitvec[iarg]);
      iarg += 2;
    } else if (splitvec[iarg] == "volFH") {
      input::isDouble(splitvec[iarg+1],volFH,splitvec[iarg]);
      iarg += 2;
    } else if (splitvec[iarg] == "gamma") {
      input::isDouble(splitvec[iarg+1],gamma,splitvec[iarg]);
      iarg += 2;
    } else if (splitvec[iarg] == "temp") {
      input::isDouble(splitvec[iarg+1],temp,splitvec[iarg]);
      iarg += 2;
    } else if (splitvec[iarg] == "chi") {
      input::isDouble(splitvec[iarg+1],chi,splitvec[iarg]);
      iarg += 2;
    } else if (splitvec[iarg] == "chi_LP") {
      input::isDouble(splitvec[iarg+1],chi_LP,splitvec[iarg]);
      iarg += 2;
    } else if (splitvec[iarg] == "chi_LL") {
      input::isDouble(splitvec[iarg+1],chi_LL,splitvec[iarg]);
      iarg += 2;
    } else if (splitvec[iarg] == "nucmax") {
      input::isDouble(splitvec[iarg+1],nucmax,splitvec[iarg]);
      iarg += 2;
    } else if (splitvec[iarg] == "nucwidth") {
      input::isDouble(splitvec[iarg+1],nucwidth,splitvec[iarg]);
      iarg += 2;
    } else {
      std::cout << splitvec[iarg] << std::endl;
      throw std::runtime_error("Error: invalid argument for build_solution.");
    }
  }


}

void SolutionParams::printall()
{
  std::cout << std::endl << std::endl;
  std::cout << "Solution parameters in this simulation: " << std::endl << std::endl;

  std::cout << "mobility: " << mobility << "." << std::endl;
  std::cout << "volFH: " << volFH << "." << std::endl;
  std::cout << "gamma: " << gamma << "." << std::endl;
  std::cout << "temp: " << temp << "." << std::endl;
  std::cout << "chi: " << chi << "." << std::endl;
  std::cout << "chi_LP: " << chi_LP << "." << std::endl;
  std::cout << "chi_LL: " << chi_LL << "." << std::endl;
  std::cout << "nucmax: " << nucmax << "." << std::endl;
  std::cout << "nucwidth: " << nucwidth << "." << std::endl;

}
