#include <iostream>

#include "solutionparams.hpp"


SolutionParams::SolutionParams(std::vector<std::string> splitvec)
  : mobility(default_mobility),volFH(default_volFH),gamma(default_gamma),
    chi(default_chi),temp(default_temp),chi_LP(default_chi_LP),nucmax(default_nucmax),nucwidth(default_nucwidth)
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
      throw std::runtime_error("Error: invalid argument for build_solution.");
    }
  }


}

void SolutionParams::printall()
{
  std::cout << std::endl << std::endl;
  std::cout << "Solution parameters in this simulation: " << std::endl << std::endl;

  std::cout << "mobility: " << mobility << ". (Default is " << default_mobility << ".)"
	    << std::endl;
  std::cout << "volFH: " << volFH << ". (Default is " << default_volFH << ".)"
	    << std::endl;
  std::cout << "gamma: " << gamma << ". (Default is " << default_gamma << ".)"
	    << std::endl;
  std::cout << "temp: " << temp << ". (Default is " << default_temp << ".)"
	    << std::endl;
  std::cout << "chi: " << chi << ". (Default is " << default_chi
	    << ".)"   << std::endl << std::endl;
  std::cout << "chi_LP: " << chi_LP << ". (Default is " << default_chi_LP
	    << ".)"   << std::endl << std::endl;
  std::cout << "chi_LL: " << chi_LL << ". (Default is " << default_chi_LL
	    << ".)"   << std::endl << std::endl;
  std::cout << "nucmax: " << nucmax << ". (Default is " << default_nucmax
	    << ".)"   << std::endl << std::endl;
  std::cout << "nucwidth: " << nucwidth << ". (Default is " << default_nucwidth
	    << ".)"   << std::endl << std::endl;


  
}
