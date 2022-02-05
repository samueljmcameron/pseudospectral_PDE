#include <iostream>

#include "solutionparams.hpp"


SolutionParams::SolutionParams(std::vector<std::string> splitvec)
  : mobility(default_mobility),volFH(default_volFH),gamma(default_gamma),
    chi(default_chi),temp(default_temp)
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


  
}
