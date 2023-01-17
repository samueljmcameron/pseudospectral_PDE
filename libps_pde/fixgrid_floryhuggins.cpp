#include "fixgrid_floryhuggins.hpp"
#include <cmath>
#include <stdexcept>


using namespace psPDE;

void FixGridFloryHuggins::readCoeffs(const std::vector<std::string> &v_line)
{


  temp = volFH = chi = 1;

  int iarg = 0;

  while (iarg < v_line.size()) {

    if (v_line[iarg] == "temp") {
      temp = std::stod(v_line[iarg+1]);
      iarg += 2;

    } else if (v_line[iarg] == "volFH") {
      volFH = std::stod(v_line[iarg+1]);
      iarg += 2;
    } else if (v_line[iarg] == "chi") {
      chi = std::stod(v_line[iarg+1]);
      iarg += 2;
    } else {
      throw std::runtime_error("Error: invalid fixgrid/floryhuggins command");
    }
  }
  return;
}


void FixGridFloryHuggins::compute(Grid &grid)
{

  fftw_MPI_3Darray<double> &phi = *(grid.phi);
  fftw_MPI_3Darray<double> &nonlinear = *(grid.nonlinear);
  
  
  for (int i = 0; i < nonlinear.Nz(); i++) 
    for (int j = 0; j < nonlinear.Ny(); j++)
      for (int k = 0; k < nonlinear.Nx(); k++)
	nonlinear(i,j,k) 
	  += temp/volFH*(log(phi(i,j,k)/(1-phi(i,j,k)))+chi*(1-2*phi(i,j,k)));

  return;

}
