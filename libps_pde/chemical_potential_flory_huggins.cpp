#include "chemical_potential_flory_huggins.hpp"
#include <cmath>

using namespace psPDE;

ChemPotFloryHuggins::ChemPotFloryHuggins()
  : ChemPot()
{
 
}


void ChemPotFloryHuggins::compute()
{

  //  std::vector<double> dlink(3);
  
  for (int i = 0; i < nonlinear.Nz(); i++) 
    for (int j = 0; j < nonlinear.Ny(); j++) 
      for (int k = 0; k < nonlinear.Nx(); k++)
	nonlinear(i,j,k)
	  = temp/volFH*(log(phi(i,j,k)/(1-phi(i,j,k)))+chi*(1-2*phi(i,j,k)));
  return;
}
