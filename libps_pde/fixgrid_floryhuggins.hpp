#ifndef PSPDE_FIXGRID_FLORYHUGGINS_HPP
#define PSPDE_FIXGRID_FLORYHUGGINS_HPP

#include <cmath>
#include <vector>
#include <string>

#include "grid.hpp"

namespace psPDE{

class FixGridFloryHuggins {
public:
  void readCoeffs(const std::vector<std::string> &);
  void compute(Grid &);
  int num_less_zero,num_great_zero;
private:
  double temp,volFH,chi;

};

}
#endif
