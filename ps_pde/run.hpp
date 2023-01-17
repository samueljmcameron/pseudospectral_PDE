#ifndef RUN_HPP
#define RUN_HPP

#include <vector>

#include "grid.hpp"
#include "conjugate_volfrac.hpp"
#include "fixgrid_floryhuggins.hpp"

void run(psPDE::Grid &,psPDE::ConjugateVolFrac &,
	 psPDE::FixGridFloryHuggins &, double, int );

#endif
