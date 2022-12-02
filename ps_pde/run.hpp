#ifndef RUN_HPP
#define RUN_HPP

#include <vector>
#include "globalparams.hpp"
#include "solutionparams.hpp"

void run(psPDE::GlobalParams gp, psPDE::SolutionParams solparams,
	 std::vector<std::vector<double>> & X_is,
	 std::vector<double> &,
	 std::vector<double> & , std::vector<double> &);

#endif
