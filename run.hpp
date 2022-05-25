#ifndef RUN_HPP
#define RUN_HPP

#include <vector>
#include "globalparams.hpp"
#include "solutionparams.hpp"

void run(GlobalParams gp, SolutionParams solparams,
	 std::vector<std::vector<double>> & X_is);

#endif
