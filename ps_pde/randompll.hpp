#ifndef PSPDE_RANDOMPLL_HPP
#define PSPDE_RANDOMPLL_HPP

#include <random>
#include <vector>


#include <mpi.h>

namespace psPDE {
class RandomPll {
public:
  RandomPll(const MPI_Comm &,const int, const int,const int);
  ~RandomPll();
  int get_processor_seed(); 
private:
  const int seed;
  const int id;
  std::mt19937 gen;

  std::uniform_int_distribution<int> integer_dist;

  std::vector<int> processor_seeds;
};
};
#endif
