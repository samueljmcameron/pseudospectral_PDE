#ifndef RANDOMPLL_HPP
#define RANDOMPLL_HPP

#include <random>
#include <vector>


#include <mpi.h>


class RandomPll {
public:
  RandomPll(MPI_Comm,const int, const int,const int);
  ~RandomPll();
  int get_processor_seed(); 
private:
  const int seed;
  const int id;
  MPI_Comm comm;
  std::mt19937 gen;

  std::uniform_int_distribution<int> integer_dist;

  std::vector<int> processor_seeds;
};

#endif
