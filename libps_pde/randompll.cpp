#include "randompll.hpp"

using namespace psPDE;

RandomPll::RandomPll(MPI_Comm comm, const int id, const int seed, const int numprocs)
  : seed(seed),id(id),comm(comm), integer_dist(), processor_seeds(numprocs)
{

  gen.seed(seed);
  
  if (id == 0) {
    for (int i = 0; i < numprocs; i++) {
      processor_seeds[i] = integer_dist(gen);
    }
  }
  MPI_Bcast(processor_seeds.data(),processor_seeds.size(),MPI_INT,0,comm);
  
}


RandomPll::~RandomPll()
{
  return;
}

int RandomPll::get_processor_seed()
{

  return processor_seeds.at(id);
};
