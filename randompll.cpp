#include "randompll.hpp"

RandomPll::RandomPll(MPI_Comm comm, const int id, const int seed, const int numprocs)
  : comm(comm),id(id),seed(seed), gen(seed), integer_dist(), processor_seeds(numprocs)
{

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
