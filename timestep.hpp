#ifndef TIMESTEP_HPP
#define TIMESTEP_HPP

#include <mpi.h>
#include "integrator.hpp"
#include "conjplane.hpp"

class TimeStep
{

  ConjPlane *cpcl;
  const MPI_Comm comm;
  const int id, mpi_size;

public:
 TimeStep(const MPI_Comm,const int , const int,
	  const int ,const int);
  ~TimeStep();
  void update(double ,Integrator &);
};

#endif
