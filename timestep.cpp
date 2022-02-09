#include "timestep.hpp"

TimeStep::TimeStep(const MPI_Comm comm,const int mpi_size, const int id,
		   const int local_fftNz,const int local_fftNy)
  : comm(comm),mpi_size(mpi_size),id(id)
{
  
  std::vector<int> all_local_fftNzs(mpi_size);
  for (int i = 0; i < mpi_size; i++) {
    
    if (i == id) {
      all_local_fftNzs[i] = local_fftNz;
      
    }
    
    MPI_Bcast(&all_local_fftNzs[i],1,MPI_INT,i,comm);

  }


  
  cpcl = new ConjPlane(mpi_size,id,comm,all_local_fftNzs,local_fftNy);


}

TimeStep::~TimeStep()
{
  delete cpcl;
}

void TimeStep::update(double t,Integrator &integrator)
{

  int local_fftNz = integrator.ft_phi.axis_size(0);
  int local_fftNy = integrator.ft_phi.axis_size(1);
  int local_fftNx = integrator.ft_phi.axis_size(2);  
  int local0start = integrator.ft_phi.get_local0start(); 
  
  // update the bulk of the system which has no constraints
  
  for (int nz = 0; nz < local_fftNz; nz ++) {
    
    for (int ny = 0; ny < local_fftNy; ny ++) {
      
      for (int nx = 1; nx < local_fftNx-1; nx++) {
	integrator.integrate(nz,ny,nx);
      }
    }
  }
  
  
  
  // update the two planes at nx = 0 and nx = Nx-1
  // which are constrained to ensure that phi remains real
  
  for (int nx = 0; nx < local_fftNx; nx += local_fftNx-1) {
    
    if (mpi_size == 1) {
      
      cpcl->single(integrator,nx);
      cpcl->line_single(integrator,nx);
      
    } else {
      
      if (cpcl->is_equal() && id == 0) {
	cpcl->first(integrator,nx);
	cpcl->line_first(integrator,nx);
      } else if (id < cpcl->get_divider()) {
	cpcl->lefthalf(integrator,nx);
	cpcl->line_lefthalf(integrator,nx);
      } else if (id == cpcl->get_divider()) {
	if (mpi_size % 2 == 0) {
	  cpcl->middle_even(integrator,nx);
	  cpcl->line_middle_even(integrator,nx);
	} else {
	  cpcl->middle_odd(integrator,nx);
	  cpcl->line_middle_odd(integrator,nx);
	}
      } else if (!cpcl->is_equal() && id == mpi_size-1) {
	cpcl->last(integrator,nx);
	cpcl->line_last(integrator,nx);
      } else {
	cpcl->righthalf(integrator,nx);
	cpcl->line_righthalf(integrator,nx);
      }
    }
  }
  
  
  // update the eight points where ft_phi must be real
  if (id == 0) {
    int nx = local_fftNx-1;
    integrator.integrate_real(0,local_fftNy/2,0);
    integrator.integrate_real(0,local_fftNy/2,nx);
    integrator.integrate_real(0,0,nx);

    integrator.ft_phi(0,0,0) = integrator.ft_phi(0,0,0)
      /(1.0*integrator.ft_phi.grid.get_Nx()*integrator.ft_phi.grid.get_Ny()
	*integrator.ft_phi.grid.get_Nz());


  }
  if (local0start <= integrator.ft_phi.grid.get_Nz()/2
      && local0start + local_fftNz -1 >= integrator.ft_phi.grid.get_Nz()/2) {
    int nz = integrator.ft_phi.grid.get_Nz()/2 - local0start;
    int nx = local_fftNx-1;
    integrator.integrate_real(nz,0,0);
    integrator.integrate_real(nz,local_fftNy/2,0);
    integrator.integrate_real(nz,local_fftNy/2,nx);
    integrator.integrate_real(nz,0,nx);
  }
    




  
  return;
}
