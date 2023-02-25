#include <stdexcept>

#include "iovtk.hpp"
#include "input.hpp"
#include "run.hpp"

void run(psPDE::Domain &domain,psPDE::Grid &grid,psPDE::ConjugateVolFrac &conjvfrac,
	 psPDE::FixGridFloryHuggins &fxgridFH, double dt, int Nsteps)
{

  
  double t = 0;
  int step = 0;

  const double dx = domain.period[0]/grid.boxgrid[0];
  const double dy = domain.period[1]/grid.boxgrid[1];
  const double dz = domain.period[2]/grid.boxgrid[2];
  
  std::string prefix = "vtkfiles/unstable_p%";

  input::replacePercentages(prefix,domain.me);

  std::string fname = prefix + std::string("_") + std::to_string(0) + std::string(".vti");
  std::string cname = prefix + std::string(".pvd");

  if (domain.me == 0)
    std::cout << "Saving on step : 0" << std::endl;
  psPDE::ioVTK::writeVTKcollectionHeader(cname);
  psPDE::ioVTK::writeVTKImageData(fname,{grid.phi.get()},domain.boxlo,
				  {dx,dy,dz});

  psPDE::ioVTK::writeVTKcollectionMiddle(cname,fname,t);


  int errflag = 0;
  int totalerr;

  for (int step = 1; step <= Nsteps; step++) {


    grid.nonlinear->setZero();
    // flory huggins contribution
    fxgridFH.compute(grid);


    
    // FFT phi(r,t) 

    fftw_execute(grid.forward_phi);
    fftw_execute(grid.forward_nonlinear);


    // compute phi(q,t+dt) 

    conjvfrac.update();

    // IFFT back to get phi(r,t+dt) 
    fftw_execute(grid.backward_phi);

    if (std::isnan((*grid.phi)(0,0,0) )) 
      errflag = 1;

    MPI_Allreduce(&errflag,&totalerr,1,MPI_INT,MPI_SUM,grid.comm);

    if (totalerr) {
      psPDE::ioVTK::writeVTKcollectionFooter(cname);
      throw std::runtime_error("NAN encountered in phi.");

    }
    
    t += dt;


    if (step % 1000 == 0) {

      if (domain.me == 0)
	std::cout << "Saving on step : " << step << std::endl;

      
      std::string fname = prefix + std::string("_") + std::to_string(step)
	+ std::string(".vti");
      psPDE::ioVTK::writeVTKImageData(fname,{grid.phi.get()},domain.boxlo,
				      {dx,dy,dz});
      
      psPDE::ioVTK::writeVTKcollectionMiddle(cname,fname,t);

    }

  }

  psPDE::ioVTK::writeVTKcollectionFooter(cname);

  return;
}
