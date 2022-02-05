#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include <random>
#include <complex>

#include "fftw_mpi_3darray.hpp"
#include "griddata.hpp"
#include "solutionparams.hpp"

class Integrator
{
private:
  int local0start;

  double qx,qy,qz,q2;

  const GridData fourier;

  const int seed;


  const double mobility;
  const double gamma;
  const double temp;
  const double chi;
  const double volFH;
  
  const double normalization;
  double dt,sqrtdt;
  std::complex<double> noise;

  double complexprefactor,realprefactor;
  
  std::mt19937 gen;

  std::uniform_real_distribution<double> real_dist;

  
public:
  void ode(std::complex<double> &, std::complex<double>,
	   std::complex<double>, double);
  Integrator(MPI_Comm,const GridData&,const int,
	     const SolutionParams&,const double);
  ~Integrator();
  fftw_MPI_3Darray<std::complex<double>> ft_phi;
  fftw_MPI_3Darray<std::complex<double>> ft_nonlinear;


  void initialize(fftw_MPI_3Darray<double> &,
		  const double , const double );
  
  void nonlinear(fftw_MPI_3Darray<double>&,
		 const fftw_MPI_3Darray<double>&);
    

  double get_dt() { return dt; };
  
  void update(int , int , int);
  void update_real(int , int , int);
};

#endif
