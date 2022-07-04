#ifndef PSPDE_INTEGRATOR_HPP
#define PSPDE_INTEGRATOR_HPP

#include <random>
#include <complex>

#include "fftw_mpi_3darray.hpp"
#include "griddata.hpp"
#include "solutionparams.hpp"

namespace psPDE {
class Integrator
{
private:
  int local0start;

  double qx,qy,qz,q2;


  const int seed;


  const double mobility;
  const double gamma;
  const double temp;
  const double chi;
  const double volFH;
  const double chi_LP;
  const double nucmax;
  const double nucwidth;
  
  const double normalization;
  double dt,sqrtdt;
  std::complex<double> noise;

  double complexprefactor,realprefactor;
  
  std::mt19937 gen;

  std::uniform_real_distribution<double> real_dist;

  
  void ode(std::complex<double> &, std::complex<double>,
	   std::complex<double>, double);
  
public:

  Integrator(MPI_Comm,const GridData&,const int,
	     const SolutionParams&,const double);
  ~Integrator();
  fftw_MPI_3Darray<std::complex<double>> ft_phi;
  fftw_MPI_3Darray<std::complex<double>> ft_nonlinear;


  void initialize(fftw_MPI_3Darray<double> &,
		  const double , const double );
  
  std::vector<std::vector<double>> nonlinear(fftw_MPI_3Darray<double>&,
					     const fftw_MPI_3Darray<double>&,
					     const std::vector<std::vector<double>> &,
					     double &);
    

  double get_dt() { return dt; };
  
  void integrate(int , int , int);
  void integrate_real(int , int , int);


  double linker_phi(double , double , double ,double,double,double,
		    const std::vector<std::vector<double>> & );

  void linker_derivative(std::vector<double> &,double, double ,
			 double ,double,double, double,
			 const std::vector<double> &);

  

};
};
#endif
