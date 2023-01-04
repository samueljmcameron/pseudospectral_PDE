#ifndef PSPDE_DOMAIN_HPP
#define PSPDE_DOMAIN_HPP

#include "atom.hpp"
#include <array>
#include <Eigen/Core>

namespace psPDE {
class Domain
{
public:
  Domain(double,double,double,double,double,double,int,int);
  void pbc (Atom &) const;
  int set_image() const;
  const std::array<double,3> period,boxlo,boxhi;
  
  const int me,nprocs;

  void map(Eigen::Ref<Eigen::Vector3d>,
  	   const Eigen::Ref<const Eigen::Vector3d> &,int ) const;

  std::array<double,3> sublo,subhi;

  Domain fourier_transpose(int Nx,int Ny, int Nz)
  {

    double twopi = 3.141592*2;
    
    return Domain(twopi/period[0]*Nx,twopi/period[2]*Nz,twopi/period[1]*Ny,
		  0,0,0, me,nprocs);
  }

  Domain fourier_mod_transpose(int Nx, int Ny, int Nz) const
  {
    double twopi = 3.141592*2;

    Nx = Nx/2+1;

    double Lx = period[0]/2;

    return Domain(twopi/period[0]*Nx/2.0,twopi/period[2]*Nz,twopi/period[1]*Ny,
		  0,0,0, me,nprocs);
    

  }
  
private:
  void unmap(Eigen::Ref<Eigen::Vector3d>,
  	     const Eigen::Ref<const Eigen::Vector3d> &,int image) const;

};
}
#endif
