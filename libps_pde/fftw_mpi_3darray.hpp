#ifndef PSPDE_FFTW_MPI_3DARRAY_HPP
#define PSPDE_FFTW_MPI_3DARRAY_HPP

#include <fftw3-mpi.h>
#include <array>
#include <string>
#include <complex>

#include "domain.hpp"

namespace psPDE {

template <typename T>
class fftw_MPI_3Darray 
{
private:
  T *arr;
  ptrdiff_t alloc_local, local_0_start;

  ptrdiff_t size;
  const std::array<ptrdiff_t,3> boxgrid;
  std::array<ptrdiff_t,3> sizeax; // = {Nz,Ny,Nx} because for some stupid reason 12 months
  //                                  ago I thought x index should vary fastest
  
  
  std::string array_name;
  int spacer;

  MPI_Comm comm;

  void set_subdomain(Domain &) const;
  
public:

  fftw_MPI_3Darray(const MPI_Comm &,std::string,
		   ptrdiff_t, ptrdiff_t, ptrdiff_t, Domain & );
  fftw_MPI_3Darray(const fftw_MPI_3Darray<T> &,std::string name = "");


  fftw_MPI_3Darray<std::complex<double>> make_fourier_transpose(Domain &) const;
  fftw_MPI_3Darray<double> make_fourier_mod(Domain &) const;
  
  ~fftw_MPI_3Darray();


  const double dx,dy,dz;

  
  T* data() {
    return arr;
  };

  ptrdiff_t totalsize() const {
    return size;
  };


  ptrdiff_t Nz() const
  {
    return sizeax[0];
  }

  ptrdiff_t Ny() const
  {
    return sizeax[1];
  }

  ptrdiff_t Nx() const
  {
    return sizeax[2];
  }

  
  

  ptrdiff_t get_local0start() const {
    return local_0_start;
  };

  std::string get_name() const {
    return array_name;
  }
    
  T& operator()(ptrdiff_t , ptrdiff_t , ptrdiff_t );
  T operator()(ptrdiff_t , ptrdiff_t , ptrdiff_t ) const;

  T& operator()(ptrdiff_t );
  T operator()(ptrdiff_t ) const;

  fftw_MPI_3Darray<T>& operator=(fftw_MPI_3Darray<T> other);
  fftw_MPI_3Darray<T>& operator/=(T rhs);


  
  void abs(fftw_MPI_3Darray<double>&) const;
  void mod(fftw_MPI_3Darray<double>&) const;

  void setZero();

  void running_mod(fftw_MPI_3Darray<double>&) const;

  friend void swap(fftw_MPI_3Darray<T>& first, fftw_MPI_3Darray<T>& second)
  {
    using std::swap;

    swap(first.arr,second.arr);
    swap(first.alloc_local,second.alloc_local);
    swap(first.local_0_start,second.local_0_start);
    swap(first.size,second.size);
    swap(first.sizeax,second.sizeax);
    swap(first.array_name,second.array_name);
    swap(first.spacer,second.spacer);
    swap(first.comm,second.comm);

    return;

  }

  


};




};

#endif
