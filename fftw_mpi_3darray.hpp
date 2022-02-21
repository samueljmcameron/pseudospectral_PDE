#ifndef FFTW_MPI_3DARRAY_HPP
#define FFTW_MPI_3DARRAY_HPP

#include <iostream>
#include <fftw3-mpi.h>
#include <vector>
#include <string>

#include "griddata.hpp"

template <typename T>
class fftw_MPI_3Darray 
{
private:
  T *arr;
  ptrdiff_t alloc_local, local_0_start;

  ptrdiff_t _size;
  ptrdiff_t *_sizeax;
  std::string array_name;
  int spacer;

  MPI_Comm comm;

  std::ostream& _write_ostream(std::ostream &out, std::vector<std::string> outer,
			       std::vector<std::string> middle,
			       std::vector<std::string> inner,std::string delim,
			       std::string innerdelim) const;
  void copy(const fftw_MPI_3Darray<T> &);
  
public:

  const GridData grid;
  
  fftw_MPI_3Darray(MPI_Comm ,std::string ,const int , const int ,
		   const int,const double,const double,const double );

  fftw_MPI_3Darray(MPI_Comm ,std::string ,const GridData& );

  
  fftw_MPI_3Darray(const fftw_MPI_3Darray<T> &);

  
  ~fftw_MPI_3Darray();

  T* data() {
    return arr;
  };

  ptrdiff_t size() const {
    return _size;
  };

  ptrdiff_t axis_size(int i) const
  {
    return _sizeax[i];
  }
  

  ptrdiff_t get_local0start() const {
    return local_0_start;
  };

  std::string get_name() const {
    return array_name;
  }
    
  T& operator()(ptrdiff_t , ptrdiff_t , ptrdiff_t );
  T operator()(ptrdiff_t , ptrdiff_t , ptrdiff_t ) const;

  //  template <typename T0>
  fftw_MPI_3Darray<T>& operator=(fftw_MPI_3Darray<T> other);

  fftw_MPI_3Darray<T>& operator/=(T rhs);
  void abs(fftw_MPI_3Darray<double>&) const;
  void mod(fftw_MPI_3Darray<double>&) const;

  void running_mod(fftw_MPI_3Darray<double>&) const;
  void split(fftw_MPI_3Darray<double> &,
	     fftw_MPI_3Darray<double> &) const;

  std::ostream& numpy_save(std::ostream & out) const;
  
  template <typename T0>
  friend std::ostream& operator<< (std::ostream& out,
				   const fftw_MPI_3Darray<T0> & in);



  friend void swap(fftw_MPI_3Darray<T>& first, fftw_MPI_3Darray<T>& second)
  {
    using std::swap;

    swap(first.arr,second.arr);
    swap(first.alloc_local,second.alloc_local);
    swap(first.local_0_start,second.local_0_start);
    swap(first._size,second._size);
    swap(first._sizeax,second._sizeax);
    swap(first.array_name,second.array_name);
    swap(first.spacer,second.spacer);
    swap(first.comm,second.comm);

    return;

  }

  


};





#endif
