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
  int _globalNz;
  ptrdiff_t _size;
  ptrdiff_t *_sizeax;
  std::string array_name;
  int spacer;

  std::ostream& _write_ostream(std::ostream &out, std::vector<std::string> outer,
			       std::vector<std::string> middle,
			       std::vector<std::string> inner,std::string delim,
			       std::string innerdelim) const;
  
  
public:
  fftw_MPI_3Darray(MPI_Comm ,std::string ,const int , const int ,
		   const int );

  fftw_MPI_3Darray(MPI_Comm ,std::string ,const GridData& );

  
  ~fftw_MPI_3Darray();

  T* data() {
    return arr;
  };

  ptrdiff_t size() {
    return _size;
  };

  ptrdiff_t axis_size(int i)
  {
    return _sizeax[i];
  }
  
  ptrdiff_t get_globalNz() {
    return _globalNz;
  };

  ptrdiff_t get_local0start() {
    return local_0_start;
  };

  std::string get_name() {
    return array_name;
  }
    
  T& operator()(ptrdiff_t , ptrdiff_t , ptrdiff_t );
  T operator()(ptrdiff_t , ptrdiff_t , ptrdiff_t ) const;

  void copy(fftw_MPI_3Darray &);

  std::ostream& numpy_save(std::ostream & out);
  
  template <typename T0>
  friend std::ostream& operator<< (std::ostream& out,
				   const fftw_MPI_3Darray<T0> & in);


};





#endif
