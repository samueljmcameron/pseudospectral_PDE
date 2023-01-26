#include <typeinfo>
#include "fftw_mpi_3darray.hpp"

using namespace psPDE;


template <typename T>
fftw_MPI_3Darray<T>::fftw_MPI_3Darray(const MPI_Comm &comm,std::string name,
				      ptrdiff_t iNx, ptrdiff_t iNy, ptrdiff_t iNz)
/*
  Constructor for a 3D array with axis sizes (iNz,iNy,iNx) (the x dimension
  varies the quickest). The array is not contiguous in memory for different
  y and z indices. Values are either all doubles, or all std::complex.

  Parameters
  ----------
  comm : mpi communicator
      Typically MPI_COMM_WORLD, but could be other I suppose.
  name : string
      The name of the array (useful when needing to save data).
      
      
*/
{


  ptrdiff_t local_n0;

  alloc_local = fftw_mpi_local_size_3d(iNz, iNy , iNx/2 + 1,
				       comm,&local_n0,&local_0_start);
  

  sizeax[0] = local_n0;
  sizeax[1] = iNy;
  
  
  if (typeid(T) == typeid(double)) {
    sizeax[2] = iNx;
    arr = (T*) fftw_alloc_real(2*alloc_local);
    spacer = 2*(iNx/2+1);    
    size = spacer*alloc_local;
    
  } else if (typeid(T) == typeid(std::complex<double>)) {
    sizeax[2] = iNx/2+1;
    spacer=iNx/2+1;
    arr = (T*) fftw_alloc_complex(alloc_local);
    size = alloc_local;
  } else
    throw std::runtime_error("fftw_MPI_3Darray can only have type double, "
			     "or std::complex<double>.");
  
  array_name = name;
  
};




template <typename T>
fftw_MPI_3Darray<T>::fftw_MPI_3Darray(const fftw_MPI_3Darray<T> & base,
				      std::string name)
  : alloc_local(base.alloc_local),local_0_start(base.local_0_start),
    size(base.size),array_name(base.array_name), spacer(base.spacer)
/*
  Copy array, but if name (other than "") is provided then only make an
  array of the same size with the new name, but don't copy the elements in the
  array.

  Parameters
  ----------
  base : fftw_MPI_3Darray
      The array to either copy from or set size structures from
  name : std::string (optional)
      The name of the new array, default is "" which makes this a true copy constructor.
      
*/


{


  sizeax[0] = base.sizeax[0];
  sizeax[1] = base.sizeax[1];
  sizeax[2] = base.sizeax[2];

  ptrdiff_t tmpsize;
  
  if (typeid(T) == typeid(double)) {

    arr = (T*) fftw_alloc_real(2*alloc_local);

    tmpsize = 2*alloc_local;
    
  } else if (typeid(T) == typeid(std::complex<double>)) {
    arr = (T*) fftw_alloc_complex(alloc_local);
    tmpsize = alloc_local;
  } else
    throw std::runtime_error("fftw_MPI_3Darray can only have type double, "
			     "or std::complex<double>.");

  if (name != "") array_name = name;
  else std::copy(base.arr,base.arr+tmpsize,arr);
    

  
}








template <typename T>
fftw_MPI_3Darray<T>::~fftw_MPI_3Darray() {
  fftw_free(arr);
};

template <typename T>
fftw_MPI_3Darray<T>& fftw_MPI_3Darray<T>::operator/=(T rhs)
/* Division of all elements in the array by the same value. */
{
  for (int i = 0; i < sizeax[0]; i++) {
    for (int j = 0; j < sizeax[1]; j++) {
      for (int k = 0; k < sizeax[2]; k++) {
	arr[k + (i*sizeax[1] + j ) * spacer] /= rhs;
      }
    }
  }
  return *this;
};


template <typename T>
T& fftw_MPI_3Darray<T>::operator()(ptrdiff_t i,
				   ptrdiff_t j, ptrdiff_t k)
/* read/write access array elements via (i,j,k). */
{

  return arr[k + (i*sizeax[1] + j ) * spacer];
}


template <typename T>
T fftw_MPI_3Darray<T>::operator()(ptrdiff_t i,
				  ptrdiff_t j, ptrdiff_t k) const
/* read-only access (but don't change) array elements via (i,j,k). */
{
  return arr[k + (i*sizeax[1] + j ) * spacer];

}



template <typename T>
T& fftw_MPI_3Darray<T>::operator()(ptrdiff_t flat)
/* read/write access array elements via flattened array. */
{

  int k = flat % sizeax[2];
  int j = (flat/sizeax[2]) % sizeax[1];
  int i = (flat/sizeax[2])/sizeax[1];

  return arr[k + (i*sizeax[1] + j ) * spacer];
}


template <typename T>
T fftw_MPI_3Darray<T>::operator()(ptrdiff_t flat) const
/* read-only access (but don't change) array elements via flattened array. */
{

  int k = flat % sizeax[2];
  int j = (flat/sizeax[2]) % sizeax[1];
  int i = (flat/sizeax[2])/sizeax[1];

  return arr[k + (i*sizeax[1] + j ) * spacer];

}

template <typename T>
void fftw_MPI_3Darray<T>::setZero()
{

  for (int i = 0; i < sizeax[0]; i++) {
    for (int j = 0; j < sizeax[1]; j++) {
      for (int k = 0; k < sizeax[2]; k++) {
	arr[k + (i*sizeax[1] + j) * spacer] = 0.0;
      }
    }
  }
  return;

}


template <typename T>
void fftw_MPI_3Darray<T>::abs(fftw_MPI_3Darray<double>& modulus) const
{


    
  if (Nz() != modulus.Nz() || Ny() != modulus.Ny() || Nx() != modulus.Nx())
      throw std::runtime_error("Cannot take abs of fftw_MPI_3Darray (wrong output shape).");



  for (int i = 0; i < sizeax[0]; i++) {
    for (int j = 0; j < sizeax[1]; j++) {
      for (int k = 0; k < sizeax[2]; k++) {
	modulus(i,j,k) = std::abs(arr[k + (i*sizeax[1] + j ) * spacer]);
      }
    }
  }
  return;

}

template <typename T>
void fftw_MPI_3Darray<T>::mod(fftw_MPI_3Darray<double>& modulus) const
{


  if (Nz() != modulus.Nz() || Ny() != modulus.Ny() || Nx() != modulus.Nx())  
    throw std::runtime_error("Cannot take abs of fftw_MPI_3Darray (wrong output shape).");



  for (int i = 0; i < sizeax[0]; i++) {
    for (int j = 0; j < sizeax[1]; j++) {
      for (int k = 0; k < sizeax[2]; k++) {
	modulus(i,j,k) = std::abs(arr[k + (i*sizeax[1] + j ) * spacer])
	  *std::abs(arr[k + (i*sizeax[1] + j ) * spacer]);
      }
    }
  }
  return;

}


template <typename T>
void fftw_MPI_3Darray<T>::running_mod(fftw_MPI_3Darray<double>& modulus) const
{

  if (Nz() != modulus.Nz() || Ny() != modulus.Ny() || Nx() != modulus.Nx())  
    throw std::runtime_error("Cannot take abs of fftw_MPI_3Darray (wrong output shape).");


  for (int i = 0; i < sizeax[0]; i++) {
    for (int j = 0; j < sizeax[1]; j++) {
      for (int k = 0; k < sizeax[2]; k++) {
	modulus(i,j,k) += std::abs(arr[k + (i*sizeax[1] + j ) * spacer])
	  *std::abs(arr[k + (i*sizeax[1] + j ) * spacer]);
      }
    }
  }
  return;

}

template <typename T>
void fftw_MPI_3Darray<T>::reverseFlat(int gridindex, int &i, int &j, int &k) const
{

  k = gridindex % sizeax[2];
  j = (gridindex / sizeax[2]) % sizeax[1];
  i = (gridindex / sizeax[2]) / sizeax[1];

}


template <typename T>
fftw_MPI_3Darray<T>& fftw_MPI_3Darray<T>::operator=(fftw_MPI_3Darray<T> other)
{
  swap(*this,other);
  
  return *this;
}



template class fftw_MPI_3Darray<double>;
template class fftw_MPI_3Darray<std::complex<double>>;

// Template below doesn't work because it forces return of double [2] array.
//template class fftw_MPI_3Darray<fftw_complex>; 
