#include <complex>
#include <typeinfo>
#include "fftw_mpi_3darray.hpp"

using namespace psPDE;
template <typename T>
fftw_MPI_3Darray<T>::fftw_MPI_3Darray(MPI_Comm comm,std::string name,
				      const int Nz,const int Ny,
				      const int Nx,const double Lz,
				      const double Ly, const double Lx)
  : comm(comm), grid(Nz,Ny,Nx,Lz,Ly,Lx)
/*
  Constructor for a 3D array with axis sizes (Nz,Ny,Nx) (the x dimension
  varies the quickest). The array is not contiguous in memory for different
  y and z indices. Values are either all doubles, or all std::complex.

  Parameters
  ----------
  comm : mpi communicator
      Typically MPI_COMM_WORLD, but could be other I suppose.
  name : string
      The name of the array (useful when needing to save data).
  Nz : const int
      The axis size in the z dimension (first index).
  Ny : const int
      The axis size in the y dimension (second index).
  Nx : const int
      The axis size in the x direction (third index).

*/
{


  ptrdiff_t local_n0;


  alloc_local = fftw_mpi_local_size_3d(Nz, Ny , Nx/2 + 1, comm,
				       &local_n0,&local_0_start);
  
  
  _sizeax = new ptrdiff_t[3];
  _sizeax[0] = local_n0;
  _sizeax[1] = Ny;
  
  
  if (typeid(T) == typeid(double)) {
    _sizeax[2] = Nx;
    arr = (T*) fftw_alloc_real(2*alloc_local);
    spacer = 2*(Nx/2+1);    
    _size = spacer*alloc_local;
    
  } else if (typeid(T) == typeid(std::complex<double>)) {
    _sizeax[2] = Nx/2+1;
    spacer=Nx/2+1;
    arr = (T*) fftw_alloc_complex(alloc_local);
    _size = alloc_local;
  } else
    throw std::runtime_error("fftw_MPI_3Darray can only have type double, "
			     "or std::complex<double>.");
  
  array_name = name;

};


template <typename T>
fftw_MPI_3Darray<T>::fftw_MPI_3Darray(MPI_Comm comm,std::string name,
				      const GridData& grid)
  : comm(comm), grid(grid)
/*
  Constructor for a 3D array with axis sizes (Nz,Ny,Nx) (the x dimension
  varies the quickest). The array is not contiguous in memory for different
  y and z indices. Values are either all doubles, or all std::complex.

  Parameters
  ----------
  comm : mpi communicator
      Typically MPI_COMM_WORLD, but could be other I suppose.
  name : string
      The name of the array (useful when needing to save data).
  Nz : const int
      The axis size in the z dimension (first index).
  Ny : const int
      The axis size in the y dimension (second index).
  Nx : const int
      The axis size in the x direction (third index).

*/
{


  ptrdiff_t local_n0;

  alloc_local = fftw_mpi_local_size_3d(grid.get_Nz(), grid.get_Ny() , grid.get_Nx()/2 + 1, comm,
				       &local_n0,&local_0_start);
  
  

  _sizeax = new ptrdiff_t[3];
  _sizeax[0] = local_n0;
  _sizeax[1] = grid.get_Ny();
  
  
  if (typeid(T) == typeid(double)) {
    _sizeax[2] = grid.get_Nx();
    arr = (T*) fftw_alloc_real(2*alloc_local);
    spacer = 2*(grid.get_Nx()/2+1);    
    _size = spacer*alloc_local;
    
  } else if (typeid(T) == typeid(std::complex<double>)) {
    _sizeax[2] = grid.get_Nx()/2+1;
    spacer=grid.get_Nx()/2+1;
    arr = (T*) fftw_alloc_complex(alloc_local);
    _size = alloc_local;
  } else
    throw std::runtime_error("fftw_MPI_3Darray can only have type double, "
			     "or std::complex<double>.");
  
  array_name = name;

};

template <typename T>
fftw_MPI_3Darray<T>::fftw_MPI_3Darray(const fftw_MPI_3Darray<T> & base)
  : alloc_local(base.alloc_local),local_0_start(base.local_0_start),
    _size(base._size),array_name(base.array_name),
    spacer(base.spacer),comm(base.comm),grid(base.grid)
{

  _sizeax = new ptrdiff_t[3];
  _sizeax[0] = base._sizeax[0];
  _sizeax[1] = base._sizeax[1];
  _sizeax[2] = base._sizeax[2];
  
  
  if (typeid(T) == typeid(double)) {

    arr = (T*) fftw_alloc_real(2*alloc_local);
    
  } else if (typeid(T) == typeid(std::complex<double>)) {


    arr = (T*) fftw_alloc_complex(alloc_local);

  } else
    throw std::runtime_error("fftw_MPI_3Darray can only have type double, "
			     "or std::complex<double>.");
  
  copy(base);
}


template <typename T>
fftw_MPI_3Darray<T>::~fftw_MPI_3Darray() {
  delete [] _sizeax;
  fftw_free(arr);
};

template <typename T>
fftw_MPI_3Darray<T>& fftw_MPI_3Darray<T>::operator/=(T rhs)
{
  for (int i = 0; i < _sizeax[0]; i++) {
    for (int j = 0; j < _sizeax[1]; j++) {
      for (int k = 0; k < _sizeax[2]; k++) {
	arr[k + (i*_sizeax[1] + j ) * spacer] /= rhs;
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

  return arr[k + (i*_sizeax[1] + j ) * spacer];
}


template <typename T>
T fftw_MPI_3Darray<T>::operator()(ptrdiff_t i,
				  ptrdiff_t j, ptrdiff_t k) const
/* read-only access (but don't change) array elements via (i,j,k). */
{
  return arr[k + (i*_sizeax[1] + j ) * spacer];

}



template <typename T>
void fftw_MPI_3Darray<T>::copy(const fftw_MPI_3Darray<T> &input) 
/* explicitly copy elements of input to arr. */
{

  for (int i = 0; i < 3; i++) {
    
    if (_sizeax[i] != input.axis_size(i))
      throw std::runtime_error("Cannot copy fftw_MPI_3Darray (wrong shape).");
  }
  for (int i = 0; i < _sizeax[0]; i++) {
    for (int j = 0; j < _sizeax[1]; j++) {
      for (int k = 0; k < _sizeax[2]; k++) {
	arr[k + (i*_sizeax[1] + j ) * spacer] = input(i,j,k);
      }
    }
  }
  return;
};



template <typename T>
std::ostream& operator<< (std::ostream& out,
			  const fftw_MPI_3Darray<T> &in)
{
  std::string delim = ", ";
  std::string innerdelim = ", ";
  std::vector<std::string> outer(2), middle(2), inner(2);
  
  outer[0]= middle[0] = inner[0] = "[ ";
  outer[1]= middle[1] = inner[1] = "] ";


  return in._write_ostream(out,outer,middle,inner,
			   delim,innerdelim);
}



template <typename T>
std::ostream& fftw_MPI_3Darray<T>::numpy_save(std::ostream &out) const
{
  
  std::string delim = "\t";
  std::string innerdelim = "\n";

  std::vector<std::string> outer(2), middle(2), inner(2);
  outer[0] = "";
  middle[0] = "#start of new first index\n";
  inner[0] = "#start of new second index\n";
  
  outer[1] = middle[1] = inner[1] = "";


  return _write_ostream(out,outer,middle,inner,
			delim,innerdelim);

  
};

template <typename T>
void fftw_MPI_3Darray<T>::abs(fftw_MPI_3Darray<double>& modulus) const
{

  for (int i = 0; i < 3; i++) {
    
    if (axis_size(i) != modulus.axis_size(i))
      throw std::runtime_error("Cannot take abs of fftw_MPI_3Darray (wrong output shape).");
  }


  for (int i = 0; i < _sizeax[0]; i++) {
    for (int j = 0; j < _sizeax[1]; j++) {
      for (int k = 0; k < _sizeax[2]; k++) {
	modulus(i,j,k) = std::abs(arr[k + (i*_sizeax[1] + j ) * spacer]);
      }
    }
  }
  return;

}

template <typename T>
void fftw_MPI_3Darray<T>::mod(fftw_MPI_3Darray<double>& modulus) const
{

  for (int i = 0; i < 3; i++) {
    
    if (axis_size(i) != modulus.axis_size(i))
      throw std::runtime_error("Cannot take abs of fftw_MPI_3Darray (wrong output shape).");
  }


  for (int i = 0; i < _sizeax[0]; i++) {
    for (int j = 0; j < _sizeax[1]; j++) {
      for (int k = 0; k < _sizeax[2]; k++) {
	modulus(i,j,k) = std::abs(arr[k + (i*_sizeax[1] + j ) * spacer])
	  *std::abs(arr[k + (i*_sizeax[1] + j ) * spacer]);
      }
    }
  }
  return;

}


template <typename T>
void fftw_MPI_3Darray<T>::running_mod(fftw_MPI_3Darray<double>& modulus) const
{

  for (int i = 0; i < 3; i++) {
    
    if (axis_size(i) != modulus.axis_size(i))
      throw std::runtime_error("Cannot take abs of fftw_MPI_3Darray (wrong output shape).");
  }


  for (int i = 0; i < _sizeax[0]; i++) {
    for (int j = 0; j < _sizeax[1]; j++) {
      for (int k = 0; k < _sizeax[2]; k++) {
	modulus(i,j,k) += std::abs(arr[k + (i*_sizeax[1] + j ) * spacer])
	  *std::abs(arr[k + (i*_sizeax[1] + j ) * spacer]);
      }
    }
  }
  return;

}



template <typename T>
fftw_MPI_3Darray<T>& fftw_MPI_3Darray<T>::operator=(fftw_MPI_3Darray<T> other)
{
  swap(*this,other);
  
  return *this;
}

template <typename T>
std::ostream& fftw_MPI_3Darray<T>::_write_ostream(std::ostream &out,
						  std::vector<std::string> outer,
						  std::vector<std::string> middle,
						  std::vector<std::string> inner,
						  std::string delim,
						  std::string innerdelim) const
{
  
  int L0 = _sizeax[0];
  int M0 = _sizeax[1];
  int N0 = _sizeax[2];
  
  out << outer[0];
  for (int i = 0; i < L0-1; i++) {
    out << middle[0];
    for (int j = 0; j < M0-1; j++) {
      out << inner[0];
      for (int k = 0; k < N0-1; k++) {
	out << arr[k + (i*M0 + j ) * spacer] << delim;
      }
      out << arr[N0-1 + (i*M0 + j ) * spacer] <<  inner[1]
	  << innerdelim;
    }
    
    int jdum = M0 -1;
    out << inner[0];
    for (int k = 0; k < N0-1; k++) {
      out << arr[k + (i*M0 + jdum ) * spacer] << delim;
    }
    out << arr[N0-1 + (i*M0 + jdum ) * spacer]
	<< inner[1] << middle[1] << innerdelim ;
  }
  
  
  
  int i = L0-1;
  out << outer[0];
  for (int j = 0; j < M0-1; j++) {
    out << middle[0];
    for (int k = 0; k < N0-1; k++) {
      out << arr[k + (i*M0 + j ) * spacer] << delim;
    }
    out << arr[N0-1 + (i*M0 + j ) * spacer] <<  inner[1]
	<< innerdelim;
  }

  out << inner[0];
  for (int k = 0; k < N0-1; k++) {
    out << arr[k + (i*M0 + M0-1 ) * spacer] << delim;
  }
  out << arr[N0-1 + (i*M0 + M0-1 ) * spacer]
      << inner[1] << middle[1] << outer[1];
  
  
  return out;
  
}


template class fftw_MPI_3Darray<double>;
template class fftw_MPI_3Darray<std::complex<double>>;
//template std::ostream& operator<< (std::ostream& ,
//				   const fftw_MPI_3Darray<double>&);
//template std::ostream& operator<< (std::ostream& ,
//				   const fftw_MPI_3Darray<std::complex<double>>&);


//template std::ostream& operator= (fftw_MPI_3Darray<double>);
//template std::ostream& operator= (fftw_MPI_3Darray<std::complex<double>>);

// Template below doesn't work because it forces return of double [2] array.
//template class fftw_MPI_3Darray<fftw_complex>; 
