#include <complex>
#include <stdexcept>

#include "smatrix.hpp"

using namespace psPDE;

template <typename T>
sMatrix<T>::sMatrix()
/*
  Constructor for a 2D array without setting axis sizes.
  Values are either all doubles, or all std::complex.

*/
{
  allocd = false;

};



template <typename T>
sMatrix<T>::sMatrix(const int Ny, const int Nx)
/*
  Constructor for a 2D array with axis sizes (Ny,Nx) (the x dimension
  varies the quickest). The array is contiguous in memory.
  Values are either all doubles, or all std::complex.

  Parameters
  ----------
  Ny : const int
      The axis size in the y dimension (first index).
  Nx : const int
      The axis size in the x dimension (second index).
*/
{

  allocd = true;

  _sizeax = new int[2];
  _sizeax[0] = Ny;
  _sizeax[1] = Nx;

  _size = Ny*Nx;
  
  arr = new T[Nx*Ny];


};


template <typename T>
sMatrix<T>::~sMatrix() noexcept(false) {
  if (allocd) {
    delete [] _sizeax;
    delete [] arr;
  } else {
    throw std::bad_alloc();
  }

}


template <typename T>
void sMatrix<T>::size(const int Ny, const int Nx)
{

  if (allocd) {
    throw std::bad_alloc();
  }
  allocd = true;

  _sizeax = new int[2];
  _sizeax[0] = Ny;
  _sizeax[1] = Nx;

  _size = Ny*Nx;
  
  arr = new T[Nx*Ny];

  return ;
}

template <typename T>
int sMatrix<T>::size_including(int firstindex)
{
  return _sizeax[1]*firstindex;
}


template <typename T>
T& sMatrix<T>::operator()(int i, int j)
/* read/write access array elements via (i,j,k). */
{

  return arr[i*_sizeax[1] + j ];
}


template <typename T>
T sMatrix<T>::operator()(int i,int j) const
/* read-only access (but don't change) array elements via (i,j,k). */
{
  return arr[i*_sizeax[1] + j ];

}




template <typename T>
void sMatrix<T>::copy(sMatrix &input)
/* explicitly copy elements of input to arr. */
{
  
  for (int i = 0; i < 2; i++) {
    
    if (this->axis_size(i) != input.axis_size(i))
      throw std::runtime_error("Cannot copy sMatrix array (wrong shape).");
  }
  for (int i = 0; i < _sizeax[0]; i++) {
    for (int j = 0; j < _sizeax[1]; j++) {
      arr[i*_sizeax[1] + j ] = input(i,j);
    }
  }
  return;
};



template <typename T>
std::ostream& operator<< (std::ostream& out,
			  const sMatrix<T> &in)
{
  std::string delim = ", ";
  std::string innerdelim = ", ";
  std::vector<std::string> outer(2), middle(2), inner(2);
  
  outer[0]= inner[0] = "[ ";
  outer[1]= inner[1] = "] ";


  return in._write_ostream(out,outer,inner,
			   delim,innerdelim);
}



template <typename T>
std::ostream& sMatrix<T>::numpy_save(std::ostream &out)
{
  
  std::string delim = "\t";
  std::string innerdelim = "\n";

  std::vector<std::string> outer(2), inner(2);
  outer[0] = "";
  inner[0] = "#start of new first index\n";
  
  outer[1] = inner[1] = "";


  return _write_ostream(out,outer,inner,
			delim,innerdelim);

  
};





template <typename T>
std::ostream& sMatrix<T>::_write_ostream(std::ostream &out,
					 std::vector<std::string> outer,
					 std::vector<std::string> inner,
					 std::string delim,
					 std::string innerdelim) const
{
  
  int L0 = _sizeax[0];
  int M0 = _sizeax[1];
  
  out << outer[0];
  for (int i = 0; i < L0; i++) {
    out << inner[0];
    for (int j = 0; j < M0-1; j++) {
      out << arr[i*M0 + j ] << delim;
    }
    out << arr[i*M0 + M0-1] <<  inner[1]
	<< innerdelim;
  }

  out << outer[1];

  
  return out;
  
}


template class sMatrix<double>;
template class sMatrix<std::complex<double>>;
//template std::ostream& operator<< (std::ostream& ,
//				   const sMatrix<double>&);
//template std::ostream& operator<< (std::ostream& ,
//				   const sMatrix<std::complex<double>>&);


