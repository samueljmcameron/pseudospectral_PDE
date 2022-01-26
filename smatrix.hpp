#ifndef SMATRIX_HPP
#define SMATRIX_HPP

#include <iostream>
#include <vector>
#include <string>



template <typename T>
class sMatrix
{
private:
  T *arr;


  bool allocd;
  int _size;
  int *_sizeax;
  std::string array_name;
  int spacer;

  std::ostream& _write_ostream(std::ostream &out, std::vector<std::string> outer,
			       std::vector<std::string> inner,std::string delim,
			       std::string innerdelim) const;
  
  
public:

  sMatrix();


  sMatrix(const int Nx, const int Ny);

  ~sMatrix();

  T* data() {
      return arr;
  };

  int size() {
    return _size;
  };

  void size(const int , const int);
  
  int axis_size(int i)
  {
    return _sizeax[i];
  }

  int size_including(int);
  
    
  T& operator()(int,int );
  T operator()(int,int) const;

  void copy(sMatrix &);

  std::ostream& numpy_save(std::ostream & );
  
  template <typename T0>
  friend std::ostream& operator<< (std::ostream& out,
				   const sMatrix<T0> & in);


};





#endif
