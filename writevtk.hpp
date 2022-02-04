#ifndef WRITEVTK_HPP
#define WRITEVTK_HPP

#include <string>
#include <vector>



#include "fftw_mpi_3darray.hpp"
#include "griddata.hpp"
//#include "globalinfompi.hpp"



namespace writeVTK {
  typedef fftw_MPI_3Darray<double> ft_dub;
  //  template <typename T>
  //  void writeVTKImageData(std::string, const std::vector<T*>&, GridData);

  //  template <ptrdiff_t NZ, ptrdiff_t NY, ptrdiff_t NX>
  void writeVTKImageData(std::string,const std::vector<ft_dub*>&,
			 GridData);
  void readVTKImageData(std::vector<ft_dub*>&, std::string , GridData);
  
  //  void writeVTK_P_ImageData(std::string,
  //			    const std::vector<ft_dub*>& ,
  //			    GridData, GlobalInfoMPI);
}



#endif
