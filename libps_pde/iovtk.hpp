#ifndef PSPDE_IOVTK_HPP
#define PSPDE_IOVTK_HPP

#include <string>
#include <vector>



#include "fftw_mpi_3darray.hpp"
#include "griddata.hpp"
//#include "globalinfompi.hpp"


namespace psPDE {
namespace ioVTK {
  typedef fftw_MPI_3Darray<double> ft_dub;
  //  template <typename T>
  //  void writeVTKImageData(std::string, const std::vector<T*>&, GridData);

  //  template <ptrdiff_t NZ, ptrdiff_t NY, ptrdiff_t NX>
  void writeVTKImageData(std::string,const std::vector<ft_dub*>,
			 const GridData&);
  void readVTKImageData(std::vector<ft_dub*>, std::string);
  
  void writeVTKcollectionHeader(const std::string);
  void writeVTKcollectionMiddle(const std::string,
				const std::string, const double);
  void writeVTKcollectionFooter(const std::string);
  void restartVTKcollection(const std::string, const MPI_Comm );
  

  //  void writeVTK_P_ImageData(std::string,
  //			    const std::vector<ft_dub*>& ,
  //			    GridData, GlobalInfoMPI);

  
}
};


#endif
