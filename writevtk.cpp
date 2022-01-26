#include <stdexcept>
#include <fstream>
#include <complex>

#include "writevtk.hpp"




/* Simple function to write a vtk file with binary data. */
void writeVTK::writeVTKImageData(std::string fname,
				 const std::vector<ft_dub*>& scalar_outputs,
				 GridData grid)
/*============================================================================*/
/*
  Write scalar image data to a vtk (and paraview) compatible file
  (extension .vti).

  Parameters
  ----------

  fname : string
      Name of file to save with extension (either ".vti" or ".pvti").

  scalar_outputs : vector of pointers
      Each component of the vector points to a unique 3D array of data
      of type fftw_MPI_3Darray<double>.

  grid : struct of grid data
      Contains the following attributes: Ox, Oy, Oz for the
      start points of the grid (i.e. if you want x component to start
      at -4, then Ox=-4), and dx, dy, dz for the grid spacings.
*/
/*============================================================================*/

{

  // determine byte length of input data
  const int zstart = scalar_outputs[0]->get_local0start();

  
  const int Nz0 = scalar_outputs[0]->axis_size(0);
  const int Ny0 = scalar_outputs[0]->axis_size(1);
  const int Nx0 = scalar_outputs[0]->axis_size(2);

  
  unsigned int bytelength = Nx0*Ny0*Nz0*sizeof(double);


  
  for (const auto &elem : scalar_outputs)
    if ( elem->axis_size(0) != Nz0 || elem->axis_size(1) != Ny0
	 || elem->axis_size(2) != Nx0)
      throw std::runtime_error("All scalars must be on same grid.");
  
  
  auto myfile = std::fstream(fname, std::ios::out | std::ios::binary);
  myfile << "<?xml version=\"1.0\"?>" << std::endl
	 << "<VTKFile type=\"ImageData\"  version=\"1.0\""
	 << " byte_order=\"LittleEndian\">" << std::endl
	 << "<ImageData WholeExtent=\"0 " << Nx0-1 << " 0 "
	 << Ny0-1 << " " << zstart <<  "  " << Nz0+zstart-1 << "\" "
	 <<"Origin=\"" << grid.Ox << " " << grid.Oy << " " << grid.Oz << " \" "
	 << "Spacing=\"" << grid.dx << " " << grid.dy << " " << grid.dz << " "
	 << "\">" << std::endl
	 << "<Piece Extent=\"0 " << Nx0-1 << " 0 "
	 << Ny0-1 << " " << zstart << " " << Nz0+zstart-1 << "\">" << std::endl
	 << "<PointData Scalars=\"scalars\">" << std::endl;

  for (int i = 0; i < scalar_outputs.size(); i++) {
    
    int offset = i*(Nx0*Ny0*Nz0+sizeof(bytelength));
    
    myfile << "<DataArray Name=\"" << scalar_outputs[i]->get_name()
	   << "\" type=\"Float64\" format=\"appended\" "
	   << "offset=\"" << offset << "\"/>" << std::endl;
    
  }

  myfile << "</PointData>" << std::endl
	 << "</Piece>" << std::endl
	 << "</ImageData>" << std::endl
	 << "<AppendedData encoding=\"raw\">" << std::endl << "_";


  for (const auto &elem : scalar_outputs)    {
    //  for (int isc = 0; isc < scalar_outputs.size(); isc++) {
    myfile.write((char*)&bytelength,sizeof(bytelength));
    // since real fftw arrays aren't contiguous, need to write each row separately.
    for (int i = 0; i < Nz0; i++) {
      for (int j = 0; j < Ny0; j++) {
	myfile.write((char*)&(*elem)(i,j,0),sizeof(double)*Nx0);
	//
      }
    }
  }

  myfile << std::endl << "</AppendedData>" << std::endl
	 << "</VTKFile>" << std::endl;    
  myfile.close();
  
}

/*
void writeVTK::writeVTK_P_ImageData(std::string fname,
				 const std::vector<ft_dub*>& scalar_outputs,
				 GridData grid, GlobalInfoMPI info)
{

  // determine byte length of input data
  const int zstart = scalar_outputs[0]->get_local0start();


  const int globalNz = scalar_outputs[0]->get_globalNz();
  const int Nz0 = scalar_outputs[0]->axis_size(0);
  const int Ny0 = scalar_outputs[0]->axis_size(1);
  const int Nx0 = scalar_outputs[0]->axis_size(2);
  
  unsigned int bytelength = Nx0*Ny0*Nz0*sizeof(double);


  
  for (const auto &elem : scalar_outputs)
    if ( elem->axis_size(0) != Nz0 || elem->axis_size(1) != Ny0
	 || elem->axis_size(2) != Nx0)
      throw std::runtime_error("All scalars must be on same grid.");
  
  
  auto myfile = std::fstream(fname, std::ios::out | std::ios::binary);
  myfile << "<?xml version=\"1.0\"?>" << std::endl
	 << "<VTKFile type=\"PImageData\"  version=\"1.0\""
	 << " byte_order=\"LittleEndian\">" << std::endl
	 << "<PImageData WholeExtent=\"0 " << Nx0-1 << " 0 "
	 << Ny0-1 << " 0 "  << globalNz-1 << "\" "
	 <<"Origin=\"" << grid.Ox << " " << grid.Oy << " " << grid.Oz << " \" "
	 << "Spacing=\"" << grid.dx << " " << grid.dy << " " << grid.dz << "\" "
	 << "GhostLevel=\"0\">" << std::endl
	 << "<PPointData Scalars=\"scalars\">" << std::endl;

  for (int i = 0; i < scalar_outputs.size(); i++) {
    
    myfile << "<PDataArray Name=\"" << scalar_outputs[i]->get_name()
	   << "\" type=\"Float64\" />" << std::endl;
    
  }

  myfile << "</PPointData>" << std::endl;


  for (int id = 0; id < info.total_processors; id ++ ) {
  
    myfile << "<Piece Extent=\"0 " << Nx0-1 << " 0 "
	   << Ny0-1 << " " << info.zstarts[id] << " "
	   << info.Nz0s[id]+info.zstarts[id]-1 
	   << "\" Source=\"" << info.fnames[id] << "\"/>" << std::endl;
  }
  myfile << "</PImageData>" << std::endl
	 << "</VTKFile>" << std::endl;    
  myfile.close();
  
}

*/
