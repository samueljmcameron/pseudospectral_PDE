#include <stdexcept>
#include <fstream>

#include <algorithm>

#include "iovtk.hpp"


using namespace psPDE;

/* Simple function to write a vtk file with binary data. */
void ioVTK::writeVTKImageData(std::string fname,
			      const std::vector<ft_dub*> scalar_outputs,
			      const std::array<double,3> &boxlo)
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
      Contains the following attributes: get_Ox(), get_Oy(), get_Oz() for the
      start points of the grid (i.e. if you want x component to start
      at -4, then Ox=-4), and get_dx(), get_dy(), get_dz() for the grid spacings.
*/
/*============================================================================*/

{

  // determine byte length of input data
  const int zstart = scalar_outputs[0]->get_local0start();

  
  const int Nz0 = scalar_outputs[0]->Nz();
  const int Ny0 = scalar_outputs[0]->Ny();
  const int Nx0 = scalar_outputs[0]->Nx();


  const double dz = scalar_outputs[0]->dx;
  const double dy = scalar_outputs[0]->dx;
  const double dx = scalar_outputs[0]->dx;

  
  unsigned int bytelength = Nx0*Ny0*Nz0*sizeof(double);
  
  for (const auto &elem : scalar_outputs)
    if ( elem->axis_size(0) != Nz0 || elem->axis_size(1) != Ny0
	 || elem->axis_size(2) != Nx0)
      throw std::runtime_error("All scalars must be on same grid.");
  
  
  auto myfile = std::fstream(fname, std::ios::out | std::ios::binary);

  if (not myfile.is_open()) {
    throw std::runtime_error(std::string("Cannot open file ") + fname);
  }

  
  myfile << "<?xml version=\"1.0\"?>" << std::endl
	 << "<VTKFile type=\"ImageData\"  version=\"1.0\""
	 << " byte_order=\"LittleEndian\">" << std::endl
	 << "<ImageData WholeExtent=\"0 " << Nx0-1 << " 0 "
	 << Ny0-1 << " " << zstart <<  "  " << Nz0+zstart-1 << "\" "
	 <<"Origin=\"" << boxlo[0] << " " << boxlo[1] << " " << boxlo[2] << " \" "
	 << "Spacing=\"" << dx << " " << dy << " " << dz << " "
	 << "\">" << std::endl
	 << "<Piece Extent=\"0 " << Nx0-1 << " 0 "
	 << Ny0-1 << " " << zstart << " " << Nz0+zstart-1 << "\">" << std::endl
	 << "<PointData Scalars=\"scalars\">" << std::endl;

  for (unsigned i = 0; i < scalar_outputs.size(); i++) {
    
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
    //  for (unsigned isc = 0; isc < scalar_outputs.size(); isc++) {
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

void ioVTK::writeVTKcollectionHeader(const std::string fname)
{
  auto myfile = std::ofstream(fname);
  if (not myfile.is_open()) {
    throw std::runtime_error(std::string("Cannot open file ") + fname);
  }
  
  myfile << "<?xml version=\"1.0\"?>" << std::endl
	 << "<VTKFile type=\"Collection\"  version=\"1.0\""
	 << " byte_order=\"LittleEndian\">" << std::endl
	 << "<Collection>" << std::endl;


  myfile.close();
  

}

void ioVTK::writeVTKcollectionMiddle(const std::string collectionfname,
				     const std::string filename,
				     const double time)
{

  //size_t num_slashes = std::count(collectionfname.begin(), collectionfname.end(), '/');

  size_t firstslash = filename.find_last_of("\\/");
  std::string file_no_path = filename;
  if (firstslash != std::string::npos) {
    file_no_path = filename.substr(firstslash+1);
  }
  
  auto myfile = std::fstream(collectionfname,std::ios_base::app);
  if (not myfile.is_open()) {
    throw std::runtime_error(std::string("Cannot open file ") + collectionfname);
  }
  
  myfile << "<DataSet timestep=\"" << time << "\" group=\"\" part=\"0\""
	 << " file=\"" << file_no_path << "\"/>" << std::endl;

  myfile.close();
  
}


void ioVTK::writeVTKcollectionFooter(const std::string fname)
{
  auto myfile = std::fstream(fname, std::ios_base::app);
  if (not myfile.is_open()) {
    throw std::runtime_error(std::string("Cannot open file ") + fname);
  }
  
  myfile << "</Collection>" << std::endl
	 << "</VTKFile>";
  
  myfile.close();
}

void ioVTK::restartVTKcollection(const std::string fname, const MPI_Comm comm)
/* restart the collection by copying the current file (up to the collection
   end) into a new file called "restart" + oldfname. */
   
{

  auto oldfile = std::ifstream(fname);
  if (not oldfile.is_open()) {
    throw std::runtime_error(std::string("Cannot open file ") + fname);
  }


  std::vector<std::string> lines;
  std::string stopline = "";
  std::getline(oldfile,stopline);
  lines.push_back(stopline);

  while (std::getline(oldfile,stopline) && stopline != "</Collection>") {
    lines.push_back(stopline);
  }
  oldfile.close();

  if (stopline != "</Collection>") {
    throw std::runtime_error("Invalid restart file." );
  }

  // if one process throws error, then don't want to erase data in other processes.
  MPI_Barrier(comm);

  auto newfile = std::ofstream(fname);
  if (not newfile.is_open()) {
    throw std::runtime_error(std::string("Cannot open file ") + fname);
  }

  for (const auto &line : lines)
    newfile << line << std::endl;;
  


  newfile.close();

  return;
}


/* Simple function to write a vtk file with binary data. */
void ioVTK::readVTKImageData(std::vector<ft_dub*> scalar_outputs,
			     std::string fname)
/*============================================================================*/
/*
  Write scalar image data to a vtk (and paraview) compatible file
  (extension .vti).

  Parameters
  ----------

  scalar_outputs : vector of pointers
      Each component of the vector points to a unique 3D array of data
      of type fftw_MPI_3Darray<double>.

  fname : string
      Name of file to read from with extension (either ".vti" or ".pvti").

  grid : struct of grid data
      Contains the following attributes: get_Ox(), get_Oy(), get_Oz() for the
      start points of the grid (i.e. if you want x component to start
      at -4, then get_Ox()=-4), and get_dx(), get_dy(), get_dz() for the grid spacings.
*/
/*============================================================================*/

{
  
  const int Nz0 = scalar_outputs[0]->Nz();
  const int Ny0 = scalar_outputs[0]->Ny();
  const int Nx0 = scalar_outputs[0]->Nx();

  
  unsigned int bytelength = Nx0*Ny0*Nz0*sizeof(double);

  std::string stopline = "";

  
  for (const auto &elem : scalar_outputs)
    if ( elem->axis_size(0) != Nz0 || elem->axis_size(1) != Ny0
	 || elem->axis_size(2) != Nx0)
      throw std::runtime_error("All scalars must be on same grid.");
  

  auto myfile = std::fstream(fname, std::ios::in | std::ios::binary);

  if (not myfile.is_open()) 
    throw std::runtime_error("Image data reading file doesn't exist."); 

  while (stopline != "<AppendedData encoding=\"raw\">") {
    std::getline(myfile,stopline);
  }

  char memblock[2];
  
  myfile.read(memblock,1);
  stopline = memblock[0];



  
  for (const auto &elem : scalar_outputs)    {

    myfile.read((char*)&bytelength,sizeof(bytelength));
    //    bytelength = *(unsigned int*) itmp;

    // since real fftw arrays aren't contiguous, need to write each row separately.
    for (int i = 0; i < Nz0; i++) {
      for (int j = 0; j < Ny0; j++) {
	myfile.read((char*)&(*elem)(i,j,0),sizeof(double)*Nx0);
	//
      }
    }
  }

  myfile.close();
  
}


/*
void ioVTK::writeVTK_P_ImageData(std::string fname,
				 const std::vector<ft_dub*> scalar_outputs,
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
	 <<"Origin=\"" << grid.get_Ox() << " " << grid.get_Oy() << " " << grid.get_Oz() << " \" "
	 << "Spacing=\"" << grid.get_dx() << " " << grid.get_dy() << " " << grid.get_dz() << "\" "
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
