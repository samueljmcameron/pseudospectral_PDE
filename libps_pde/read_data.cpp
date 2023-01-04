#include <iostream>
#include <algorithm>


#include "read_data.hpp"
#include "input.hpp"

using namespace psPDE;

ReadData::ReadData(int id, int mpi_size, const std::string & fname)
  : id(id),mpi_size(mpi_size),split_input(false)
{

  std::string multiprocname = format_fname(fname);

  if (multiprocname != fname)
    split_input = true;
  
  datafile.open(multiprocname);

}


std::string ReadData::format_fname(std::string fname)
{

  const std::string label = "%";
  std::string newstr = std::to_string(id);
  auto tmppos = fname.find(label);
  if (tmppos != fname.npos) {
    fname.erase(tmppos,label.size()-1);
    fname.replace(tmppos,1,newstr);
  }

  return fname;

}


int ReadData::read_file(std::unique_ptr<Atom> &atoms,
			const MPI_Comm & comm)
{

  if (split_input)
    errflag = read_file_template<true>(atoms);
  else
    errflag =  read_file_template<false>(atoms);

  if (errflag != SUCCESS) errflag = 1;
  else errflag = 0;


  int total_errflag;
  
  MPI_Allreduce(&errflag,&total_errflag,1,MPI_INT,MPI_SUM,comm);

  if (total_errflag) return FORMAT_ERROR;
  else return SUCCESS;
  
}


template <bool MULTIPROC>
int ReadData::read_file_template(std::unique_ptr<Atom> &atoms)
/*
  For reading from a single data file, and splitting data onto
  multiple processors
*/
{
  std::vector<std::string> split_vec;

  if (datafile) {


    // first find the first non-empty, non-commented line of the file,
    //  which must start with "natoms", then instantiate the Atom
    //  unique_ptr with the number of atoms per processor
    while (std::getline(datafile,line)) {
      
      if (line == "" || line[0] == '#') continue;
      
      split_vec = input::split_line(line);

      
      split_vec.erase(std::find(split_vec.begin(),split_vec.end(),"#"),
		      split_vec.end());

      
      if (split_vec.size() != 2 || split_vec[0] != "natoms")
	return FORMAT_ERROR;

      try {
	input::isInt(split_vec[1],natoms,split_vec[0]);
      } catch (std::runtime_error &err) {
	return FORMAT_ERROR;
      }

      if (!MULTIPROC) {
	natoms_this_processor = (natoms/mpi_size );

	if (id < natoms % mpi_size )
	  natoms_this_processor += 1;

	
      } else {
	
	natoms_this_processor = natoms;

      }

      atoms = std::make_unique<Atom>(natoms_this_processor);

      
      break;
      
    }


    // next, read in the line which says what properties these
    //  atoms will have
    while (std::getline(datafile,line)) {
      if (line == "" || line[0] == '#') continue;

      atom_properties = input::split_line(line);

      atom_properties.erase(std::find(atom_properties.begin(),atom_properties.end(),"#"),
			    atom_properties.end());

      

      break;
    }


    // finally, fill in the atom properties
    return read_properties_template<MULTIPROC>(atoms.get());

    
  } else atoms = std::make_unique<Atom>(0);


  return SUCCESS;
}

template <bool MULTIPROC>
int ReadData::read_properties_template(Atom *atoms)
/*
  Expects datafile to have already read in the number of atoms.

  For reading from a single data file, and splitting data onto
  multiple processors

*/
{
  std::string line;
  int count = natoms; // set to natoms just for check at end

  if (!MULTIPROC) count = 0; 

  int iatom = 0;



  std::vector<std::string> split_vec;


  // skip all lines which are blank or uncommented
  while(std::getline(datafile,line)) {
    
    if (line == "" || line[0] == '#') continue;
    
    split_vec = input::split_line(line);
    split_vec.erase(std::find(split_vec.begin(),split_vec.end(),"#"),
		    split_vec.end());

    

    // the plus two is for the fact that positions take up three values
    if (split_vec.size() != atom_properties.size() + 2)
      return FORMAT_ERROR;

    for (auto dbl : split_vec)
      try {
	double tmp;
	input::isDouble(dbl,tmp,"atom property");
      } catch (std::runtime_error &err) {
	return FORMAT_ERROR;
    }

    // splitting atoms up by processor (instead of location) - this
    //  is a design choice based on the way atoms are forward/reverse
    //  communicated between processors.

    if (!MULTIPROC) {
      if (count % mpi_size == id) {
	int errflag = get_properties(split_vec,atoms,iatom);
	if (errflag != SUCCESS) return errflag;
	iatom += 1;
      }

      count += 1;
    } else {
      int errflag = get_properties(split_vec,atoms,iatom);
      if (errflag != SUCCESS) return errflag;
      
      iatom +=1;
    }
  }

  if (iatom != natoms_this_processor || count != natoms)
    return FORMAT_ERROR;


  
  return SUCCESS;
}


int ReadData::get_properties(const std::vector<std::string> & split_vec,
			     Atom *atoms,int iatom)
{
  int item = 0;
  if (iatom >= atoms->nowned) return FORMAT_ERROR;
  for (auto &property : atom_properties)
    if (property == "positions") {
      atoms->xs(0,iatom) = std::stod(split_vec[item++]);
      atoms->xs(1,iatom) = std::stod(split_vec[item++]);
      atoms->xs(2,iatom) = std::stod(split_vec[item++]);
    } else if (property == "nucmaxs") {
      atoms->nucmaxs(iatom) = std::stod(split_vec[item++]);
    } else if (property == "viscosities") {
      atoms->viscosities(iatom) = std::stod(split_vec[item++]);
    } else if (property == "radii") {
      atoms->radii(iatom) = std::stod(split_vec[item++]);
    } else {
      return FORMAT_ERROR;
    }
  return SUCCESS;
}
