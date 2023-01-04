#ifndef PSPDE_READ_DATA_HPP
#define PSPDE_READ_DATA_HPP

#include <string>
#include <fstream>
#include <memory>
#include <mpi.h>

#include "atom.hpp"

namespace psPDE {

class ReadData {
public:

  static inline int SUCCESS = 0;
  static inline int FORMAT_ERROR = 1;
  ReadData(int,int,const std::string &);

  int read_file(std::unique_ptr<Atom> &,const MPI_Comm &);


private:

  std::ifstream datafile;

  int id, mpi_size;

  bool split_input;

  std::string line;

  std::vector<std::string> atom_properties;

  int nproperties, natoms_this_processor, natoms;

  int errflag;
  
  std::string format_fname(std::string);

  
  template <bool MULTIPROC>
  int read_file_template(std::unique_ptr<Atom> &);

  template <bool MULTIPROC>
  int read_properties_template(Atom *);
  int get_properties(const std::vector<std::string> &,Atom *,int);

};


}

#endif
