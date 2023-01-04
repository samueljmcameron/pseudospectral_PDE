#include <fstream>

#include "read_input.hpp"


ReadInput::ReadInput(const std::string & fname)
  : gridset(false),atomset(false),griddump(false),atomdump(false),
    
{

  datafile.open(fname)
  
}




int ReadInput::read_all(std::map<std::string,std::string> const& varMap)
{

  std::vector<std::string> split_vec;
  while(std::getline(datafile,line)) {
    
    if (line == "" || line[0] == '#') continue;

    split_vec = input::split_line(line);
    
      
    split_vec.erase(std::find(split_vec.begin(),split_vec.end(),"#"),
		    split_vec.end());

    for (auto &word : split_vec)
      input::convertVariable(word,varMap);

    if (split_vec[0] == "grid_setup") {
      if (gridset) return FORMAT_ERROR;
      set_grid(split_vec);
      gridset=true;
    } else if (split_vec[0] == "atom_setup") {
      if (atomset) return FORMAT_ERROR;
      set_atom(split_vec);
      atomset=true;
    } else if (split_vec[0] == "grid_dump") {
      set_grid_dump(split_vec);
    } else if (split_vec[0] == "atom_dump") {
      set_atom_dump(split_vec);
    } else if (split_vec[0] == "thermo") {
      set_thermo(split_vec);
    } else if (split_vec[0] == "dt") {
      set_dt(split_vec);
    } else if (split_vec[0] == "seed") {
      set_seed(split_vec);
    } else if (split_vec[0] == "chemical_potential") {
      set_chemical_potential(split_vec);
    } else if (split_vec[0] == "run") {
      set_run(split_vec);
    }
  }

}


void ReadInput::set_grid()
{
  

}
