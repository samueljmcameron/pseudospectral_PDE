#include <iostream>
#include <fstream>

#include "input.hpp"

namespace psPDE {

namespace input {
  std::vector<std::string> split_line(std::string& line)
  {
    
    std::vector<std::string> out;
    
    std::string tmp;
    
    
    decltype(line.size()) index;
    
    index = 0;
    while(index < line.size()) {
      
      for (; index != line.size()
	     && !isspace(line[index]); ++index)
	tmp += line[index];
      
      if (tmp != "") {
	out.push_back(tmp);
      }
      tmp = "";
      
      index += 1;
    }
    return out;
  }
  
  bool isInt(std::string &val,int &out,std::string varname)
  {
    std::string::size_type sz;
    try { 
      out = std::stoi(val,&sz);
    } catch (std::invalid_argument &inv) {
      throw std::runtime_error(std::string("Error in format of ") + varname
			       + std::string(", raising ") + inv.what()
			       + std::string(" exception."));
    }
    
    return (sz == val.size());
    
  }
  
  bool isDouble(std::string &val, double &out,
		std::string varname)
  {
    std::string::size_type sz;
    try { 
      out = std::stod(val,&sz);
    } catch (std::invalid_argument &inv) {
      throw std::runtime_error(std::string("Error in format of ") + varname
			       + std::string(", raising ") + inv.what()
			       + std::string(" exception."));
      
    }
    
    return (sz == val.size());
    
  }
  
  
  void convertVariable(std::string &raw,
		       std::map<std::string, std::string> const& varMap)
  {
    
    
    
    std::string::size_type vstart,vend;
    
    vend = 0;
    while (vend != raw.size()) {
      
      vstart = raw.find("${");
      
      if (vstart != std::string::npos) {
	
	vend = raw.find("}",vstart+2);
	
	if (vend != std::string::npos) {
	  std::string tmp = raw.substr(vstart+2,vend-vstart-2);
	  if (tmp == "")
	    throw std::runtime_error("No content ('${}') in input file.");
	  bool found_key = false;
	  
	  for (const auto &xm: varMap) {
	    if (xm.first == tmp) {
	      found_key = true;
	      raw.erase(vstart,vend - vstart+1);
	      raw.insert(vstart,xm.second);
	    }
	  }
	  if (!found_key) {
	    std::string errorMessage
	      = std::string("No such command line variable ") + tmp;
	    throw std::runtime_error(errorMessage);
	  }
	  
	  
	  
	} else {
	  throw::std::runtime_error("Missing '}' in input script.");
	}
	
      } else {
	vend = raw.size();
      }
    }
    
    return;
  }


  // read in nucleation sites from another file
  void put_in_vectors(std::vector<std::vector<double>> & X_is,
		      const std::vector<int> & nucs_to_keep,
		      std::string filename, MPI_Comm comm,int mpi_id,
		      double starttime)
  {
    
    std::string finalline;
    
    // first check if file exists
    int signal = 0;
    if (mpi_id == 0) {
      
      try {
	finalline = getMatchingLine(filename,starttime);
      }
      catch (const std::runtime_error & error) {
	signal = -1;
      }
      
    }
    
    MPI_Bcast(&signal,1,MPI_INT,0,comm);
    
    if (signal == -1) {
      throw std::runtime_error("File " + std::string(filename) + " does not exist");    
    }
    
    std::vector<std::string> splitvec = split_line(finalline);
    
    // a thermo file with N nucleation sites in it should have 2*N*3+2 entries per line,
    //  (the + 2 is due to time and free energy, the *3 is due to three components per
    //  vector, and the *2 is due to each nucleation site also having a free energy
    //  derivative vector).
    
    
    int maxsites;

    if (mpi_id == 0)
     maxsites = (splitvec.size()-2)/6;
    MPI_Bcast(&maxsites,1,MPI_INT,0,comm);
    
    if (nucs_to_keep.size() > maxsites)
      throw std::runtime_error("Nucleation sites to keep is inconsistent with number of "
			       + std::string("nucleation sites in") + filename);

    double Xx,Xy,Xz;
    
    // add nucleation sites to list of X_is
    for (const auto &nuc : nucs_to_keep) {
      
      if (nuc >= maxsites || nuc < 0)
	throw std::runtime_error("Invalid nucleation site.");
      

      
      if (mpi_id == 0) {
	
	isDouble(splitvec.at(1+nuc*3),Xx,"");
	isDouble(splitvec.at(1+nuc*3+1),Xy,"");
	isDouble(splitvec.at(1+nuc*3+2),Xz,"");
      }
      
      MPI_Bcast(&Xx,1,MPI_DOUBLE,0,comm);
      MPI_Bcast(&Xy,1,MPI_DOUBLE,0,comm);
      MPI_Bcast(&Xz,1,MPI_DOUBLE,0,comm);
      
      X_is.push_back({Xx,Xy,Xz});
    }
    
    return;
    
  }

  // read in nucleation sites from another file (and keep all of them)
  void put_in_vectors(std::vector<std::vector<double>> & X_is,
		      std::string filename, MPI_Comm comm,int mpi_id,
		      double starttime)
  {
    
    std::string finalline;
    
    // first check if file exists
    int signal = 0;
    if (mpi_id == 0) {
      
      try {
	finalline = getMatchingLine(filename,starttime);
      }
      catch (const std::runtime_error & error) {
	signal = -1;
      }
      
    }
    
    MPI_Bcast(&signal,1,MPI_INT,0,comm);
    
    if (signal == -1) {
      throw std::runtime_error("File " + std::string(filename) + " does not exist");    
    }


    
    // a thermo file with N nucleation sites in it should have 2*N*3+2 entries per line,
    //  (the + 2 is due to time and free energy, the *3 is due to three components per
    //  vector, and the *2 is due to each nucleation site also having a free energy
    //  derivative vector).
    
    int maxsites;
    std::vector<std::string> splitvec;
    if (mpi_id == 0) {
      splitvec = split_line(finalline);
      maxsites = (splitvec.size()-2)/6;
    }
    MPI_Bcast(&maxsites,1,MPI_INT,0,comm);
    
    double Xx,Xy,Xz;

    for (int nuc = 0; nuc < maxsites; nuc ++) {
      

      if (mpi_id == 0) {
	
	isDouble(splitvec.at(1+nuc*3),Xx,"");
	isDouble(splitvec.at(1+nuc*3+1),Xy,"");
	isDouble(splitvec.at(1+nuc*3+2),Xz,"");

      }
      
      MPI_Bcast(&Xx,1,MPI_DOUBLE,0,comm);
      MPI_Bcast(&Xy,1,MPI_DOUBLE,0,comm);
      MPI_Bcast(&Xz,1,MPI_DOUBLE,0,comm);
      
      X_is.push_back({Xx,Xy,Xz});
    }
    
    return;
    
  }
  // return the line whose first entry is closest (either equal to or less than) starttime
  std::string getMatchingLine(std::string filename,double starttime)
  {
    std::ifstream myfile (filename);
    std::string previousline = "";
    std::string line;
    std::vector<std::string> split;
    std::string::size_type sz;

    // choosing a very negative value to start with
    double t = -1e14;
    double tol = 1e-6;

    if (myfile) {
      
      while(std::getline(myfile,line)) {

	split = split_line(line);
	try {
	  t = std::stod(split[0],&sz);
	} catch (std::invalid_argument &inv) {
	  continue;
	}

	if (std::abs(t-starttime) < tol) {
	  previousline = line;
	  break;
	} else if (t > starttime) {
	  break;
	}
	previousline = line;

	
	
      }
      
    } else {
      throw std::runtime_error("File " + std::string(filename) + " does not exist");
    }


    if (previousline == "")
      throw std::runtime_error("File " + std::string(filename) + " does not have t <= "
			       + std::to_string(starttime));
    myfile.close();
    return previousline;
  }


    
  std::string getLastLine(std::string filename)
  {
    
    std::ifstream myfile (filename);
    std::string lastline;
    
    if (myfile) {
      
      // assumes that last line is just a newline char, hence the -2 below (vs -1)
      myfile.seekg(-2,myfile.end);
      
      char ch;
      
      myfile.get(ch);
      
      while (ch != '\n' && myfile.tellg() >1) {
	
	myfile.seekg(-2,myfile.cur);
	myfile.get(ch);
      }
      
      
      
      
      std::getline(myfile,lastline);
      
    } else {
      throw std::runtime_error("File " + std::string(filename) + " does not exist");
    }
    
    myfile.close();
    return lastline;
  }
  
}
}
