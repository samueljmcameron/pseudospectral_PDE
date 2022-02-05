#include <iostream>


#include "input.hpp"

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
	  
	  for (const auto &[key,value]: varMap) {
	    if (key == tmp) {
	      found_key = true;
	      raw.erase(vstart,vend - vstart+1);
	      raw.insert(vstart,value);
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
}
