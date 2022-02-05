#ifndef INPUT_HPP
#define INPUT_HPP

#include <string>
#include <vector>
#include <map>
namespace input {
  std::vector<std::string> split_line(std::string&);
  bool isInt(std::string&,int&,std::string);
  bool isDouble(std::string&,double&,std::string);
  void convertVariable(std::string &,
		       std::map<std::string, std::string> const&);
  
}
#endif
