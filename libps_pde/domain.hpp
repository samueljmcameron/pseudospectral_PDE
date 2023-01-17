#ifndef PSPDE_DOMAIN_HPP
#define PSPDE_DOMAIN_HPP

#include <vector>
#include <array>
#include <string>

namespace psPDE {
class Domain
{
public:
  Domain(int,int,double,double,double,double,double,double);
  Domain(int , int ,const std::vector<std::string> & );

  const int me,nprocs;
  std::array<double,3> period,boxlo,boxhi;
  std::array<double,3> sublo,subhi;
  
private:

};
}
#endif
