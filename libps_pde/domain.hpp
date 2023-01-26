#ifndef PSPDE_DOMAIN_HPP
#define PSPDE_DOMAIN_HPP

#include <vector>
#include <array>
#include <string>

namespace psPDE {
class Grid;
class Domain
{
public:
  Domain(int,int,double,double,double,double,double,double);
  Domain(int , int ,const std::vector<std::string> & );

  const int me,nprocs;
  std::array<double,3> period,boxlo,boxhi;


  void partition(const Grid *);


  double dqx() { return 2*3.141592/period[0];};
  double dqy() { return 2*3.141592/period[2];};  
  double dqz() { return 2*3.141592/period[1];};

  
  
private:

};
}
#endif
