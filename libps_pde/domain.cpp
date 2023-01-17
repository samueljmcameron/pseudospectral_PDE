
#include <iostream>

#include "domain.hpp"


using namespace psPDE;

Domain::Domain(int id,int mpi_size,double _xprd, double _yprd,double _zprd,
	       double _xlo, double _ylo, double _zlo)
  : period{_xprd,_yprd,_zprd},boxlo{_xlo,_ylo,_zlo},
    boxhi{_xlo+_xprd,_ylo+_yprd,_zlo+_zprd},me(id),nprocs(mpi_size)
{

  // initially just set subdomains to be same as domains

  sublo = boxlo;
  subhi = boxhi;
  
};

Domain::Domain(int id, int mpi_size,const std::vector<std::string> &v_line)
  : me(id),nprocs(mpi_size)
{



  int iargs = 0;

  // default values

  period = {1.0,1.0,1.0};
  boxlo = {-0.5,-0.5,-0.5};

  bool period_flag = false;
  bool origin_flag = false;
  
  while (iargs < v_line.size()) {
    
    if (v_line[iargs] == "boxdims") {
      period_flag = true;
      iargs ++;
      period[0] = std::stod(v_line[iargs++]);
      period[1] = std::stod(v_line[iargs++]);
      period[2] = std::stod(v_line[iargs++]);
    } else if (v_line[iargs] == "boxorigin") {
      origin_flag = true;
      iargs ++;
      boxlo[0] = std::stod(v_line[iargs++]);
      boxlo[1] = std::stod(v_line[iargs++]);
      boxlo[2] = std::stod(v_line[iargs++]);
      
    } else {
      throw std::invalid_argument("Bad arguments for domain_setup.");
      
    }
  }

  if (!period_flag  && id == 0)
    std::cerr << "Warning, domain size not set, reverting to default values."
	      << std::endl;

  
  if (!origin_flag  && id == 0)
    std::cerr << "Warning, domain origin not set, reverting to default values."
	      << std::endl;

  
  for (int i = 0; i < 3; i++)  
    boxhi[i] = boxlo[i] + period[i];

  sublo = boxlo;
  subhi = boxhi;
  
};



