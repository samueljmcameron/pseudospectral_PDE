#include <iostream>

#include "atom.hpp"
#include "utility.hpp"



using namespace psPDE;


int Atom::unpack_reverse(int n,const std::vector<int> & sendlist, double *buf)
{

  int m = 0;
  int j;
  for (int i = 0; i < n; i++) {
    j = sendlist[i];
    Fs(0,j) += buf[m++];
    Fs(1,j) += buf[m++];
    Fs(2,j) += buf[m++];
  }

  return m;
  
}

int Atom::pack_comm(int n,const std::vector<int> & sendlist, double *buf,
		    const std::array<double,3> & pbc)
{

  int m = 0;
  int j;
  for (int i = 0; i < n; i++) {
    j = sendlist[i];
    buf[m++] = xs(0,j) + pbc[0];
    buf[m++] = xs(1,j) + pbc[1];
    buf[m++] = xs(2,j) + pbc[2];
  }

  return m;
  
}

int Atom::pack_comm_in_z(int n,const std::vector<int> & sendlist, double * buf,
			 const std::vector<double> &pbcz )
{
 int m = 0;
  int j;
  for (int i = 0; i < n; i++) {
    j = sendlist[i];
    buf[m++] = xs(0,j);
    buf[m++] = xs(1,j);
    buf[m++] = xs(2,j) + pbcz[i];
  }
  return m;

}

int Atom::pack_border_in_z(int nsend, const std::vector<int> & sendlist,double *buf,
			   const std::vector<double> & pbcz,
			   const std::vector<int> & labelz)
{

  int j;

  int m = 0;

  for (int i = 0; i < nsend; i++) {
    j = sendlist[i];
    buf[m++] = xs(0,j);
    buf[m++] = xs(1,j);
    buf[m++] = xs(2,j) + pbcz[i];
    buf[m++] = ubuf(tags[j]).d;
    buf[m++] = ubuf(types[j]).d;
    buf[m++] = ubuf(labelz[i]).d;
    if (solvent_interactions) {
      buf[m++] = phi[j];
      buf[m++] = nucwidth[j];
      buf[m++] = epsilonstrength[j];
    }
  }

  return m;
}


int Atom::pack_border(int nsend, const std::vector<int> & sendlist,double *buf,
		      const std::array<double,3> & pbc)
{

  int j;

  int m = 0;

  for (int i = 0; i < nsend; i++) {
    j = sendlist[i];
    buf[m++] = xs(0,j) + pbc[0];
    buf[m++] = xs(1,j) + pbc[1];
    buf[m++] = xs(2,j) + pbc[2];
    buf[m++] = ubuf(tags[j]).d;
    buf[m++] = ubuf(types[j]).d;
    buf[m++] = ubuf(Atom::GHOST).d;
    if (solvent_interactions) {
      buf[m++] = phi[j];
      buf[m++] = nucwidth[j];
      buf[m++] = epsilonstrength[j];
    }
    
  }

  return m;
}



void Atom::unpack_border(int nrecv,int first, const double *buf)
{
  int m = 0;
  int last = first + nrecv;

  xs.conservativeResize(Eigen::NoChange,last);
  Fs.conservativeResize(Eigen::NoChange,last);
  tags.resize(last);
  types.resize(last);
  labels.resize(last);
  if (solvent_interactions) {
    phi.resize(last);
    nucwidth.resize(last);
    epsilonstrength.resize(last);
  }
  
  
  for (int i = first; i < last; i++) {

    xs(0,i) = buf[m++];
    xs(1,i) = buf[m++];
    xs(2,i) = buf[m++];
    tags[i] = (int) ubuf(buf[m++]).i;
    types[i] = (int) ubuf(buf[m++]).i;
    labels[i] = (int) ubuf(buf[m++]).i;
    if (solvent_interactions) {
      phi[i] = buf[m++];
      nucwidth[i] = buf[m++];
      epsilonstrength[i] = buf[m++];
    }
    
  }

  
  return;
}


