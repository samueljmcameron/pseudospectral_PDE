
#include "domain.hpp"

#define IMGMASK 1023
#define IMGMAX 512
#define IMGBITS 10
#define IMG2BITS 20


using namespace psPDE;

Domain::Domain(double _xprd, double _yprd,double _zprd,
	       double _xlo, double _ylo, double _zlo,
	       int id,int mpi_size)
  : period{_xprd,_yprd,_zprd},boxlo{_xlo,_ylo,_zlo},
    boxhi{_xlo+_xprd,_ylo+_yprd,_zlo+_zprd},me(id),nprocs(mpi_size)
{

  // initially just set subdomains to be same as domains

  sublo = boxlo;
  subhi = boxhi;
  
};


void Domain::pbc(Atom & atoms) const
{

  int idim,otherdims;
  for (int i = 0; i < atoms.nowned; i++) {
    if (atoms.xs(0,i) < boxlo[0]) {
      atoms.xs(0,i) += period[0];
      idim = atoms.images[i] & IMGMASK;
      otherdims = atoms.images[i] ^ idim;
      idim --;
      idim &= IMGMASK;
      atoms.images[i] = otherdims | idim;
    }
    if (atoms.xs(0,i) >= boxhi[0]) {
      atoms.xs(0,i) -= period[0];
      atoms.xs(0,i) = (atoms.xs(0,i) > boxlo[0] ? atoms.xs(0,i) : boxlo[0]);
      idim = atoms.images[i] & IMGMASK;
      otherdims = atoms.images[i] ^ idim;
      idim ++;
      idim &= IMGMASK;
      atoms.images[i] = otherdims | idim;
    }
    if (atoms.xs(1,i) < boxlo[1]) {
      atoms.xs(1,i) += period[1];
      idim = (atoms.images[i] >> IMGBITS) & IMGMASK;
      otherdims = atoms.images[i] ^ (idim << IMGBITS);
      idim--;
      idim &= IMGMASK;
      atoms.images[i] = otherdims | (idim << IMGBITS);
    }
    if (atoms.xs(1,i) >= boxhi[1]) {
      atoms.xs(1,i) -= period[1];
      atoms.xs(1,i) = (atoms.xs(1,i) > boxlo[1] ? atoms.xs(1,i) : boxlo[1]);
      idim = (atoms.images[i]  >> IMGBITS) & IMGMASK;
      otherdims = atoms.images[i] ^ (idim << IMGBITS);
      idim ++;
      idim &= IMGMASK;
      atoms.images[i] = otherdims | (idim << IMGBITS);
    }
    if (atoms.xs(2,i) < boxlo[2]) {
      atoms.xs(2,i) += period[2];
      idim = atoms.images[i] >> IMG2BITS;
      otherdims = atoms.images[i] ^ (idim << IMG2BITS);
      idim--;
      idim &= IMGMASK;
      atoms.images[i] = otherdims | (idim << IMG2BITS);
    }
    if (atoms.xs(2,i) >= boxhi[2]) {
      atoms.xs(2,i) -= period[2];
      atoms.xs(2,i) = (atoms.xs(2,i) > boxlo[2] ? atoms.xs(2,i) : boxlo[2]);
      idim = atoms.images[i]  >> IMG2BITS;
      otherdims = atoms.images[i] ^ (idim << IMG2BITS);
      idim ++;
      idim &= IMGMASK;
      atoms.images[i] = otherdims | (idim << IMG2BITS);
    }
    // store atom in unwrapped coords (needed for polymer updates).
    unmap(atoms.uxs.col(i),atoms.xs.col(i),atoms.images[i]);
  }

  return;
}


void Domain::unmap(Eigen::Ref<Eigen::Vector3d> unwrapped,
		   const Eigen::Ref<const Eigen::Vector3d> & wrapped,
		   int image) const
{
  int xbox = (image & IMGMASK) - IMGMAX;
  int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
  int zbox = (image >> IMG2BITS) - IMGMAX;

  unwrapped(0) = wrapped(0) + xbox*period[0];
  unwrapped(1) = wrapped(1) + ybox*period[1];
  unwrapped(2) = wrapped(2) + zbox*period[2];

}


void Domain::map(Eigen::Ref<Eigen::Vector3d> wrapped,
		 const Eigen::Ref<const Eigen::Vector3d> & unwrapped,
		 int image) const
{
  int xbox = (image & IMGMASK) - IMGMAX;
  int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
  int zbox = (image >> IMG2BITS) - IMGMAX;

  wrapped(0) = unwrapped(0) - xbox*period[0];
  wrapped(1) = unwrapped(1) - ybox*period[1];
  wrapped(2) = unwrapped(2) - zbox*period[2];

}


int Domain::set_image() const
{
  int i1 = IMGMAX << IMG2BITS;
  int i2 = IMGMAX << IMGBITS;
  int i3 = IMGMAX;

  return i1 | i2 | i3;

}
