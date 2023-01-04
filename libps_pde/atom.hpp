#ifndef PSPDE_ATOM_HPP
#define PSPDE_ATOM_HPP

#include <Eigen/Core>
#include <vector>
#include <array>

namespace psPDE {
class Atom {
public:
  static inline int OWNED = 0;
  static inline int LOCAL = 1;
  static inline int GHOST = 2;
  Atom(int Natoms, bool sol_int=false)
    : size_forward(3), size_reverse(3),size_border(6),
      solvent_interactions(sol_int) {

    xs.resize(3,Natoms); 
    uxs.resize(3,Natoms); 
    Fs.resize(3,Natoms);
    tags.resize(Natoms);
    types.resize(Natoms);
    labels.resize(Natoms);
    images.resize(Natoms);

    if (solvent_interactions) {
      phi.resize(Natoms);
      nucwidth.resize(Natoms);
      epsilonstrength.resize(Natoms);
      size_border += 3;
    }

    nowned = Natoms;
    ngathered = 0;
    nlocal = 0;
    nghost = 0;

    
  };
  ~Atom() {};
  
  int nowned; // atoms which are owned by this processor (as they are in a polymer)
  int nlocal; // atoms which are local to this processor
  int nghost; // atoms which are not local, but within cutoff of the processor
  int ngathered; // sum of local and ghost atoms (but not owned)


  // atom positions within periodic box and outside of periodic box, respectively
  Eigen::Matrix3Xd xs, uxs,Fs;


  // tag is atom id (unique to every atom)
  // type is a specific atom type (could be e.g. different polymers)
  // images is a flag for converting to/from PBC
  // label indicates if the atom is owned by the processor (and so time-stepped on this processor),
  //  local to this processor (so if partitioned into physical space via domain, is the atom in
  //  the domain), or a ghost atom (close enough in physical space to the processor).
  std::vector<int> tags,types,images,labels;


  // these are set if 
  std::vector<double> phi,nucwidth,epsilonstrength;



  int unpack_reverse(int, const std::vector<int> &, double *);

  int pack_comm(int,const std::vector<int> &, double *,
		const std::array<double,3> &);

  int pack_comm_in_z(int,const std::vector<int> &, double *,
		     const std::vector<double> &);
  
  int pack_border(int,const std::vector<int> &, double *,
		  const std::array<double,3> &);

  
  int pack_border_in_z(int,const std::vector<int> &, double *,
		       const std::vector<double> &,
		       const std::vector<int> &);

  
  void unpack_border(int, int , const double *);

  /*
  int pack_comm(int,const std::vector<int> &, double *,
		const std::array<double,3> &);

  void unpack_comm(int, int , const double *);
  */



  bool solvent_interactions;

  int get_size_fwd() const { return size_forward;};
  int get_size_rev() const { return size_reverse;};
  int get_size_bdr() const { return size_border;};
  

  
  
private:

  // size of packing in forward comm, reverse comm, and border comm
  int size_forward, size_reverse, size_border;

};
}

#endif
