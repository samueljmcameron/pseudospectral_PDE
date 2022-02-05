#ifndef GRIDDATA_HPP
#define GRIDDATA_HPP

class GridData {
  
public:
  const int Nz, Ny, Nx;
  const double dz, dy, dx;
  const double Oz, Oy, Ox;


  GridData(const int Nz, const int Ny, const int Nx,
	   const double Lz, const double Ly, const double Lx)
    : Nz(Nz), Ny(Ny), Nx(Nx), dz(Lz/Nz), dy(Ly/Ny), dx(Lx/Nx),
      Oz(-Lz/2),Oy(-Ly/2),Ox(-Lx/2) {};

};

#endif
