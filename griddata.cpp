

#include "griddata.hpp"
#include <cmath>

GridData::GridData(const int Nz, const int Ny, const int Nx,
		   const double Lz, const double Ly, const double Lx)
  : Nz(Nz), Ny(Ny), Nx(Nx), Lz(Lz), Ly(Ly), Lx(Lx)
{

  dz = Lz/Nz;
  dy = Ly/Ny;
  dx = Lx/Nx;
  Oz = -Lz/2;
  Oy = -Ly/2;
  Ox = -Lx/2;
  
}

GridData::GridData(const int Nz, const int Ny, const int Nx,
		   const double Lz, const double Ly, const double Lx,
		   const double Oz,const double Oy,const double Ox)
  : Nz(Nz), Ny(Ny), Nx(Nx), Lz(Lz), Ly(Ly), Lx(Lx),Oz(Oz),Oy(Oy),Ox(Ox)
{

  dz = Lz/Nz;
  dy = Ly/Ny;
  dx = Lx/Nx;
  
}
  
GridData GridData::fft_grid() const
{
  
  const double fftLz = 2*M_PI/Lz*Nz;
  const double fftLy = 2*M_PI/Ly*Ny;
  const double fftLx = 2*M_PI/Lz*Nx;
  
  return GridData(Nz,Ny,Nx,fftLz,fftLy,fftLx,0.0,0.0,0.0);
  
}

void GridData::transpose_yz()
{
  int Ntmp = Nz;
  Nz = Ny;
  Ny = Ntmp;
  
  double tmp = Lz;
  Lz = Ly;
  Ly = tmp;

  tmp = dz;
  dz = dy;
  dy = tmp;
  
  tmp = Oz;
  Oz = Oy;
  Oy = tmp;  
  
  return;
  
}

