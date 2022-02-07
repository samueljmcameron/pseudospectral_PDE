

#include "griddata.hpp"
#include <cmath>

GridData::GridData(const int Nz, const int Ny, const int Nx,
		   const double Lz, const double Ly, const double Lx)
  : Nz(Nz), Ny(Ny), Nx(Nx), Lz(Lz), Ly(Ly), Lx(Lx)
{
  update_all();
  
}

GridData::GridData(const int N, const double L)
  : Nz(N), Ny(N), Nx(N), Lz(L), Ly(L), Lx(L)
{
  update_all();
}
  
GridData::GridData(const GridData& grid)
  : Nz(grid.Nz), Ny(grid.Ny), Nx(grid.Nx), Lz(grid.Lz), Ly(grid.Ly), Lx(grid.Lx)
{
  update_all();
}
  
void GridData::update_origins()
{
  Oz = -Lz/2;
  Oy = -Ly/2;
  Ox = -Lx/2;
  
}

void GridData::update_spacings()
{
  dz = Lz/Nz;
  dy = Ly/Ny;
  dx = Lx/Nx;
  
}

void GridData::update_all()
{
  update_origins();
  update_spacings();
  
}



GridData GridData::fft_grid()
{
  
  const double fftLz = 2*M_PI/Lz*Nz;
  const double fftLy = 2*M_PI/Ly*Ny;
  const double fftLx = 2*M_PI/Lz*Nx;
  
  return GridData(Nz,Ny,Nx,fftLz,fftLy,fftLx);
  
}

void GridData::transpose_yz()
{
  int Ntmp = Nz;
  Nz = Ny;
  Ny = Ntmp;
  
  double Ltmp = Lz;
  Lz = Ly;
  Ly = Ltmp;
  
  update_all();
  return;
  
}
  
  
void GridData::reset(const int tNz, const int tNy, const int tNx,
		     const double tLz, const double tLy, const double tLx)
{
  Nz = tNz;
  Ny = tNy;
  Nx = tNx;
  Lz = tLz;
  Ly = tLy;
  Lx = tLx;
  
  update_all();
  
  return;
}

