#ifndef GRIDDATA_HPP
#define GRIDDATA_HPP

class GridData {
  
  int Nz, Ny, Nx;
  double Lz, Ly, Lx;
  double Oz,Oy,Ox;
  double dz,dy,dx;
  
public:
  
  
  
  
  GridData(const int , const int , const int ,
	   const double , const double , const double );

  GridData(const int , const int , const int ,
	   const double , const double , const double,
	   const double, const double, const double);


  GridData fft_grid() const;

  void transpose_yz();


  int get_Nx() const
  {
    return Nx;
  };

  int get_Ny() const
  {
    return Ny;
  };

  int get_Nz() const
  {
    return Nz;
  };

  double get_Ox() const
  {
    return Ox;
  };
  double get_Oy() const
  {
    return Oy;
  };
  double get_Oz() const 
  {
    return Oz;
  };


  double get_dx() const
  {
    return dx;
  };
  double get_dy() const
  {
    return dy;
  };
  double get_dz() const
  {
    return dz;
  };


  double get_Lx() const
  {
    return Lx;
  };
  double get_Ly() const
  {
    return Ly;
  };
  double get_Lz() const
  {
    return Lz;
  };

  GridData get_positiveNx_grid() const
  {
    return GridData(Nz,Ny,Nx/2+1,Lz,Ly,Lx/2+dx,Oz,Oy,0);
  };
  

  
};

#endif
