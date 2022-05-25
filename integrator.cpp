#include "integrator.hpp"




Integrator::Integrator(MPI_Comm comm,const GridData& fourier, const int seed,
		       const SolutionParams& solparams,const double dt)
  : ft_phi(comm,"ft_phi",fourier), ft_nonlinear(comm,"ft_nl",fourier),
    seed(seed),
    real_dist(-0.5,0.5),mobility(solparams.mobility),gamma(solparams.gamma),
    temp(solparams.temp),chi(solparams.chi), volFH(solparams.volFH),
    chi_LP(solparams.chi_LP),nucmax(solparams.nucmax),nucwidth(solparams.nucwidth),
    dt(dt),normalization(sqrt(1.0/(fourier.get_Nx()*fourier.get_Ny()*fourier.get_Nz())))
{

  gen.seed(seed);
  double tmp = fourier.get_dx()*fourier.get_dy()*fourier.get_dz()/(2*M_PI*2*M_PI*2*M_PI);
  local0start = ft_phi.get_local0start();
  complexprefactor = sqrt(12*temp*mobility*tmp);
  realprefactor = sqrt(24*temp*mobility*tmp);
  sqrtdt = sqrt(dt);

  //  std::cout << (fourier.get_dz()*fourier.get_dz()*fourier.get_Nz()/2*fourier.get_Nz()/2+fourier.get_dy()*fourier.get_dy()*fourier.get_Ny()/2*fourier.get_Ny()/2+fourier.get_dx()*fourier.get_dx()*fourier.get_Nx()/2*fourier.get_Nx()/2) << std::endl;
  std::cout << "RN : " << real_dist(gen) << std::endl;
  
}

Integrator::~Integrator()
{
}

void Integrator::initialize(fftw_MPI_3Darray<double> & phi,
			  const double initial_volfrac, const double variance)
{
  
  for (int i = 0; i < phi.axis_size(0); i++) {
    for (int j = 0; j < phi.axis_size(1); j++) {
      for (int k = 0; k < phi.axis_size(2); k++) {
	phi(i,j,k) = initial_volfrac + variance*real_dist(gen);
      }
    }
  }
  return;
}


std::vector<std::vector<double>> Integrator::nonlinear(fftw_MPI_3Darray<double>& nonlinear,
						       const fftw_MPI_3Darray<double>& phi,
						       const std::vector<std::vector<double>> & X_is,
						       double & free_energy)
{

  int real0start = phi.get_local0start();
  double dx = phi.grid.get_dx();
  double dy = phi.grid.get_dy();
  double dz = phi.grid.get_dz();
  double Ox = phi.grid.get_Ox();
  double Oy = phi.grid.get_Oy();
  double Oz = phi.grid.get_Oz();

  double Lx = phi.grid.get_Lx();
  double Ly = phi.grid.get_Ly();
  double Lz = phi.grid.get_Lz();

  double x,y,z;
  double philink;

  free_energy = 0;
  double allfree_energy;

  // sum will store the free energy derivative integral
  std::vector<std::vector<double>> sum ;
  // allsum is sum after MPI_Allreduce
  std::vector<std::vector<double>> allsum ;

  for (unsigned i = 0; i < X_is.size(); i++) {
    sum.push_back({0,0,0});
    allsum.push_back({0,0,0});
  }


  std::vector<double> dlink(3);
  
  for (int i = 0; i < nonlinear.axis_size(0); i++) {
    z = dz*(i+  local0start) + Oz;

    
    for (int j = 0; j < nonlinear.axis_size(1); j++) {
      
      y = dy*j + Oy;

      
      for (int k = 0; k < nonlinear.axis_size(2); k++) {
	x = dx*k + Ox;


	philink = linker_phi(x,y,z,Lx,Ly,Lz,X_is);
	nonlinear(i,j,k)
	  = temp/volFH*(log(phi(i,j,k)/(1-phi(i,j,k)))+chi*(1-2*phi(i,j,k))
			+philink*(phi(i,j,k)-nucmax));
	//	  = 2*temp/volFH*phi(i,j,k);
	//= -quad*phi(i,j,k)+quartic*phi(i,j,k)*phi(i,j,k)*phi(i,j,k);


	for (unsigned index = 0; index < sum.size() ; index ++) {
	  linker_derivative(dlink,x,y,z,Lx,Ly,Lz,X_is[index]);
	  for (int coord = 0; coord < 3; coord++)
	    sum[index][coord] += (dlink[coord]*(phi(i,j,k)-nucmax)
				  *(phi(i,j,k)-nucmax)*dx*dy*dz);
	}
	free_energy += philink*(phi(i,j,k)-nucmax)*(phi(i,j,k)-nucmax)*dx*dy*dz;

      }
    }
  }

  for (unsigned index = 0; index < sum.size() ; index ++) {

    for (int coord = 0; coord < 3; coord ++) {
      MPI_Allreduce(&(sum[index][coord]),&(allsum[index][coord]),
		    1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      //      std::cout << "dfdx = " << allsum[index][coord] << std::endl;
    }
  }
  MPI_Allreduce(&free_energy,&allfree_energy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  free_energy = allfree_energy;

  return allsum;

}

void Integrator::integrate(int i, int j, int k)
{

  if (i + local0start > ft_phi.grid.get_Nz()/2) {
    qz = ft_phi.grid.get_dz()*(ft_phi.grid.get_Nz()- i - local0start);
  } else {
    qz = ft_phi.grid.get_dz()*(i + local0start);
  }

  if (j > ft_phi.grid.get_Ny()/2) {
    qy = ft_phi.grid.get_dy()*(ft_phi.grid.get_Ny()-j);
  } else {
    qy = ft_phi.grid.get_dy()*j;
  }  

  qx = ft_phi.grid.get_dx()*k;
  q2 = qx*qx + qy*qy + qz*qz;
  noise.real(complexprefactor*sqrt(q2)*sqrtdt*real_dist(gen));
  noise.imag(complexprefactor*sqrt(q2)*sqrtdt*real_dist(gen));
  
  ode(ft_phi(i,j,k),ft_nonlinear(i,j,k),noise,q2);
  
  
  //  ft_phi(i,j,k) = 1.0*(i+local0start) + 1i*(1.0*j);
  return;
}


void Integrator::integrate_real(int i, int j, int k)
{

  if (i + local0start > ft_phi.grid.get_Nz()/2) {
    qz = ft_phi.grid.get_dz()*(ft_phi.grid.get_Nz()- i - local0start);
  } else {
    qz = ft_phi.grid.get_dz()*(i + local0start);
  }

  if (j > ft_phi.grid.get_Ny()/2) {
    qy = ft_phi.grid.get_dy()*(ft_phi.grid.get_Ny()-j);
  } else {
    qy = ft_phi.grid.get_dy()*j;
  }
  
  qx = ft_phi.grid.get_dx()*k;
  q2 = qx*qx + qy*qy + qz*qz;
  
  noise = realprefactor*sqrt(q2)*sqrtdt*real_dist(gen);

  
  ode(ft_phi(i,j,k),ft_nonlinear(i,j,k),noise,q2);
  
  //  ft_phi(i,j,k) = (1.0*(i+local0start) + 1i*(1.0*j))/2.0;
  return;
}

void Integrator::ode(std::complex<double> & y, std::complex<double> ynl,
		   std::complex<double> rnd, double q2)
{
  y = (y-mobility*q2*dt*(ynl+gamma*q2*y))*normalization*normalization +  rnd;
  return;
}





/* function expects that x-X_i etc are appropriately wrapped. */
double Integrator::linker_phi(double x, double y, double z,double Lx,
			      double Ly, double Lz,
			      const std::vector<std::vector<double>> & X_is)
{


  
  double tmp = 0;

  double xtmp,ytmp,ztmp;

  
  for (const auto & X_i : X_is) {

    xtmp = x;
    ytmp = y;
    ztmp = z;

    while (xtmp-X_i[0] > Lx/2)
      xtmp -= Lx;
    while (xtmp-X_i[0] < -Lx/2)
      xtmp += Lx;
    
    while (ytmp-X_i[1] > Ly/2)
      ytmp -= Ly;
    while (ytmp-X_i[1] < -Ly/2)
      ytmp += Ly;

    while (ztmp-X_i[2] > Lz/2)
      ztmp -= Lz;
    while (ztmp-X_i[2] < -Lz/2)
      ztmp += Lz;
    
    tmp += exp(-((xtmp-X_i[0])*(xtmp-X_i[0]) + (ytmp-X_i[1])*(ytmp-X_i[1])
		 + (ztmp-X_i[2])*(ztmp-X_i[2]))/(2*nucwidth*nucwidth));
  }


  
  return chi_LP*tmp;
}

/* function expects that x-X_i etc are appropriately wrapped. */
void Integrator::linker_derivative(std::vector<double> & out,
				   double x, double y, double z,
				   double Lx, double Ly, double Lz,
				   const std::vector<double> & X_i)
{



  double xtmp = x;
  double ytmp = y;
  double ztmp = z;
  while (xtmp-X_i[0] > Lx/2)
    xtmp -= Lx;
  while (xtmp-X_i[0] < -Lx/2)
    xtmp += Lx;
  
  while (ytmp-X_i[1] > Ly/2)
    ytmp -= Ly;
  while (ytmp-X_i[1] < -Ly/2)
    ytmp += Ly;
  
  while (ztmp-X_i[2] > Lz/2)
    ztmp -= Lz;
  while (ztmp-X_i[2] < -Lz/2)
    ztmp += Lz;
  
    
  double tmp = exp(-((xtmp-X_i[0])*(xtmp-X_i[0]) + (ytmp-X_i[1])*(ytmp-X_i[1])
		     + (ztmp-X_i[2])*(ztmp-X_i[2]))/(2*nucwidth*nucwidth));

  
  out[0] = chi_LP*(xtmp-X_i[0])/(nucwidth*nucwidth)*tmp;
  out[1] = chi_LP*(ytmp-X_i[1])/(nucwidth*nucwidth)*tmp;
  out[2] = chi_LP*(ztmp-X_i[2])/(nucwidth*nucwidth)*tmp;

  
  return;
  
}


/*
void Integrator::ode(std::complex<double> & y, std::complex<double> ynl,
		   std::complex<double> rnd, double q2)
{
  y = (((y-mobility*q2*ynl*dt)*normalization*normalization + rnd )/(1+mobility*gamma*q2*q2*dt));
  return;
}
*/
