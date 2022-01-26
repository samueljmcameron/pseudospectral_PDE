#include "conjplane.hpp"


/* This class writes ALMOST all data in a plane parallel to whatever
   axis (call it y) is split into different MPI processes
   and perpendicular to an axis (call it x) which doesn't contain
   data in the <0 plane. The only data it doesn't write is along
   the third axis (call it z) at z = 0 and z = Nz/2.
*/
ConjPlane::ConjPlane(int numprocs,int mpi_id,
		     MPI_Comm communicator,std::vector<int> all_cNys,
		     int Nz)
  :   block0(), block1(), block2(), block3(),
      line0(), line1(), line2(), line3(), nreq(4),
      mpi_size(numprocs),id(mpi_id),comm(communicator), cNz(Nz)
{



  // I'm assuming that the complex array blocks will always be split in a certain
  // configuration (see logic in if statement on next line).
  for (int i = 1; i < mpi_size; i++ ) {
    if (all_cNys.at(mpi_size-i) == 0 ) {
      std::cout << "EXITING!!! number of processors chosen such that id = " << mpi_size-i
		<< " has array size zero." << std::endl;
      MPI_Abort(comm,1);
    } else if (all_cNys.at(0) < all_cNys.at(i)) {
      std::cout << "EXITING!!! Block size of process (id = " << i
		<< ") is larger than block size of initial process id = 0." << std::endl;
      MPI_Abort(comm,1);
    } else if (all_cNys.at(i-1) != all_cNys.at(i) && i != mpi_size-1) {
      std::cout << "EXITING!!! Block size of process (id = " << i
		<< ") is not the same size as the previous block." << std::endl;
      MPI_Abort(comm,1);
    }
    
  }


  int x = all_cNys.at(0);
  int y = all_cNys.at(mpi_size-1);


  
  if (id == 0 ) {
    loopstart = 1;
  } else {
    loopstart = 0;
  }

  block_requests = new MPI_Request[nreq];
  line_requests = new MPI_Request[nreq];

  if (x == y) {

    equal_flag = true;
    divider = mpi_size/2;
    
    block0.size(1,Nz/2);
    block1.size(x-1,Nz/2);
    block2.size(1,Nz/2);
    block3.size(x-1,Nz/2);

    line0.size(1,1);
    line1.size(x-1,1);
    line2.size(1,1);
    line3.size(x-1,1);
    
  } else {
    
    equal_flag = false;
    divider = (mpi_size-1)/2;

    
    if (id == 0 || id == mpi_size-1) {
      block0.size(y,Nz/2);
      block2.size(y,Nz/2);

      line0.size(y,1);
      line2.size(y,1);
      
    } else {
      block0.size(y+1,Nz/2);
      block2.size(y+1,Nz/2);
      
      line0.size(y+1,1);
      line2.size(y+1,1);
      
    }
    
    block1.size(x-y-1,Nz/2);
    block3.size(x-y-1,Nz/2);

    line1.size(x-y-1,1);
    line3.size(x-y-1,1);
    
  }    

}

ConjPlane::~ConjPlane()
{
  delete[] block_requests;
  delete[] line_requests;
}



void ConjPlane::single(TimeStep & timestep,int nx)
/* Only to be used when there is only one processor. */
{
  
  int cNy = timestep.ft_phi.axis_size(0);
  
  
  // write line at ny = 0  
  int ny = 0;
  for (int j = 1; j < cNz/2; j++) {
    
    timestep.update(ny,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    timestep.ft_phi(ny,cNz-j,nx) = std::conj(timestep.ft_phi(ny,j,nx));
					     
  }
  
  
  // write bottom left square to top right square where possible
  for (int i = 1; i < cNy/2; i++) {
    
    for (int j = 1; j < cNz/2; j++) {
      
      timestep.update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      timestep.ft_phi(cNy-i,cNz-j,nx) = std::conj(timestep.ft_phi(i,j,nx));
      
    }
  }
  // write bottom right square to top left square where possible
  for (int i = cNy/2+1; i < cNy; i++) {
    for (int j = 1; j < cNz/2; j++ ) {
      timestep.update(i,j,nx);
      timestep.ft_phi(cNy-i,cNz-j,nx) = std::conj(timestep.ft_phi(i,j,nx));
    }
  }
  
  
  // don't forget Ny/2 term!
  ny = cNy/2;
  
  
  for (int j = 1; j < cNz/2; j++) {
    
    timestep.update(ny,j,nx);//1.0*(ny+complex_local) + 1i*(1.0*j);
    timestep.ft_phi(ny,cNz-j,nx) = std::conj(timestep.ft_phi(ny,j,nx));
    
  }
  
  return;
  
}

void ConjPlane::line_single(TimeStep & timestep,int nx)
/* Only to be used when there is only one processor. */
{

  int cNy = timestep.ft_phi.axis_size(0);



  
  int nz = 0;
  // write bottom left square to bottom right square where possible
  for (int i = 1; i < cNy/2; i++) {
    timestep.update(i,nz,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    timestep.ft_phi(cNy-i,nz,nx) = std::conj(timestep.ft_phi(i,nz,nx));
      
  }
  nz = cNz/2;
  
  // write top right square to top left square where possible
  for (int i = cNy/2+1; i < cNy; i++) {
    timestep.update(i,nz,nx);
    timestep.ft_phi(cNy-i,nz,nx) = std::conj(timestep.ft_phi(i,nz,nx));
  }

  return;
  
}

  


void ConjPlane::first(TimeStep & timestep,int nx)
/* only to be used when all processors have same
   number of data points on them. */
{
  
  int cNy = timestep.ft_phi.axis_size(0);

  // write the line at ny = 0;

  int ny = 0;
  for (int j = 1; j < cNz/2; j++) {
    
    timestep.update(ny,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    timestep.ft_phi(ny,cNz-j,nx) = std::conj(timestep.ft_phi(ny,j,nx));
					     
  }
  
  // send from far left corner of bottom left square

  for (int i = 1; i < cNy; i++) {

    for (int j = 1; j < cNz/2; j++) {

      timestep.update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      //(ft_phi(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);
      
      
      block1(block1.axis_size(0)-i,block1.axis_size(1)-j) = std::conj(timestep.ft_phi(i,j,nx));
      
    }
  }
  
  MPI_Isend(block1.data(),block1.size()*2,MPI_DOUBLE,mpi_size-1,1,
	    comm,&block_requests[0]);
  
  // receive from bottom right square 
  MPI_Recv(block3.data(), block3.size()*2,MPI_DOUBLE,mpi_size-1,3,
	   comm,MPI_STATUS_IGNORE);
  
  for (int i = 1; i < cNy; i++) {
    for (int j = cNz/2+1; j < cNz; j++) {
      
      timestep.ft_phi(i,j,nx) = block3(i-1,j-block3.axis_size(1));
      
      
    }
  }
  
  
  MPI_Wait(&block_requests[0],MPI_STATUS_IGNORE);
  
  return;

}

void ConjPlane::last(TimeStep & timestep,int nx)
/* only to be used when the last processor has
   less data points than the rest. */

{
  
  int cNy = timestep.ft_phi.axis_size(0);

  
  // send from far right corner of bottom right square
  for (int i = 0; i < block2.axis_size(0); i++) {

    for (int j = 1; j < cNz/2; j++) {
      
      
      timestep.update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      //(ft_phi(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);
      
      
      block2(block2.axis_size(0)-i-1,block2.axis_size(1)-j) = std::conj(timestep.ft_phi(i,j,nx));
      
    }
  }
  
  MPI_Isend(block2.data(),block2.size()*2,MPI_DOUBLE,mpi_size-id-1,2,
	    comm,&block_requests[2]);
  
  // receive from bottom left square 1 <= y < Ny/2, z <= 1 < Nz/2
  MPI_Recv(block0.data(), block0.size()*2,MPI_DOUBLE,mpi_size-id-1,0,
	   comm,MPI_STATUS_IGNORE);
  
  for (int i = 0; i < block0.axis_size(0); i++) {
    for (int j = cNz/2+1; j < cNz; j++) {
      
      timestep.ft_phi(i,j,nx) = block0(i,j-block0.axis_size(1));
      
      
    }
  }
  
  
  MPI_Wait(&block_requests[2],MPI_STATUS_IGNORE);
  
  return;

}


void ConjPlane::lefthalf(TimeStep & timestep,int nx)
{
  if (equal_flag)
    lefthalf_equal(timestep,nx);
  else
    lefthalf_unequal(timestep,nx);
  return;
}

void ConjPlane::righthalf(TimeStep & timestep,int nx)
{
  if (equal_flag)
    righthalf_equal(timestep,nx);
  else
    righthalf_unequal(timestep,nx);
  return;
}


void ConjPlane::middle_odd(TimeStep & timestep,int nx)
{
  if (equal_flag)
    middle_odd_equal(timestep,nx);
  else
    middle_odd_unequal(timestep,nx);
  return;
}

void ConjPlane::middle_even(TimeStep & timestep,int nx)
{
  if (equal_flag)
    middle_even_equal(timestep,nx);
  else
    middle_even_unequal(timestep,nx);
  return;
}


void ConjPlane::line_first(TimeStep &timestep, int nx)
/* only to be used when all processors have same
   number of data points on them. */

{
  
  int cNy = timestep.ft_phi.axis_size(0);
  int nz = 0;
  
  // send from far left corner of bottom left square
  for (int i = 1; i < cNy; i++) {

    
    timestep.update(i,nz,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    //(ft_phi(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);
      
      
    line1(line1.axis_size(0)-i,0) = std::conj(timestep.ft_phi(i,nz,nx));
      
    
  }
  
  MPI_Isend(line1.data(),line1.size()*2,MPI_DOUBLE,mpi_size-1,1,
	    comm,&line_requests[0]);
  
  // receive from bottom right square 
  MPI_Recv(line3.data(), line3.size()*2,MPI_DOUBLE,mpi_size-1,3,
	   comm,MPI_STATUS_IGNORE);
  
  nz = cNz/2;
  for (int i = 1; i < cNy; i++) {
    timestep.ft_phi(i,nz,nx) = line3(i-1,0);
  }
  
  
  MPI_Wait(&line_requests[0],MPI_STATUS_IGNORE);
  
  return;

}


void ConjPlane::line_last(TimeStep & timestep,int nx)
/* only to be used when the last processor has
   less data points than the rest. */
  
{
  
  int cNy = timestep.ft_phi.axis_size(0);
  int nz = cNz/2;
  
  // send from far right corner of bottom right square
  for (int i = 0; i < line2.axis_size(0); i++) {
    timestep.update(i,nz,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    //(ft_phi(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);

    line2(line2.axis_size(0)-i-1,0) = std::conj(timestep.ft_phi(i,nz,nx));

  }
  
  MPI_Isend(line2.data(),line2.size()*2,MPI_DOUBLE,mpi_size-id-1,2,
	    comm,&line_requests[2]);
  
  // receive from bottom left square 1 <= y < Ny/2, z <= 1 < Nz/2
  MPI_Recv(line0.data(), line0.size()*2,MPI_DOUBLE,mpi_size-id-1,0,
	   comm,MPI_STATUS_IGNORE);

  nz = 0;
  for (int i = 0; i < line0.axis_size(0); i++) {
    
    timestep.ft_phi(i,nz,nx) = line0(i,0);

  }
  
  MPI_Wait(&line_requests[2],MPI_STATUS_IGNORE);
  
  return;

}


void ConjPlane::line_lefthalf(TimeStep & timestep,int nx)
{
  if (equal_flag)
    line_lefthalf_equal(timestep,nx);
  else
    line_lefthalf_unequal(timestep,nx);
  return;
}

void ConjPlane::line_righthalf(TimeStep & timestep,int nx)
{
  if (equal_flag)
    line_righthalf_equal(timestep,nx);
  else
    line_righthalf_unequal(timestep,nx);
  return;
}


void ConjPlane::line_middle_odd(TimeStep & timestep,int nx)
{
  if (equal_flag)
    line_middle_odd_equal(timestep,nx);
  else
    line_middle_odd_unequal(timestep,nx);
  return;
}

void ConjPlane::line_middle_even(TimeStep & timestep,int nx)
{
  if (equal_flag)
    line_middle_even_equal(timestep,nx);
  else
    line_middle_even_unequal(timestep,nx);
  return;
}

















void ConjPlane::lefthalf_equal(TimeStep & timestep,int nx)
/* only to be used when all processors have same
   number of data points on them. */
{

  int cNy = timestep.ft_phi.axis_size(0);

  
  // single line at start of bottom left plane square

  int ny = 0;
  for (int j = 1; j < cNz/2; j++) {
    
    timestep.update(ny,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    //(ft_phi(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);
    
    block0(ny,block0.axis_size(1)-j) = std::conj(timestep.ft_phi(ny,j,nx));
    
  }
  

  
  MPI_Isend(block0.data(),block0.size()*2,MPI_DOUBLE,mpi_size-id,0,
	    comm,&block_requests[0]);
  
  for (int i = 1; i < cNy ; i++) {
    
    for (int j = 1; j < cNz/2; j++ ) {
      
      timestep.update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      block1(block1.axis_size(0)-i,block1.axis_size(1)-j)
	= std::conj(timestep.ft_phi(i,j,nx));
      
    }
  }
  
  MPI_Isend(block1.data(),block1.size()*2,MPI_DOUBLE,mpi_size-id-1,1,
	    comm,&block_requests[1]);
  
  
  
  // receive from bottom right square to fill in top left square
  // (1 <= y < Ny/2, Nz/2 < z < Nz)
  
  
  MPI_Irecv(block2.data(), block2.size()*2,MPI_DOUBLE,mpi_size-id,2,
	    comm,&block_requests[2]);
  
  MPI_Irecv(block3.data(), block3.size()*2,MPI_DOUBLE,mpi_size-id-1,3,
	    comm,&block_requests[3]);
  
  
  // now fill phi with block3 and block4 received data in an efficient manner
  
  int b2flag=0;
  int b3flag = 0;
  
  
  // next block of code ensures that whichever block data is received first
  // is filled in first
  while (!(b2flag || b3flag)) {
    MPI_Test(&block_requests[2],&b2flag,MPI_STATUS_IGNORE);
    MPI_Test(&block_requests[3],&b3flag,MPI_STATUS_IGNORE);
  }
  
  
  if (b2flag) { // if block2 is received first

    for (int j = cNz/2+1; j < cNz; j++) {
	
	timestep.ft_phi(ny,j,nx) = block2(ny,j-block2.axis_size(1));
	
    }
    
    MPI_Wait(&block_requests[3],MPI_STATUS_IGNORE);
    
    for (int i = 1; i < cNy; i++) {
      for (int j = cNz/2+1; j < cNz; j++) {
	
	timestep.ft_phi(i,j,nx) = block3(i-1,j-block3.axis_size(1));
	
      }
    }
  } else if (b3flag) { // if block3 is received first

    for (int i = 1; i < cNy; i++) {
      for (int j = cNz/2+1; j < cNz; j++) {
	
	timestep.ft_phi(i,j,nx) = block3(i-1,j-block3.axis_size(1));
	
      }
    }
    
    MPI_Wait(&block_requests[2],MPI_STATUS_IGNORE);
    
    
    
    for (int j = cNz/2+1; j < cNz; j++) {
      
      timestep.ft_phi(ny,j,nx) = block2(ny,j-block2.axis_size(1));
      
    }


  }
  
  MPI_Waitall(nreq-2,block_requests,MPI_STATUSES_IGNORE);    

  return;
}


void ConjPlane::righthalf_equal(TimeStep & timestep,int nx)
/* only to be used when all processors have same
   number of data points on them. */
{

  int cNy = timestep.ft_phi.axis_size(0);

  
  // send from far right corner of bottom right square
  
  int ny = 0;

  for (int j = 1; j < cNz/2; j++) {
      
    
    timestep.update(ny,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      //(ft_phi(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);
    
    
    block2(0,block2.axis_size(1)-j) = std::conj(timestep.ft_phi(ny,j,nx));
    
  }


  MPI_Isend(block2.data(),block2.size()*2,MPI_DOUBLE,mpi_size-id,2,
	      comm,&block_requests[2]);
    
    
  for (int i = 1; i < cNy; i++) {
    for (int j = 1; j < cNz/2; j++) {
      
      timestep.update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      
      block3(block3.axis_size(0)-i,block3.axis_size(1)-j)
	= std::conj(timestep.ft_phi(i,j,nx));
    }
    
  }
  
  
  MPI_Isend(block3.data(),block3.size()*2,MPI_DOUBLE,mpi_size-id-1,3,
	    comm,&block_requests[3]);
  
  
  // receive both from bottom left square, then fill in top right square
  
  MPI_Irecv(block0.data(), block0.size()*2,MPI_DOUBLE,mpi_size-id,0,
	    comm,&block_requests[0]);
  
  
  MPI_Irecv(block1.data(), block1.size()*2,MPI_DOUBLE,mpi_size-id-1,1,
	    comm,&block_requests[1]);
  
  
  int b0flag=0;
  int b1flag = 0;

  
  // next block of code ensures that whichever block data is received first
  // is filled in first
  while (!(b0flag || b1flag)) {
    MPI_Test(&block_requests[0],&b0flag,MPI_STATUS_IGNORE);
    MPI_Test(&block_requests[1],&b1flag,MPI_STATUS_IGNORE);
  }
  
  if (b0flag) { // if block0 is received first
    
    for (int j = cNz/2+1; j < cNz; j++) {
      
      timestep.ft_phi(ny,j,nx) = block0(ny,j-block0.axis_size(1));
      
    }
    
    
    MPI_Wait(&block_requests[1],MPI_STATUS_IGNORE);
    
    for (int i = 1; i < cNy; i++) {
      for (int j = cNz/2+1; j < cNz; j++) {
	
	timestep.ft_phi(i,j,nx) = block1(i-1,j-block1.axis_size(1));
	
      }
    }
    
    
  } else if (b1flag) { // if block1 is received first
    
    for (int i = 1; i < cNy; i++) {
      for (int j = cNz/2+1; j < cNz; j++) {
	
	timestep.ft_phi(i,j,nx) = block1(i-1,j-block1.axis_size(1));
	
      }
    }
    
    MPI_Wait(&block_requests[0],MPI_STATUS_IGNORE);
    

    for (int j = cNz/2+1; j < cNz; j++) {
      
      timestep.ft_phi(ny,j,nx) = block0(ny,j-block0.axis_size(1));
      
    }
    
    
  }
  
  MPI_Waitall(nreq-2,&block_requests[2],MPI_STATUSES_IGNORE);


  
  return;
}
  


void ConjPlane::middle_odd_equal(TimeStep & timestep,int nx)
/* only to be used when all processors have same
   number of data points on them. */
{

  int cNy = timestep.ft_phi.axis_size(0);



  // send one to next processor

  int ny = 0;
  for (int j = 1; j < cNz/2; j++) {
    
    timestep.update(ny,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    //(ft_phi(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);
    
    block0(ny,block0.axis_size(1)-j) = std::conj(timestep.ft_phi(ny,j,nx));
    
  }
  

  
  MPI_Isend(block0.data(),block0.size()*2,MPI_DOUBLE,mpi_size-id,0,
	    comm,&block_requests[0]);


  // receive one from next processor

  MPI_Recv(block2.data(), block2.size()*2,MPI_DOUBLE,mpi_size-id,2,
	   comm,MPI_STATUS_IGNORE);

  for (int j = cNz/2+1; j < cNz; j++) {
    
    timestep.ft_phi(ny,j,nx) = block2(ny,j-block2.axis_size(1));
    
  }


  
  // write bottom left square to top right square where possible
  for (int i = 1; i < (block1.axis_size(0)-1)/2+1; i++) {

    for (int j = 1; j < cNz/2; j++) {
      
      
      timestep.update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      timestep.ft_phi(cNy-i,cNz-j,nx) = std::conj(timestep.ft_phi(i,j,nx));
      
    }
  }
  // write bottom right square to top left square where possible
  for (int i = (block1.axis_size(0)-1)/2 + 2; i < cNy; i++) {
    for (int j = 1; j < cNz/2; j++ ) {
      timestep.update(i,j,nx);
      timestep.ft_phi(cNy-i,cNz-j,nx) = std::conj(timestep.ft_phi(i,j,nx));
    }
  }

  
  // don't forget Ny/2 term!
  ny = (block1.axis_size(0)-1)/2+1;
  
  
  for (int j = 1; j < cNz/2; j++) {
    
    timestep.update(ny,j,nx);//1.0*(ny+complex_local) + 1i*(1.0*j);
    timestep.ft_phi(ny,cNz-j,nx) = std::conj(timestep.ft_phi(ny,j,nx));
    
  }

  MPI_Wait(&block_requests[0],MPI_STATUS_IGNORE);

  return;
  
}


void ConjPlane::middle_even_equal(TimeStep & timestep,int nx)
/* only to be used when all processors have same
   number of data points on them. */
{




  // don't forget Ny/2 term!
  int ny = 0;
  
  
  for (int j = 1; j < cNz/2; j++) {
    
    timestep.update(ny,j,nx);//1.0*(ny+complex_local) + 1i*(1.0*j);
    timestep.ft_phi(ny,cNz-j,nx) = std::conj(timestep.ft_phi(ny,j,nx));
    
  }

  int cNy = timestep.ft_phi.axis_size(0);

  // fill and send bottom right square  
    
  for (int i = 1; i < cNy; i++) {
    for (int j = 1; j < cNz/2; j++) {
      
      timestep.update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      
      block3(block3.axis_size(0)-i,block3.axis_size(1)-j)
	= std::conj(timestep.ft_phi(i,j,nx));
    }
    
  }
  
  
  MPI_Isend(block3.data(),block3.size()*2,MPI_DOUBLE,mpi_size-id-1,3,
	    comm,&block_requests[3]);
  
  
  // receive from bottom left square, then fill in top right square
  
  MPI_Recv(block1.data(), block1.size()*2,MPI_DOUBLE,mpi_size-id-1,1,
	    comm,MPI_STATUS_IGNORE);
  

  for (int i = 1; i < cNy; i++) {
    for (int j = cNz/2+1; j < cNz; j++) {
      
      timestep.ft_phi(i,j,nx) = block1(i-1,j-block1.axis_size(1));
      
    }
  }

  MPI_Wait(&block_requests[3],MPI_STATUS_IGNORE);
  
  return;
  
}



void ConjPlane::lefthalf_unequal(TimeStep & timestep,int nx)
/* only to be used when the last processor has
   less data points than the rest. */
{

  int cNy = timestep.ft_phi.axis_size(0);

  if (id == 0) {
    // write the line at ny = 0;
    
    int ny = 0;
    for (int j = 1; j < cNz/2; j++) {
    
      timestep.update(ny,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      timestep.ft_phi(ny,cNz-j,nx) = std::conj(timestep.ft_phi(ny,j,nx));
      
    }
  }

  
  
  // bottom left plane square (1 <= y < Ny/2, 1 <= z < Nz/2)
  
  for (int i = loopstart; i < block0.axis_size(0)+loopstart; i++) {
    
    for (int j = 1; j < cNz/2; j++) {
      
      
      timestep.update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      //(ft_phi(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);

      block0(block0.axis_size(0)-i-(1-loopstart),block0.axis_size(1)-j) = std::conj(timestep.ft_phi(i,j,nx));
      
    }

  }
  
  MPI_Isend(block0.data(),block0.size()*2,MPI_DOUBLE,mpi_size-id-1,0,
	    comm,&block_requests[0]);
  
  for (int i = block0.axis_size(0)+loopstart; i < cNy ; i++) {
    
    for (int j = 1; j < cNz/2; j++ ) {
      
      timestep.update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      block1(block1.axis_size(0)-(i-block0.axis_size(0)+1-loopstart),block1.axis_size(1)-j)
	= std::conj(timestep.ft_phi(i,j,nx));
      
    }
  }
  
  MPI_Isend(block1.data(),block1.size()*2,MPI_DOUBLE,mpi_size-id-2,1,
	    comm,&block_requests[1]);
  
  
  
  // receive from bottom right square to fill in top left square
  // (1 <= y < Ny/2, Nz/2 < z < Nz)
  
  
  MPI_Irecv(block2.data(), block2.size()*2,MPI_DOUBLE,mpi_size-id-1,2,
	    comm,&block_requests[2]);
  
  MPI_Irecv(block3.data(), block3.size()*2,MPI_DOUBLE,mpi_size-id-2,3,
	    comm,&block_requests[3]);
  
  
  // now fill phi with block3 and block4 received data in an efficient manner
  
  int b2flag=0;
  int b3flag = 0;
  
  
  // next block of code ensures that whichever block data is received first
  // is filled in first
  while (!(b2flag || b3flag)) {
    MPI_Test(&block_requests[2],&b2flag,MPI_STATUS_IGNORE);
    MPI_Test(&block_requests[3],&b3flag,MPI_STATUS_IGNORE);
  }
  
  
  if (b2flag) { // if block2 is received first
    for (int i = loopstart; i < block2.axis_size(0)+loopstart; i++) {
      for (int j = cNz/2+1; j < cNz; j++) {
	
	timestep.ft_phi(i,j,nx) = block2(i-loopstart,j-block2.axis_size(1));
	
      }
    }
    
    MPI_Wait(&block_requests[3],MPI_STATUS_IGNORE);
    
    for (int i = block2.axis_size(0)+loopstart; i < cNy; i++) {
      for (int j = cNz/2+1; j < cNz; j++) {
	
	timestep.ft_phi(i,j,nx) = block3(i-(block2.axis_size(0)+loopstart),j-block3.axis_size(1));
	
      }
    }
  } else if (b3flag) { // if block3 is received first

    for (int i = block2.axis_size(0)+loopstart; i < cNy; i++) {
      for (int j = cNz/2+1; j < cNz; j++) {
	
	timestep.ft_phi(i,j,nx) = block3(i-(block2.axis_size(0)+loopstart),j-block3.axis_size(1));
	
      }
    }
    
    MPI_Wait(&block_requests[2],MPI_STATUS_IGNORE);
    
    
    for (int i = loopstart; i < block2.axis_size(0)+loopstart; i++) {
      for (int j = cNz/2+1; j < cNz; j++) {
	
	timestep.ft_phi(i,j,nx) = block2(i-loopstart,j-block2.axis_size(1));
	
      }
    }      

  }
  
  MPI_Waitall(nreq-2,block_requests,MPI_STATUSES_IGNORE);    

  return;
}


void ConjPlane::righthalf_unequal(TimeStep & timestep,int nx)
/* only to be used when the last processor has
   less data points than the rest. */

{

  int cNy = timestep.ft_phi.axis_size(0);

  
  // send from far right corner of bottom right square
  
  for (int i = 0; i < block2.axis_size(0); i++) {

    for (int j = 1; j < cNz/2; j++) {
      
      
      timestep.update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      //(ft_phi(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);
      
      
      block2(block2.axis_size(0)-i-1,block2.axis_size(1)-j) = std::conj(timestep.ft_phi(i,j,nx));
      
    }
  }
  
  MPI_Isend(block2.data(),block2.size()*2,MPI_DOUBLE,mpi_size-id-1,2,
	      comm,&block_requests[2]);
    
    
  for (int i = block2.axis_size(0); i < cNy; i++) {
    for (int j = 1; j < cNz/2; j++) {
      
      timestep.update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      
      block3(block3.axis_size(0)-(i-block2.axis_size(0)+1),block3.axis_size(1)-j)
	= std::conj(timestep.ft_phi(i,j,nx));
    }
    
  }
  
  
  MPI_Isend(block3.data(),block3.size()*2,MPI_DOUBLE,mpi_size-id-2,3,
	    comm,&block_requests[3]);
  
  
  // receive both from bottom left square, then fill in top right square
  
  MPI_Irecv(block0.data(), block0.size()*2,MPI_DOUBLE,mpi_size-id-1,0,
	    comm,&block_requests[0]);
  
  
  MPI_Irecv(block1.data(), block1.size()*2,MPI_DOUBLE,mpi_size-id-2,1,
	    comm,&block_requests[1]);
  
  
  int b0flag=0;
  int b1flag = 0;

  
  // next block of code ensures that whichever block data is received first
  // is filled in first
  while (!(b0flag || b1flag)) {
    MPI_Test(&block_requests[0],&b0flag,MPI_STATUS_IGNORE);
    MPI_Test(&block_requests[1],&b1flag,MPI_STATUS_IGNORE);
  }
  
  if (b0flag) { // if block0 is received first
    for (int i = 0; i < block0.axis_size(0); i++) {
      for (int j = cNz/2+1; j < cNz; j++) {
	
	timestep.ft_phi(i,j,nx) = block0(i,j-block0.axis_size(1));
	
      }
    }

    MPI_Wait(&block_requests[1],MPI_STATUS_IGNORE);
    
    for (int i = block0.axis_size(0); i < cNy; i++) {
      for (int j = cNz/2+1; j < cNz; j++) {
	
	timestep.ft_phi(i,j,nx) = block1(i-block0.axis_size(0),j-block1.axis_size(1));
	
      }
    }
    
    
  } else if (b1flag) { // if block1 is received first
    
    for (int i = block0.axis_size(0); i < cNy; i++) {
      for (int j = cNz/2+1; j < cNz; j++) {
	
	timestep.ft_phi(i,j,nx) = block1(i-block0.axis_size(0),j-block1.axis_size(1));
	
      }
    }
    
    MPI_Wait(&block_requests[0],MPI_STATUS_IGNORE);
    
    for (int i = 0; i < block0.axis_size(0); i++) {
      for (int j = cNz/2+1; j < cNz; j++) {
	
	timestep.ft_phi(i,j,nx) = block0(i,j-block0.axis_size(1));
	
      }
    }      
    
  }
  
  MPI_Waitall(nreq-2,&block_requests[2],MPI_STATUSES_IGNORE);    
  return;
}
  


void ConjPlane::middle_odd_unequal(TimeStep & timestep,int nx)
/* only to be used when the last processor has
   less data points than the rest. */

{

  int cNy = timestep.ft_phi.axis_size(0);

  
  // write bottom left square to top right square where possible
  for (int i = 0; i < (block0.axis_size(0)-1)/2; i++) {

    for (int j = 1; j < cNz/2; j++) {
      
      
      timestep.update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      timestep.ft_phi(cNy-i-block1.axis_size(0)-1,cNz-j,nx) = std::conj(timestep.ft_phi(i,j,nx));
      
    }
    /*
    ny = (globalNy- i - complex_local);
    for (int j = cNz/2 + 1; j < cNz; j++) {
      
      // since conjugatation for top, change expression for ft_phi
      ft_phi(i,j,nx) = std::conj(timestep.update(ny,Nz-j,nx,ft_phi));//1.0*(Ny-i-complex_local)-1i*(1.0*(Nz-j));
      ft_phi(cNy-i-block1.axis_size(0)-1,cNz-j,nx) = std::conj(ft_phi(i,j,nx));
      
      }*/
  }
  // write bottom right square to top left square where possible
  for (int i = (block0.axis_size(0)-1)/2 + 1; i < block0.axis_size(0); i++) {
    for (int j = 1; j < cNz/2; j++ ) {
      timestep.update(i,j,nx);
      timestep.ft_phi(cNy-i-block1.axis_size(0)-1,cNz-j,nx) = std::conj(timestep.ft_phi(i,j,nx));
    }
  }

  
  // don't forget Ny/2 term!
  int ny = (block0.axis_size(0)-1)/2;
  
  
  for (int j = 1; j < cNz/2; j++) {
    
    timestep.update(ny,j,nx);//1.0*(ny+complex_local) + 1i*(1.0*j);
    timestep.ft_phi(ny,cNz-j,nx) = std::conj(timestep.ft_phi(ny,j,nx));
    
  }
  
  
  // write to block and send to previous processor
  for (int i = cNy-block3.axis_size(0); i < cNy; i++) {
    
    for (int j = 1; j < cNz/2; j++) {
      
      timestep.update(i,j,nx);//1.0*(i+complex_local) + 1i*(1.0*j);
      block3(cNy-i-1,block3.axis_size(1)-j) = std::conj(timestep.ft_phi(i,j,nx));
    }
    
  }
  
  MPI_Isend(block3.data(),block3.size()*2,MPI_DOUBLE,mpi_size-id-2,3,
	    comm,&block_requests[3]);
  
  
  // and receive extra data from the previous processor
  MPI_Recv(block1.data(), block1.size()*2,MPI_DOUBLE,mpi_size-id-2,1,
	   comm,MPI_STATUS_IGNORE);
  
  for (int i = cNy-block1.axis_size(0); i < cNy; i++) {
    
    for (int j = cNz/2+1; j < cNz; j++) {
      timestep.ft_phi(i,j,nx) = block1(i-(cNy-block1.axis_size(0)),j-block1.axis_size(1));
    }
  }
  
  
  MPI_Wait(&block_requests[3],MPI_STATUS_IGNORE);

  return;
  
}


void ConjPlane::middle_even_unequal(TimeStep & timestep,int nx)
/* only to be used when the last processor has
   less data points than the rest. */

{

  int cNy = timestep.ft_phi.axis_size(0);



  // write bottom left to next processor

  for (int i = 0; i < block0.axis_size(0); i++) {
    
    for (int j = 1; j < cNz/2; j++) {
      
      
      timestep.update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      //(ft_phi(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);

      block0(block0.axis_size(0)-i-1,block0.axis_size(1)-j) = std::conj(timestep.ft_phi(i,j,nx));
      
    }

  }

  MPI_Isend(block0.data(),block0.size()*2,MPI_DOUBLE,mpi_size-id-1,0,
	    comm,&block_requests[0]);

  // write bottom left to top right in same processor where possible
  for (int i = block0.axis_size(0); i < block0.axis_size(0)+(block1.axis_size(0)-1)/2; i++) {

    for (int j = 1; j < cNz/2; j++) {

      timestep.update(i,j,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      timestep.ft_phi(cNy-(i-block0.axis_size(0))-1,cNz-j,nx) = std::conj(timestep.ft_phi(i,j,nx));

    }
  }

  // write bottom right to top left in the same processor where possible
  for (int i = block0.axis_size(0)+(block1.axis_size(0)-1)/2+1; i < cNy; i++) {
    for (int j = 1; j < cNz/2; j++ ) {
      timestep.update(i,j,nx);
      timestep.ft_phi(cNy-(i-block0.axis_size(0))-1,cNz-j,nx) = std::conj(timestep.ft_phi(i,j,nx));
    }
  }
  
  // don't forget Ny/2 term (on same processor again)!
  int ny = block0.axis_size(0)+(block1.axis_size(0)-1)/2;
  
  
  for (int j = 1; j < cNz/2; j++) {
    
    timestep.update(ny,j,nx);//1.0*(ny+complex_local) + 1i*(1.0*j);
    timestep.ft_phi(ny,cNz-j,nx) = std::conj(timestep.ft_phi(ny,j,nx));
    
  }
  
  
  
  // and receive extra data from the next processor to put in top 
  MPI_Recv(block2.data(), block2.size()*2,MPI_DOUBLE,mpi_size-id-1,2,
	   comm,MPI_STATUS_IGNORE);
  
  for (int i = 0; i < block2.axis_size(0); i++) {
    
    for (int j = cNz/2+1; j < cNz; j++) {
      timestep.ft_phi(i,j,nx) = block2(i,j-block2.axis_size(1));
      //      std::cout << "ft_phi(" << i+timestep.ft_phi.get_local0start() << "," << j
      //		<< ") = " << timestep.ft_phi(i,j,nx) << std::endl;
    }
  }
  
  
  MPI_Wait(&block_requests[0],MPI_STATUS_IGNORE);

  return;
  
}




void ConjPlane::line_lefthalf_equal(TimeStep & timestep,int nx)
{

  int cNy = timestep.ft_phi.axis_size(0);


  int ny = 0;
  int nz = 0;
    
  timestep.update(ny,nz,nx);//1.0*(i+complex_local)+1i*(1.0*j);
  //(ft_phi(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);
    
  line0(0,0) = std::conj(timestep.ft_phi(ny,nz,nx));


  
  MPI_Isend(line0.data(),line0.size()*2,MPI_DOUBLE,mpi_size-id,0,
	    comm,&line_requests[0]);
  
  for (int i = 1; i < cNy ; i++) {
    
    timestep.update(i,nz,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    line1(line1.axis_size(0)-i,0)
      = std::conj(timestep.ft_phi(i,nz,nx));
  }
  
  MPI_Isend(line1.data(),line1.size()*2,MPI_DOUBLE,mpi_size-id-1,1,
	    comm,&line_requests[1]);
  
  
  
  // receive from bottom right square to fill in top left square
  // (1 <= y < Ny/2, Nz/2 < z < Nz)
  
  
  MPI_Irecv(line2.data(), line2.size()*2,MPI_DOUBLE,mpi_size-id,2,
	    comm,&line_requests[2]);
  
  MPI_Irecv(line3.data(), line3.size()*2,MPI_DOUBLE,mpi_size-id-1,3,
	    comm,&line_requests[3]);
  

  // now fill phi with line3 and line4 received data in an efficient manner
  
  int b2flag=0;
  int b3flag = 0;
  nz = cNz/2;  
  
  // next line of code ensures that whichever line data is received first
  // is filled in first
  while (!(b2flag || b3flag)) {
    MPI_Test(&line_requests[2],&b2flag,MPI_STATUS_IGNORE);
    MPI_Test(&line_requests[3],&b3flag,MPI_STATUS_IGNORE);
  }
  
  
  if (b2flag) { // if line2 is received first

	
    timestep.ft_phi(ny,nz,nx) = line2(ny,0);

    
    MPI_Wait(&line_requests[3],MPI_STATUS_IGNORE);
    
    for (int i = 1; i < cNy; i++) {
	timestep.ft_phi(i,nz,nx) = line3(i-1,0);
    }

  } else if (b3flag) { // if line3 is received first

    for (int i = 1; i < cNy; i++) {
      timestep.ft_phi(i,nz,nx) = line3(i-1,0);
    }
    
    MPI_Wait(&line_requests[2],MPI_STATUS_IGNORE);

      
    timestep.ft_phi(ny,nz,nx) = line2(ny,0);
  }
  
  MPI_Waitall(nreq-2,line_requests,MPI_STATUSES_IGNORE);    

  return;
}


void ConjPlane::line_righthalf_equal(TimeStep & timestep,int nx)
{

  int cNy = timestep.ft_phi.axis_size(0);
  int nz = cNz/2;


  int ny = 0;

  timestep.update(ny,nz,nx);//1.0*(i+complex_local)+1i*(1.0*j);
  //(ft_phi(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);
    
    
  line2(0,0) = std::conj(timestep.ft_phi(ny,nz,nx));



  MPI_Isend(line2.data(),line2.size()*2,MPI_DOUBLE,mpi_size-id,2,
	      comm,&line_requests[2]);
    
    
  for (int i = 1; i < cNy; i++) {
    
    timestep.update(i,nz,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      
    line3(line3.axis_size(0)-i,0)
      = std::conj(timestep.ft_phi(i,nz,nx));
  }
  
  
  MPI_Isend(line3.data(),line3.size()*2,MPI_DOUBLE,mpi_size-id-1,3,
	    comm,&line_requests[3]);
  
  
  // receive both from bottom left square, then fill in top right square
  
  MPI_Irecv(line0.data(), line0.size()*2,MPI_DOUBLE,mpi_size-id,0,
	    comm,&line_requests[0]);
  
  
  MPI_Irecv(line1.data(), line1.size()*2,MPI_DOUBLE,mpi_size-id-1,1,
	    comm,&line_requests[1]);
  

  nz = 0;
  
  int b0flag=0;
  int b1flag = 0;

  
  // next line of code ensures that whichever line data is received first
  // is filled in first
  while (!(b0flag || b1flag)) {
    MPI_Test(&line_requests[0],&b0flag,MPI_STATUS_IGNORE);
    MPI_Test(&line_requests[1],&b1flag,MPI_STATUS_IGNORE);
  }
  
  if (b0flag) { // if line0 is received first
    
    timestep.ft_phi(ny,nz,nx) = line0(ny,0);
    
    
    MPI_Wait(&line_requests[1],MPI_STATUS_IGNORE);
    
    for (int i = 1; i < cNy; i++) {
      timestep.ft_phi(i,nz,nx) = line1(i-1,0);
    }
    
    
  } else if (b1flag) { // if line1 is received first
    
    for (int i = 1; i < cNy; i++) {
      timestep.ft_phi(i,nz,nx) = line1(i-1,nz);
    }
    
    MPI_Wait(&line_requests[0],MPI_STATUS_IGNORE);
    


      
    timestep.ft_phi(ny,nz,nx) = line0(ny,0);
  }
  
  MPI_Waitall(nreq-2,&line_requests[2],MPI_STATUSES_IGNORE);
  
  return;
}
  


void ConjPlane::line_middle_odd_equal(TimeStep & timestep,int nx)
{

  int cNy = timestep.ft_phi.axis_size(0);
  int nz = 0;
  int ny = 0;


  // send one to next processor

  timestep.update(ny,nz,nx);//1.0*(i+complex_local)+1i*(1.0*j);
  //(ft_phi(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);
    
  line0(ny,0) = std::conj(timestep.ft_phi(ny,nz,nx));

  MPI_Isend(line0.data(),line0.size()*2,MPI_DOUBLE,mpi_size-id,0,
	    comm,&line_requests[0]);


  // receive one from next processor


  nz = cNz/2;
  MPI_Recv(line2.data(), line2.size()*2,MPI_DOUBLE,mpi_size-id,2,
	   comm,MPI_STATUS_IGNORE);
    
  timestep.ft_phi(ny,nz,nx) = line2(ny,0);



  nz = 0;
  // write bottom left square to bottom right square where possible
  for (int i = 1; i < (line1.axis_size(0)-1)/2+1; i++) {
  
    timestep.update(i,nz,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    timestep.ft_phi(cNy-i,nz,nx) = std::conj(timestep.ft_phi(i,nz,nx));
      
  }

  nz = cNz/2;
  // write top right square to top left square where possible
  for (int i = (line1.axis_size(0)-1)/2 + 2; i < cNy; i++) {

    timestep.update(i,nz,nx);
    timestep.ft_phi(cNy-i,nz,nx) = std::conj(timestep.ft_phi(i,nz,nx));

  }
  MPI_Wait(&line_requests[0],MPI_STATUS_IGNORE);

  return;
  
}


void ConjPlane::line_middle_even_equal(TimeStep & timestep,int nx)
/* only to be used when all processors have same
   number of data points on them. */
{



  // don't forget Ny/2 term!
  int ny = 0;

  int nz = cNz/2;

  int cNy = timestep.ft_phi.axis_size(0);

  // fill and send top right square 
    
  for (int i = 1; i < cNy; i++) {
    timestep.update(i,nz,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      
    line3(line3.axis_size(0)-i,0)
	= std::conj(timestep.ft_phi(i,nz,nx));

  }
  
  
  MPI_Isend(line3.data(),line3.size()*2,MPI_DOUBLE,mpi_size-id-1,3,
	    comm,&line_requests[3]);
  
  
  // receive from bottom left square, then fill in bottom right square
  
  MPI_Recv(line1.data(), line1.size()*2,MPI_DOUBLE,mpi_size-id-1,1,
	   comm,MPI_STATUS_IGNORE);
  
  nz = 0;
  for (int i = 1; i < cNy; i++) {
    timestep.ft_phi(i,nz,nx) = line1(i-1,0);
  }
    
  MPI_Wait(&line_requests[3],MPI_STATUS_IGNORE);
  
  return;
  
}



void ConjPlane::line_lefthalf_unequal(TimeStep & timestep,int nx)
{

  int cNy = timestep.ft_phi.axis_size(0);

  int nz = 0;
  
  // write and send bottom left line (1 <= y < Ny/2, z =0)
  
  for (int i = loopstart; i < line0.axis_size(0)+loopstart; i++) {

    timestep.update(i,nz,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    //(ft_phi(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);

    line0(line0.axis_size(0)-i-(1-loopstart),0) = std::conj(timestep.ft_phi(i,nz,nx));
      
  }
  
  MPI_Isend(line0.data(),line0.size()*2,MPI_DOUBLE,mpi_size-id-1,0,
	    comm,&line_requests[0]);
  
  for (int i = line0.axis_size(0)+loopstart; i < cNy ; i++) {

    timestep.update(i,nz,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    line1(line1.axis_size(0)-(i-line0.axis_size(0)+1-loopstart),0)
      = std::conj(timestep.ft_phi(i,nz,nx));
  }
  
  MPI_Isend(line1.data(),line1.size()*2,MPI_DOUBLE,mpi_size-id-2,1,
	    comm,&line_requests[1]);
  
  
  
  // receive from top right side line to fill in top left line
  // (1 <= y < Ny/2, z = Nz/2)
  

  nz = cNz/2;
  
  MPI_Irecv(line2.data(), line2.size()*2,MPI_DOUBLE,mpi_size-id-1,2,
	    comm,&line_requests[2]);
  
  MPI_Irecv(line3.data(), line3.size()*2,MPI_DOUBLE,mpi_size-id-2,3,
	    comm,&line_requests[3]);
  
  
  // now fill phi with line3 and line4 received data in an efficient manner
  
  int b2flag=0;
  int b3flag = 0;
  
  
  // next line of code ensures that whichever line data is received first
  // is filled in first
  while (!(b2flag || b3flag)) {
    MPI_Test(&line_requests[2],&b2flag,MPI_STATUS_IGNORE);
    MPI_Test(&line_requests[3],&b3flag,MPI_STATUS_IGNORE);
  }
  
  
  if (b2flag) { // if line2 is received first
    for (int i = loopstart; i < line2.axis_size(0)+loopstart; i++) {
      timestep.ft_phi(i,nz,nx) = line2(i-loopstart,0);
      
    }

    
    MPI_Wait(&line_requests[3],MPI_STATUS_IGNORE);
    
    for (int i = line2.axis_size(0)+loopstart; i < cNy; i++) {
      timestep.ft_phi(i,nz,nx) = line3(i-(line2.axis_size(0)+loopstart),0);
    }
  } else if (b3flag) { // if line3 is received first

    for (int i = line2.axis_size(0)+loopstart; i < cNy; i++) {
	timestep.ft_phi(i,nz,nx) = line3(i-(line2.axis_size(0)+loopstart),0);
    }
    
    MPI_Wait(&line_requests[2],MPI_STATUS_IGNORE);
    
    
    for (int i = loopstart; i < line2.axis_size(0)+loopstart; i++) {
      timestep.ft_phi(i,nz,nx) = line2(i-loopstart,0);
	
    }

  }
  
  MPI_Waitall(nreq-2,line_requests,MPI_STATUSES_IGNORE);    

  
  return;
}


void ConjPlane::line_righthalf_unequal(TimeStep & timestep,int nx)
{

  int cNy = timestep.ft_phi.axis_size(0);
  int nz = cNz/2;
  
  // write and send top right line (1 <= y < Ny/2, z = cNz/2)
  
  for (int i = 0; i < line2.axis_size(0); i++) {
    timestep.update(i,nz,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    //(ft_phi(i,j,nx) - mobility*q2*ft_nonlinear(i,j,nx)*dt + noise*dt )/(1+mobility*gamma*q2*q2);

    line2(line2.axis_size(0)-i-1,0) = std::conj(timestep.ft_phi(i,nz,nx));
    
  }
  
  MPI_Isend(line2.data(),line2.size()*2,MPI_DOUBLE,mpi_size-id-1,2,
	      comm,&line_requests[2]);
    
    
  for (int i = line2.axis_size(0); i < cNy; i++) {
      
    timestep.update(i,nz,nx);//1.0*(i+complex_local)+1i*(1.0*j);
      
    line3(line3.axis_size(0)-(i-line2.axis_size(0)+1),0)
      = std::conj(timestep.ft_phi(i,nz,nx));
  }

  MPI_Isend(line3.data(),line3.size()*2,MPI_DOUBLE,mpi_size-id-2,3,
	    comm,&line_requests[3]);
  
  
  // receive from bottom left line to fill bottom right line
  // (1 <= y < Ny/2, z = 0)

  nz = 0;
  
  MPI_Irecv(line0.data(), line0.size()*2,MPI_DOUBLE,mpi_size-id-1,0,
	    comm,&line_requests[0]);
  
  
  MPI_Irecv(line1.data(), line1.size()*2,MPI_DOUBLE,mpi_size-id-2,1,
	    comm,&line_requests[1]);
  
  
  int b0flag=0;
  int b1flag = 0;

  
  // next line of code ensures that whichever line data is received first
  // is filled in first
  while (!(b0flag || b1flag)) {
    MPI_Test(&line_requests[0],&b0flag,MPI_STATUS_IGNORE);
    MPI_Test(&line_requests[1],&b1flag,MPI_STATUS_IGNORE);
  }
  
  if (b0flag) { // if line0 is received first
    for (int i = 0; i < line0.axis_size(0); i++) {
      timestep.ft_phi(i,nz,nx) = line0(i,0);

    }

    MPI_Wait(&line_requests[1],MPI_STATUS_IGNORE);
    
    for (int i = line0.axis_size(0); i < cNy; i++) {
	
      timestep.ft_phi(i,nz,nx) = line1(i-line0.axis_size(0),0);

    }
    
    
  } else if (b1flag) { // if line1 is received first
    
    for (int i = line0.axis_size(0); i < cNy; i++) {
	timestep.ft_phi(i,nz,nx) = line1(i-line0.axis_size(0),0);

    }
    
    MPI_Wait(&line_requests[0],MPI_STATUS_IGNORE);
    
    for (int i = 0; i < line0.axis_size(0); i++) {
      timestep.ft_phi(i,nz,nx) = line0(i,0);

    }      
    
  }
  
  MPI_Waitall(nreq-2,&line_requests[2],MPI_STATUSES_IGNORE);    
  return;
}
  


void ConjPlane::line_middle_odd_unequal(TimeStep & timestep,int nx)
{

  int cNy = timestep.ft_phi.axis_size(0);
  int nz = 0;
  
  // write bottom left square to bottom right square where possible
  for (int i = 0; i < (line0.axis_size(0)-1)/2; i++) {
    
    timestep.update(i,nz,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    timestep.ft_phi(cNy-i-line1.axis_size(0)-1,nz,nx) = std::conj(timestep.ft_phi(i,nz,nx));
      
  }

  // write top right square to top left square where possible

  nz = cNz/2;
  for (int i = (line0.axis_size(0)-1)/2 + 1; i < line0.axis_size(0); i++) {
    timestep.update(i,nz,nx);
    timestep.ft_phi(cNy-i-line1.axis_size(0)-1,nz,nx) = std::conj(timestep.ft_phi(i,nz,nx));
  }

  
  // write to line and send to previous processor

  nz = cNz/2;
  
  for (int i = cNy-line3.axis_size(0); i < cNy; i++) {
    timestep.update(i,nz,nx);//1.0*(i+complex_local) + 1i*(1.0*j);
    line3(cNy-i-1,0) = std::conj(timestep.ft_phi(i,nz,nx));
  }
  
  MPI_Isend(line3.data(),line3.size()*2,MPI_DOUBLE,mpi_size-id-2,3,
	    comm,&line_requests[3]);
  
  
  // and receive extra data from the previous processor
  MPI_Recv(line1.data(), line1.size()*2,MPI_DOUBLE,mpi_size-id-2,1,
	   comm,MPI_STATUS_IGNORE);

  nz = 0;
  for (int i = cNy-line1.axis_size(0); i < cNy; i++) {
    timestep.ft_phi(i,nz,nx) = line1(i-(cNy-line1.axis_size(0)),0);
  }
  
  
  MPI_Wait(&line_requests[3],MPI_STATUS_IGNORE);

  return;
  
}

void ConjPlane::line_middle_even_unequal(TimeStep & timestep,int nx)
{

  int cNy = timestep.ft_phi.axis_size(0);
  int nz = 0;

  // write bottom left and send to next processor
  for (int i = 0; i < line0.axis_size(0); i++) {
      
    timestep.update(i,nz,nx);//1.0*(i+complex_local)+1i*(1.0*j);


    line0(line0.axis_size(0)-i-1,0) = std::conj(timestep.ft_phi(i,nz,nx));
      
  }

  // 

  MPI_Isend(line0.data(),line0.size()*2,MPI_DOUBLE,mpi_size-id-1,0,
	    comm,&line_requests[0]);

  // write bottom left to bottom right in same processor where possible
  for (int i = line0.axis_size(0); i < line0.axis_size(0)+(line1.axis_size(0)-1)/2; i++) {
    timestep.update(i,nz,nx);//1.0*(i+complex_local)+1i*(1.0*j);
    timestep.ft_phi(cNy-(i-line0.axis_size(0))-1,nz,nx) = std::conj(timestep.ft_phi(i,nz,nx));
  }

  nz = cNz/2;
  // write top right to top left in the same processor where possible
  for (int i = line0.axis_size(0)+(line1.axis_size(0)-1)/2+1; i < cNy; i++) {
    timestep.update(i,nz,nx);
    timestep.ft_phi(cNy-(i-line0.axis_size(0))-1,nz,nx) = std::conj(timestep.ft_phi(i,nz,nx));
  }

  // and receive extra data from the next processor to put in left

  nz = cNz/2;
  MPI_Recv(line2.data(), line2.size()*2,MPI_DOUBLE,mpi_size-id-1,2,
	   comm,MPI_STATUS_IGNORE);
  
  for (int i = 0; i < line2.axis_size(0); i++) {
    timestep.ft_phi(i,nz,nx) = line2(i,0);
  }

  MPI_Wait(&line_requests[0],MPI_STATUS_IGNORE);


  return;
  
}




