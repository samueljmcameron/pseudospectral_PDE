

mpigxx -std=c++11 main.cpp fftw_mpi_3darray.cpp smatrix.cpp conjplane.cpp integrator.cpp randompll.cpp iovtk.cpp solutionparams.cpp griddata.cpp input.cpp timestep.cpp globalparams.cpp run.cpp -I/mnt/storage/software/libraries/intel/fftw-3.3.6-mpi/include -L/mnt/storage/software/libraries/intel/fftw-3.3.6-mpi/lib  -lfftw3_mpi -lfftw3 -lm -g -o main.mpi
#-I/mnt/storage/software/libraries/intel/fftw-3.3.6-mpi/include 
#-L/mnt/storage/software/libraries/intel/fftw-3.3.6-mpi/lib 

#mpic++ -std=c++1z main.cpp fftw_mpi_3darray.cpp smatrix.cpp conjplaneunequal.cpp conjlinesunequal.cpp timestep.cpp -lfftw3_mpi -lfftw3 -lm -g -o main
