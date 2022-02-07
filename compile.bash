

mpic++ -std=c++1z main.cpp fftw_mpi_3darray.cpp smatrix.cpp conjplane.cpp integrator.cpp randompll.cpp iovtk.cpp solutionparams.cpp griddata.cpp input.cpp timestep.cpp -lfftw3_mpi -lfftw3 -lm -g -o main

#mpic++ -std=c++1z main.cpp fftw_mpi_3darray.cpp smatrix.cpp conjplaneunequal.cpp conjlinesunequal.cpp timestep.cpp -lfftw3_mpi -lfftw3 -lm -g -o main
