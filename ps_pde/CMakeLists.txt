
add_executable(ps_pde-bin main.cpp run.cpp)

target_link_libraries(ps_pde-bin PRIVATE MPI::MPI_CXX PkgConfig::FFTW ${FFTW_MPI_LIB} ps_pde)


include(GNUInstallDirs)

install(TARGETS ps_pde-bin DESTINATION bin)
