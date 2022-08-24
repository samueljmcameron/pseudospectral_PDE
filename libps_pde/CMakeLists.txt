file(GLOB LIB_SOURCES [^.]*.cpp)


add_library(ps_pde STATIC ${LIB_SOURCES})



target_include_directories(ps_pde PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(ps_pde PUBLIC MPI::MPI_CXX PkgConfig::FFTW ${FFTW_MPI_LIB})


set_target_properties(ps_pde PROPERTIES PUBLIC_HEADER "conjplane.hpp;fftw_mpi_3darray.hpp;globalparams.hpp;griddata.hpp;input.hpp;integrator.hpp;iovtk.hpp;randompll.hpp;smatrix.hpp;solutionparams.hpp;timestep.hpp")


include(GNUInstallDirs)

install(TARGETS ps_pde LIBRARY DESTINATION lib PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ps_pde)


# generate pc file for pkg-config

set(target1 ps_pde)
configure_file(ps_pde.pc.in ${CMAKE_CURRENT_BINARY_DIR}/pkgconfig/ps_pde.pc @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pkgconfig/ps_pde.pc DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig)
