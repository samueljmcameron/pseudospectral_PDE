
add_executable(ps_pde-bin main.cpp run.cpp)


target_include_directories(ps_pde-bin PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(ps_pde-bin PRIVATE MPI::MPI_CXX PkgConfig::FFTW ${FFTW_MPI_LIB} ps_pde)

