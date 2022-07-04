file(GLOB LIB_SOURCES [^.]*.cpp)


add_library(ps_pde ${LIB_SOURCES})




target_include_directories(ps_pde INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(ps_pde PRIVATE MPI::MPI_CXX PkgConfig::FFTW ${FFTW_MPI_LIB})

set(LIB_HEADERS conjplane.hpp fftw_mpi_3darray.hpp globalparams.hpp griddata.hpp
    input.hpp integrator.hpp iovtk.hpp randompll.hpp smatrix.hpp solutionparams.hpp timestep.hpp)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/includes/ps_pde)

message(STATUS "LIB_HEADERS: ${LIB_HEADERS}")

foreach(_HEADER ${LIB_HEADERS})
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/includes/ps_pde/${_HEADER} COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${_HEADER} ${CMAKE_CURRENT_BINARY_DIR}/includes/ps_pde/${_HEADER} DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_HEADER})
    add_custom_target(${_HEADER} DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/includes/ps_pde/${_HEADER})
    add_dependencies(ps_pde ${_HEADER})
endforeach()
