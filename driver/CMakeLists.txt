enable_language(Fortran)
find_package(MPI REQUIRED)

set(TIOGA_EXE_SOURCES
  modules.f90
  gridutils.f90
  poisson_utils.f90
)

if (TIOGA_ENABLE_ARBORX)
  list(APPEND TIOGA_EXE_SOURCES
      kokkos.f90
      kokkos.C
  )
endif()

add_library(tiogadriver ${TIOGA_EXE_SOURCES})

if (TIOGA_ENABLE_ARBORX)
  target_link_libraries(tiogadriver PRIVATE ArborX::ArborX)
endif()

target_link_libraries(tiogadriver PUBLIC tioga)
target_link_libraries(tiogadriver PUBLIC MPI::MPI_Fortran)

add_executable(tioga.exe testTioga.F90)
# Using #ifdef requires preprocessing
target_compile_options(tioga.exe PRIVATE "-cpp")

target_link_libraries(tioga.exe tiogadriver)

add_executable(tioga_read.exe testTioga_dataRead.f90)
target_link_libraries(tioga_read.exe tiogadriver)

add_executable(poisson_mms.exe poisson_mms.f90)
target_link_libraries(poisson_mms.exe tiogadriver)

set(GMSH_SOURCES
  gmsh_io.C
)

add_library(gmsh_lib ${GMSH_SOURCES})
add_executable(tioga_gmsh.exe gmsh_test.C)
target_link_libraries(tioga_gmsh.exe gmsh_lib)

add_executable(tioga_abuttingmesh.exe testTioga_abuttingMesh.f90)
target_link_libraries(tioga_abuttingmesh.exe tiogadriver gmsh_lib)

#add_executable(gpu_test.exe)
#target_sources(gpu_test.exe PRIVATE
#  gpu_test.C)
#target_link_libraries(gpu_test.exe tioga)
#if (TIOGA_ENABLE_CUDA)
#  get_target_property(GPUTEST_CPP_SRC gpu_test.exe SOURCES)
#  list(FILTER GPUTEST_CPP_SRC INCLUDE REGEX "\\.C")
#  set_source_files_properties(${GPUTEST_CPP_SRC} PROPERTIES LANGUAGE CUDA)
#endif()
#
#if (TIOGA_ENABLE_HIP)
#  separate_arguments(TIOGA_HIP_ARCH_FLAGS)
#  target_compile_options(gpu_test.exe PRIVATE
#    $<$<COMPILE_LANGUAGE:CXX>:${TIOGA_HIP_ARCH_FLAGS}>)
#  target_compile_options(gpu_test.exe PRIVATE
#    $<$<COMPILE_LANGUAGE:CXX>:-x hip>)
#endif()

install(TARGETS tiogadriver tioga.exe tioga_read.exe poisson_mms.exe
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
