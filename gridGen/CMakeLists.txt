enable_language(Fortran)

set(GRIDGEN_SOURCES
  modules.f
  buildGrid.f
  clust.f
  partSphere.f90
  findfact.f90)

add_executable(buildGrid ${GRIDGEN_SOURCES})

target_link_libraries(buildGrid PRIVATE tioga)
target_link_libraries(buildGrid PRIVATE MPI::MPI_Fortran)

install(TARGETS buildGrid
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
