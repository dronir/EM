ADD_LIBRARY(simulation gth_hemisphere.f90 sampler.f90 cnt_grid3d.f90 medium.f90 trace.f90)
TARGET_LINK_LIBRARIES(simulation core emrandom rfield xr ${NETCDF_LIBRARIES})
install(TARGETS simulation LIBRARY DESTINATION lib
                           ARCHIVE DESTINATION lib/uhoem)
