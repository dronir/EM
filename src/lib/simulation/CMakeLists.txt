ADD_LIBRARY(simulation gth_hemisphere.f90 sampler.f90 cnt_grid3d.f90 medium.f90 trace.f90)
TARGET_LINK_LIBRARIES(simulation core emrandom rfield ${NETCDF_LIBRARIES})
install(TARGETS simulation LIBRARY DESTINATION lib
                           ARCHIVE DESTINATION lib/uhoem)

ADD_LIBRARY(trace_xr trace_xr.f90)
TARGET_LINK_LIBRARIES(trace_xr core emrandom rfield xr simulation ${NETCDF_LIBRARIES})
install(TARGETS trace_xr LIBRARY DESTINATION lib
                           ARCHIVE DESTINATION lib/uhoem)
