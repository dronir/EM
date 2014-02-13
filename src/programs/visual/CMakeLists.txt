ADD_EXECUTABLE(vScatter vScatter.f90)
ADD_EXECUTABLE(fit_2hg fit_2hg.f90)

TARGET_LINK_LIBRARIES(vScatter core emrandom xr simulation scatter rfield ${NETCDF_LIBRARIES} ${XRAYLIB_LIBRARIES} ${FFTW3_LIBRARIES} ${ACML_LIBRARIES} simulation core gomp)

TARGET_LINK_LIBRARIES(fit_2hg core emrandom xr simulation scatter rfield ${NETCDF_LIBRARIES} ${XRAYLIB_LIBRARIES} ${FFTW3_LIBRARIES} ${ACML_LIBRARIES} simulation core gomp)

install(TARGETS vScatter RUNTIME DESTINATION bin)
install(TARGETS fit_2hg RUNTIME DESTINATION bin)
