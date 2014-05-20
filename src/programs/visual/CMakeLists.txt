ADD_EXECUTABLE(vScatter vScatter.f90)
ADD_EXECUTABLE(geomScatter geomScatter.f90)
ADD_EXECUTABLE(fit_2hg fit_2hg.f90)
ADD_EXECUTABLE(hemiScatter hemiScatter.f90)

TARGET_LINK_LIBRARIES(hemiScatter core emrandom simulation scatter rfield ${NETCDF_LIBRARIES} ${FFTW3_LIBRARIES} ${ACML_LIBRARIES} simulation core gomp)

TARGET_LINK_LIBRARIES(geomScatter core emrandom simulation scatter rfield ${NETCDF_LIBRARIES} ${FFTW3_LIBRARIES} ${ACML_LIBRARIES} simulation core gomp)

TARGET_LINK_LIBRARIES(vScatter core emrandom simulation scatter rfield ${NETCDF_LIBRARIES} ${FFTW3_LIBRARIES} ${ACML_LIBRARIES} simulation core gomp)

TARGET_LINK_LIBRARIES(fit_2hg core emrandom simulation scatter rfield ${NETCDF_LIBRARIES} ${FFTW3_LIBRARIES} ${ACML_LIBRARIES} simulation core gomp)

install(TARGETS geomScatter RUNTIME DESTINATION bin)
install(TARGETS vScatter RUNTIME DESTINATION bin)
install(TARGETS fit_2hg RUNTIME DESTINATION bin)
install(TARGETS hemiScatter RUNTIME DESTINATION bin)
