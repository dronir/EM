ADD_EXECUTABLE(medgen medGen.f90)
TARGET_LINK_LIBRARIES(medgen core emrandom xr simulation rfield ${XRAYLIB_LIBRARIES} ${ACML_LIBRARIES} ${FFTW3_LIBRARIES}  ${NETCDF_LIBRARIES} ${CFITSIO_LIBRARY} gomp)

#ADD_EXECUTABLE(medarf medApplyRandomField.f90)
#ADD_DEPENDENCIES(medarf core emrandom rfield)
#TARGET_LINK_LIBRARIES(medarf core cfitsio emrandom core simulation rfield netcdf cfitsio ${FFTW3_LIBRARIES})

install(TARGETS medgen RUNTIME DESTINATION bin)