ADD_EXECUTABLE(xrfip xrfip.f90 xrfip_m_analytical.f90 xrfip_m_firstorder.f90 xrfip_m_mcpeel.f90 xrfip_m_mcforce.f90)
TARGET_LINK_LIBRARIES(xrfip core emrandom xr simulation rfield ${XRAYLIB_LIBRARIES} ${ACML_LIBRARIES} ${FFTW3_LIBRARIES} ${NETCDF_LIBRARIES} ${CFITSIO_LIBRARY})

ADD_EXECUTABLE(xrfpm xrfpm.f90 xrfpm_m_analytical.f90 xrfpm_m_firstorder.f90 xrfpm_m_mcpeel.f90)
TARGET_LINK_LIBRARIES(xrfpm core emrandom xr simulation rfield ${XRAYLIB_LIBRARIES} ${ACML_LIBRARIES} ${FFTW3_LIBRARIES}  ${NETCDF_LIBRARIES} ${CFITSIO_LIBRARY})

#ADD_EXECUTABLE(xrfipsd xrfipsd.f90)
#TARGET_LINK_LIBRARIES(xrfipsd core emrandom xr simulation rfield ${XRAYLIB_LIBRARIES} ${ACML_LIBRARIES} ${FFTW3_LIBRARIES}  ${NETCDF_LIBRARIES} ${CFITSIO_LIBRARY} gomp)

#ADD_EXECUTABLE(xrfiphs xrfiphs.f90)
#TARGET_LINK_LIBRARIES(xrfiphs core emrandom xr simulation rfield ${XRAYLIB_LIBRARIES} ${ACML_LIBRARIES} ${FFTW3_LIBRARIES}  ${NETCDF_LIBRARIES} ${CFITSIO_LIBRARY} gomp)

#ADD_EXECUTABLE(xrfpmhs xrfpmhs.f90)
#TARGET_LINK_LIBRARIES(xrfpmhs core emrandom xr simulation rfield ${ACML_LIBRARIES} ${XRAYLIB_LIBRARIES} ${FFTW3_LIBRARIES}  ${NETCDF_LIBRARIES} ${CFITSIO_LIBRARY} gomp)

install(TARGETS xrfip xrfpm RUNTIME DESTINATION bin)
