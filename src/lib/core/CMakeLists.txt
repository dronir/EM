ADD_LIBRARY(core core.f90 util.F90 geometry.f90 sphere.f90 spectrum.f90)
TARGET_LINK_LIBRARIES(core ${CFITSIO_LIBRARY})
install(TARGETS core LIBRARY DESTINATION lib ARCHIVE DESTINATION lib/uhoem)
