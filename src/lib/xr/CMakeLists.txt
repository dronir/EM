ADD_LIBRARY(xr material.f90 bfdf.f90)
TARGET_LINK_LIBRARIES(xr core emrandom)
install(TARGETS xr LIBRARY DESTINATION lib
                   ARCHIVE DESTINATION lib/uhoem)
