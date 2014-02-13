ADD_LIBRARY(emrandom stdtypes.f90 igdev.f90 mtprng.f90 random.F90)
TARGET_LINK_LIBRARIES(emrandom core)
install(TARGETS emrandom LIBRARY DESTINATION lib
                         ARCHIVE DESTINATION lib/uhoem)