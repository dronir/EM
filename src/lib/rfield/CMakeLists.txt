ADD_LIBRARY(rfield rfield.F90)
TARGET_LINK_LIBRARIES(rfield core emrandom)
install(TARGETS rfield LIBRARY DESTINATION lib
                       ARCHIVE DESTINATION lib/uhoem)
