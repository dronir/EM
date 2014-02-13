ADD_LIBRARY(scatter brdf.f90)
TARGET_LINK_LIBRARIES(scatter core)
install(TARGETS scatter LIBRARY DESTINATION lib
                        ARCHIVE DESTINATION lib/uhoem)
