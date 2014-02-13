ADD_EXECUTABLE(t_geometry t_geometry.f90)
TARGET_LINK_LIBRARIES(t_geometry core)
install(TARGETS t_geometry RUNTIME DESTINATION bin)
