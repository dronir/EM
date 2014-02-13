build_interface_module(typesizes typeSizes.f90)
build_interface_module(netcdf    netcdf4.f90)

add_dependencies(netcdf_interface typesizes_interface)
