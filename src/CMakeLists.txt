cmake_minimum_required(VERSION 2.6)
include(CheckLibraryExists)

project(PSR_EM Fortran C)

set(CMAKE_BUILD_TYPE RELEASE)

if(NOT UNIX)
  message(FATAL_ERROR "Sorry, only Unix systems supported at the moment.")
endif()

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/lib)

set( ${PROJECT_NAME}_MAJOR_VERSION 1 )
set( ${PROJECT_NAME}_MINOR_VERSION 2 )
set( ${PROJECT_NAME}_PATCH_LEVE    0 )

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

if(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  set(F_MOD_FLAG "-J${CMAKE_Fortran_MODULE_DIRECTORY}")
else(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  set(F_MOD_FLAG -module ${CMAKE_Fortran_MODULE_DIRECTORY})
endif()

##--- MACROS ---
##

macro(build_interface_module name source)
  add_custom_command(OUTPUT ${CMAKE_Fortran_MODULE_DIRECTORY}/${name}.mod
    COMMAND ${CMAKE_Fortran_COMPILER} ARGS "-c" ${F_MOD_FLAG} ${CMAKE_CURRENT_SOURCE_DIR}/${source}
    DEPENDS "${source}"
    COMMENT "Compiling ${name} external interface module"
  )

  add_custom_target(${name}_interface ALL
    DEPENDS ${CMAKE_Fortran_MODULE_DIRECTORY}/${name}.mod
  )
endmacro(build_interface_module)

##--- OPTIONS ---
##
option(WITH_ACML          "Enable ACML for random number generation" OFF )
option(WITH_RFIELD        "Enable random fields, requires fftw3"      ON )
option(WITH_CFITSIO       "Enable fits output"                       OFF )
option(WITH_OPENMP        "Enable OpenMP"                            OFF )
option(BUILD_SHARED_LIBS  "Build shared libraries"                   OFF )
option(STATIC_BUILD       "Generate static binaries"                 OFF )

## Should we link with static or shared libraries?
##
if(STATIC_BUILD)
  message(STATUS "Static build")
  if(B_LINK_SHARED)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(FFTW3_LIBRARIES "NOTFOUND" CACHE FILEPATH "" FORCE)
    set(CFITSIO_LIBRARY "NOTFOUND" CACHE FILEPATH "" FORCE)
    set(NETCDF_LIBRARY  "NOTFOUND" CACHE FILEPATH "" FORCE)
    set(NETCDFF_LIBRARY "NOTFOUND" CACHE FILEPATH "" FORCE)
    set(XRAYLIB_LIBRARY "NOTFOUND" CACHE FILEPATH "" FORCE)
    set(ACML_LIBRARY    "NOTFOUND" CACHE FILEPATH "" FORCE)
    set(ACML_MP_LIBRARY "NOTFOUND" CACHE FILEPATH "" FORCE)
    if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
      set(ACML_MV_LIBRARY "NOTFOUND" CACHE FILEPATH "" FORCE)
    endif()
    set(B_LINK_SHARED 0 CACHE INTERNAL "")
  endif()
else()
  if(NOT B_LINK_SHARED)
    set(FFTW3_LIBRARIES "NOTFOUND" CACHE FILEPATH "" FORCE)
    set(CFITSIO_LIBRARY "NOTFOUND" CACHE FILEPATH "" FORCE)
    set(NETCDF_LIBRARY  "NOTFOUND" CACHE FILEPATH "" FORCE)
    set(NETCDFF_LIBRARY "NOTFOUND" CACHE FILEPATH "" FORCE)
    set(XRAYLIB_LIBRARY "NOTFOUND" CACHE FILEPATH "" FORCE)
    set(ACML_LIBRARY    "NOTFOUND" CACHE FILEPATH "" FORCE)
    set(ACML_MP_LIBRARY "NOTFOUND" CACHE FILEPATH "" FORCE)
    if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
      set(ACML_MV_LIBRARY "NOTFOUND" CACHE FILEPATH "" FORCE)
    endif()
    set(B_LINK_SHARED 1 CACHE INTERNAL "")
  endif()
endif()

#--- FIND FFTW3 ---
#
find_library( FFTW3_LIBRARIES  NAMES fftw3   PATHS /usr/lib /usr/local/lib)


#--- FIND FITSIO ---
#
if(WITH_CFITSIO)
  set(CFITSIO_LIBRARY "NOTFOUND" CACHE FILEPATH "" FORCE)
  find_library( CFITSIO_LIBRARY  NAMES cfitsio PATHS /usr/lib )
  add_definitions("-DWITH_CFITSIO")
else()
  set(CFITSIO_LIBRARY "" CACHE FILEPATH "" FORCE)
endif()


#--- FIND NETCDF---
#
find_library( NETCDF_LIBRARY   NAMES netcdf        PATHS /usr/lib /usr/local/lib
                                                       "$ENV{NETCDF}/gfortran/lib")
find_library( NETCDFF_LIBRARY  NAMES netcdff       PATHS /usr/lib /usr/local/lib
                                                       "$ENV{NETCDF}/gfortran/lib")
find_path   ( NETCDF_MODULE          netcdf.mod    PATHS "/usr/include" "/usr/local/include"
                                                       "$ENV{NETCDF}/gfortran/include")
set(NETCDF_LIBRARIES ${NETCDF_LIBRARY} ${NETCDFF_LIBRARY})

if(NETCDF_MODULE)
  include_directories(${NETCDF_MODULE})
endif()


#--- FIND XRAYLIB ---
#
find_library( XRAYLIB_LIBRARY NAMES xrl     PATHS /usr/lib 
                                                  /usr/local/lib 
                                                  "$ENV{XRAYLIB_DIR}/lib" 
                                                  "$ENV{HOME}/.xraylib/lib")

find_path   ( XRAYLIB_MODULE  xraylib   PATHS /usr/include 
                                              /usr/local/include 
                                              "$ENV{XRAYLIB_DIR}/include" 
                                              "$ENV{XRAYLIB_DIR}/lib" 
                                              "$ENV{HOME}/.xraylib/include")

if(XRAYLIB_MODULE AND XRAYLIB_LIBRARY)
  include_directories(${XRAYLIB_MODULE})
  set(XRAYLIB_LIBRARIES ${XRAYLIB_LIBRARY})
else(XRAYLIB_MODULE AND XRAYLIB_LIBRARY)
  message(ERROR "xraylib not found")
endif(XRAYLIB_MODULE AND XRAYLIB_LIBRARY)


#--- ACML ---
#
if(WITH_ACML)
  add_definitions("-DWITH_ACML")

  if(${CMAKE_COMPILER_IS_GNUG77})
    if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
      find_library( ACML_LIBRARY    NAMES acml    PATHS /usr/lib /opt/acml4.2.0/gfortran64/lib )
      find_library( ACML_MP_LIBRARY NAMES acml_mp PATHS /usr/lib /opt/acml4.2.0/gfortran64_mp/lib )
      find_library( ACML_MV_LIBRARY NAMES acml_mv PATHS /usr/lib /opt/acml4.2.0/gfortran64_mp/lib )
    else(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
      find_library( ACML_LIBRARY    NAMES acml    PATHS /usr/lib /opt/acml4.2.0/gfortran32/lib )
      find_library( ACML_MP_LIBRARY NAMES acml_mp PATHS /usr/lib /opt/acml4.2.0/gfortran32_mp/lib )
    endif(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
  endif(${CMAKE_COMPILER_IS_GNUG77})
endif()



##--- SET COMPILATION FLAGS ---
##
if(${CMAKE_COMPILER_IS_GNUG77})
  add_definitions("-ffree-line-length-none")

  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O2" CACHE STRING "")
   
  if(WITH_OPENMP)
    add_definitions("-fopenmp")
  endif(WITH_OPENMP)

endif()

## ACML LINKING
##
if(WITH_ACML)
  if(WITH_OPENMP)
    set(ACML_LIBRARIES ${ACML_MP_LIBRARY} ${ACML_MV_LIBRARY})
  else()
    set(ACML_LIBRARIES ${ACML_LIBRARY} ${ACML_MV_LIBRARY})
  endif()
endif()

link_libraries(${NETCDF_LIBRARY})

##--- BUILD EXTERNAL MODULES ---
##
add_subdirectory(extern)

##--- BUILD LIBRARIES ---
##
add_subdirectory(lib)

##--- BUILD TESTS ---
##
add_subdirectory(tests)

##--- BUILD PROGRAMS ---
##
add_subdirectory(programs)

install(CODE "file(WRITE ${CMAKE_INSTALL_PREFIX}/xrinit.sh \"export XRROOT=${CMAKE_INSTALL_PREFIX}\nexport PATH=$PATH:${CMAKE_INSTALL_PREFIX}/bin\nexport LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CMAKE_INSTALL_PREFIX}/lib\n\")")

## Packaging
##

#set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "University of Helsinki Observatory PSR EM software.")
#set(CPACK_PACKAGE_VENDOR "Olli Wilkman <olli.wilkman@iki.fi>")
#set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/../README.md)

#set(CPACK_GENERATOR TGZ)
#set(CPACK_SOURCE_GENERATOR TGZ)

#set(CPACK_SOURCE_IGNORE_FILES "~$")

#include(CPack)
