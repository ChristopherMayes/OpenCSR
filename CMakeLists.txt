# Minimum cmake version
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# Project name
project(OpenCSR LANGUAGES Fortran)

cmake_policy(SET CMP0074 NEW)

# Uncomment for VERBOSE Makefiles. Useful for debugging.
#set(CMAKE_VERBOSE_MAKEFILE ON)

# Useful debug for Find commands
#set(CMAKE_FIND_DEBUG_MODE ON)

# Uncomment for install messages
#set(CMAKE_INSTALL_MESSAGE ALWAYS)

# we default to Release build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
    # Uncomment for Debug Build
    #set(CMAKE_BUILD_TYPE Debug)
endif()

# Folder for helper modules (Find*.cmake)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# this is where we will place the Fortran module files
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)
# this is where they will be installed
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)


#############################
# OpenCSR Library
#############################

# Source files for OpenCSR
list(APPEND _sources
    src/csr3d_mod.f90
    src/elliptic_integral_mod.f90    
    src/fft_mod.f90
)

# opencsr library
add_library(opencsr ${_sources})

# Not needed, will use internal .f90 for elliptic integrals
# Link to GSL
#find_package(GSL REQUIRED)
#target_link_libraries(opencsr PUBLIC GSL::gsl GSL::gslcblas)

# Link to FFTW
find_package(FFTW REQUIRED)
target_link_libraries(opencsr PUBLIC FFTW::Float)

install(TARGETS opencsr LIBRARY)

#############################
# OpenCSR Test Binary
#############################

# Test executable
add_executable(csr3d-test tests/test_csr3d.f90)

target_link_libraries(csr3d-test opencsr)

install(TARGETS csr3d-test RUNTIME)


# Install Dirs configuration for Unix and Windows
if(UNIX)
  # Use /lib /bin /include /share folder structure
  include(GNUInstallDirs)
else()
  if (WIN32)
    set(${CMAKE_INSTALL_LIBDIR} "lib")
    set(${CMAKE_INSTALL_DATADIR} "share")
    set(${CMAKE_INSTALL_INCLUDEDIR} "include")
    set(${CMAKE_INSTALL_BINDIR} "bin")
    message(STATUS "Setting installation destination on Windows to: ${CMAKE_INSTALL_PREFIX}")
  else()
    message(FATAL_ERROR "System not UNIX nor WIN32 - not implemented yet")
  endif()
endif()
