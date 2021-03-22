# Minimum cmake version
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# Project name
project(OpenCSR LANGUAGES Fortran)

cmake_policy(SET CMP0074 NEW)

# Folder for helper modules (Find*.cmake)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Uncomment for VERBOSE Makefiles. Useful for debugging.
#set(CMAKE_VERBOSE_MAKEFILE ON)

# Useful debug for Find commands
#set(CMAKE_FIND_DEBUG_MODE ON)

# Uncomment for install messages
#set(CMAKE_INSTALL_MESSAGE ALWAYS)

# Uncomment for Debug Build
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

# Source files for OpenCSR
list(APPEND _sources
    src/csr3d_mod.f90
)

# opencsr library
add_library(opencsr ${_sources})

install(TARGETS opencsr LIBRARY)


# Link to GSL
find_package(GSL REQUIRED)
target_link_libraries(opencsr PUBLIC GSL::gsl GSL::gslcblas)


# Link to FFTW
find_package(FFTW REQUIRED)
target_link_libraries(opencsr PUBLIC FFTW::Float)

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
