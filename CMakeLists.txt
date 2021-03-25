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

# Temporary path for FFTW in case we need to build it
set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage)
message(STATUS "${PROJECT_NAME} staged install: ${STAGED_INSTALL_PREFIX}")

# this is where we will place the Fortran module files
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)
# this is where they will be installed
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)

#############################
# Fortran Flags
#############################
# More info on dialect flags here: https://gcc.gnu.org/onlinedocs/gfortran/Fortran-Dialect-Options.html#Fortran-Dialect-Options
# More info on code generation flags here: https://gcc.gnu.org/onlinedocs/gfortran/Code-Gen-Options.html#Code-Gen-Options
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-Wconversion -ffree-form -ffree-line-length-none -cpp -std=f2008 -fimplicit-none")
    set(bounds "-fcheck=bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")


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


#############################
# External Libraries
#############################
# Force FFTW to look for STATIC libraries that we need
set(FFTW_USE_STATIC_LIBS ON)

# External Libraries
add_subdirectory(ext_lib)

# Link to FFTW
target_link_libraries(opencsr PUBLIC ${FFTW_DOUBLE_LIB})

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
