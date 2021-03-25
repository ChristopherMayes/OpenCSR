find_package(FFTW QUIET)

if(FFTW_FOUND)
  message(STATUS "Found FFTW3: ${FFTW_DOUBLE_LIB}")
  target_include_directories(opencsr PRIVATE ${FFTW_INCLUDE_DIRS})
else()
  message(STATUS "Suitable FFTW3 could not be located. Downloading and building!")
  include(ExternalProject)

  if(UNIX)
      ExternalProject_Add(fftw3_external
        URL
          http://www.fftw.org/fftw-3.3.8.tar.gz
        URL_HASH
          MD5=8aac833c943d8e90d51b697b27d4384d
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}
        CONFIGURE_COMMAND ""
        BUILD_COMMAND
            ${CMAKE_CURRENT_BINARY_DIR}/src/fftw3_external/configure --with-pic --prefix=${STAGED_INSTALL_PREFIX} &&
            make -j8 &&
            make install
        INSTALL_COMMAND ""
        PREFIX=${CMAKE_CURRENT_BINARY_DIR}
        )

      include(GNUInstallDirs)

      set(
        FFTW_ROOT ${STAGED_INSTALL_PREFIX}
        CACHE PATH "Path to internally built FFTW3Config.cmake"
        FORCE
      )

      set(
        FFTW_DOUBLE_LIB ${STAGED_INSTALL_PREFIX}/lib/libfftw3.a
        CACHE PATH "Path to FFTW"
        FORCE
      )

      # Libraries
      add_library(fftw3 STATIC IMPORTED)
      set_property(TARGET fftw3 PROPERTY IMPORTED_LOCATION ${STAGED_INSTALL_PREFIX}/lib/libfftw3.a)
  endif()
  if(WIN32)
      ExternalProject_Add(fftw3_external
        URL
          http://www.fftw.org/fftw-3.3.8.tar.gz
        URL_HASH
          MD5=8aac833c943d8e90d51b697b27d4384d
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        COMMAND cmake ${CMAKE_CURRENT_BINARY_DIR}/src/fftw3_external/ -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTS=OFF -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
        COMMAND cmake --build . --config Release --target install
        INSTALL_COMMAND ""
        )

      include(GNUInstallDirs)


      set(
        FFTW_ROOT ${STAGED_INSTALL_PREFIX}
        CACHE PATH "Path to internally built FFTW3Config.cmake"
        FORCE
      )

      set(
        FFTW_DOUBLE_LIB ${STAGED_INSTALL_PREFIX}/lib/fftw3.lib
        CACHE PATH "Path to FFTW"
        FORCE
      )

      # Libraries
      add_library(fftw3 STATIC IMPORTED)
      set_property(TARGET fftw3 PROPERTY IMPORTED_LOCATION ${STAGED_INSTALL_PREFIX}/lib/fftw3.lib)
  endif()

  add_dependencies(fftw3 fftw3_external)
  add_dependencies(opencsr fftw3_external)
  target_include_directories(opencsr PRIVATE ${STAGED_INSTALL_PREFIX}/include)
endif()
