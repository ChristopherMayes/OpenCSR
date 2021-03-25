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
          http://www.fftw.org/fftw-3.3.9.tar.gz
        URL_HASH
          MD5=50145bb68a8510b5d77605f11cadf8dc
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}
        CONFIGURE_COMMAND ""
        BUILD_COMMAND
            ${CMAKE_CURRENT_BINARY_DIR}/src/fftw3_external/configure --enable-threads --enable-openmp --with-pic --prefix=${STAGED_INSTALL_PREFIX} &&
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

      set(
        FFTW_DOUBLE_THREADS_LIB ${STAGED_INSTALL_PREFIX}/lib/libfftw3_threads.a
        CACHE PATH "Path to FFTW Threads"
        FORCE
      )

      set(
        FFTW_DOUBLE_OPENMP_LIB ${STAGED_INSTALL_PREFIX}/lib/libfftw3_omp.a
        CACHE PATH "Path to FFTW OpenMP"
        FORCE
      )

      # Libraries
      add_library(FFTW::Double INTERFACE IMPORTED)
      set_target_properties(FFTW::Double
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
              INTERFACE_LINK_LIBRARIES "${FFTW_DOUBLE_LIB}"
      )
      add_library(FFTW::DoubleThreads INTERFACE IMPORTED)
      set_target_properties(FFTW::DoubleThreads
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
              INTERFACE_LINK_LIBRARIES "${FFTW_DOUBLE_THREADS_LIB}"
      )
      add_library(FFTW::DoubleOpenMP INTERFACE IMPORTED)
      set_target_properties(FFTW::DoubleOpenMP
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
              INTERFACE_LINK_LIBRARIES "${FFTW_DOUBLE_OPENMP_LIB}"
      )
  endif()
  if(WIN32)
      ExternalProject_Add(fftw3_external
        URL
          http://www.fftw.org/fftw-3.3.9.tar.gz
        URL_HASH
          MD5=50145bb68a8510b5d77605f11cadf8dc
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

      set(
        FFTW_DOUBLE_THREADS_LIB ${STAGED_INSTALL_PREFIX}/lib/libfftw3_threads.lib
        CACHE PATH "Path to FFTW Threads"
        FORCE
      )

      set(
        FFTW_DOUBLE_THREADS_LIB ${STAGED_INSTALL_PREFIX}/lib/libfftw3_omp.lib
        CACHE PATH "Path to FFTW OpenMP"
        FORCE
      )

      # Libraries
      add_library(FFTW::Double INTERFACE IMPORTED)
      set_target_properties(FFTW::Double
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
              INTERFACE_LINK_LIBRARIES "${FFTW_DOUBLE_LIB}"
      )
      add_library(FFTW::DoubleThreads INTERFACE IMPORTED)
      set_target_properties(FFTW::DoubleThreads
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
              INTERFACE_LINK_LIBRARIES "${FFTW_DOUBLE_THREADS_LIB}"
      )
      add_library(FFTW::DoubleOpenMP INTERFACE IMPORTED)
      set_target_properties(FFTW::DoubleOpenMP
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
              INTERFACE_LINK_LIBRARIES "${FFTW_DOUBLE_OPENMP_LIB}"
      )
  endif()

  add_dependencies(opencsr fftw3_external)
  target_include_directories(opencsr PRIVATE ${STAGED_INSTALL_PREFIX}/include)
endif()
