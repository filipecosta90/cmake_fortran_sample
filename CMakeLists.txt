cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(HelloWorld)
enable_language(Fortran)


# Define the executable name
SET(FOOEXE project)

# we default to Release build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

############################################################
# Fourier transform
find_package(FFTW3 COMPONENTS double single threads REQUIRED)
############################################################
include_directories(${FFTW3_INCLUDE_DIR})

# the executable is built from 3 source files
add_executable(
  ${FOOEXE}
  main.F90
  )


TARGET_LINK_LIBRARIES(${FOOEXE} ${FFTW3_LIBRARIES})
