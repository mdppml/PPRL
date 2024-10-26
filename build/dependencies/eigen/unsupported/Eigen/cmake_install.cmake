# Install script for directory: /Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Library/Developer/CommandLineTools/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/AdolcForward"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/AlignedVector3"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/ArpackSupport"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/AutoDiff"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/BVH"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/EulerAngles"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/FFT"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/IterativeSolvers"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/KroneckerProduct"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/LevenbergMarquardt"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/MatrixFunctions"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/MoreVectorization"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/MPRealSupport"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/NonLinearOptimization"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/NumericalDiff"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/OpenGLSupport"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/Polynomials"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/Skyline"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/SparseExtra"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/SpecialFunctions"
    "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/Splines"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "/Users/sesame/Desktop/CECILIA/dependencies/eigen/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/sesame/Desktop/CECILIA/build/dependencies/eigen/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

