cmake_minimum_required (VERSION 2.6)

project (stefan_enthalpy)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -mtune=native -funroll-loops -ffast-math -ftree-vectorize -ftree-vectorizer-verbose=0")


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)


#set(CMAKE_COLOR_MAKEFILE ON)
#set(CMAKE_VERBOSE_MAKEFILE ON)


#FIND_PACKAGE(OpenMP)
#if (OPENMP_FOUND)
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
#  message("OpenMP is found\n")
#endif (OPENMP_FOUND)


set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/bin)
MAKE_DIRECTORY(${CMAKE_SOURCE_DIR}/bin)

set(BINDIR ${CMAKE_SOURCE_DIR}/bin)

MAKE_DIRECTORY(${BINDIR})

if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall")
    set(BinaryFileName ${PROJECT_NAME}Debug)
endif()

if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
    set(BinaryFileName ${PROJECT_NAME})
endif()


set(sources
  ./src/stefan/stefan.cpp
  ./src/stefan/material_info.h
  ./src/stefan/mesh.h
  ./src/stefan/solver.h
  ./src/stefan/settings_parser.h
)

set(TinyXmlSrc
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinystr.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinyxml.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinyxmlerror.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinyxmlparser.cpp
)

add_executable (${BinaryFileName} ${sources}  ${TinyXmlSrc})
#target_link_libraries (${BinaryFileName} ${GSL_LIBRARIES})

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${CMAKE_SOURCE_DIR}/3rdparty/gsl -L${CMAKE_SOURCE_DIR}/3rdparty/gsl/ -lgsl -lgslcblas")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${CMAKE_SOURCE_DIR}/3rdparty/gsl/ -lgsl -lgslcblas")


