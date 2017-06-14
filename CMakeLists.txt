cmake_minimum_required (VERSION 2.6)

project (stefan_enthalpy)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mtune=native -funroll-loops -ffast-math -ftree-vectorize -ftree-vectorizer-verbose=0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=native -funroll-loops -ffast-math -ftree-vectorize -ftree-vectorizer-verbose=0")

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_BUILD_TYPE)

string(TOUPPER ${CMAKE_SYSTEM} BuildPlatform)
message("Generating makefiles for: " ${BuildPlatform} "\n")
message("Build type: " ${CMAKE_BUILD_TYPE} "\n")

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)

#FIND_PACKAGE(OpenMP)
#if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  message("OpenMP is found\n")
#endif (OPENMP_FOUND)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/bin)
MAKE_DIRECTORY(${CMAKE_SOURCE_DIR}/bin)

set(BINDIR ${CMAKE_SOURCE_DIR}/bin)

MAKE_DIRECTORY(${BINDIR})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BINDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${BINDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${BINDIR})

IF (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    SET(BinaryFileName ${PROJECT_NAME}Debug)
ENDIF()

IF (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    SET(BinaryFileName ${PROJECT_NAME})
ENDIF()

set(math
  ./src/math/VectorMaths.h
  ./src/math/Vector3.h
  ./src/math/Vector2.h
  ./src/math/Util.h
  ./src/math/Tensor3.h
  ./src/math/Tensor2.h
  ./src/math/Spaces.h
  ./src/math/Random.h
  ./src/math/AABB3.h
  ./src/math/AABB2.h
)

set(utils
  ./src/utils/Utils.h
  ./src/utils/Base64.h
)

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
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinystr.h
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinyxml.h
)

add_executable (${BinaryFileName} ${sources} ${math} ${utils} ${TinyXmlSrc})
target_link_libraries (${BinaryFileName})

