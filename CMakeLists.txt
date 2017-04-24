cmake_minimum_required (VERSION 2.6)
project (stefan_enthalpy)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

set(${PROJECT_NAME}_sources stefan.cpp)

add_executable (${PROJECT_NAME} ${${PROJECT_NAME}_sources})
target_link_libraries (${PROJECT_NAME})

