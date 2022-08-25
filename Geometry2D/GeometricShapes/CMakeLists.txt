cmake_minimum_required(VERSION 3.23)
project(Geometry CXX)

set(CMAKE_C_STANDARD 20)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/geometry.cpp)
set(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/geometry.hpp)

add_library(GeomLib ${HEADERS} ${SOURCES})

add_executable(Geometry ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

target_link_libraries(Geometry PUBLIC GeomLib)
target_include_directories(Geometry PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
