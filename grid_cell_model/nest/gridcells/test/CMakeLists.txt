########################################################################
# CMake build script for grid cell NEST module tests.
project(gridcells_test CXX)
cmake_minimum_required(VERSION 2.6.2)

set (gtest_SOURCE_DIR "${PROJECT_SOURCE_DIR}/googletest/gtest-1.7.0")
set (gridcells_SOURCE_DIR "${PROJECT_SOURCE_DIR}/..")

include_directories(
  ${PROJECT_SOURCE_DIR}
  ${gridcells_SOURCE_DIR}
  ${gtest_SOURCE_DIR}/include
  ${gtest_SOURCE_DIR})

link_directories(${gtest_SOURCE_DIR})

add_executable(iaf_gridcells_unittest iaf_gridcells_unittest.cpp)
target_link_libraries(iaf_gridcells_unittest gtest)
