cmake_minimum_required(VERSION 2.6)
project(GridCells)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(NUMPY_BOOST_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/numpy_boost/include)

#add_subdirectory(python/testmodule)
add_subdirectory(grid_cell_model/analysis/src)
