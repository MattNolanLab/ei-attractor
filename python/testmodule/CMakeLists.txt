find_package(Boost REQUIRED COMPONENTS Python)
find_package(PythonLibs REQUIRED)
find_package(Numpy REQUIRED)

include_directories(${Boost_INCLUDE_DIR})
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${NUMPY_INCLUDE_DIRS})
include_directories(${NUMPY_BOOST_INCLUDE_DIR})

link_directories(
    ${Boost_LIBRARY_DIRS}
    ${PYTHON_LIBRARY_DIRS}
)

SET(CMAKE_CXX_FLAGS "-Wall -g ${CMAKE_CXX_FLAGS}")

add_library(testmodule SHARED testmodule.cpp)
add_executable(test_numpy_boost test_numpy_boost.cpp)

if (APPLE)
    set_target_properties(
        testmodule PROPERTIES
        PREFIX      ""
        OUTPUT_NAME "_test"
        SUFFIX      ".so"
    )
endif()

target_link_libraries(testmodule ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
target_link_libraries(test_numpy_boost ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

install(
    TARGETS testmodule
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})

