find_package(PythonLibrary REQUIRED)
find_package(Numpy REQUIRED)
find_package(Blitz REQUIRED)

include_directories(${PYTHON_INCLUDE_PATH})
include_directories(${NUMPY_INCLUDE_DIRS})
include_directories(${Blitz_INCLUDE_DIR})

link_directories(
    ${PYTHON_LIBRARY_DIRS}
    ${Boost_LIBRARY_DIRS}
    ${Blitz_LIBRARY}
)

#SET(CMAKE_CXX_FLAGS "-Wall -g -O3 -DBZ_DEBUG ${CMAKE_CXX_FLAGS}")
SET(CMAKE_CXX_FLAGS "-Wall -O3 ${CMAKE_CXX_FLAGS}")

add_library(_signalmodule SHARED _signalmodule.cpp)
add_library(_spikesmodule SHARED _spikesmodule.cpp)

set_target_properties(
    _signalmodule _spikesmodule PROPERTIES
    PREFIX      ""
    OUTPUT_NAME "_signal"
)

set_target_properties(
    _signalmodule PROPERTIES
    OUTPUT_NAME "_signal"
)

set_target_properties(
    _spikesmodule PROPERTIES
    OUTPUT_NAME "_spikes"
)
if (APPLE)
    set_target_properties(
        _signalmodule _spikesmodule PROPERTIES
        SUFFIX      ".so"
    )
endif()

target_link_libraries(_signalmodule ${PYTHON_LIBRARIES} ${Blitz_LIBRARIES})
target_link_libraries(_spikesmodule ${PYTHON_LIBRARIES} ${Blitz_LIBRARIES})

install(
    TARGETS _signalmodule _spikesmodule
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/..)

