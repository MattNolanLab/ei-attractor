find_package(PythonLibs REQUIRED)
find_package(Numpy REQUIRED)
find_package(Blitz REQUIRED)

include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${NUMPY_INCLUDE_DIRS})
include_directories(${Blitz_INCLUDE_DIR})

link_directories(
    ${PYTHON_LIBRARY_DIRS}
    ${Boost_LIBRARY_DIRS}
    ${Blitz_LIBRARY}
)

SET(CMAKE_CXX_FLAGS "-Wall -g -DBZ_DEBUG ${CMAKE_CXX_FLAGS}")

add_library(_signalmodule SHARED _signalmodule.cpp)

if (APPLE)
    set_target_properties(
        _signalmodule PROPERTIES
        PREFIX      ""
        OUTPUT_NAME "_signal"
        SUFFIX      ".so"
    )
endif()

target_link_libraries(_signalmodule ${PYTHON_LIBRARIES} ${Blitz_LIBRARIES})

install(
    TARGETS _signalmodule
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/..)

