find_package(Boost REQUIRED COMPONENTS Python)
find_package(PythonLibs REQUIRED)

include_directories(${Boost_INCLUDE_DIR})
include_directories(${PYTHON_INCLUDE_DIRS})

add_library(testmodule SHARED testmodule.c)

if (APPLE)
    set_target_properties(
        testmodule PROPERTIES
        PREFIX      ""
        OUTPUT_NAME "_test"
        SUFFIX      ".so"
    )
endif()

target_link_libraries(testmodule ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

install(
    TARGETS testmodule
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})

