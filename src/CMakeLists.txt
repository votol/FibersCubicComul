cmake_minimum_required(VERSION 3.0)

get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)

set (SOURCES main.cpp
             NetCdfWriter.cpp
             AbOA.cpp
             APolA.cpp
             AbComul.cpp
             eqManager.cpp
             ParametersInit.cpp
             ModelBuilder.cpp
             OutputWrapper.cpp
             OutputTime.cpp
             OutputInterface.h
             NetCdfWriter.h
             AbOA.h
             APolA.h
             AbComul.h
             eqManager.h
             ParameterDefines.h
             ParametersInit.h
             ModelBuilder.h
             OutputWrapper.h
             OutputTime.h
             )

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(TAR ${PROJECT_NAME}_c++)

add_custom_command(
    OUTPUT  schema.h
    COMMAND ${PARENT_DIR}/Workflow/scripts/generateSchema.py -d ${PARENT_DIR}/description.yaml -f ${CMAKE_CURRENT_BINARY_DIR}/schema.h
    DEPENDS ${PARENT_DIR}/description.yaml
)

add_custom_command(
        OUTPUT ${EXECUTABLE_OUTPUT_PATH}/CalcProj
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/src/CalcProj
                ${EXECUTABLE_OUTPUT_PATH}/CalcProj
        DEPENDS ${CMAKE_SOURCE_DIR}/src/CalcProj
        )

get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
include_directories (${CMAKE_CURRENT_BINARY_DIR})
include_directories (~/include)
include_directories (../CLDE/include)


add_subdirectory(../CLDE/src ${CMAKE_CURRENT_BINARY_DIR}/clde)
add_executable(${TAR}
    ${SOURCES}
    ${CMAKE_CURRENT_BINARY_DIR}/schema.h
    ${EXECUTABLE_OUTPUT_PATH}/CalcProj
)

set_property(TARGET ${TAR} PROPERTY CXX_STANDARD 17)
target_link_libraries (${TAR} ~/lib/libyaml-cpp.so
                                        netcdf_c++4
                                        netcdf
                                        clde)



