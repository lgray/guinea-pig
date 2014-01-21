cmake_minimum_required(VERSION 2.6)

project(guinea-pig)

# add source dir as a place for CMake modules (e.g. FindGSL.cmake)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

include(guinea_definitions)

include(compilers)

# With this true we automatically include current src/binary directory
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

add_subdirectory(src)

enable_testing()
include(CTest)

add_subdirectory(testing)
