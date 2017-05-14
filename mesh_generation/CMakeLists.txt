cmake_minimum_required(VERSION 2.6)

PROJECT(MESHES)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(CMAKE_SOURCE_DIR)

add_executable(mesh_gen main.cc meshes.cc $<TARGET_OBJECTS:parser>)

configure_file(object.data ${CMAKE_CURRENT_BINARY_DIR}/object.data COPYONLY)
