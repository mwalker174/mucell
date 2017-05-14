cmake_minimum_required(VERSION 2.6)

PROJECT(PIO)

include_directories(CMAKE_SOURCE_DIR)

add_library(pio OBJECT PIO.cc)
add_executable(test_pio main.cc $<TARGET_OBJECTS:pio> $<TARGET_OBJECTS:helpers> $<TARGET_OBJECTS:parser>)

