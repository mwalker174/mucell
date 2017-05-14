cmake_minimum_required(VERSION 2.6)

PROJECT(IonicModels)

include_directories(CMAKE_SOURCE_DIR)

add_subdirectory(stoch/stoch_spark)

add_library(tt06 OBJECT TT06_RRG.cc stoch.cc ModelFactory.cc)
add_executable(tt06_test main.cc TT06_RRG.cc)

