cmake_minimum_required(VERSION 2.6)

PROJECT(PARSER)

include_directories(CMAKE_SOURCE_DIR)

#if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})  
#else()
#  message(SEND_ERROR "Please set Boost folders manually")      
#endif()

add_library(parser OBJECT spirit_wrapper.cc)

