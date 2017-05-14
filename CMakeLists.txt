cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_COMPILER "mpicxx")

PROJECT(RYR)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(CMAKE_SOURCE_DIR)

#list(APPEND CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")
list(APPEND CMAKE_CXX_FLAGS "-fopenmp ${CMAKE_CXX_FLAGS}")

find_package(Boost 1.36.0 COMPONENTS regex filesystem system)
#find_package(VTK REQUIRED NO_MODULE)

add_subdirectory(models)
add_subdirectory(parser)
add_subdirectory(helpers)
add_subdirectory(pio)
#add_subdirectory(unit_tests)
add_subdirectory(mesh_generation)
#add_subdirectory(types_and_constants)

SET(SCOTCH_PATH "" CACHE PATH "Path to Scotch header")
SET(SCOTCH_LIB_PATH "" CACHE PATH "Path to Scotch lib")

include_directories("${SCOTCH_PATH}" "/home-3/mwalke49@jhu.edu/scratch/apps/hypre-2.10.0b/build/include")
LINK_DIRECTORIES("${SCOTCH_LIB_PATH}" "/home-3/mwalke49@jhu.edu/scratch/apps/hypre-2.10.0b/build/lib")

#add_library(test Mesh.cc Matrix.cc mpi_io_stream.cc Problem.cc $<TARGET_OBJECTS:parser>)
add_executable(main main.cc Mesh.cc Matrix.cc Cell.cc mpi_io_stream.cc Problem.cc
                                                                                $<TARGET_OBJECTS:parser> $<TARGET_OBJECTS:tt06>
                                                                                $<TARGET_OBJECTS:stochlib> $<TARGET_OBJECTS:helpers>
                                                                                $<TARGET_OBJECTS:pio>)
target_link_libraries(     main
													 ptscotch
													 scotch
													 scotchmetis
													 ptscotcherr
													 ptscotcherrexit
													 HYPRE
													 )
													 

configure_file(object.data ${CMAKE_CURRENT_BINARY_DIR}/object.data COPYONLY)
