project (PuassonGradient)
cmake_minimum_required(VERSION 3.1)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

#set (CMAKE_CXX_STANDARD 11)
set(CMAKE_C_COMPILER "mpic++")
find_package(MPI REQUIRED)
#find_package(Armadillo REQUIRED)
#include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})

set(SEQ_SOURCES main_seq.cpp mesh.h mesh.cpp mat.h iter.h iter.cpp)
set(PAR_SOURCES main_par.cpp mesh.h mesh.cpp mat.h iter.h iter.cpp)

add_executable(sequential ${SEQ_SOURCES})
add_executable(parallel ${PAR_SOURCES})
target_link_libraries(parallel ${MPI_CXX_LIBRARIES}) #${ARMADILLO_LIBRARIES})
target_link_libraries(sequential ${MPI_CXX_LIBRARIES}) #${ARMADILLO_LIBRARIES})
#target_link_libraries(sequential) #${ARMADILLO_LIBRARIES})
