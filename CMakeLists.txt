## CMakeList.txt : Top-level CMake project file, do global configuration
## and include sub-projects here.
##
#
#project(pp)
#
#cmake_minimum_required (VERSION 2.6)
#set(CMAKE_CXX_STANDARD 11)
##set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -idirafter /usr/include")
##IF(NOT CMAKE_BUILD_TYPE)
##    #SET(CMAKE_BUILD_TYPE "DEBUG")
##    SET(CMAKE_BUILD_TYPE "RELEASE")
##    #SET(CMAKE_BUILD_TYPE "RELWITHDEBINFO")
##    #SET(CMAKE_BUILD_TYPE "MINSIZEREL")
##ENDIF()
##include_directories ("training")
##include_directories ("testing")
#include_directories("inc")
#include_directories("src")
#file(GLOB SOURCES "src/*.cpp")
##file(GLOB SOURCES "inc/*.h")
## Find Boost
#find_package( Boost 1.49.0 REQUIRED COMPONENTS program_options system filesystem)
#include_directories(${Boost_INCLUDE_DIRS})
#
#find_package(dlib REQUIRED)
##add_library(${SOURCES})
##add_executable(pp ${SOURCES} )
#add_executable(pp src/driver_train.cpp src/CBS.cpp)
#target_link_libraries(pp dlib::dlib)
#target_link_libraries(pp ${Boost_LIBRARIES})
#target_include_directories(pp PRIVATE "src")
#target_include_directories(pp PRIVATE "inc")
##project(pp)
##cmake_minimum_required (VERSION 2.6)
##
##include_directories("src")
##include_directories("inc")
##add_executable(pp training/driver_train.cpp)


project(pp)

cmake_minimum_required (VERSION 2.6)
set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -idirafter /usr/include")
#IF(NOT CMAKE_BUILD_TYPE)
#    #SET(CMAKE_BUILD_TYPE "DEBUG")
#    SET(CMAKE_BUILD_TYPE "RELEASE")
#    #SET(CMAKE_BUILD_TYPE "RELWITHDEBINFO")
#    #SET(CMAKE_BUILD_TYPE "MINSIZEREL")
#ENDIF()

include_directories("inc")
file(GLOB SOURCES "src/*.cpp")

# Find Boost
find_package( Boost 1.49.0 REQUIRED COMPONENTS program_options system filesystem)
include_directories(${Boost_INCLUDE_DIRS})

find_package(dlib REQUIRED)

add_executable(pp ${SOURCES} )
target_link_libraries(pp dlib::dlib)
target_link_libraries(pp ${Boost_LIBRARIES})