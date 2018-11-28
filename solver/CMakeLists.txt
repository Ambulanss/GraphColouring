cmake_minimum_required(VERSION 3.0)

# Project name
project(solver)

#Set output folder
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories("${PROJECT_SOURCE_DIR}")
include_directories(include)

#either this
set(SOURCES "src/main.cpp")
# or this
#file(GLOB SOURCES "src/*.cpp")

add_executable(solver ${SOURCES})