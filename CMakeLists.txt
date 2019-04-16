cmake_minimum_required(VERSION 3.13)
project(LU)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fmax-errors=1")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -gdwarf-3")


set(Source src/lu.cpp src/resol.cpp src/sistema.cpp)
include_directories(include libs/include test/include)

add_library(LU OBJECT ${Source})
add_executable(Main src/main.cpp $<TARGET_OBJECTS:LU>)

file(GLOB Tests "test/src/*.cpp")
add_executable(Test ${Tests} $<TARGET_OBJECTS:LU>)
