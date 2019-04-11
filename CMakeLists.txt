cmake_minimum_required(VERSION 3.13)
project(LU)

set(CMAKE_CXX_STANDARD 11)


set(Source src/lu.cpp src/resol.cpp src/sistema.cpp)
include_directories(include libs/include)

add_library(LU OBJECT ${Source})
add_executable(Main src/main.cpp LU)

file(GLOB Tests "test/*.cpp")
add_executable(Test ${Tests} LU)
