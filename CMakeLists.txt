cmake_minimum_required(VERSION 3.13)
project(mandelbulb C)

set(CMAKE_C_STANDARD 11)

add_executable(mandelbulb main.c qbmp/qdbmp.h qbmp/qdbmp.c)
target_link_libraries(mandelbulb m)