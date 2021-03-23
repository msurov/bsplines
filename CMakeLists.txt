cmake_minimum_required(VERSION 3.2)

project(bspline CXX)
set(CMAKE_CXX_STANDARD 14)

add_library(bspline STATIC
    inc/bspline/bspline.h
    src/bspline.cpp
)
target_include_directories(bspline PUBLIC ${PROJECT_SOURCE_DIR}/inc)

find_package(GTest)
if (GTEST_FOUND)
    enable_testing()
    add_subdirectory(tests)
endif()
