cmake_minimum_required(VERSION 3.2)

project(bspline CXX)
set(CMAKE_CXX_STANDARD 17)
find_package(Eigen3 3.3 REQUIRED)
add_library(bspline STATIC
    inc/bspline/splines.h
    inc/bspline/bisect.h
    inc/bspline/eigen.h
    src/splines.cpp
)
target_include_directories(bspline PUBLIC ${PROJECT_SOURCE_DIR}/inc)
target_link_libraries(bspline Eigen3::Eigen)

if (CMAKE_TESTING_ENABLED)
    add_subdirectory(tests)
endif()
