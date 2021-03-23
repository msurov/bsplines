cmake_minimum_required(VERSION 3.2)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIR})
find_package(Threads REQUIRED)

# test bspline
add_executable(test_bspline test_bspline.cpp)
target_link_libraries(test_bspline bspline Threads::Threads ${GTEST_LIBRARIES})
add_test(NAME test_bspline COMMAND test_bspline)
