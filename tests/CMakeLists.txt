cmake_minimum_required(VERSION 3.2)

find_package(GTest REQUIRED)
find_package(Threads REQUIRED)

add_executable(test_bspline test_bspline.cpp)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
target_link_libraries(test_bspline bspline Threads::Threads ${GTEST_LIBRARIES})
add_test(NAME test_bspline COMMAND test_bspline)
