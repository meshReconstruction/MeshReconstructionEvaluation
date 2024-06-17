cmake_minimum_required(VERSION 3.22)
project(edgeDetection)

set(CMAKE_CXX_STANDARD 23)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

add_executable(edgeDetection main.cpp)

target_link_libraries(edgeDetection PUBLIC CGAL::CGAL_Basic_viewer)