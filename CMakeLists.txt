cmake_minimum_required(VERSION 3.16)
project(structurePlotter CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(Boost REQUIRED)

find_path(CAIRO_INCLUDE_DIR cairo-pdf.h PATH_SUFFIXES cairo)
find_library(CAIRO_LIBRARY cairo)
if(NOT CAIRO_INCLUDE_DIR OR NOT CAIRO_LIBRARY)
    message(FATAL_ERROR "Could not find cairo (cairo-pdf.h / libcairo)")
endif()

add_executable(structurePlotter src/main.cpp)

target_include_directories(structurePlotter PRIVATE ${CAIRO_INCLUDE_DIR})
target_link_libraries(structurePlotter PRIVATE Eigen3::Eigen Boost::boost ${CAIRO_LIBRARY})
