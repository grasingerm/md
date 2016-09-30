cmake_minimum_required(VERSION 2.8)

project(md-hw2)

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -g -std=c++14 -pedantic -pedantic-errors")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -std=c++14 -pedantic -pedantic-errors -O3")

set(CMAKE_BINARY_DIR build)
set(LIBRARY_OUTPUT_PATH lib)
set(EXECUTABLE_OUTPUT_PATH bin)

include_directories(include)

file(GLOB libsources "src/*.cpp" "include/*.cpp")

add_library(mmd SHARED ${libsources})
add_executable(prob2i prob/prob2i.cpp)
add_executable(prob2ii prob/prob2ii.cpp)
add_executable(prob3 prob/prob3.cpp)
add_executable(bonus prob/bonus.cpp)

if (UNIX)
    target_link_libraries(mmd m)
endif (UNIX)

target_link_libraries(mmd armadillo)
target_link_libraries(prob2i armadillo mmd)
target_link_libraries(prob2ii armadillo mmd)
target_link_libraries(prob3 armadillo mmd)
target_link_libraries(bonus armadillo mmd)
