cmake_minimum_required(VERSION 3.0)

project(md-hw2)

set(CC clang++-3.8)

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -g -std=c++14 -pedantic -pedantic-errors")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -std=c++14 -pedantic -pedantic-errors -O3")

set(CMAKE_BINARY_DIR build)
set(EXECUTABLE_OUTPUT_PATH bin)

include_directories(include)

file(GLOB libsources "src/*.cpp" "include/*.cpp")

add_library(mmd SHARED ${libsources})
add_executable(prob2i prob/prob2i.cpp)

if (UNIX)
    target_link_libraries(mmd m)
endif (UNIX)

target_link_libraries(mmd armadillo)
target_link_libraries(prob2i armadillo mmd)
