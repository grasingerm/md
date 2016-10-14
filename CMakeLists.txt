cmake_minimum_required(VERSION 2.8)

project(md)

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -g -std=c++14 -pedantic -pedantic-errors")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -std=c++14 -pedantic -pedantic-errors -O3")

set(CMAKE_BINARY_DIR build)
set(LIBRARY_OUTPUT_PATH lib)
set(EXECUTABLE_OUTPUT_PATH bin)

include_directories(include)

file(GLOB libsources "src/*.cpp" "include/*.cpp")

add_library(mmd SHARED ${libsources})
add_executable(spring examples/spring/spring.cpp)
add_executable(poly_springs examples/spring/poly_springs.cpp)
add_executable(nondimensional_lj examples/lj1/nondimensional_lj.cpp)
add_executable(quenched_lj examples/lj1/quenched_lj.cpp)
add_executable(liquid_lj examples/lj2/liquid_lj.cpp)

if (UNIX)
    target_link_libraries(mmd m)
endif (UNIX)

target_link_libraries(mmd armadillo)
target_link_libraries(spring armadillo mmd)
target_link_libraries(poly_springs armadillo mmd)
target_link_libraries(nondimensional_lj armadillo mmd)
target_link_libraries(quenched_lj armadillo mmd)
target_link_libraries(liquid_lj armadillo mmd)
