cmake_minimum_required(VERSION 2.8)

project(mmd)

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -g -std=c++14 -pedantic -pedantic-errors")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -std=c++14 -pedantic -pedantic-errors -O3")

# use OpenMP if compiler supports it
find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_BINARY_DIR build)
set(LIBRARY_OUTPUT_PATH lib)
set(EXECUTABLE_OUTPUT_PATH bin)

include_directories(include)

file(GLOB libsources "src/*.cpp" "include/*.cpp")

add_library(mmd SHARED ${libsources})
add_executable(spring examples/spring.cpp)
add_executable(poly_springs examples/poly_springs.cpp)
add_executable(nondimensional_lj examples/nondimensional_lj.cpp)
add_executable(quenched_lj examples/quenched_lj.cpp)

if (UNIX)
    target_link_libraries(mmd m)
endif (UNIX)

target_link_libraries(mmd armadillo)
target_link_libraries(spring armadillo mmd)
target_link_libraries(poly_springs armadillo mmd)
target_link_libraries(nondimensional_lj armadillo mmd)
target_link_libraries(quenched_lj armadillo mmd)
