# CMake configuration file

# minimum version requirement
cmake_minimum_required(VERSION 3.10)

# set project info
project(ljmd VERSION 1.0 LANGUAGES C)

# add executable and its properties
add_executable(ljmd src/ljmd.c)
target_include_directories(ljmd PRIVATE include)
target_compile_definitions(ljmd PRIVATE -DLJMD_VERSION=${CMAKE_PROJECT_VERSION})
set_target_properties(ljmd PROPERTIES OUTPUT_NAME ljmd-serial.x)

# add library
add_library(utils
        src/force_compute.c
        src/kinetic.c
        src/output.c
        src/utils.c
        src/verlet.c
)

# add compilation options
add_compile_options(-Wall -g -O3 -ffast-math -fomit-frame-pointer)
target_include_directories(utils PRIVATE include)

# add lik libraries
target_link_libraries(ljmd PRIVATE utils)
target_link_libraries(utils PRIVATE m)
