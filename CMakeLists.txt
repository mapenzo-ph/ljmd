# CMake configuration file

# minimum version requirement
cmake_minimum_required(VERSION 3.10)

# set project info
project(ljmd VERSION 1.0 LANGUAGES C)

# add executable
add_executable(ljmd src/ljmd.c)

# add compilation options
add_compile_options(-Wall -g -O3 -ffast-math -fomit-frame-pointer)

# add lik libraries
target_link_libraries(ljmd PRIVATE m)

# add compile-time definition for version
target_compile_definitions(ljmd PRIVATE -DLJMD_VERSION=${CMAKE_PROJECT_VERSION})

# set executable name
set_target_properties(ljmd PROPERTIES OUTPUT_NAME ljmd-serial.x)
