# CMake configuration file

# minimum version requirement
cmake_minimum_required(VERSION 3.10)

# set project info
project(ljmd VERSION 0.1 LANGUAGES C)

# add executable
add_executable(ljmd src/ljmd.c)

#configure_file()
#include_directories( ljmd PUBLIC include)
include_directories(${ljmd_SOURCE_DIR}/include)



# add compilation options
add_compile_options(-Wall -g -O3 -ffast-math -fomit-frame-pointer)


add_library(utils src/utils.c src/force_compute.c src/kinetic.c src/output.c src/verlet.c )

# add link libraries
target_link_libraries(ljmd PRIVATE m utils)
target_link_libraries(utils PRIVATE m)
target_link_directories(ljmd PRIVATE src)

# add compile-time definition for version
target_compile_definitions(ljmd PRIVATE -DLJMD_VERSION=${CMAKE_PROJECT_VERSION})

# set executable name
set_target_properties(ljmd PROPERTIES OUTPUT_NAME ljmd-serial.x)
