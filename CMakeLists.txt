# CMake configuration file

# minimum version requirement
cmake_minimum_required(VERSION 3.10)

# set project info
project(ljmd VERSION 1.0 LANGUAGES C)

# add executable and its properties
add_executable(ljmd src/ljmd.c src/force_compute.c src/kinetic.c src/output.c src/utils.c src/verlet.c)

target_include_directories(ljmd PRIVATE include)
target_compile_definitions(ljmd PRIVATE -DLJMD_VERSION=${CMAKE_PROJECT_VERSION})
set_target_properties(ljmd PROPERTIES OUTPUT_NAME ljmd-serial.x)


# add compilation options
add_compile_options(-Wall -g -O3 -ffast-math -fomit-frame-pointer)

# add lik libraries
target_link_libraries(ljmd PRIVATE m)


if(ENABLE_TESTING)
    enable_testing()
    file(WRITE ${CMAKE_BINARY_DIR}/test_write.sh "cp -r ../test . && cd test && ../ljmd-serial.x < simple_test.inp ") 
    
    add_test(NAME Runs COMMAND bash test_write.sh)

endif()