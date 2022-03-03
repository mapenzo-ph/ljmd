# CMake configuration file

# minimum version requirement
cmake_minimum_required(VERSION 3.10)

# set project info
project(ljmd VERSION 1.0 LANGUAGES C)

# add executable and its properties
add_executable(ljmd src/ljmd.c)

target_include_directories(ljmd PRIVATE include)
target_compile_definitions(ljmd PRIVATE -DLJMD_VERSION=${CMAKE_PROJECT_VERSION})
set_target_properties(ljmd PROPERTIES OUTPUT_NAME ljmd.x)

add_library(ljmd_lib src/force_compute.c src/kinetic.c src/output.c src/utils.c src/verlet.c)

target_include_directories(ljmd_lib PRIVATE include)

# add general compilation options
add_compile_options(-Wall -g -O3 -ffast-math -fomit-frame-pointer)

# add OpenMP support
if(USE_OMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(ljmd_lib PRIVATE OpenMP::OpenMP_C)
endif()

#add MPI support
if(USE_MPI)
    find_package(MPI REQUIRED)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
    target_link_libraries(ljmd PRIVATE ${MPI_C_LIBRARIES} )
    target_link_libraries(ljmd_lib PRIVATE ${MPI_C_LIBRARIES} )
    target_compile_definitions(ljmd PRIVATE -DUSE_MPI)
    target_compile_definitions(ljmd_lib PRIVATE -DUSE_MPI)
endif()

# add link libraries for all executables
target_link_libraries(ljmd PRIVATE ljmd_lib)
target_link_libraries(ljmd_lib PRIVATE m)

if(ENABLE_TESTING)
    enable_testing()
    enable_language(CXX)
    
    # test functions
    set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
    include(GTest)

    # test force
    add_executable(test_force test/test_force.cpp)
    target_include_directories(test_force PRIVATE include)
    target_link_libraries(test_force PRIVATE GTest::GTestMain GTest::GTest ljmd_lib)
    if(USE_OMP)
        target_link_libraries(test_force PRIVATE OpenMP::OpenMP_C)
    endif()
    add_test(NAME Force COMMAND test_force)

    # test kinetic energy
    add_executable(test_ekin test/test_ekin.cpp)
    target_include_directories(test_ekin PRIVATE include)
    target_link_libraries(test_ekin PRIVATE GTest::GTestMain GTest::GTest ljmd_lib)
    if(USE_OMP)
        target_link_libraries(test_force PRIVATE OpenMP::OpenMP_C)
    endif()
    add_test(NAME EKin COMMAND test_ekin)

    # test verlet
    add_executable(test_verlet test/test_verlet.cpp)
    target_include_directories(test_verlet PRIVATE include)
    target_link_libraries(test_verlet PRIVATE GTest::GTestMain GTest::GTest ljmd_lib)
    if(USE_OMP)
        target_link_libraries(test_verlet PRIVATE OpenMP::OpenMP_C)
    endif()
    add_test(NAME Verlet COMMAND test_verlet)

    # test if application runs
    file(WRITE ${CMAKE_BINARY_DIR}/test_write.sh "cp -r ../test . && cd test && ../ljmd.x < simple_test.inp ") 
    add_test(NAME Runs COMMAND bash test_write.sh)

endif()
