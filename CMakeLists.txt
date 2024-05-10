cmake_minimum_required(VERSION 3.5)
project(LatticeBoltzmanMethod C CXX)

# Set the compilers and the compiler flags
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_C_FLAGS "-O3 -mavx2 -mfma")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-O3 -mavx2 -mfma")
set(CMAKE_CXX_STANDARD 11)

# Include directories for JSON and CSV parsers (if still needed)
include_directories(
    "${CMAKE_SOURCE_DIR}/lib/rapidjson" 
    "${CMAKE_SOURCE_DIR}/lib/fast-cpp-csv-parser"
    "${CMAKE_SOURCE_DIR}/include/2d"
    "${CMAKE_SOURCE_DIR}/include/3d"
    "${CMAKE_SOURCE_DIR}/include/3d/stream"
    "${CMAKE_SOURCE_DIR}/include/shared"
    )

# Find the Threads package
find_package(Threads)


# -------------------------------------------- 3D ----------------------------------#
# Define the source files for the 3D executables
file(GLOB LBM_SOURCES ${CMAKE_SOURCE_DIR}/src/3d/*.c)
file(GLOB STREAM_SOURCES ${CMAKE_SOURCE_DIR}/src/3d/stream/*.c)
add_library(${CLIB} STATIC ${LBM_SOURCES})

# Add the executables
add_executable(cmdline ${CMAKE_SOURCE_DIR}/scripts/3d/main.cpp ${CMAKE_SOURCE_DIR}/scripts/3d/output.cpp ${LBM_SOURCES} ${STREAM_SOURCES})
add_executable(timing ${CMAKE_SOURCE_DIR}/scripts/3d/timing.cpp ${CMAKE_SOURCE_DIR}/scripts/3d/output.cpp ${LBM_SOURCES} ${STREAM_SOURCES})

# Link libraries to executables if any common libraries are to be used
target_link_libraries(cmdline ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(timing ${CMAKE_THREAD_LIBS_INIT})


# -------------------------------------------- 2D ----------------------------------#
#Define the source files for the 2D executables
file(GLOB 2D_SOURCES ${CMAKE_SOURCE_DIR}/src/2d/lbm2d_baseline.c)
file(GLOB 2D_STRUCT_SOURCES ${CMAKE_SOURCE_DIR}/src/2d/struct/*.c)

#add the executables
add_executable(2D_validate_baseline ${CMAKE_SOURCE_DIR}/scripts/2d/validate_baseline.cpp ${2D_SOURCES})
add_executable(2D_timing ${CMAKE_SOURCE_DIR}/scripts/2d/timing.cpp ${2D_SOURCES})
add_executable(2D_testing ${CMAKE_SOURCE_DIR}/scripts/2d/testing.cpp ${2D_SOURCES})
#add_executable(2D_struct_testing  ${CMAKE_SOURCE_DIR}/scripts/2d/testing_struct.cpp ${2D_STRUCT_SOURCES} )



