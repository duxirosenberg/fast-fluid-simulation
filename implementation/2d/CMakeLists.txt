cmake_minimum_required(VERSION 3.10)
project(LBM_Simulation)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall -std=c++11")

# No need to include the directory now, because CMakeLists.txt is in the 2d directory
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/2d)

set(SOURCES_C lbm2d.c)
set(SOURCES_CPP_MAIN main.cpp)
set(SOURCES_CPP_VALIDATE validate_baseline.cpp)

add_executable(lbm_simulation ${SOURCES_C} ${SOURCES_CPP_MAIN})
add_executable(validate_baseline ${SOURCES_C} ${SOURCES_CPP_VALIDATE})

set_target_properties(lbm_simulation validate_baseline PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)
