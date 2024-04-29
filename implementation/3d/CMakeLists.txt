cmake_minimum_required (VERSION 3.5)
set(CMAKE_C_COMPILER "gcc")
SET(CMAKE_C_FLAGS  "-O3")
set(CMAKE_CXX_COMPILER g++)
SET(CMAKE_CXX_FLAGS  "-O3")
set(CMAKE_CXX_STANDARD 11)
set(CLIB clib)

project(cmdline_LBM C CXX)
include_directories("${CMAKE_SOURCE_DIR}/clib" "${CMAKE_SOURCE_DIR}/clib/headers" "${CMAKE_SOURCE_DIR}/lib" "${CMAKE_SOURCE_DIR}/lib/headers" "${CMAKE_SOURCE_DIR}/cmdline/rapidjson" "${CMAKE_SOURCE_DIR}/cmdline/fast-cpp-csv-parser")
find_package (Threads)

file(GLOB LBM_SOURCES ${CMAKE_SOURCE_DIR}/clib/source/*.c)

add_library(${CLIB} STATIC ${LBM_SOURCES})
add_executable(timing  ${CMAKE_SOURCE_DIR}/lib/source/output.cpp ${CMAKE_SOURCE_DIR}/cmdline/timing.cpp ${LBM_SOURCES})
add_executable(cmdline ${CMAKE_SOURCE_DIR}/lib/source/output.cpp ${CMAKE_SOURCE_DIR}/cmdline/main.cpp ${LBM_SOURCES})


# target_include_directories(??)
target_sources(cmdline
PRIVATE
  ${CMAKE_SOURCE_DIR}/lib/source/output.cpp
  ${CMAKE_SOURCE_DIR}/cmdline/rapidjson/document.h
  ${CMAKE_SOURCE_DIR}/cmdline/rapidjson/writer.h
  ${CMAKE_SOURCE_DIR}/cmdline/rapidjson/stringbuffer.h
  ${CMAKE_SOURCE_DIR}/cmdline/fast-cpp-csv-parser/csv.h
)

target_link_libraries(timing ${CMAKE_THREAD_LIBS_INIT} ${CLIB})
target_link_libraries(cmdline ${CMAKE_THREAD_LIBS_INIT} ${CLIB})

