cmake_minimum_required (VERSION 3.5)
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_STANDARD 11)
set(CLIB clib)

project(cmdline_LBM C CXX)
include_directories("${CMAKE_SOURCE_DIR}/clib" "${CMAKE_SOURCE_DIR}/clib/headers" "${CMAKE_SOURCE_DIR}/lib" "${CMAKE_SOURCE_DIR}/lib/headers" "${CMAKE_SOURCE_DIR}/cmdline/rapidjson" "${CMAKE_SOURCE_DIR}/cmdline/fast-cpp-csv-parser")
find_package (Threads)


add_library(${CLIB} STATIC ${CMAKE_SOURCE_DIR}/clib/source/LBM.c )
add_executable(main_cmdline ${CMAKE_SOURCE_DIR}/lib/source/output.cpp ${CMAKE_SOURCE_DIR}/cmdline/main.cpp ${CMAKE_SOURCE_DIR}/clib/source/LBM.c)
target_sources(main_cmdline
PRIVATE
  ${CMAKE_SOURCE_DIR}/lib/source/output.cpp
  ${CMAKE_SOURCE_DIR}/cmdline/rapidjson/document.h
  ${CMAKE_SOURCE_DIR}/cmdline/rapidjson/writer.h
  ${CMAKE_SOURCE_DIR}/cmdline/rapidjson/stringbuffer.h
  ${CMAKE_SOURCE_DIR}/cmdline/fast-cpp-csv-parser/csv.h
)
target_link_libraries(main_cmdline ${CMAKE_THREAD_LIBS_INIT} ${CLIB})
