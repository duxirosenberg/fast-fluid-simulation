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


add_library( interfaceOld INTERFACE )
target_compile_definitions(interfaceOld INTERFACE)

add_library( interfaceStruct INTERFACE )
target_compile_definitions(interfaceStruct INTERFACE -DLBM_STRUCT )


add_library( interfaceBase INTERFACE )
target_compile_definitions(interfaceBase INTERFACE -DLBM_STRUCT -DLBM_BASELINE)



add_library(${CLIB} STATIC ${CMAKE_SOURCE_DIR}/clib/source/LBM.c )
add_executable(main_cmdlineOld  ${CMAKE_SOURCE_DIR}/lib/source/output.cpp ${CMAKE_SOURCE_DIR}/cmdline/main.cpp ${CMAKE_SOURCE_DIR}/clib/source/measurement.c ${CMAKE_SOURCE_DIR}/clib/source/LBM.c )
add_executable(main_cmdline ${CMAKE_SOURCE_DIR}/lib/source/output.cpp ${CMAKE_SOURCE_DIR}/cmdline/main.cpp ${CMAKE_SOURCE_DIR}/clib/source/measurement.c ${CMAKE_SOURCE_DIR}/clib/source/LBMstruct.c )
add_executable(main_cmdlineBaseline ${CMAKE_SOURCE_DIR}/lib/source/output.cpp ${CMAKE_SOURCE_DIR}/cmdline/main.cpp ${CMAKE_SOURCE_DIR}/clib/source/measurement.c ${CMAKE_SOURCE_DIR}/clib/source/LBMbaseline.c )


# target_include_directories(??)
target_sources(main_cmdline
PRIVATE
  ${CMAKE_SOURCE_DIR}/lib/source/output.cpp
  ${CMAKE_SOURCE_DIR}/cmdline/rapidjson/document.h
  ${CMAKE_SOURCE_DIR}/cmdline/rapidjson/writer.h
  ${CMAKE_SOURCE_DIR}/cmdline/rapidjson/stringbuffer.h
  ${CMAKE_SOURCE_DIR}/cmdline/fast-cpp-csv-parser/csv.h
)

target_link_libraries(main_cmdlineOld ${CMAKE_THREAD_LIBS_INIT} ${CLIB} interfaceOld)
target_link_libraries(main_cmdline ${CMAKE_THREAD_LIBS_INIT} ${CLIB} interfaceStruct)
target_link_libraries(main_cmdlineBaseline ${CMAKE_THREAD_LIBS_INIT} ${CLIB} interfaceBase)



# add_library(${CLIB} STATIC ${CMAKE_SOURCE_DIR}/clib/source/LBM.c )
# add_executable(main_cmdline ${CMAKE_SOURCE_DIR}/lib/source/output.cpp ${CMAKE_SOURCE_DIR}/cmdline/main.cpp ${CMAKE_SOURCE_DIR}/clib/source/LBM.c)
# target_sources(main_cmdline
# PRIVATE
#   ${CMAKE_SOURCE_DIR}/lib/source/output.cpp
#   ${CMAKE_SOURCE_DIR}/cmdline/rapidjson/document.h
#   ${CMAKE_SOURCE_DIR}/cmdline/rapidjson/writer.h
#   ${CMAKE_SOURCE_DIR}/cmdline/rapidjson/stringbuffer.h
#   ${CMAKE_SOURCE_DIR}/cmdline/fast-cpp-csv-parser/csv.h
# )
# target_link_libraries(main_cmdline ${CMAKE_THREAD_LIBS_INIT} ${CLIB})

