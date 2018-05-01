cmake_minimum_required(VERSION 3.10)
project(MapTran Fortran)
enable_testing()

set(CTEST_OUTPUT_ON_FAILURE ON)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules/)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_compile_options(-g -O0)
else()
  add_compile_options(-O3)
endif()

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  add_compile_options(-check all -traceback -warn -debug extended)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
  add_compile_options(-O3 -march=native -ffpe-trap=invalid,zero,overflow -g -fbacktrace 
                    -Wall -Wextra -Wpedantic -Warray-bounds)# -Wfatal-errors)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI)
  add_compile_options(-Mallocatable=03)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Flang)
  add_compile_options(-Mallocatable=03)
  link_libraries(-static-flang-libs)
endif()
    


add_library(maptran fortran/maptran.f90 fortran/vallado.f90 fortran/assert.f90 )
# --- testing
add_executable(testmaptran fortran/Test.f90)
target_link_libraries(testmaptran maptran)

add_test(NAME Fortran COMMAND testmaptran)

find_package(PythonInterp)
#if(PYTHON_VERSION_STRING STRGREATER_EQUAL 3.6)
if(PYTHONINTERP_FOUND)
  add_test(NAME Python COMMAND pytest -v
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
endif()

find_package(Octave)
if (OCTAVE_MAJOR_VERSION GREATER_EQUAL 4)
  add_test(NAME MatlabOctave COMMAND octave-cli -q Test.m
           WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
endif()

