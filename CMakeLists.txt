cmake_minimum_required(VERSION 3.7)
project(MapTran Fortran)

add_compile_options(-O3 -march=native -ffpe-trap=invalid,zero,overflow -g -fbacktrace 
                    -Wall -Wextra -Wpedantic -Warray-bounds)# -Wfatal-errors)
                    
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

add_library(maptran fortran/maptran.f90 fortran/vallado.f90 fortran/assert.f90 )
# --- testing
enable_testing()

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
if (OCTAVE_EXECUTABLE)
  add_test(NAME MatlabOctave COMMAND octave-cli -q Test.m
           WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
endif()

