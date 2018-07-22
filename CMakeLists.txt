cmake_minimum_required(VERSION 3.7)
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
  list(APPEND FFLAGS -check all -traceback -warn -debug extended)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
  if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 8)
    list(APPEND FFLAGS -std=f2018)
  endif()

  list(APPEND FFLAGS -march=native -ffpe-trap=invalid,zero,overflow -fbacktrace 
                    -Wall -Wextra -Wpedantic -Warray-bounds)# -Wfatal-errors)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI)
  list(APPEND FFLAGS -Mallocatable=03)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Flang)
  list(APPEND FFLAGS -Mallocatable=03)
  link_libraries(-static-flang-libs)
endif()
    

# --- assert library
add_library(assert fortran/assert.f90)
target_compile_options(assert PRIVATE ${FFLAGS})

set(f08 GNU Intel)
if(CMAKE_Fortran_COMPILER_ID IN_LIST f08)
 target_sources(assert PRIVATE fortran/error2008.f90)
else()  # non-Fortran 2008 compilers
  target_sources(assert PRIVATE fortran/error2003.f90)
endif()
# --- other libraries
add_library(maptran fortran/maptran.f90 fortran/vallado.f90)
target_link_libraries(maptran PRIVATE assert)
target_compile_options(maptran PRIVATE ${FFLAGS})
# --- testing
add_executable(testmaptran fortran/Test.f90)
target_link_libraries(testmaptran PRIVATE assert maptran)
target_compile_options(testmaptran PRIVATE ${FFLAGS})
add_test(NAME FortranMap COMMAND testmaptran)
set_tests_properties(Fortran PROPERTIES TIMEOUT 30)

find_package(PythonInterp)
if(PYTHON_VERSION_STRING STRGREATER_EQUAL 3.6)
  add_test(NAME Python COMMAND pytest -sv
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  set_tests_properties(Python PROPERTIES TIMEOUT 30)
endif()

find_package(Octave)
if (OCTAVE_MAJOR_VERSION GREATER_EQUAL 4)
  add_test(NAME MatlabOctave COMMAND octave-cli -q Test.m
           WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
  set_tests_properties(MatlabOctave PROPERTIES TIMEOUT 30)
endif()

