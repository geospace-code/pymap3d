#!/bin/bash
# boilerplate to test with popular Fortran compilers, helping check for quirks

# --- colors https://stackoverflow.com/a/20983251
red=`tput setaf 1`
reset=`tput sgr0`

fcomp=(gfortran gfortran-5 gfortran-6 gfortran-7 gfortran-8 ifort pgf95 flang   nagfor)
ccomp=(gcc      gcc-5      gcc-6      gcc-7      gcc-8      icc   pgcc  clang   gcc)
pcomp=(g++      g++-5      g++-6      g++-7      g++-8      icpc  pgc++ clang++ g++)
# --- loops
for i in {1..8}
do

(

export FC=${fcomp[$i]}
export CC=${ccomp[$i]}
export CXX=${pcomp[$i]}
[[ $FC == ifort ]] && . ~/intel.sh

cd "${0%/*}"  # change to directory of this script

echo
es=$(command -v $FC)
ec=$(command -v $CC)
ep=$(command -v $CXX)

if [[ $es && $ec && $ep ]]
then
  echo  
  echo "testing with"
  echo $FC  $CC  $CXX
  echo "press Enter to proceed."
  read skip
  [[ $skip == "s" ]] && continue
else
  echo
  echo "${red}*** skipping $FC $CC $CXX *** ${reset}"
  continue
fi

touch ../bin/junk
rm -r ../bin/*
cd ../bin

cmake ..
cmake --build .
ret=$?
[[ $ret == 0 ]] && ctest -V || exit $ret

)
  
done

cd "${0%/*}"  # change to directory of this script
rm -r ../bin/*
