#!/bin/bash
# Shell script to compile qd library for Multifan-CL easier on Linux and MacOS 
#
#  Author Yukio Takeuchi
#  2019-01-09 Now can be run both Linux and MacOS
#  
SECONDS=0
LIB="$(pwd)/libs"
if [ ! -d "$LIB"  ]; then
  mkdir "$LIB"
fi

DIRORG="$(pwd)"
if [ "$(uname)" == 'Darwin' ]; then
  OS='Mac'
  CXX=g++-8
  CC=gcc-8
  CXXFLAGS='-O3' 
  FCFLAGS="-O3" 
elif [ "$(expr substr $(uname -s) 1 5)" == 'Linux' ]; then
  OS='Linux'
  CXX=g++
elif [ "$(expr substr $(uname -s) 1 10)" == 'MINGW32_NT' ]; then                                                                                           
  OS='Cygwin'
else
  echo "Your platform ($(uname -a)) is not supported."
  exit 1
fi
##

OPTARRAY=("default" "O3" "O3fma" "native")
#OPTARRAY=("default")
#OPTARRAY=("O3")
#OPTARRAY=("native")
#OPTARRAY=("O3fma" "O3")
#OPTARRAY=("O3fma")
#VERSION=2.3.20
VERSION=2.3.22

tar xvf qd-$VERSION.tar.gz
cd qd-$VERSION

for OPT in ${OPTARRAY[@]} ; 
do
  echo "OPT is $OPT"
  if [ -d "$LIB/qd-${VERSION}_$OPT" ]; then
    rm -r "$LIB/qd-${VERSION}_$OPT" 
  fi
  if [ $OPT="dynamic" ]; then
    MAKEARGS="DYNAMIC_ARCH=1"
  else
    MAKEARCH=""
  fi 
  if [ -f "Nakefile" ]; then
    make distclean
  fi 
  PREFIX="${LIB}/qd-${VERSION}_$OPT"
  if [ ! -d $PREFIX ]; then
    mkdir "$PREFIX"
    if [ ! -d "$PREFIX/include/qd/" ] ; then
      mkdir -p "$PREFIX/include/qd/"
    fi
    cp -p "$DIRORG/qd.h" "$PREFIX/include/qd/qd.h"
    if [ $? -ne 0 ]; then
      echo "Could not copy qd.h to $PREFIX/include/qd/qd.h" ; exit
    fi
  fi 
  PREFIX=`cd $PREFIX && echo "$(pwd)"`
  echo "$(pwd)"
  echo "PREFIX=$PREFIX"
  
  if [ $OPT = "default"  ]; then 
    echo "OPT=$OPT"
    ./configure --prefix=$PREFIX CXX=$CXX
  elif [ $OPT = "O3"  ]; then 
    echo "OPT=$OPT"
    ./configure --prefix=$PREFIX CXXFLAGS='-O3'  FCFLAGS="-O3" CFLAGS="-O3"  CXX=$CXX 
  elif [ $OPT = "O3fma" ]; then
    echo "OPT=$OPT"
    ./configure --prefix=$PREFIX CXXFLAGS='-O3'  FCFLAGS="-O3" CFLAGS='-O3' CXX=$CXX  --enable-fma="c99"
  elif [ $OPT = "native" ]; then
    echo "OPT=$OPT"
    ./configure --prefix=$PREFIX CXXFLAGS='-O3 -march=native'  FCFLAGS='-O3 -march=native' CFLAGS='-O3 -march=native'  CXX=$CXX  --enable-fma="c99"
  else
    echo "No valid OPT! OPT=$OPT"
    exit
  fi 
  echo "CC=$CC"
  echo "CXX=$CXX"
  echo "FCFLAGS=$FCFLAGS"
  echo "CXXFLAGS=$CXXFLAGS"
  #exit
  make clean
  make 
  if [ $? -ne 0 ]; then
    echo "Make failed"
    echo "PREFIX=$PREFIX" 
    exit
  fi
  make check 
  make time
  make install
  
  #unset CXXFLAGS
  #unset FCFLAGS
  #unset CFLAGS
  echo "CXX=$CXX"
  echo "FCFLAGS=$FCFLAGS"
  echo "CXXFLAGS=$CXXFLAGS" 
done
time=$SECONDS ; echo "$time sec ellapsed"
exit
