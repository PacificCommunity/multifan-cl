#!/bin/bash
# Shell script to compile qd library for Multifan-CL easier on Linux and MacOS 
#
#  Author Yukio Takeuchi
#  2019-01-09 Now can be run both Linux and MacOS
#  
SECONDS=0

LIB="$(pwd)/libs"
if [ ! -d "$LIB" ]; then
  mkdir -p "$LIB"
fi

if [ "$(uname)" == 'Darwin' ]; then
  OS='Mac'
  export CC=gcc-8
  VERSION=0.3.5
elif [ "$(expr substr $(uname -s) 1 5)" == 'Linux' ]; then
  OS='Linux'
  export CC=gcc
  VERSION=0.2.20
elif [ "$(expr substr $(uname -s) 1 10)" == 'MINGW32_NT' ]; then                                                                                           
  OS='Cygwin'
  export CC=gcc
else
  echo "Your platform ($(uname -a)) is not supported."
  exit 1
fi
##
OPTARRAY=("default" "generic" "dynamic" )
#OPTARRAY=("dynamic")
#OPTARRAY=("generic")
#OPTARRAY=("dynamic" "generic")
#`VERSION=0.2.20
#VERSION=0.3.3
tar xvf OpenBLAS-$VERSION.tar.gz
cd OpenBLAS-$VERSION

for OPT in ${OPTARRAY[@]} ; 
do
  echo "OPT is $OPT"
  if [ -d "$LIB/OpenBLAS-${VERSION}_${OPT}" ]; then
    rm -r "$LIB/OpenBLAS-${VERSION}_${OPT}" 
  fi
  MAKEARGS="  NUM_THREADS=2  NO_SHARED=1"
  if [ $OPT = "dynamic" ]; then
    MAKEARGS+=" DYNAMIC_ARCH=1 "
  elif [ $OPT = "generic" ]; then
    MAKEARGS+=" TARGET=GENERIC "
  else
    MAKEARGS+=""
  fi 
  make clean
  #RET=`echo $(make $MAKEARGS)`
  make $MAKEARGS
  RET=$?
  echo "RET=$RET"
  echo "RET"
  #time=$SECONDS ; echo "$time sec ellapsed" 
  #exit
  if [ $RET -ne 0 ]; then
    echo "Build failed"
    echo "MAKEARGS are ${MAKEARGS}"
  fi
  #RET=`echo $(make PREFIX=${LIB}/OpenBLAS-${VERSION}_$OPT $MAKEARGS install)`
  make PREFIX=${LIB}/OpenBLAS-${VERSION}_${OPT} ${MAKEARGS} install
  RET=$?
  if [ $RET -ne 0 ]; then
    echo "Build failed"
    echo "MAKEARGS are ${MAKEARGS}"
    echo "PREFIX is ${PREFIX}"
  fi
done
time=$SECONDS ; echo "$time sec ellapsed"
exit

