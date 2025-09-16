#!/bin/bash
# Shell script to make compilation of Multifan-CL easier on Linux and MacOS 
#
#  Author Yukio Takeuchi
#  2019-01-09 Now can be run both Linux and MacOS
#  
SECONDS=0


if [ "$(which realpath)" == "" ]; then
  echo "realpath command is needed. But could not find it"
  exit 1
fi
DEBUG=1  # for release (optimized) build. This line should be DEBUG=0
#DEBUG=0  # for debug build. Yhis line should be DEBUG=1
## ADMBSRCDIR  full or relative path of root of ofp-sam-admb repository
## MFCLSRCDIR  full or relative path of root of ofp-sam-mfcl repository
## You need to adjust ADMBSRCDIR and MFCLSRCDIR 
#ADMBSRCDIR="$HOME/repository/ofp-sam-admb/"
ADMBSRCDIR="./ofp-sam-admb/"
#MFCLSRCDIR="$HOME/repository/ofp-sam-mfcl"
MFCLSRCDIR="./ofp-sam-mfcl"
LIB="$(pwd)/libs"
DIRORG="$(pwd)"
if [ ! -d "$LIB"  ]; then
  echo "$LIB does not exist"
  exit 1
fi

if [ "$(uname)" == 'Darwin' ]; then
  OS='Mac'
  CXX=g++-8
  CC=gcc-8
  CXXFLAGS='-O3' 
  FCFLAGS="-O3" 
  OPENBLASVER="0.3.5"
  QDVER="2.3.22"
  ADMBMAKEFILEORG="newmake-threaded_mac.mak"
  MFCLMAKEFILEORG="mfcl_thrd_mac64.mak"
   if [ "$DEBUG" -ne 1 ]; then 
    MFCLEXE="mfclo64"
    MFCLOPT="OPT=TRUE"
  else
    MFCLEXE="mfclsdbg64"
    MFCLOPT=""
  fi 
  ### Count number of cpu's available
  NCPU=$(sysctl -n hw.physicalcpu)
#  echo "This script is not yet ported for mac os" ; exit 1
elif [ "$(expr substr $(uname -s) 1 5)" == 'Linux' ]; then
  OS='Linux'
  CXX=g++
  NAME=$(cat /etc/os-release | grep "^NAME" | tr '=' ' ' |gawk '{gsub(/\"/,"") ;print $2}')
  NAME+=$(cat /etc/os-release | grep "^VERSION_ID" | tr '=' ' ' |gawk '{gsub(/\"/,"") ;print $2}')
  OPENBLASVER="0.2.20"
  QDVER="2.3.22"
  ADMBMAKEFILEORG="newmake-threaded_NOUM02.mak"
  if [ "$DEBUG" -ne 1 ]; then 
    MFCLMAKEFILEORG="mfcl_thrd_linux64_opt.mak"
    MFCLEXE="mfclo64"
    MFCLOPT=""
  else
    MFCLMAKEFILEORG="mfcl_thrd_linux64_debug.mak"
    MFCLEXE="mfclsdbg64"
    MFCLOPT=""
  fi 
  ### Count number of cpu's available
  NCPU=$(grep -c ^processor /proc/cpuinfo)
elif [ "$(expr substr $(uname -s) 1 10)" == 'MINGW32_NT' ]; then
  OS='Cygwin'
else
  echo "Your platform ($(uname -a)) is not supported."
  exit 1
fi


ADMBSRCDIR="$(realpath $ADMBSRCDIR)"
if [ ! -d "$ADMBSRCDIR" ]; then
  echo "$ADMBSRCDIR does not exist"
    exit
fi
MFCLSRCDIR="$(realpath $MFCLSRCDIR)"

MFCLSRCDIR="${MFCLSRCDIR}"/src
echo "ADMBSRCDIR is $ADMBSRCDIR"
echo "MFCLSRCDIR is $MFCLSRCDIR"

QDOPTS=( "default" "O3" "O3fma" "native" )
#QDOPTS=( "default" )
OPENBLASOPTS=( "generic" "default" "dynamic" )
#OPENBLASOPTS=( "generic" )
#OPTS=( "O2" "O3" )
OPTS=( "O2" )
if [ $DEBUG -eq 1 ]; then
  OPTS=( "O0" )
fi
if [ ${#OPTS[@]} -eq 0 ]; then
  echo "OPTS is empty. Stop"
  exit
fi

#TARGETDIR0="$HOME/bin/devvsn11/${NAME}/`date '+%Y%m%d'`_"
TARGETDIR0="${DIRORG}/devvsn11/${NAME}/$(date '+%Y%m%d')_"
echo "TARGETDIR0=$TARGETDIR0"
echo "HERE81"
#exit
#for OPT in ${OPTS[@]} ;
#do	
#  echo "OPT=$OPT"

for QDOPT in "${QDOPTS[@]}" ;
do
  echo "QDOPT is $QDOPT"
  #QDNEW="${DIRORG}/lib/qd-${QDVER}_${QDOPT}/"
  QDNEW="${LIB}/qd-${QDVER}_${QDOPT}/"
  echo "QDNEW=$QDNEW"
  ##########
  #ADMBHOMENEW="${DIRORG}/lib/ADMB_qd-${QDVER}_${QDOPT}/"
  ADMBHOMENEW="${LIB}/ADMB_qd-${QDVER}_${QDOPT}/"
      ###########
  if [ -d $ADMBSRCDIR ]; then
    cd "$ADMBSRCDIR" || exit
    # ONELINER='s/\/home\/mfcl\/MFCL\/libs\/qd-2.3.17\//\/home\/yukiot\/buildScripts\/lib\/qd-2.3.22_'"$QDOPT"'\//g'
    if [ ! -f "$ADMBMAKEFILEORG" ]; then
      echo "$ADMBMAKEFILEORG does not exist"
      exit
    fi
    echo "ADMBMAKEFILEORG=$ADMBMAKEFILEORG" 
    
    #QDORG=$(cat "$ADMBMAKEFILEORG" | gawk '{$1=$1;print}' | grep "^QD_HOME=" | tr "=" " " | gawk '{print $2}')
    #QDORG=`cat "$ADMBMAKEFILEORG" | gawk '{$1=$1;print}' | grep "^QD_HOME="` 
    QDORG=$(cat "$ADMBMAKEFILEORG" | gawk '{gsub(/^[ \t]+/,"",$2);gsub(/[ \t]+$/,"",$2)}1' | grep -e "^QD_HOME=" -e "^QD_HOME:=" | gawk -F'=' '{print $2}' )
    #QDORG=$(cat "$ADMBMAKEFILEORG" | gawk '{gsub(/^[ \t]+/,"",$2);gsub(/[ \t]+$/,"",$2)}1' | grep "^QD_HOME=" | gawk -F'=' '{print $2}' ) 
    echo "QDORG=$QDORG"
    if [ "$QDORG" == "" ]; then
      echo "could not find QD_HOME" ;exit 
    fi
    echo "Original QD_HOME=$QDORG"
    echo "New QD_HOME=$QDNEW"
    #exit
    ADMBHOMEORG=$(cat "$ADMBMAKEFILEORG" | gawk '{$1=$1;print}' | grep -e "^ADMB_HOME=" -e "^ADMB_HOME:=" | gawk -F'=' '{print $2}' )
    echo "ADMBHOMEORG=$ADMBHOMEORG"
    #cat newmake-threaded_NOUM02.mak | sed -e $ONELINER > newmake-threaded_build.mak
    #cat "$ADMBMAKEFILEORG" | gawk "{gsub(  "\""$QDORG"\"", "\""$QDNEW"\"" ) ; print \$0 }" >newmake-threaded_build.mak
    cat "$ADMBMAKEFILEORG" | gawk "{gsub(  "\""$QDORG"\"", "\""$QDNEW"\"" ) ; print \$0 }" | \
	    gawk "{gsub(  "\""$ADMBHOMEORG"\"", "\""$ADMBHOMENEW"\"" ) ; print \$0 }"  >newmake-threaded_build.mak
    echo "HERE143"
    if [ $? -ne 0 ]; then
      echo "Failed on line145" ; exit
    fi
    #exit
    make -f newmake-threaded_build.mak clean
    make -f newmake-threaded_build.mak  -j$NCPU
    RET=$?
    if [ $RET -ne 0 ]; then
      echo "Build  of ADMB  failed"
      echo "OPT=$OPT" ;  echo "QDOPT is $QDOPT" ; echo "OPENBLASOPT=$OPENBLASOPT" ;echo "TARGETDIR=$TARGETDIR"
      exit
    fi
  else
    echo "$ADMBSRCDIR does not exist"
    echo "OPT=$OPT" ;echo "QDOPT is $QDOPT" ; echo "OPENBLASOPT=$OPENBLASOPT" ; echo "TARGETDIR=$TARGETDIR"
    exit
  fi
  for OPENBLASOPT in "${OPENBLASOPTS[@]}" ;
  do
    echo "OPENBLASOPT=$OPENBLASOPT"
    for OPT in "${OPTS[@]}" ;
    do
      echo "OPT=$OPT"
      if [ -d "$MFCLSRCDIR" ]; then
        TARGETDIR="${TARGETDIR0}${OPT}_QD_${QDOPT}_${OPENBLASOPT}_OpenBLAS"
        echo "TARGETDIR=$TARGETDIR"
        if [ ! -d "$TARGETDIR" ] ; then
          mkdir -p "$TARGETDIR"
          if [ $? -ne 0 ]; then
            echo "Could not create $TARGETDIR"
          else
            echo "$TARGETDIR was created"
          fi
        else
          echo "$TARGETDIR exists"
        fi
        cd "$MFCLSRCDIR" || exit
        #ONELINER='s/\/home\/mfcl\/MFCL\/libs\/qd-2.3.17/\/home\/yukiot\/buildScripts\/lib\/qd-2.3.22_'"$QDOPT"'/g'
        #cat mfcl_thrd_linux64_opt.mak | sed -e $ONELINER > tmp.mak
        
        #ONELINER2='s/\/home\/mfcl\/MFCL\/libs\/OpenBLAS-0.2.19/\/home\/yukiot\/buildScripts\/lib\/OpenBLAS-0.2.20_'"${OPENBLASOPT}"'/g'
        #cat tmp.mak | sed -e $ONELINER2 > mfcl_thrd_linux64_build_opt.mak
        ## xargs work to trim space (see https://qiita.com/kazuhidet/items/09eb9dad8b44003d555f )
        #QDORG=$(cat $MFCLMAKEFILEORG |  gawk '{$1=$1;print}' | grep "^QD_HOME=" | tr "=" " " | gawk '{print $2}')
        if [ -f "$MFCLMAKEFILEORG" ]; then
          QDORG=$(cat $MFCLMAKEFILEORG |  gawk '{gsub(/^[ \t]+/,"",$2);gsub(/[ \t]+$/,"",$2)}1' | grep -e "^QD_HOME=" -e "^QD_HOME:=" | gawk -F'=' '{print $2}' )
        else
          echo "$MFCLMAKEFILEORG does not exist";exit
        fi 
        echo "HERE174"
        #cat $MFCLMAKEFILEORG |  gawk '{gsub(/^[ \t]+/,"",$2);gsub(/[ \t]+$/,"",$2)}1' 
        echo "176 $MFCLMAKEFILEORG"
        echo "QDORG=$QDORG"
        #exit
        #OPENBLASORG="/home/mfcl/MFCL/libs/OpenBLAS-0.2.19/"
        OPENBLASORG=$(cat $MFCLMAKEFILEORG |  gawk '{$1=$1;print}' | grep -e "^OPENBLAS_HOME=" -e "^OPENBLAS_HOME:=" | gawk -F'=' '{print $2}' )
        OPENBLASNEW="${LIB}/OpenBLAS-${OPENBLASVER}_${OPENBLASOPT}/"
        echo "Original QD_HOME=$QDORG"
        echo "New QD_HOME=$QDNEW"
        echo "Original OPENBLAS_HOME=$OPENBLASORG"
        echo "New OPENBLAS_HOME=$OPENBLASNEW"
        ADMBHOMEORG=$(cat "$MFCLMAKEFILEORG" | gawk '{$1=$1;print}' | grep -e "^ADMB_HOME=" -e "^ADMB_HOME:=" | gawk -F'=' '{print $2}')
     	  #cat mfcl_thrd_linux64_opt.mak |   gawk '{gsub( "$OPENBLASORG", "$OPENBLASNEW" ) ; print $_ }' > tmp2.mak
        cat $MFCLMAKEFILEORG |   gawk "{gsub( "\""$QDORG"\"", "\""$QDNEW"\"" ) ; print \$0 }" | \
          gawk "{gsub( "\""$OPENBLASORG"\"", "\""$OPENBLASNEW"\"" ) ; print \$0 }" | \
	      gawk "{gsub(  "\""$ADMBHOMEORG"\"", "\""$ADMBHOMENEW"\"" ) ; print \$0 }" | \
               gawk "{gsub(  "\""src/\.libs"\"", "\""lib"\"" ) ; print \$0 }" > tmp.mak
	       #    tr  "src/\.libs", "lib"   > tmp.mak
        #exit
        if [ "$OPT" = "O3" ]; then 
          #cat tmp2.mak | sed -e $ONELINER2 > tmp3.mak
          cat tmp.mak | sed -e 's/\-O2/\-O3/g' > mfcl_thrd_build.mak 
        else
	        cat tmp.mak  > mfcl_thrd_build.mak
        fi
        echo "MFCLOPT=$MFCLOPT"
        if [ -f $MFCLEXE ]; then
          rm $MFCLEXE
        fi
        # 	exit
        if [ "$OS" == 'Linux' ]; then
	  make -f mfcl_thrd_build.mak clean
          make -f mfcl_thrd_build.mak  -j$NCPU
        else
          make -f mfcl_thrd_build.mak $MFCLOPT clean
          make -f mfcl_thrd_build.mak $MFCLOPT -j$NCPU
        fi
        RET=$? ;echo RET=$RET
        if [ $RET -ne 0 ]; then
          echo "RET!=0; Build failed"
          exit
        fi
        if [ -f $MFCLEXE ]; then 
          cp -p $MFCLEXE "$TARGETDIR"
          if [ $? -ne 0 ]; then
            make -f mfcl_thrd_build.mak $MFCLOPT clean
            echo "$MFCLEXE could not be copied to $TARGETDIR"
            exit
          else
            make -f mfcl_thrd_build.mak  $MFCLOPT clean
          fi
        else
          make -f mfcl_thrd_build.mak $MFCLOPT clean
          echo "Build failed"
          exit
        fi
      fi 
    done ## End for-loop of $OPT
  done  ## End for-loop of $OPENBLASOPT
done  ## End for-loop of $QDOPT
time=$SECONDS ; echo "$time sec ellapsed"

