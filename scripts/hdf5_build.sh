#!/bin/bash

# helper to build HDF5
# automatic download anb build szip, zlib and HDF5 libraries
# 2025-05, Stefano Zaghi, stefano.zaghi@gmail.com

# global var defaults
HERE=$(pwd)           # paths prefix
IS_PATHS_ABS=no       # use absolute paths, relative is assumed by default (HERE=pwd)
LSRC_PATH=libhdf5_src # path of sources
HDF5_PATH=libhdf5/    # path of built libraries
USE_CMAKE=yes         # use CMAKE to build HDF5 (default), if set to no use autotools

print_usage () {
   echo
   echo "`basename $0`"
   echo "build HDF5 library"
   echo "usage: `basename $0` cmd [opts]"
   echo
   echo "cmds:"
   echo "   -get"
   echo "      `basename $0` -get"
   echo "      get latest release of HDF5 (within szip/zlib)"
   echo "   -build"
   echo "      `basename $0` -build"
   echo "      build HDF5 (within szip/zlib)"
   echo "   -h/--help"
   echo "      `basename $0` -h"
   echo "      print this help message"
   echo "      `basename $0` --help"
   echo "      print this help message"
   echo "opts:"
   echo "   -lsrc LSRC_PATH"
   echo "      specify path where sources are placed"
   echo "      default LSRC_PATH=libhdf5_src (relative to where script is run)"
   echo "       + `basename $0` -get -lsrc my_cumbersome_src_path"
   echo "         get latest release of HDF5 (within szip/zlib) and place in my_cumbersome_src_path"
   echo "       + `basename $0` -build -lsrc my_cumbersome_src_path"
   echo "         build HDF5 (within szip/zlib) and place sources in my_cumbersome_src_path"
   echo
   echo "   -hdf5 HDF5_PATH"
   echo "      specify path where built libraries are placed"
   echo "      default HDF5_PATH=libhdf5 (relative to where script is run)"
   echo "       + `basename $0` -build -hdf5 my_cumbersome_path"
   echo "         build HDF5 (within szip/zlib) and place sources in my_cumbersome_path"
   echo
   echo "   -abs-paths"
   echo "      enable absolute paths, default are relative to where script is run"
   echo "       + `basename $0` -get -lsrc /my_absolute_cumbersome_src_path -abs-paths"
   echo "         get latest release of HDF5 (within szip/zlib) and place in /my_absolute_cumbersome_src_path"
   echo "       + `basename $0` -build -lsrc /my_absolute_cumbersome_src_path -abs-paths"
   echo "         build HDF5 (within szip/zlib) and place sources in /my_absolute_cumbersome_src_path"
   echo
   echo "   -use-autot"
   echo "      disable CMAKE and use autotools to build HDF5"
   echo "       + `basename $0` -build -use-autot"
   echo "         build HDF5 (within szip/zlib) using autotools instead of CMAKE"
   echo
   echo "examples:"
   echo "   + `basename $0` -get"
   echo "     get latest release of HDF5 (within szip/zlib) and place in libhdf5_src"
   echo
   echo "   + `basename $0` -build"
   echo "     build HDF5 (within szip/zlib) and install in libhdf5 (relative to where script is run)"
   echo
   echo "   + `basename $0` -get -lsrc libhdf5_github_src"
   echo "     get latest release of HDF5 (within szip/zlib) and place in libhdf5_github_src"
   echo
   echo "   + `basename $0` -build -hdf5 libhdf5_github -lsrc libhdf5_github_src"
   echo "     build HDF5 (within szip/zlib) and install in libhdf5_github and place sources"
   echo "     in libhdf5_github_src (relative to where script is run)"
   echo
   echo "   + `basename $0` -build -hdf5 /opt/libhdf5/gnu -abs-paths"
   echo "     build HDF5 (within szip/zlib) and install in /opt/libhdf5/gnu"
   echo
   echo "   + `basename $0` -build -hdf5 /opt/libhdf5/gnu -abs-paths -use-autot"
   echo "     build HDF5 (within szip/zlib) and install in /opt/libhdf5/gnu using autotools instead of CMAKE"
   echo
   echo "   + `basename $0` -h"
   echo "     print this help message"
   echo
   echo "   + `basename $0` --help"
   echo "     print this help message"
}

get_libaec () {
   mkdir -p $HERE/$LSRC_PATH
   echo "get LIBAEC sources"
   wget https://github.com/MathisRosenhauer/libaec/releases/download/v1.1.3/libaec-1.1.3.tar.gz
   tar xf libaec-1.1.3.tar.gz
   rm -f libaec-1.1.3.tar.gz
   mv libaec-* $HERE/$LSRC_PATH/libaec
}

get_szip () {
   mkdir -p $HERE/$LSRC_PATH
   echo "get SZIP sources"
   wget https://github.com/erdc/szip/archive/refs/heads/master.zip
   unzip master.zip
   rm -f master.zip
   mv szip-master $HERE/$LSRC_PATH/szip
}

get_zlib () {
   mkdir -p $HERE/$LSRC_PATH
   echo "get ZLIB sources"
   wget https://github.com/madler/zlib/releases/download/v1.3.1/zlib-1.3.1.tar.gz
   tar xf zlib-1.3.1.tar.gz
   rm -f zlib-1.3.1.tar.gz
   mv zlib-* $HERE/$LSRC_PATH/zlib
}

get_hdf5 () {
   mkdir -p $HERE/$LSRC_PATH
   echo "get HDF5 sources"
   wget https://github.com/HDFGroup/hdf5/releases/latest/download/hdf5.tar.gz
   tar xf hdf5.tar.gz
   rm -f hdf5.tar.gz
   mv hdf5-* $HERE/$LSRC_PATH/hdf5
}

build_libaec () {
   echo "build LIBAEC"
   mkdir -p $HERE/$HDF5_PATH/libaec
   cd $HERE/$LSRC_PATH/libaec
   CC=gcc CXX=g++ FC=gfortran CFLAGS='-O3' CXXFLAGS='-O3' FCFLAGS='-O3' ./configure --prefix=$HERE/$HDF5_PATH/libaec
   make
   make check
   make install
   cd -
}

build_szip () {
   echo "build SZIP"
   mkdir -p $HERE/$HDF5_PATH/szip
   cd $HERE/$LSRC_PATH/szip
   CC=gcc CXX=g++ FC=gfortran CFLAGS='-O3' CXXFLAGS='-O3' FCFLAGS='-O3' ./configure --prefix=$HERE/$HDF5_PATH/szip
   make
   make check
   make install
   cd -
}

build_zlib () {
   echo "build ZLIB"
   mkdir -p $HERE/$HDF5_PATH/zlib
   cd $HERE/$LSRC_PATH/zlib
   CC=gcc CFLAGS='-O3 -fPIC' ./configure --prefix=$HERE/$HDF5_PATH/zlib
   make
   make check
   make install
   cd -
}

build_hdf5_by_autot () {
   echo "build HDF5 by autotools"
   mkdir -p $HERE/$HDF5_PATH
   cd $HERE/$LSRC_PATH/hdf5
   CFLAGS="-fPIC" FCFLAGS="-fPIC" ./configure \
      --prefix=$HERE/$HDF5_PATH               \
      --enable-fortran                        \
      --enable-optimization=high              \
      --enable-shared=no                      \
      --enable-static=yes                     \
      --enable-parallel                       \
      --with-zlib=$HERE/$HDF5_PATH/zlib       \
      --with-szip=$HERE/$HDF5_PATH/szip
   make -j 4
   make install
   cd -
}

build_hdf5_by_cmake () {
   echo "build HDF5 by cmake"
   mkdir -p $HERE/$HDF5_PATH
   mkdir -p $HERE/$LSRC_PATH/hdf5/build
   cd $HERE/$LSRC_PATH/hdf5/build
   CFLAGS="-fPIC" FCFLAGS="-fPIC" cmake ../                       \
      -DCMAKE_INSTALL_PREFIX:PATH=$HERE/$HDF5_PATH                \
      -DCMAKE_BUILD_TYPE:STRING=Release                           \
      -DBUILD_SHARED_LIBS:BOOL=OFF                                \
      -DBUILD_STATIC_LIBS:BOOL=ON                                 \
      -DHDF5_ENABLE_PARALLEL:BOOL=ON                              \
      -DHDF5_BUILD_FORTRAN:BOOL=ON                                \
      -DZLIB_LIBRARY:FILEPATH=$HERE/$HDF5_PATH/zlib/lib/libz.a    \
      -DZLIB_INCLUDE_DIR:PATH=$HERE/$HDF5_PATH/zlib/include       \
      -DZLIB_USE_EXTERNAL:BOOL=OFF                                \
      -DSZIP_LIBRARY:FILEPATH=$HERE/$HDF5_PATH/szip/lib/libsz.a   \
      -DSZIP_INCLUDE_DIR:PATH=$HERE/$HDF5_PATH/szip/include       \
      -Dlibaec_LIBRARY:FILEPATH=$HERE/$HDF5_PATH/szip/lib/libsz.a \
      -Dlibaec_INCLUDE_DIR:PATH=$HERE/$HDF5_PATH/szip/include     \
      -DSZIP_USE_EXTERNAL:BOOL=OFF
   cmake --build . --config Release
   make install
   cd -
}

# parsing command line
if [ $# -eq 0 ] ; then
   print_usage; exit 0
fi
DO_BUILD=no
DO_GET=no
while [ $# -gt 0 ]; do
   case "$1" in
      "-get")
         DO_GET=yes
         ;;
      "-build")
         DO_GET=yes
         DO_BUILD=yes
         ;;
      "-hdf5")
         shift
         HDF5_PATH=$1
         ;;
      "-lsrc")
         shift
         LSRC_PATH=$1
         ;;
      "-abs-paths")
         IS_PATHS_ABS=yes
         ;;
      "-use-autot")
         USE_CMAKE=no
         ;;
      "-h")
         print_usage; exit 0
         ;;
      "--help")
         print_usage; exit 0
         ;;
      *)
         echo; echo "unknown command $1"; print_usage; exit 1
         ;;
   esac
   shift
done

if [ "$IS_PATHS_ABS" == "yes" ] ; then
   HERE=''
fi
if [ "$DO_GET" == "yes" ] ; then
   get_szip
   get_zlib
   get_hdf5
fi
if [ "$DO_BUILD" == "yes" ] ; then
   build_szip
   build_zlib
  if [ "$USE_CMAKE" == "yes" ] ; then
     build_hdf5_by_cmake
   else
     build_hdf5_by_autot
  fi
fi

# while [ $# -gt 0 ]; do
#    case "$1" in
#       "-build")
#          shift
#          if [ $# -gt 0 ] ; then
#             while [ $# -gt 0 ]; do
#                case "$1" in
#                   "-hdf5")
#                      shift
#                      HDF5_PATH=$1
#                      ;;
#                   "-lsrc")
#                      shift
#                      LSRC_PATH=$1
#                      ;;
#                   *)
#                      echo; echo "unknown command $1"; print_usage; exit 1
#                      ;;
#                esac
#                shift
#             done
#          fi
#          # get_libaec
#          get_szip
#          get_zlib
#          get_hdf5
#          # build_libaec
#          build_szip
#          build_zlib
#          # build_hdf5_by_autot
#          build_hdf5_by_cmake
#          exit 0
#          ;;
#       "-get")
#          shift
#          if [ $# -gt 0 ] ; then
#             while [ $# -gt 0 ]; do
#                case "$1" in
#                   "-lsrc")
#                      shift
#                      LSRC_PATH=$1
#                      ;;
#                   *)
#                      echo; echo "unknown command $1"; print_usage; exit 1
#                      ;;
#                esac
#                shift
#             done
#          fi
#          get_szip
#          get_zlib
#          get_hdf5
#          exit 0
#          ;;
#       "-h")
#          print_usage; exit 0
#          ;;
#       "--help")
#          print_usage; exit 0
#          ;;
#       *)
#          echo; echo "unknown command $1"; print_usage; exit 1
#          ;;
#    esac
#    shift
# done
