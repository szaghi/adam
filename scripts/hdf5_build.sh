#!/bin/bash

HDF5_PATH=lib/hdf5/gnu

print_usage () {
   echo
   echo "`basename $0`"
   echo "build HDF5 library"
   echo "usage: `basename $0` cmd"
   echo
   echo "   `basename $0` -build [HDF5_PATH]"
   echo "   build HDF5 and install in HDF5_PATH (default lib/hdf5/gnu)"
   echo
   echo "   `basename $0` -get"
   echo "   get latest release of HDF5 sources"
   echo
   echo "   `basename $0` -h"
   echo "   print this help message"
   echo
   echo "   `basename $0` --help"
   echo "   print this help message"
}

get_hdf5 () {
   echo "get HDF5 sources"
   wget https://github.com/HDFGroup/hdf5/releases/latest/download/hdf5.tar.gz
   tar xf hdf5.tar.gz
   rm -f hdf5.tar.gz
   mv hdf5-* hdf5-src
}

build_hdf5 () {
   here=$(pwd)
   mkdir -p $HDF5_PATH
   cd hdf5-src
   mkdir build
   cd build
   CC=mpicc CXX=mpicxx FC=mpif90 cmake ../ \
     -DCMAKE_INSTALL_PREFIX=$here/$HDF5_PATH \
     -DBUILD_SHARED_LIBS:BOOL=ON \
     -DBUILD_STATIC_LIBS:BOOL=ON \
     -DHDF5_BUILD_FORTRAN:BOOL=ON \
     -DHDF5_ENABLE_PARALLEL:BOOL=ON \
     -DHDF5_ENABLE_NONSTANDARD_FEATURE_FLOAT16:BOOL=OFF \
     -DCMAKE_ANSI_CFLAGS:STRING=-fPIC \
     -DCMAKE_ANSI_FCFLAGS:STRING=-fPIC \
     -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=OFF
   cmake --build .
   make install
   cd ../
   cd ../
   rm -rf hdf5-src
}

if [ $# -eq 0 ] ; then
   print_usage; exit 0
fi
while [ $# -gt 0 ]; do
   case "$1" in
      "-build")
         shift
         if [ $# -gt 0 ] ; then
            HDF5_PATH=$1
         fi
         get_hdf5
         build_hdf5
         exit 0
         ;;
      "-get")
         get_hdf5
         exit 0
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
