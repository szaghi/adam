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
   CFLAGS="-fPIC" FCFLAGS="-fPIC" ./configure --prefix=$here/$HDF5_PATH --enable-shared --enable-parallel --enable-fortran --disable-libtool-lock FC=mpif90
   make
   make install
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
