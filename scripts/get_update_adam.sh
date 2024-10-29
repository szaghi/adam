#!/bin/bash

print_usage () {
   echo
   echo "`basename $0`"
   echo "get and update ADAM sources"
   echo "usage: `basename $0` cmd"
   echo
   echo "   `basename $0` -get"
   echo "   get and update ADAM sources"
   echo
   echo "   `basename $0` -update"
   echo "   only update ADAM sources"
   echo
   echo "   `basename $0` -h"
   echo "   print this help message"
   echo
   echo "   `basename $0` --help"
   echo "   print this help message"
}

get_adam () {
   echo "get ADAM sources"
   git clone git@github.com:szaghi/adam.git
   cd adam
   update_adam
}

update_adam () {
   echo "update ADAM sources"
   git pull
   for dir in src/third_party/*; do
      if [ -d "$dir" ]; then
         git submodule update --init $dir
      fi
   done
}

if [ $# -eq 1 ] ; then
   case "$1" in
      "-get")
         get_adam; exit 0
         ;;
      "-update")
         update_adam; exit 0
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
else
   echo "error: run `basename $0` with exactly one command"; print_usage; exit 2
fi
