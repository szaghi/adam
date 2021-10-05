#!/bin/bash

# M100 compilation
module load gnu/8.4.0
g++ -std=c++17 -frounding-math -O2 -I/cineca/prod/opt/libraries/boost/1.76.0/spectrum_mpi--10.4.0--binary/include -I../../../lib/CGAL-5.0.2-INSTALL/include/ -I /cineca/prod/opt/libraries/mpfr/4.0.2/gnu--8.4.0/include -c cgal_c_wrappers.cpp

# Local compilation
#g++ -std=c++17 -frounding-math -O2 -I/opt/cgal/5.2.1/include/ -c cgal_c_wrappers.cpp
