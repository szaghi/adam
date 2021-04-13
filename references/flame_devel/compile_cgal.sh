#!/bin/bash
module load gnu/8.4.0
g++ -std=c++17 -frounding-math -O2 -I/m100_scratch/userinternal/fsalvado/CGAL/INSTALL_CGAL-5.0.2/include/ -I /cineca/prod/opt/libraries/mpfr/4.0.2/gnu--8.4.0/include -c cgal_c_wrappers.cpp
