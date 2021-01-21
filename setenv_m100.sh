# GNU
module purge ; module load wmlce/1.6.2 gnu/9.3.0 spectrum_mpi/10.3.1--binary
../FoBiS/src/main/python/FoBiS.py build -verbose -f fobos_m100 -mode tests-gnu-debug-mpi -t adam_test_adam_object_mpi.F90

# NVFORTRAN
module purge ; module load wmlce/1.6.2 hpc-sdk spectrum_mpi/10.3.1--binary
../FoBiS/src/main/python/FoBiS.py build -verbose -f fobos_m100 -mode tests-nvf-debug-mpi-cuda -t adam_test_adam_object_mpi_cuda.F90

# XLF
# How to compiler hdf5 1.12.0 with xlf
export CC=xlc_r # o simili...
export CXX=xlc++_r
export FC=xlf2008_r
./configure CFLAGS=-qpic --enable-cxx --enable-fortran --enable-fortran2003 --with-szlib=/cineca/prod/opt/libraries/szip/2.1.1/gnu--8.4.0 --with-zlib=/cineca/prod/opt/libraries/zlib/1.2.11/gnu--8.4.0/ --prefix=/m100_work/cin_staff/fsalvado/HDF5_XLF/install_hdf5-1.12.0

module purge ; module load wmlce/1.6.2 xl/16.1.1--binary spectrum_mpi/10.3.1--binary
export OMPI_FC=xlf2008_r
../FoBiS/src/main/python/FoBiS.py build -verbose -f fobos -mode tests-xlf-debug-mpi -t adam_test_laplace_convect1D.F90
