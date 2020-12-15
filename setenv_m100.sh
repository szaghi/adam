# GNU
module purge ; module load wmlce/1.6.2 gnu/9.3.0 spectrum_mpi/10.3.1--binary
../FoBiS/src/main/python/FoBiS.py build -verbose -f fobos_m100 -mode tests-gnu-debug-mpi -t adam_test_adam_object_mpi.F90

# NVFORTRAN
module purge ; module load wmlce/1.6.2 hpc-sdk spectrum_mpi/10.3.1--binary
../FoBiS/src/main/python/FoBiS.py build -verbose -f fobos_m100 -mode tests-nvf-debug-mpi-cuda -t adam_test_adam_object_mpi_cuda.F90
