gcc -c cm.c
gfortran -c cm_wrapper.F90 
gfortran -c cm_wrapper_test.F90 
gfortran cm.o cm_wrapper.o cm_wrapper_test.o 
