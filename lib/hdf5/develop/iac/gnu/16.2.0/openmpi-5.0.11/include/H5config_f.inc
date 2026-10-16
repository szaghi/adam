! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!   Copyright by The HDF Group.                                               *
!   All rights reserved.                                                      *
!                                                                             *
!   This file is part of HDF5.  The full HDF5 copyright notice, including     *
!   terms governing use, modification, and redistribution, is contained in    *
!   the LICENSE file, which can be found at the root of the source code       *
!   distribution tree, or in https://www.hdfgroup.org/licenses.               *
!   If you do not have access to either file, you may request a copy from     *
!   help@hdfgroup.org.                                                        *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! fortran/H5config_f.inc. Generated from fortran/src/H5config_f.inc.cmake by CMake

! Define if there is parallel support
#define H5_HAVE_PARALLEL

! Define if MPI supports mpi_f08 module
#define H5_HAVE_MPI_F08

! Define if there is subfiling support
/* #undef H5_HAVE_SUBFILING_VFD */

! Define if on APPLE
/* #undef H5_HAVE_DARWIN */

! Define if the intrinsic function STORAGE_SIZE exists
#define H5_FORTRAN_HAVE_STORAGE_SIZE

! Define if the intrinsic function SIZEOF exists
#define H5_FORTRAN_HAVE_SIZEOF

! Define if the intrinsic function C_SIZEOF exists
#define H5_FORTRAN_HAVE_C_SIZEOF

! Define if allocatable character is supported
#define H5_FORTRAN_HAVE_CHAR_ALLOC

! Define if the intrinsic function C_LONG_DOUBLE exists
#define H5_FORTRAN_HAVE_C_LONG_DOUBLE

! Define if Fortran C_LONG_DOUBLE is different from C_DOUBLE
#define H5_FORTRAN_C_LONG_DOUBLE_IS_UNIQUE

! Define if Fortran C_BOOL is different from default LOGICAL
#define H5_FORTRAN_C_BOOL_IS_UNIQUE

! Define  MPI Fortran KIND of LOGICAL
#define H5_MPI_LOGICAL_KIND 4

! Define if Fortran supports ISO_FORTRAN_ENV (F08)
#define H5_HAVE_ISO_FORTRAN_ENV

! Define the size of C's double
#define H5_SIZEOF_DOUBLE 8

! Define the size of C's long double
#define H5_SIZEOF_LONG_DOUBLE 16

! Define the maximum decimal precision for reals
#define H5_PAC_FC_MAX_REAL_PRECISION 33

! If C has quad precision
#define H5_HAVE_FLOAT128

! Define if INTEGER*16 is available
#define H5_HAVE_Fortran_INTEGER_SIZEOF_16

! Maximum decimal precision for C
#define H5_PAC_C_MAX_REAL_PRECISION 33

! number of valid REAL KINDs
#define H5_H5CONFIG_F_NUM_RKIND INTEGER, PARAMETER :: num_rkinds = 4

! valid REAL KINDs (need to have a matching C counter-part)
#define H5_H5CONFIG_F_RKIND INTEGER, DIMENSION(1:num_rkinds) :: rkind = (/4,8,10,16/)

! valid REAL KINDs (need to have a matching C counter-part)
#define H5_H5CONFIG_F_RKIND_SIZEOF INTEGER, DIMENSION(1:num_rkinds) :: rkind_sizeof = (/4,8,16,16/)

! number of valid INTEGER KINDs
#define H5_H5CONFIG_F_NUM_IKIND INTEGER, PARAMETER :: num_ikinds = 5

! valid INTEGER KINDs (need to have a matching C counter-part)
#define H5_H5CONFIG_F_IKIND INTEGER, DIMENSION(1:num_ikinds) :: ikind = (/1,2,4,8,16/)

! Fortran compiler id
#define H5_Fortran_COMPILER_ID GNU

! Define if deprecated public API symbols are disabled
/* #undef H5_NO_DEPRECATED_SYMBOLS */

! For major interface/format changes
#define H5_VERS_MAJOR 2

! For minor interface/format changes
#define H5_VERS_MINOR 3

! For tweaks, bug-fixes, or development
#define H5_VERS_RELEASE 0

! macros for comparing versions
#define H5_VERSION_GE(Maj, Min, Rel)                                                   \
    (((H5_VERS_MAJOR == Maj) && (H5_VERS_MINOR == Min) && (H5_VERS_RELEASE >= Rel)) || \
     ((H5_VERS_MAJOR == Maj) && (H5_VERS_MINOR > Min)) || (H5_VERS_MAJOR > Maj))

#define H5_VERSION_LE(Maj, Min, Rel)                                                   \
    (((H5_VERS_MAJOR == Maj) && (H5_VERS_MINOR == Min) && (H5_VERS_RELEASE <= Rel)) || \
     ((H5_VERS_MAJOR == Maj) && (H5_VERS_MINOR < Min)) || (H5_VERS_MAJOR < Maj))

