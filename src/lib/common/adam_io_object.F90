!< ADAM, IO class definition.
module adam_io_object
!< ADAM, IO class definition.

use adam_field_object
use adam_grid_object
use adam_mpih_object
use motion
use penf
use stringifor
! use mpi

implicit none
private
public :: io_object

type :: io_object
   !< ADAM class definition.
   type(mpih_object)           :: mpih          !< The MPI handler.
   type(grid_object),  pointer :: grid=>null()  !< The grid.
   type(field_object), pointer :: field=>null() !< The field.
   ! auxiliary fields data, pointer to user data
   ! vector data (q) have [nv,ni,nj.nk,nb] dimensions
   ! scalar data (s) have [   ni,nj.nk,nb] dimensions
   real(R8P),    pointer :: q1_R8P(:,:,:,:,:)=>null() !< Auxiliary (1) vector cell centered fields [nv,ni,nj,nk,nb], kind R8P.
   real(R8P),    pointer :: q2_R8P(:,:,:,:,:)=>null() !< Auxiliary (2) vector cell centered fields [nv,ni,nj,nk,nb], kind R8P.
   real(R8P),    pointer :: q3_R8P(:,:,:,:,:)=>null() !< Auxiliary (3) vector cell centered fields [nv,ni,nj,nk,nb], kind R8P.
   real(R8P),    pointer :: q4_R8P(:,:,:,:,:)=>null() !< Auxiliary (4) vector cell centered fields [nv,ni,nj,nk,nb], kind R8P.
   real(R8P),    pointer :: q5_R8P(:,:,:,:,:)=>null() !< Auxiliary (5) vector cell centered fields [nv,ni,nj,nk,nb], kind R8P.
   real(R4P),    pointer :: q1_R4P(:,:,:,:,:)=>null() !< Auxiliary (1) vector cell centered fields [nv,ni,nj,nk,nb], kind R4P.
   real(R4P),    pointer :: q2_R4P(:,:,:,:,:)=>null() !< Auxiliary (2) vector cell centered fields [nv,ni,nj,nk,nb], kind R4P.
   real(R4P),    pointer :: q3_R4P(:,:,:,:,:)=>null() !< Auxiliary (3) vector cell centered fields [nv,ni,nj,nk,nb], kind R4P.
   real(R4P),    pointer :: q4_R4P(:,:,:,:,:)=>null() !< Auxiliary (4) vector cell centered fields [nv,ni,nj,nk,nb], kind R4P.
   real(R4P),    pointer :: q5_R4P(:,:,:,:,:)=>null() !< Auxiliary (5) vector cell centered fields [nv,ni,nj,nk,nb], kind R4P.
   integer(I8P), pointer :: q1_I8P(:,:,:,:,:)=>null() !< Auxiliary (1) vector cell centered fields [nv,ni,nj,nk,nb], kind I8P.
   integer(I8P), pointer :: q2_I8P(:,:,:,:,:)=>null() !< Auxiliary (2) vector cell centered fields [nv,ni,nj,nk,nb], kind I8P.
   integer(I8P), pointer :: q3_I8P(:,:,:,:,:)=>null() !< Auxiliary (3) vector cell centered fields [nv,ni,nj,nk,nb], kind I8P.
   integer(I8P), pointer :: q4_I8P(:,:,:,:,:)=>null() !< Auxiliary (4) vector cell centered fields [nv,ni,nj,nk,nb], kind I8P.
   integer(I8P), pointer :: q5_I8P(:,:,:,:,:)=>null() !< Auxiliary (5) vector cell centered fields [nv,ni,nj,nk,nb], kind I8P.
   integer(I4P), pointer :: q1_I4P(:,:,:,:,:)=>null() !< Auxiliary (1) vector cell centered fields [nv,ni,nj,nk,nb], kind I4P.
   integer(I4P), pointer :: q2_I4P(:,:,:,:,:)=>null() !< Auxiliary (2) vector cell centered fields [nv,ni,nj,nk,nb], kind I4P.
   integer(I4P), pointer :: q3_I4P(:,:,:,:,:)=>null() !< Auxiliary (3) vector cell centered fields [nv,ni,nj,nk,nb], kind I4P.
   integer(I4P), pointer :: q4_I4P(:,:,:,:,:)=>null() !< Auxiliary (4) vector cell centered fields [nv,ni,nj,nk,nb], kind I4P.
   integer(I4P), pointer :: q5_I4P(:,:,:,:,:)=>null() !< Auxiliary (5) vector cell centered fields [nv,ni,nj,nk,nb], kind I4P.
   integer(I2P), pointer :: q1_I2P(:,:,:,:,:)=>null() !< Auxiliary (1) vector cell centered fields [nv,ni,nj,nk,nb], kind I2P.
   integer(I2P), pointer :: q2_I2P(:,:,:,:,:)=>null() !< Auxiliary (2) vector cell centered fields [nv,ni,nj,nk,nb], kind I2P.
   integer(I2P), pointer :: q3_I2P(:,:,:,:,:)=>null() !< Auxiliary (3) vector cell centered fields [nv,ni,nj,nk,nb], kind I2P.
   integer(I2P), pointer :: q4_I2P(:,:,:,:,:)=>null() !< Auxiliary (4) vector cell centered fields [nv,ni,nj,nk,nb], kind I2P.
   integer(I2P), pointer :: q5_I2P(:,:,:,:,:)=>null() !< Auxiliary (5) vector cell centered fields [nv,ni,nj,nk,nb], kind I2P.
   integer(I1P), pointer :: q1_I1P(:,:,:,:,:)=>null() !< Auxiliary (1) vector cell centered fields [nv,ni,nj,nk,nb], kind I1P.
   integer(I1P), pointer :: q2_I1P(:,:,:,:,:)=>null() !< Auxiliary (2) vector cell centered fields [nv,ni,nj,nk,nb], kind I1P.
   integer(I1P), pointer :: q3_I1P(:,:,:,:,:)=>null() !< Auxiliary (3) vector cell centered fields [nv,ni,nj,nk,nb], kind I1P.
   integer(I1P), pointer :: q4_I1P(:,:,:,:,:)=>null() !< Auxiliary (4) vector cell centered fields [nv,ni,nj,nk,nb], kind I1P.
   integer(I1P), pointer :: q5_I1P(:,:,:,:,:)=>null() !< Auxiliary (5) vector cell centered fields [nv,ni,nj,nk,nb], kind I1P.
   real(R8P),    pointer :: s1_R8P(  :,:,:,:)=>null() !< Auxiliary (1) scalar cell centered fields [   ni,nj,nk,nb], kind R8P.
   real(R8P),    pointer :: s2_R8P(  :,:,:,:)=>null() !< Auxiliary (2) scalar cell centered fields [   ni,nj,nk,nb], kind R8P.
   real(R8P),    pointer :: s3_R8P(  :,:,:,:)=>null() !< Auxiliary (3) scalar cell centered fields [   ni,nj,nk,nb], kind R8P.
   real(R8P),    pointer :: s4_R8P(  :,:,:,:)=>null() !< Auxiliary (4) scalar cell centered fields [   ni,nj,nk,nb], kind R8P.
   real(R8P),    pointer :: s5_R8P(  :,:,:,:)=>null() !< Auxiliary (5) scalar cell centered fields [   ni,nj,nk,nb], kind R8P.
   real(R4P),    pointer :: s1_R4P(  :,:,:,:)=>null() !< Auxiliary (1) scalar cell centered fields [   ni,nj,nk,nb], kind R4P.
   real(R4P),    pointer :: s2_R4P(  :,:,:,:)=>null() !< Auxiliary (2) scalar cell centered fields [   ni,nj,nk,nb], kind R4P.
   real(R4P),    pointer :: s3_R4P(  :,:,:,:)=>null() !< Auxiliary (3) scalar cell centered fields [   ni,nj,nk,nb], kind R4P.
   real(R4P),    pointer :: s4_R4P(  :,:,:,:)=>null() !< Auxiliary (4) scalar cell centered fields [   ni,nj,nk,nb], kind R4P.
   real(R4P),    pointer :: s5_R4P(  :,:,:,:)=>null() !< Auxiliary (5) scalar cell centered fields [   ni,nj,nk,nb], kind R4P.
   integer(I8P), pointer :: s1_I8P(  :,:,:,:)=>null() !< Auxiliary (1) scalar cell centered fields [   ni,nj,nk,nb], kind I8P.
   integer(I8P), pointer :: s2_I8P(  :,:,:,:)=>null() !< Auxiliary (2) scalar cell centered fields [   ni,nj,nk,nb], kind I8P.
   integer(I8P), pointer :: s3_I8P(  :,:,:,:)=>null() !< Auxiliary (3) scalar cell centered fields [   ni,nj,nk,nb], kind I8P.
   integer(I8P), pointer :: s4_I8P(  :,:,:,:)=>null() !< Auxiliary (4) scalar cell centered fields [   ni,nj,nk,nb], kind I8P.
   integer(I8P), pointer :: s5_I8P(  :,:,:,:)=>null() !< Auxiliary (5) scalar cell centered fields [   ni,nj,nk,nb], kind I8P.
   integer(I4P), pointer :: s1_I4P(  :,:,:,:)=>null() !< Auxiliary (1) scalar cell centered fields [   ni,nj,nk,nb], kind I4P.
   integer(I4P), pointer :: s2_I4P(  :,:,:,:)=>null() !< Auxiliary (2) scalar cell centered fields [   ni,nj,nk,nb], kind I4P.
   integer(I4P), pointer :: s3_I4P(  :,:,:,:)=>null() !< Auxiliary (3) scalar cell centered fields [   ni,nj,nk,nb], kind I4P.
   integer(I4P), pointer :: s4_I4P(  :,:,:,:)=>null() !< Auxiliary (4) scalar cell centered fields [   ni,nj,nk,nb], kind I4P.
   integer(I4P), pointer :: s5_I4P(  :,:,:,:)=>null() !< Auxiliary (5) scalar cell centered fields [   ni,nj,nk,nb], kind I4P.
   integer(I2P), pointer :: s1_I2P(  :,:,:,:)=>null() !< Auxiliary (1) scalar cell centered fields [   ni,nj,nk,nb], kind I2P.
   integer(I2P), pointer :: s2_I2P(  :,:,:,:)=>null() !< Auxiliary (2) scalar cell centered fields [   ni,nj,nk,nb], kind I2P.
   integer(I2P), pointer :: s3_I2P(  :,:,:,:)=>null() !< Auxiliary (3) scalar cell centered fields [   ni,nj,nk,nb], kind I2P.
   integer(I2P), pointer :: s4_I2P(  :,:,:,:)=>null() !< Auxiliary (4) scalar cell centered fields [   ni,nj,nk,nb], kind I2P.
   integer(I2P), pointer :: s5_I2P(  :,:,:,:)=>null() !< Auxiliary (5) scalar cell centered fields [   ni,nj,nk,nb], kind I2P.
   integer(I1P), pointer :: s1_I1P(  :,:,:,:)=>null() !< Auxiliary (1) scalar cell centered fields [   ni,nj,nk,nb], kind I1P.
   integer(I1P), pointer :: s2_I1P(  :,:,:,:)=>null() !< Auxiliary (2) scalar cell centered fields [   ni,nj,nk,nb], kind I1P.
   integer(I1P), pointer :: s3_I1P(  :,:,:,:)=>null() !< Auxiliary (3) scalar cell centered fields [   ni,nj,nk,nb], kind I1P.
   integer(I1P), pointer :: s4_I1P(  :,:,:,:)=>null() !< Auxiliary (4) scalar cell centered fields [   ni,nj,nk,nb], kind I1P.
   integer(I1P), pointer :: s5_I1P(  :,:,:,:)=>null() !< Auxiliary (5) scalar cell centered fields [   ni,nj,nk,nb], kind I1P.
   ! auxiliary fields names
   type(string), allocatable :: q1_R8P_name(:) !< Auxiliary (1) vector fields names [nv], kind R8P.
   type(string), allocatable :: q2_R8P_name(:) !< Auxiliary (2) vector fields names [nv], kind R8P.
   type(string), allocatable :: q3_R8P_name(:) !< Auxiliary (3) vector fields names [nv], kind R8P.
   type(string), allocatable :: q4_R8P_name(:) !< Auxiliary (4) vector fields names [nv], kind R8P.
   type(string), allocatable :: q5_R8P_name(:) !< Auxiliary (5) vector fields names [nv], kind R8P.
   type(string), allocatable :: q1_R4P_name(:) !< Auxiliary (1) vector fields names [nv], kind R4P.
   type(string), allocatable :: q2_R4P_name(:) !< Auxiliary (2) vector fields names [nv], kind R4P.
   type(string), allocatable :: q3_R4P_name(:) !< Auxiliary (3) vector fields names [nv], kind R4P.
   type(string), allocatable :: q4_R4P_name(:) !< Auxiliary (4) vector fields names [nv], kind R4P.
   type(string), allocatable :: q5_R4P_name(:) !< Auxiliary (5) vector fields names [nv], kind R4P.
   type(string), allocatable :: q1_I8P_name(:) !< Auxiliary (1) vector fields names [nv], kind I8P.
   type(string), allocatable :: q2_I8P_name(:) !< Auxiliary (2) vector fields names [nv], kind I8P.
   type(string), allocatable :: q3_I8P_name(:) !< Auxiliary (3) vector fields names [nv], kind I8P.
   type(string), allocatable :: q4_I8P_name(:) !< Auxiliary (4) vector fields names [nv], kind I8P.
   type(string), allocatable :: q5_I8P_name(:) !< Auxiliary (5) vector fields names [nv], kind I8P.
   type(string), allocatable :: q1_I4P_name(:) !< Auxiliary (1) vector fields names [nv], kind I4P.
   type(string), allocatable :: q2_I4P_name(:) !< Auxiliary (2) vector fields names [nv], kind I4P.
   type(string), allocatable :: q3_I4P_name(:) !< Auxiliary (3) vector fields names [nv], kind I4P.
   type(string), allocatable :: q4_I4P_name(:) !< Auxiliary (4) vector fields names [nv], kind I4P.
   type(string), allocatable :: q5_I4P_name(:) !< Auxiliary (5) vector fields names [nv], kind I4P.
   type(string), allocatable :: q1_I2P_name(:) !< Auxiliary (1) vector fields names [nv], kind I2P.
   type(string), allocatable :: q2_I2P_name(:) !< Auxiliary (2) vector fields names [nv], kind I2P.
   type(string), allocatable :: q3_I2P_name(:) !< Auxiliary (3) vector fields names [nv], kind I2P.
   type(string), allocatable :: q4_I2P_name(:) !< Auxiliary (4) vector fields names [nv], kind I2P.
   type(string), allocatable :: q5_I2P_name(:) !< Auxiliary (5) vector fields names [nv], kind I2P.
   type(string), allocatable :: q1_I1P_name(:) !< Auxiliary (1) vector fields names [nv], kind I1P.
   type(string), allocatable :: q2_I1P_name(:) !< Auxiliary (2) vector fields names [nv], kind I1P.
   type(string), allocatable :: q3_I1P_name(:) !< Auxiliary (3) vector fields names [nv], kind I1P.
   type(string), allocatable :: q4_I1P_name(:) !< Auxiliary (4) vector fields names [nv], kind I1P.
   type(string), allocatable :: q5_I1P_name(:) !< Auxiliary (5) vector fields names [nv], kind I1P.
   type(string)              :: s1_R8P_name    !< Auxiliary (1) scalar field name, kind R8P.
   type(string)              :: s2_R8P_name    !< Auxiliary (2) scalar field name, kind R8P.
   type(string)              :: s3_R8P_name    !< Auxiliary (3) scalar field name, kind R8P.
   type(string)              :: s4_R8P_name    !< Auxiliary (4) scalar field name, kind R8P.
   type(string)              :: s5_R8P_name    !< Auxiliary (5) scalar field name, kind R8P.
   type(string)              :: s1_R4P_name    !< Auxiliary (1) scalar field name, kind R4P.
   type(string)              :: s2_R4P_name    !< Auxiliary (2) scalar field name, kind R4P.
   type(string)              :: s3_R4P_name    !< Auxiliary (3) scalar field name, kind R4P.
   type(string)              :: s4_R4P_name    !< Auxiliary (4) scalar field name, kind R4P.
   type(string)              :: s5_R4P_name    !< Auxiliary (5) scalar field name, kind R4P.
   type(string)              :: s1_I8P_name    !< Auxiliary (1) scalar field name, kind I8P.
   type(string)              :: s2_I8P_name    !< Auxiliary (2) scalar field name, kind I8P.
   type(string)              :: s3_I8P_name    !< Auxiliary (3) scalar field name, kind I8P.
   type(string)              :: s4_I8P_name    !< Auxiliary (4) scalar field name, kind I8P.
   type(string)              :: s5_I8P_name    !< Auxiliary (5) scalar field name, kind I8P.
   type(string)              :: s1_I4P_name    !< Auxiliary (1) scalar field name, kind I4P.
   type(string)              :: s2_I4P_name    !< Auxiliary (2) scalar field name, kind I4P.
   type(string)              :: s3_I4P_name    !< Auxiliary (3) scalar field name, kind I4P.
   type(string)              :: s4_I4P_name    !< Auxiliary (4) scalar field name, kind I4P.
   type(string)              :: s5_I4P_name    !< Auxiliary (5) scalar field name, kind I4P.
   type(string)              :: s1_I2P_name    !< Auxiliary (1) scalar field name, kind I2P.
   type(string)              :: s2_I2P_name    !< Auxiliary (2) scalar field name, kind I2P.
   type(string)              :: s3_I2P_name    !< Auxiliary (3) scalar field name, kind I2P.
   type(string)              :: s4_I2P_name    !< Auxiliary (4) scalar field name, kind I2P.
   type(string)              :: s5_I2P_name    !< Auxiliary (5) scalar field name, kind I2P.
   type(string)              :: s1_I1P_name    !< Auxiliary (1) scalar field name, kind I1P.
   type(string)              :: s2_I1P_name    !< Auxiliary (2) scalar field name, kind I1P.
   type(string)              :: s3_I1P_name    !< Auxiliary (3) scalar field name, kind I1P.
   type(string)              :: s4_I1P_name    !< Auxiliary (4) scalar field name, kind I1P.
   type(string)              :: s5_I1P_name    !< Auxiliary (5) scalar field name, kind I1P.
   contains
      ! public methods
      procedure, pass(self) :: initialize             !< Initialize class.
      procedure, pass(self) :: save_xh5f              !< Save in XH5F (XDMF/HDF5) format.
      generic               :: save_field =>           &
                               save_xh5f_field_4D_R8P, &
                               save_xh5f_field_4D_R4P, &
                               save_xh5f_field_4D_I8P, &
                               save_xh5f_field_4D_I4P, &
                               save_xh5f_field_4D_I2P, &
                               save_xh5f_field_4D_I1P, &
                               save_xh5f_field_5D_R8P, &
                               save_xh5f_field_5D_R4P, &
                               save_xh5f_field_5D_I8P, &
                               save_xh5f_field_5D_I4P, &
                               save_xh5f_field_5D_I2P, &
                               save_xh5f_field_5D_I1P !< Save fields by XH5F (XDMF/HDF5) file handler.
      ! private methods
      procedure, pass(self), private :: save_xh5f_field_4D_R8P !< Save fields by XH5F file handler, rank 4D, kind R8P.
      procedure, pass(self), private :: save_xh5f_field_4D_R4P !< Save fields by XH5F file handler, rank 4D, kind R4P.
      procedure, pass(self), private :: save_xh5f_field_4D_I8P !< Save fields by XH5F file handler, rank 4D, kind I8P.
      procedure, pass(self), private :: save_xh5f_field_4D_I4P !< Save fields by XH5F file handler, rank 4D, kind I4P.
      procedure, pass(self), private :: save_xh5f_field_4D_I2P !< Save fields by XH5F file handler, rank 4D, kind I2P.
      procedure, pass(self), private :: save_xh5f_field_4D_I1P !< Save fields by XH5F file handler, rank 4D, kind I1P.
      procedure, pass(self), private :: save_xh5f_field_5D_R8P !< Save fields by XH5F file handler, rank 5D, kind R8P.
      procedure, pass(self), private :: save_xh5f_field_5D_R4P !< Save fields by XH5F file handler, rank 5D, kind R4P.
      procedure, pass(self), private :: save_xh5f_field_5D_I8P !< Save fields by XH5F file handler, rank 5D, kind I8P.
      procedure, pass(self), private :: save_xh5f_field_5D_I4P !< Save fields by XH5F file handler, rank 5D, kind I4P.
      procedure, pass(self), private :: save_xh5f_field_5D_I2P !< Save fields by XH5F file handler, rank 5D, kind I2P.
      procedure, pass(self), private :: save_xh5f_field_5D_I1P !< Save fields by XH5F file handler, rank 5D, kind I1P.
endtype io_object

interface register_aux_field
   !< Register auxiliary fields data into ADAM IO class.
   module procedure register_aux_field_4D_R8P, &
                    register_aux_field_4D_R4P, &
                    register_aux_field_4D_I8P, &
                    register_aux_field_4D_I4P, &
                    register_aux_field_4D_I2P, &
                    register_aux_field_4D_I1P, &
                    register_aux_field_5D_R8P, &
                    register_aux_field_5D_R4P, &
                    register_aux_field_5D_I8P, &
                    register_aux_field_5D_I4P, &
                    register_aux_field_5D_I2P, &
                    register_aux_field_5D_I1P
endinterface register_aux_field
contains
   ! public methods
   subroutine initialize(self, grid, field,                                                                              &
                         q1_R8P,q1_R8P_name,q2_R8P,q2_R8P_name,q3_R8P,q3_R8P_name,q4_R8P,q4_R8P_name,q5_R8P,q5_R8P_name, &
                         q1_R4P,q1_R4P_name,q2_R4P,q2_R4P_name,q3_R4P,q3_R4P_name,q4_R4P,q4_R4P_name,q5_R4P,q5_R4P_name, &
                         q1_I8P,q1_I8P_name,q2_I8P,q2_I8P_name,q3_I8P,q3_I8P_name,q4_I8P,q4_I8P_name,q5_I8P,q5_I8P_name, &
                         q1_I4P,q1_I4P_name,q2_I4P,q2_I4P_name,q3_I4P,q3_I4P_name,q4_I4P,q4_I4P_name,q5_I4P,q5_I4P_name, &
                         q1_I2P,q1_I2P_name,q2_I2P,q2_I2P_name,q3_I2P,q3_I2P_name,q4_I2P,q4_I2P_name,q5_I2P,q5_I2P_name, &
                         q1_I1P,q1_I1P_name,q2_I1P,q2_I1P_name,q3_I1P,q3_I1P_name,q4_I1P,q4_I1P_name,q5_I1P,q5_I1P_name, &
                         s1_R8P,s1_R8P_name,s2_R8P,s2_R8P_name,s3_R8P,s3_R8P_name,s4_R8P,s4_R8P_name,s5_R8P,s5_R8P_name, &
                         s1_R4P,s1_R4P_name,s2_R4P,s2_R4P_name,s3_R4P,s3_R4P_name,s4_R4P,s4_R4P_name,s5_R4P,s5_R4P_name, &
                         s1_I8P,s1_I8P_name,s2_I8P,s2_I8P_name,s3_I8P,s3_I8P_name,s4_I8P,s4_I8P_name,s5_I8P,s5_I8P_name, &
                         s1_I4P,s1_I4P_name,s2_I4P,s2_I4P_name,s3_I4P,s3_I4P_name,s4_I4P,s4_I4P_name,s5_I4P,s5_I4P_name, &
                         s1_I2P,s1_I2P_name,s2_I2P,s2_I2P_name,s3_I2P,s3_I2P_name,s4_I2P,s4_I2P_name,s5_I2P,s5_I2P_name, &
                         s1_I1P,s1_I1P_name,s2_I1P,s2_I1P_name,s3_I1P,s3_I1P_name,s4_I1P,s4_I1P_name,s5_I1P,s5_I1P_name)
   !< Initialize class.
   class(io_object),   intent(inout)                 :: self              !< IO handler.
   type(grid_object),  intent(in), target            :: grid              !< The grid.
   type(field_object), intent(in), target            :: field             !< The field.
   real(R8P),          intent(in), pointer, optional :: q1_R8P(:,:,:,:,:) !< Auxiliary (1) vector cell centered fields, kind R8P.
   real(R8P),          intent(in), pointer, optional :: q2_R8P(:,:,:,:,:) !< Auxiliary (2) vector cell centered fields, kind R8P.
   real(R8P),          intent(in), pointer, optional :: q3_R8P(:,:,:,:,:) !< Auxiliary (3) vector cell centered fields, kind R8P.
   real(R8P),          intent(in), pointer, optional :: q4_R8P(:,:,:,:,:) !< Auxiliary (4) vector cell centered fields, kind R8P.
   real(R8P),          intent(in), pointer, optional :: q5_R8P(:,:,:,:,:) !< Auxiliary (5) vector cell centered fields, kind R8P.
   real(R4P),          intent(in), pointer, optional :: q1_R4P(:,:,:,:,:) !< Auxiliary (1) vector cell centered fields, kind R4P.
   real(R4P),          intent(in), pointer, optional :: q2_R4P(:,:,:,:,:) !< Auxiliary (2) vector cell centered fields, kind R4P.
   real(R4P),          intent(in), pointer, optional :: q3_R4P(:,:,:,:,:) !< Auxiliary (3) vector cell centered fields, kind R4P.
   real(R4P),          intent(in), pointer, optional :: q4_R4P(:,:,:,:,:) !< Auxiliary (4) vector cell centered fields, kind R4P.
   real(R4P),          intent(in), pointer, optional :: q5_R4P(:,:,:,:,:) !< Auxiliary (5) vector cell centered fields, kind R4P.
   integer(I8P),       intent(in), pointer, optional :: q1_I8P(:,:,:,:,:) !< Auxiliary (1) vector cell centered fields, kind I8P.
   integer(I8P),       intent(in), pointer, optional :: q2_I8P(:,:,:,:,:) !< Auxiliary (2) vector cell centered fields, kind I8P.
   integer(I8P),       intent(in), pointer, optional :: q3_I8P(:,:,:,:,:) !< Auxiliary (3) vector cell centered fields, kind I8P.
   integer(I8P),       intent(in), pointer, optional :: q4_I8P(:,:,:,:,:) !< Auxiliary (4) vector cell centered fields, kind I8P.
   integer(I8P),       intent(in), pointer, optional :: q5_I8P(:,:,:,:,:) !< Auxiliary (5) vector cell centered fields, kind I8P.
   integer(I4P),       intent(in), pointer, optional :: q1_I4P(:,:,:,:,:) !< Auxiliary (1) vector cell centered fields, kind I4P.
   integer(I4P),       intent(in), pointer, optional :: q2_I4P(:,:,:,:,:) !< Auxiliary (2) vector cell centered fields, kind I4P.
   integer(I4P),       intent(in), pointer, optional :: q3_I4P(:,:,:,:,:) !< Auxiliary (3) vector cell centered fields, kind I4P.
   integer(I4P),       intent(in), pointer, optional :: q4_I4P(:,:,:,:,:) !< Auxiliary (4) vector cell centered fields, kind I4P.
   integer(I4P),       intent(in), pointer, optional :: q5_I4P(:,:,:,:,:) !< Auxiliary (5) vector cell centered fields, kind I4P.
   integer(I2P),       intent(in), pointer, optional :: q1_I2P(:,:,:,:,:) !< Auxiliary (1) vector cell centered fields, kind I2P.
   integer(I2P),       intent(in), pointer, optional :: q2_I2P(:,:,:,:,:) !< Auxiliary (2) vector cell centered fields, kind I2P.
   integer(I2P),       intent(in), pointer, optional :: q3_I2P(:,:,:,:,:) !< Auxiliary (3) vector cell centered fields, kind I2P.
   integer(I2P),       intent(in), pointer, optional :: q4_I2P(:,:,:,:,:) !< Auxiliary (4) vector cell centered fields, kind I2P.
   integer(I2P),       intent(in), pointer, optional :: q5_I2P(:,:,:,:,:) !< Auxiliary (5) vector cell centered fields, kind I2P.
   integer(I1P),       intent(in), pointer, optional :: q1_I1P(:,:,:,:,:) !< Auxiliary (1) vector cell centered fields, kind I1P.
   integer(I1P),       intent(in), pointer, optional :: q2_I1P(:,:,:,:,:) !< Auxiliary (2) vector cell centered fields, kind I1P.
   integer(I1P),       intent(in), pointer, optional :: q3_I1P(:,:,:,:,:) !< Auxiliary (3) vector cell centered fields, kind I1P.
   integer(I1P),       intent(in), pointer, optional :: q4_I1P(:,:,:,:,:) !< Auxiliary (4) vector cell centered fields, kind I1P.
   integer(I1P),       intent(in), pointer, optional :: q5_I1P(:,:,:,:,:) !< Auxiliary (5) vector cell centered fields, kind I1P.
   real(R8P),          intent(in), pointer, optional :: s1_R8P(  :,:,:,:) !< Auxiliary (1) scalar cell centered fields, kind R8P.
   real(R8P),          intent(in), pointer, optional :: s2_R8P(  :,:,:,:) !< Auxiliary (2) scalar cell centered fields, kind R8P.
   real(R8P),          intent(in), pointer, optional :: s3_R8P(  :,:,:,:) !< Auxiliary (3) scalar cell centered fields, kind R8P.
   real(R8P),          intent(in), pointer, optional :: s4_R8P(  :,:,:,:) !< Auxiliary (4) scalar cell centered fields, kind R8P.
   real(R8P),          intent(in), pointer, optional :: s5_R8P(  :,:,:,:) !< Auxiliary (5) scalar cell centered fields, kind R8P.
   real(R4P),          intent(in), pointer, optional :: s1_R4P(  :,:,:,:) !< Auxiliary (1) scalar cell centered fields, kind R4P.
   real(R4P),          intent(in), pointer, optional :: s2_R4P(  :,:,:,:) !< Auxiliary (2) scalar cell centered fields, kind R4P.
   real(R4P),          intent(in), pointer, optional :: s3_R4P(  :,:,:,:) !< Auxiliary (3) scalar cell centered fields, kind R4P.
   real(R4P),          intent(in), pointer, optional :: s4_R4P(  :,:,:,:) !< Auxiliary (4) scalar cell centered fields, kind R4P.
   real(R4P),          intent(in), pointer, optional :: s5_R4P(  :,:,:,:) !< Auxiliary (5) scalar cell centered fields, kind R4P.
   integer(I8P),       intent(in), pointer, optional :: s1_I8P(  :,:,:,:) !< Auxiliary (1) scalar cell centered fields, kind I8P.
   integer(I8P),       intent(in), pointer, optional :: s2_I8P(  :,:,:,:) !< Auxiliary (2) scalar cell centered fields, kind I8P.
   integer(I8P),       intent(in), pointer, optional :: s3_I8P(  :,:,:,:) !< Auxiliary (3) scalar cell centered fields, kind I8P.
   integer(I8P),       intent(in), pointer, optional :: s4_I8P(  :,:,:,:) !< Auxiliary (4) scalar cell centered fields, kind I8P.
   integer(I8P),       intent(in), pointer, optional :: s5_I8P(  :,:,:,:) !< Auxiliary (5) scalar cell centered fields, kind I8P.
   integer(I4P),       intent(in), pointer, optional :: s1_I4P(  :,:,:,:) !< Auxiliary (1) scalar cell centered fields, kind I4P.
   integer(I4P),       intent(in), pointer, optional :: s2_I4P(  :,:,:,:) !< Auxiliary (2) scalar cell centered fields, kind I4P.
   integer(I4P),       intent(in), pointer, optional :: s3_I4P(  :,:,:,:) !< Auxiliary (3) scalar cell centered fields, kind I4P.
   integer(I4P),       intent(in), pointer, optional :: s4_I4P(  :,:,:,:) !< Auxiliary (4) scalar cell centered fields, kind I4P.
   integer(I4P),       intent(in), pointer, optional :: s5_I4P(  :,:,:,:) !< Auxiliary (5) scalar cell centered fields, kind I4P.
   integer(I2P),       intent(in), pointer, optional :: s1_I2P(  :,:,:,:) !< Auxiliary (1) scalar cell centered fields, kind I2P.
   integer(I2P),       intent(in), pointer, optional :: s2_I2P(  :,:,:,:) !< Auxiliary (2) scalar cell centered fields, kind I2P.
   integer(I2P),       intent(in), pointer, optional :: s3_I2P(  :,:,:,:) !< Auxiliary (3) scalar cell centered fields, kind I2P.
   integer(I2P),       intent(in), pointer, optional :: s4_I2P(  :,:,:,:) !< Auxiliary (4) scalar cell centered fields, kind I2P.
   integer(I2P),       intent(in), pointer, optional :: s5_I2P(  :,:,:,:) !< Auxiliary (5) scalar cell centered fields, kind I2P.
   integer(I1P),       intent(in), pointer, optional :: s1_I1P(  :,:,:,:) !< Auxiliary (1) scalar cell centered fields, kind I1P.
   integer(I1P),       intent(in), pointer, optional :: s2_I1P(  :,:,:,:) !< Auxiliary (2) scalar cell centered fields, kind I1P.
   integer(I1P),       intent(in), pointer, optional :: s3_I1P(  :,:,:,:) !< Auxiliary (3) scalar cell centered fields, kind I1P.
   integer(I1P),       intent(in), pointer, optional :: s4_I1P(  :,:,:,:) !< Auxiliary (4) scalar cell centered fields, kind I1P.
   integer(I1P),       intent(in), pointer, optional :: s5_I1P(  :,:,:,:) !< Auxiliary (5) scalar cell centered fields, kind I1P.
   character(*),       intent(in),          optional :: q1_R8P_name(:)    !< Auxiliary (1) vector fields names, kind R8P.
   character(*),       intent(in),          optional :: q2_R8P_name(:)    !< Auxiliary (2) vector fields names, kind R8P.
   character(*),       intent(in),          optional :: q3_R8P_name(:)    !< Auxiliary (3) vector fields names, kind R8P.
   character(*),       intent(in),          optional :: q4_R8P_name(:)    !< Auxiliary (4) vector fields names, kind R8P.
   character(*),       intent(in),          optional :: q5_R8P_name(:)    !< Auxiliary (5) vector fields names, kind R8P.
   character(*),       intent(in),          optional :: q1_R4P_name(:)    !< Auxiliary (1) vector fields names, kind R4P.
   character(*),       intent(in),          optional :: q2_R4P_name(:)    !< Auxiliary (2) vector fields names, kind R4P.
   character(*),       intent(in),          optional :: q3_R4P_name(:)    !< Auxiliary (3) vector fields names, kind R4P.
   character(*),       intent(in),          optional :: q4_R4P_name(:)    !< Auxiliary (4) vector fields names, kind R4P.
   character(*),       intent(in),          optional :: q5_R4P_name(:)    !< Auxiliary (5) vector fields names, kind R4P.
   character(*),       intent(in),          optional :: q1_I8P_name(:)    !< Auxiliary (1) vector fields names, kind I8P.
   character(*),       intent(in),          optional :: q2_I8P_name(:)    !< Auxiliary (2) vector fields names, kind I8P.
   character(*),       intent(in),          optional :: q3_I8P_name(:)    !< Auxiliary (3) vector fields names, kind I8P.
   character(*),       intent(in),          optional :: q4_I8P_name(:)    !< Auxiliary (4) vector fields names, kind I8P.
   character(*),       intent(in),          optional :: q5_I8P_name(:)    !< Auxiliary (5) vector fields names, kind I8P.
   character(*),       intent(in),          optional :: q1_I4P_name(:)    !< Auxiliary (1) vector fields names, kind I4P.
   character(*),       intent(in),          optional :: q2_I4P_name(:)    !< Auxiliary (2) vector fields names, kind I4P.
   character(*),       intent(in),          optional :: q3_I4P_name(:)    !< Auxiliary (3) vector fields names, kind I4P.
   character(*),       intent(in),          optional :: q4_I4P_name(:)    !< Auxiliary (4) vector fields names, kind I4P.
   character(*),       intent(in),          optional :: q5_I4P_name(:)    !< Auxiliary (5) vector fields names, kind I4P.
   character(*),       intent(in),          optional :: q1_I2P_name(:)    !< Auxiliary (1) vector fields names, kind I2P.
   character(*),       intent(in),          optional :: q2_I2P_name(:)    !< Auxiliary (2) vector fields names, kind I2P.
   character(*),       intent(in),          optional :: q3_I2P_name(:)    !< Auxiliary (3) vector fields names, kind I2P.
   character(*),       intent(in),          optional :: q4_I2P_name(:)    !< Auxiliary (4) vector fields names, kind I2P.
   character(*),       intent(in),          optional :: q5_I2P_name(:)    !< Auxiliary (5) vector fields names, kind I2P.
   character(*),       intent(in),          optional :: q1_I1P_name(:)    !< Auxiliary (1) vector fields names, kind I1P.
   character(*),       intent(in),          optional :: q2_I1P_name(:)    !< Auxiliary (2) vector fields names, kind I1P.
   character(*),       intent(in),          optional :: q3_I1P_name(:)    !< Auxiliary (3) vector fields names, kind I1P.
   character(*),       intent(in),          optional :: q4_I1P_name(:)    !< Auxiliary (4) vector fields names, kind I1P.
   character(*),       intent(in),          optional :: q5_I1P_name(:)    !< Auxiliary (5) vector fields names, kind I1P.
   character(*),       intent(in),          optional :: s1_R8P_name       !< Auxiliary (1) scalar field name, kind R8P.
   character(*),       intent(in),          optional :: s2_R8P_name       !< Auxiliary (2) scalar field name, kind R8P.
   character(*),       intent(in),          optional :: s3_R8P_name       !< Auxiliary (3) scalar field name, kind R8P.
   character(*),       intent(in),          optional :: s4_R8P_name       !< Auxiliary (4) scalar field name, kind R8P.
   character(*),       intent(in),          optional :: s5_R8P_name       !< Auxiliary (5) scalar field name, kind R8P.
   character(*),       intent(in),          optional :: s1_R4P_name       !< Auxiliary (1) scalar field name, kind R4P.
   character(*),       intent(in),          optional :: s2_R4P_name       !< Auxiliary (2) scalar field name, kind R4P.
   character(*),       intent(in),          optional :: s3_R4P_name       !< Auxiliary (3) scalar field name, kind R4P.
   character(*),       intent(in),          optional :: s4_R4P_name       !< Auxiliary (4) scalar field name, kind R4P.
   character(*),       intent(in),          optional :: s5_R4P_name       !< Auxiliary (5) scalar field name, kind R4P.
   character(*),       intent(in),          optional :: s1_I8P_name       !< Auxiliary (1) scalar field name, kind I8P.
   character(*),       intent(in),          optional :: s2_I8P_name       !< Auxiliary (2) scalar field name, kind I8P.
   character(*),       intent(in),          optional :: s3_I8P_name       !< Auxiliary (3) scalar field name, kind I8P.
   character(*),       intent(in),          optional :: s4_I8P_name       !< Auxiliary (4) scalar field name, kind I8P.
   character(*),       intent(in),          optional :: s5_I8P_name       !< Auxiliary (5) scalar field name, kind I8P.
   character(*),       intent(in),          optional :: s1_I4P_name       !< Auxiliary (1) scalar field name, kind I4P.
   character(*),       intent(in),          optional :: s2_I4P_name       !< Auxiliary (2) scalar field name, kind I4P.
   character(*),       intent(in),          optional :: s3_I4P_name       !< Auxiliary (3) scalar field name, kind I4P.
   character(*),       intent(in),          optional :: s4_I4P_name       !< Auxiliary (4) scalar field name, kind I4P.
   character(*),       intent(in),          optional :: s5_I4P_name       !< Auxiliary (5) scalar field name, kind I4P.
   character(*),       intent(in),          optional :: s1_I2P_name       !< Auxiliary (1) scalar field name, kind I2P.
   character(*),       intent(in),          optional :: s2_I2P_name       !< Auxiliary (2) scalar field name, kind I2P.
   character(*),       intent(in),          optional :: s3_I2P_name       !< Auxiliary (3) scalar field name, kind I2P.
   character(*),       intent(in),          optional :: s4_I2P_name       !< Auxiliary (4) scalar field name, kind I2P.
   character(*),       intent(in),          optional :: s5_I2P_name       !< Auxiliary (5) scalar field name, kind I2P.
   character(*),       intent(in),          optional :: s1_I1P_name       !< Auxiliary (1) scalar field name, kind I1P.
   character(*),       intent(in),          optional :: s2_I1P_name       !< Auxiliary (2) scalar field name, kind I1P.
   character(*),       intent(in),          optional :: s3_I1P_name       !< Auxiliary (3) scalar field name, kind I1P.
   character(*),       intent(in),          optional :: s4_I1P_name       !< Auxiliary (4) scalar field name, kind I1P.
   character(*),       intent(in),          optional :: s5_I1P_name       !< Auxiliary (5) scalar field name, kind I1P.

   call self%mpih%initialize
   call self%mpih%print_message('io_object%initialize start')
   self%grid  => grid
   self%field => field
   ! "register" user auxiliary fields data
   if (present(q1_R8P)) call register_aux_field(q_src=q1_R8P,q_name_src=q1_R8P_name,q_reg=self%q1_R8P,q_name_reg=self%q1_R8P_name)
   if (present(q2_R8P)) call register_aux_field(q_src=q2_R8P,q_name_src=q2_R8P_name,q_reg=self%q2_R8P,q_name_reg=self%q2_R8P_name)
   if (present(q3_R8P)) call register_aux_field(q_src=q3_R8P,q_name_src=q3_R8P_name,q_reg=self%q3_R8P,q_name_reg=self%q3_R8P_name)
   if (present(q4_R8P)) call register_aux_field(q_src=q4_R8P,q_name_src=q4_R8P_name,q_reg=self%q4_R8P,q_name_reg=self%q4_R8P_name)
   if (present(q5_R8P)) call register_aux_field(q_src=q5_R8P,q_name_src=q5_R8P_name,q_reg=self%q5_R8P,q_name_reg=self%q5_R8P_name)
   if (present(q1_R4P)) call register_aux_field(q_src=q1_R4P,q_name_src=q1_R4P_name,q_reg=self%q1_R4P,q_name_reg=self%q1_R4P_name)
   if (present(q2_R4P)) call register_aux_field(q_src=q2_R4P,q_name_src=q2_R4P_name,q_reg=self%q2_R4P,q_name_reg=self%q2_R4P_name)
   if (present(q3_R4P)) call register_aux_field(q_src=q3_R4P,q_name_src=q3_R4P_name,q_reg=self%q3_R4P,q_name_reg=self%q3_R4P_name)
   if (present(q4_R4P)) call register_aux_field(q_src=q4_R4P,q_name_src=q4_R4P_name,q_reg=self%q4_R4P,q_name_reg=self%q4_R4P_name)
   if (present(q5_R4P)) call register_aux_field(q_src=q5_R4P,q_name_src=q5_R4P_name,q_reg=self%q5_R4P,q_name_reg=self%q5_R4P_name)
   if (present(q1_I8P)) call register_aux_field(q_src=q1_I8P,q_name_src=q1_I8P_name,q_reg=self%q1_I8P,q_name_reg=self%q1_I8P_name)
   if (present(q2_I8P)) call register_aux_field(q_src=q2_I8P,q_name_src=q2_I8P_name,q_reg=self%q2_I8P,q_name_reg=self%q2_I8P_name)
   if (present(q3_I8P)) call register_aux_field(q_src=q3_I8P,q_name_src=q3_I8P_name,q_reg=self%q3_I8P,q_name_reg=self%q3_I8P_name)
   if (present(q4_I8P)) call register_aux_field(q_src=q4_I8P,q_name_src=q4_I8P_name,q_reg=self%q4_I8P,q_name_reg=self%q4_I8P_name)
   if (present(q5_I8P)) call register_aux_field(q_src=q5_I8P,q_name_src=q5_I8P_name,q_reg=self%q5_I8P,q_name_reg=self%q5_I8P_name)
   if (present(q1_I4P)) call register_aux_field(q_src=q1_I4P,q_name_src=q1_I4P_name,q_reg=self%q1_I4P,q_name_reg=self%q1_I4P_name)
   if (present(q2_I4P)) call register_aux_field(q_src=q2_I4P,q_name_src=q2_I4P_name,q_reg=self%q2_I4P,q_name_reg=self%q2_I4P_name)
   if (present(q3_I4P)) call register_aux_field(q_src=q3_I4P,q_name_src=q3_I4P_name,q_reg=self%q3_I4P,q_name_reg=self%q3_I4P_name)
   if (present(q4_I4P)) call register_aux_field(q_src=q4_I4P,q_name_src=q4_I4P_name,q_reg=self%q4_I4P,q_name_reg=self%q4_I4P_name)
   if (present(q5_I4P)) call register_aux_field(q_src=q5_I4P,q_name_src=q5_I4P_name,q_reg=self%q5_I4P,q_name_reg=self%q5_I4P_name)
   if (present(q1_I2P)) call register_aux_field(q_src=q1_I2P,q_name_src=q1_I2P_name,q_reg=self%q1_I2P,q_name_reg=self%q1_I2P_name)
   if (present(q2_I2P)) call register_aux_field(q_src=q2_I2P,q_name_src=q2_I2P_name,q_reg=self%q2_I2P,q_name_reg=self%q2_I2P_name)
   if (present(q3_I2P)) call register_aux_field(q_src=q3_I2P,q_name_src=q3_I2P_name,q_reg=self%q3_I2P,q_name_reg=self%q3_I2P_name)
   if (present(q4_I2P)) call register_aux_field(q_src=q4_I2P,q_name_src=q4_I2P_name,q_reg=self%q4_I2P,q_name_reg=self%q4_I2P_name)
   if (present(q5_I2P)) call register_aux_field(q_src=q5_I2P,q_name_src=q5_I2P_name,q_reg=self%q5_I2P,q_name_reg=self%q5_I2P_name)
   if (present(q1_I1P)) call register_aux_field(q_src=q1_I1P,q_name_src=q1_I1P_name,q_reg=self%q1_I1P,q_name_reg=self%q1_I1P_name)
   if (present(q2_I1P)) call register_aux_field(q_src=q2_I1P,q_name_src=q2_I1P_name,q_reg=self%q2_I1P,q_name_reg=self%q2_I1P_name)
   if (present(q3_I1P)) call register_aux_field(q_src=q3_I1P,q_name_src=q3_I1P_name,q_reg=self%q3_I1P,q_name_reg=self%q3_I1P_name)
   if (present(q4_I1P)) call register_aux_field(q_src=q4_I1P,q_name_src=q4_I1P_name,q_reg=self%q4_I1P,q_name_reg=self%q4_I1P_name)
   if (present(q5_I1P)) call register_aux_field(q_src=q5_I1P,q_name_src=q5_I1P_name,q_reg=self%q5_I1P,q_name_reg=self%q5_I1P_name)
   if (present(s1_R8P)) call register_aux_field(q_src=s1_R8P,q_name_src=s1_R8P_name,q_reg=self%s1_R8P,q_name_reg=self%s1_R8P_name)
   if (present(s2_R8P)) call register_aux_field(q_src=s2_R8P,q_name_src=s2_R8P_name,q_reg=self%s2_R8P,q_name_reg=self%s2_R8P_name)
   if (present(s3_R8P)) call register_aux_field(q_src=s3_R8P,q_name_src=s3_R8P_name,q_reg=self%s3_R8P,q_name_reg=self%s3_R8P_name)
   if (present(s4_R8P)) call register_aux_field(q_src=s4_R8P,q_name_src=s4_R8P_name,q_reg=self%s4_R8P,q_name_reg=self%s4_R8P_name)
   if (present(s5_R8P)) call register_aux_field(q_src=s5_R8P,q_name_src=s5_R8P_name,q_reg=self%s5_R8P,q_name_reg=self%s5_R8P_name)
   if (present(s1_R4P)) call register_aux_field(q_src=s1_R4P,q_name_src=s1_R4P_name,q_reg=self%s1_R4P,q_name_reg=self%s1_R4P_name)
   if (present(s2_R4P)) call register_aux_field(q_src=s2_R4P,q_name_src=s2_R4P_name,q_reg=self%s2_R4P,q_name_reg=self%s2_R4P_name)
   if (present(s3_R4P)) call register_aux_field(q_src=s3_R4P,q_name_src=s3_R4P_name,q_reg=self%s3_R4P,q_name_reg=self%s3_R4P_name)
   if (present(s4_R4P)) call register_aux_field(q_src=s4_R4P,q_name_src=s4_R4P_name,q_reg=self%s4_R4P,q_name_reg=self%s4_R4P_name)
   if (present(s5_R4P)) call register_aux_field(q_src=s5_R4P,q_name_src=s5_R4P_name,q_reg=self%s5_R4P,q_name_reg=self%s5_R4P_name)
   if (present(s1_I8P)) call register_aux_field(q_src=s1_I8P,q_name_src=s1_I8P_name,q_reg=self%s1_I8P,q_name_reg=self%s1_I8P_name)
   if (present(s2_I8P)) call register_aux_field(q_src=s2_I8P,q_name_src=s2_I8P_name,q_reg=self%s2_I8P,q_name_reg=self%s2_I8P_name)
   if (present(s3_I8P)) call register_aux_field(q_src=s3_I8P,q_name_src=s3_I8P_name,q_reg=self%s3_I8P,q_name_reg=self%s3_I8P_name)
   if (present(s4_I8P)) call register_aux_field(q_src=s4_I8P,q_name_src=s4_I8P_name,q_reg=self%s4_I8P,q_name_reg=self%s4_I8P_name)
   if (present(s5_I8P)) call register_aux_field(q_src=s5_I8P,q_name_src=s5_I8P_name,q_reg=self%s5_I8P,q_name_reg=self%s5_I8P_name)
   if (present(s1_I4P)) call register_aux_field(q_src=s1_I4P,q_name_src=s1_I4P_name,q_reg=self%s1_I4P,q_name_reg=self%s1_I4P_name)
   if (present(s2_I4P)) call register_aux_field(q_src=s2_I4P,q_name_src=s2_I4P_name,q_reg=self%s2_I4P,q_name_reg=self%s2_I4P_name)
   if (present(s3_I4P)) call register_aux_field(q_src=s3_I4P,q_name_src=s3_I4P_name,q_reg=self%s3_I4P,q_name_reg=self%s3_I4P_name)
   if (present(s4_I4P)) call register_aux_field(q_src=s4_I4P,q_name_src=s4_I4P_name,q_reg=self%s4_I4P,q_name_reg=self%s4_I4P_name)
   if (present(s5_I4P)) call register_aux_field(q_src=s5_I4P,q_name_src=s5_I4P_name,q_reg=self%s5_I4P,q_name_reg=self%s5_I4P_name)
   if (present(s1_I2P)) call register_aux_field(q_src=s1_I2P,q_name_src=s1_I2P_name,q_reg=self%s1_I2P,q_name_reg=self%s1_I2P_name)
   if (present(s2_I2P)) call register_aux_field(q_src=s2_I2P,q_name_src=s2_I2P_name,q_reg=self%s2_I2P,q_name_reg=self%s2_I2P_name)
   if (present(s3_I2P)) call register_aux_field(q_src=s3_I2P,q_name_src=s3_I2P_name,q_reg=self%s3_I2P,q_name_reg=self%s3_I2P_name)
   if (present(s4_I2P)) call register_aux_field(q_src=s4_I2P,q_name_src=s4_I2P_name,q_reg=self%s4_I2P,q_name_reg=self%s4_I2P_name)
   if (present(s5_I2P)) call register_aux_field(q_src=s5_I2P,q_name_src=s5_I2P_name,q_reg=self%s5_I2P,q_name_reg=self%s5_I2P_name)
   if (present(s1_I1P)) call register_aux_field(q_src=s1_I1P,q_name_src=s1_I1P_name,q_reg=self%s1_I1P,q_name_reg=self%s1_I1P_name)
   if (present(s2_I1P)) call register_aux_field(q_src=s2_I1P,q_name_src=s2_I1P_name,q_reg=self%s2_I1P,q_name_reg=self%s2_I1P_name)
   if (present(s3_I1P)) call register_aux_field(q_src=s3_I1P,q_name_src=s3_I1P_name,q_reg=self%s3_I1P,q_name_reg=self%s3_I1P_name)
   if (present(s4_I1P)) call register_aux_field(q_src=s4_I1P,q_name_src=s4_I1P_name,q_reg=self%s4_I1P,q_name_reg=self%s4_I1P_name)
   if (present(s5_I1P)) call register_aux_field(q_src=s5_I1P,q_name_src=s5_I1P_name,q_reg=self%s5_I1P,q_name_reg=self%s5_I1P_name)
   call self%mpih%print_message('io_object%initialize finish')
   endsubroutine initialize

   subroutine save_xh5f(self, basename, q, q_name, directory, with_ghost, with_cell_morton, t, time)
   !< Save ADAM in XH5F format.
   class(io_object), intent(inout)        :: self                    !< IO handler.
   character(*),     intent(in)           :: basename                !< Base name of output files.
   real(R8P),        intent(in)           :: q(1:,              &
                                               1-self%grid%ngc:,&
                                               1-self%grid%ngc:,&
                                               1-self%grid%ngc:,&
                                               1:)                   !< Q-vector variables [nv,ni,nj,nk,nb].
   character(*),     intent(in), optional :: q_name(:)               !< Q-vector variables names [nv].
   character(*),     intent(in), optional :: directory               !< Directory name of output files.
   logical,          intent(in), optional :: with_ghost              !< Flag to save ghost cells.
   logical,          intent(in), optional :: with_cell_morton        !< Flag to save Morton code also in cells.
   integer(I4P),     intent(in), optional :: t                       !< Time iteration.
   real(R8P),        intent(in), optional :: time                    !< Time.
   type(string), allocatable              :: q_name_(:)              !< Q-vector variables names [nv].
   character(:), allocatable              :: directory_              !< Directory name of output files, local var.
   logical                                :: with_ghost_             !< Flag to save ghost cells, local var.
   logical                                :: with_cell_morton_       !< Flag to save Morton code also in cells, local var.
   integer(I4P)                           :: t_                      !< Time iteration, local var.
   real(R8P)                              :: time_                   !< Time, local var.
   integer(I4P)                           :: ngc                     !< Ghost cells saved.
   integer(I4P)                           :: ijk(2,3)                !< Blocks extents.
   integer(I8P)                           :: nijk(3)                 !< Blocks dimensions.
   real(R8P)                              :: emin(3)                 !< Minimum abscissa of current block.
   character(:), allocatable              :: filename_hdf5           !< HDF5 file name.
   character(:), allocatable              :: filename_xdmf           !< XDMF file name.
   character(:), allocatable              :: bn                      !< Block name.
   type(xh5f_file_object)                 :: xh5f                    !< XH5F file handler.
   integer(I4P)                           :: i, b, v                 !< Counter.

   if (present(q_name)) then
      allocate(q_name_(size(q, dim=1)))
      do v=1, size(q, dim=1)
         q_name_(v) = trim(adjustl(q_name(v)))
      enddo
   else
      do v=1, size(q, dim=1)
         q_name_(v) = 'q-'//trim(strz(v,2))
      enddo
   endif
   directory_        = ''      ; if (present(directory       )) directory_        = trim(adjustl(directory))
   with_ghost_       = .false. ; if (present(with_ghost      )) with_ghost_       = with_ghost
   with_cell_morton_ = .false. ; if (present(with_cell_morton)) with_cell_morton_ = with_cell_morton
   t_                = 0_I4P   ; if (present(t               )) t_                = t
   time_             = 0._R8P  ; if (present(time            )) time_             = time
   if (with_ghost_) then
      ngc = self%grid%ngc
   else
      ngc = 0_I4P
   endif
   associate(ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk)
   ijk(:,1) = [1-ngc,ni+ngc]
   ijk(:,2) = [1-ngc,nj+ngc]
   ijk(:,3) = [1-ngc,nk+ngc]
   nijk = [ijk(2,1)-ijk(1,1)+1, &
           ijk(2,2)-ijk(1,2)+1, &
           ijk(2,3)-ijk(1,3)+1]
   endassociate
   filename_hdf5 = directory_//trim(adjustl(basename))//'-proc'//trim(strz(self%mpih%myrank,6))//'.h5'
   filename_xdmf = directory_//trim(adjustl(basename))//'.xdmf'
   call xh5f%open_file(filename_hdf5=filename_hdf5, filename_xdmf=filename_xdmf, act=FILE_PARAMETERS%FILE_ACTION_OVERWRITE)
   call xh5f%open_grid(grid_name='adam',                                 grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
   call xh5f%open_grid(grid_name='mpi_'//trim(strz(self%mpih%myrank,6)), grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION)
   associate(q1_R8P=>self%q1_R8P, q1_R8P_n=>self%q1_R8P_name, &
             q2_R8P=>self%q2_R8P, q2_R8P_n=>self%q2_R8P_name, &
             q3_R8P=>self%q3_R8P, q3_R8P_n=>self%q3_R8P_name, &
             q4_R8P=>self%q4_R8P, q4_R8P_n=>self%q4_R8P_name, &
             q5_R8P=>self%q5_R8P, q5_R8P_n=>self%q5_R8P_name, &
             q1_R4P=>self%q1_R4P, q1_R4P_n=>self%q1_R4P_name, &
             q2_R4P=>self%q2_R4P, q2_R4P_n=>self%q2_R4P_name, &
             q3_R4P=>self%q3_R4P, q3_R4P_n=>self%q3_R4P_name, &
             q4_R4P=>self%q4_R4P, q4_R4P_n=>self%q4_R4P_name, &
             q5_R4P=>self%q5_R4P, q5_R4P_n=>self%q5_R4P_name, &
             q1_I8P=>self%q1_I8P, q1_I8P_n=>self%q1_I8P_name, &
             q2_I8P=>self%q2_I8P, q2_I8P_n=>self%q2_I8P_name, &
             q3_I8P=>self%q3_I8P, q3_I8P_n=>self%q3_I8P_name, &
             q4_I8P=>self%q4_I8P, q4_I8P_n=>self%q4_I8P_name, &
             q5_I8P=>self%q5_I8P, q5_I8P_n=>self%q5_I8P_name, &
             q1_I4P=>self%q1_I4P, q1_I4P_n=>self%q1_I4P_name, &
             q2_I4P=>self%q2_I4P, q2_I4P_n=>self%q2_I4P_name, &
             q3_I4P=>self%q3_I4P, q3_I4P_n=>self%q3_I4P_name, &
             q4_I4P=>self%q4_I4P, q4_I4P_n=>self%q4_I4P_name, &
             q5_I4P=>self%q5_I4P, q5_I4P_n=>self%q5_I4P_name, &
             q1_I2P=>self%q1_I2P, q1_I2P_n=>self%q1_I2P_name, &
             q2_I2P=>self%q2_I2P, q2_I2P_n=>self%q2_I2P_name, &
             q3_I2P=>self%q3_I2P, q3_I2P_n=>self%q3_I2P_name, &
             q4_I2P=>self%q4_I2P, q4_I2P_n=>self%q4_I2P_name, &
             q5_I2P=>self%q5_I2P, q5_I2P_n=>self%q5_I2P_name, &
             q1_I1P=>self%q1_I1P, q1_I1P_n=>self%q1_I1P_name, &
             q2_I1P=>self%q2_I1P, q2_I1P_n=>self%q2_I1P_name, &
             q3_I1P=>self%q3_I1P, q3_I1P_n=>self%q3_I1P_name, &
             q4_I1P=>self%q4_I1P, q4_I1P_n=>self%q4_I1P_name, &
             q5_I1P=>self%q5_I1P, q5_I1P_n=>self%q5_I1P_name, &
             s1_R8P=>self%s1_R8P, s1_R8P_n=>self%s1_R8P_name, &
             s2_R8P=>self%s2_R8P, s2_R8P_n=>self%s2_R8P_name, &
             s3_R8P=>self%s3_R8P, s3_R8P_n=>self%s3_R8P_name, &
             s4_R8P=>self%s4_R8P, s4_R8P_n=>self%s4_R8P_name, &
             s5_R8P=>self%s5_R8P, s5_R8P_n=>self%s5_R8P_name, &
             s1_R4P=>self%s1_R4P, s1_R4P_n=>self%s1_R4P_name, &
             s2_R4P=>self%s2_R4P, s2_R4P_n=>self%s2_R4P_name, &
             s3_R4P=>self%s3_R4P, s3_R4P_n=>self%s3_R4P_name, &
             s4_R4P=>self%s4_R4P, s4_R4P_n=>self%s4_R4P_name, &
             s5_R4P=>self%s5_R4P, s5_R4P_n=>self%s5_R4P_name, &
             s1_I8P=>self%s1_I8P, s1_I8P_n=>self%s1_I8P_name, &
             s2_I8P=>self%s2_I8P, s2_I8P_n=>self%s2_I8P_name, &
             s3_I8P=>self%s3_I8P, s3_I8P_n=>self%s3_I8P_name, &
             s4_I8P=>self%s4_I8P, s4_I8P_n=>self%s4_I8P_name, &
             s5_I8P=>self%s5_I8P, s5_I8P_n=>self%s5_I8P_name, &
             s1_I4P=>self%s1_I4P, s1_I4P_n=>self%s1_I4P_name, &
             s2_I4P=>self%s2_I4P, s2_I4P_n=>self%s2_I4P_name, &
             s3_I4P=>self%s3_I4P, s3_I4P_n=>self%s3_I4P_name, &
             s4_I4P=>self%s4_I4P, s4_I4P_n=>self%s4_I4P_name, &
             s5_I4P=>self%s5_I4P, s5_I4P_n=>self%s5_I4P_name, &
             s1_I2P=>self%s1_I2P, s1_I2P_n=>self%s1_I2P_name, &
             s2_I2P=>self%s2_I2P, s2_I2P_n=>self%s2_I2P_name, &
             s3_I2P=>self%s3_I2P, s3_I2P_n=>self%s3_I2P_name, &
             s4_I2P=>self%s4_I2P, s4_I2P_n=>self%s4_I2P_name, &
             s5_I2P=>self%s5_I2P, s5_I2P_n=>self%s5_I2P_name, &
             s1_I1P=>self%s1_I1P, s1_I1P_n=>self%s1_I1P_name, &
             s2_I1P=>self%s2_I1P, s2_I1P_n=>self%s2_I1P_name, &
             s3_I1P=>self%s3_I1P, s3_I1P_n=>self%s3_I1P_name, &
             s4_I1P=>self%s4_I1P, s4_I1P_n=>self%s4_I1P_name, &
             s5_I1P=>self%s5_I1P, s5_I1P_n=>self%s5_I1P_name)
   do b=1, self%field%blocks_number
      emin = [self%field%emin(1,b)-ngc*self%field%dxyz(1,b), &
              self%field%emin(2,b)-ngc*self%field%dxyz(2,b), &
              self%field%emin(3,b)-ngc*self%field%dxyz(3,b)]
      bn = 'block_'//trim(strz(b,9))//'-proc'//trim(strz(self%mpih%myrank,6))
      call xh5f%open_block(block_type = XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN_UNIFORM, &
                           block_name = bn,                                           &
                           nijk       = nijk,                                         &
                           emin       = emin,                                         &
                           dxyz       = self%field%dxyz(:,b),                         &
                           time       = time_)
      call xh5f%save_block_field(xdmf_field_name = 'time_iteration',                                &
                                 field           = t_,                                              &
                                 field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_GRID,           &
                                 field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF, &
                                 hdf5_field_name = bn//'-time_iteration')
      ! if (with_cell_morton_) &
      !    call self%save_field(xh5f=xh5f, block_name=bn, ijk=ijk, nijk=nijk,                                            &
      !                         q=reshape([(self%field%code(b),i=1,nijk(1)*nijk(2)*nijk(3))],[nijk(1),nijk(2),nijk(3)]), &
      !                         q_name=string('morton'))
      call self%save_field(xh5f=xh5f, block_name=bn, ijk=ijk, nijk=nijk, q=q(:,:,:,:,b), q_name=q_name_)
      ! save registered auxiliary fields
      if(associated(self%q1_R8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q1_R8P(:,:,:,:,b),q_name=q1_R8P_n)
      if(associated(self%q2_R8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q2_R8P(:,:,:,:,b),q_name=q2_R8P_n)
      if(associated(self%q3_R8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q3_R8P(:,:,:,:,b),q_name=q3_R8P_n)
      if(associated(self%q4_R8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q4_R8P(:,:,:,:,b),q_name=q4_R8P_n)
      if(associated(self%q5_R8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q5_R8P(:,:,:,:,b),q_name=q5_R8P_n)
      if(associated(self%q1_R4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q1_R4P(:,:,:,:,b),q_name=q1_R4P_n)
      if(associated(self%q2_R4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q2_R4P(:,:,:,:,b),q_name=q2_R4P_n)
      if(associated(self%q3_R4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q3_R4P(:,:,:,:,b),q_name=q3_R4P_n)
      if(associated(self%q4_R4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q4_R4P(:,:,:,:,b),q_name=q4_R4P_n)
      if(associated(self%q5_R4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q5_R4P(:,:,:,:,b),q_name=q5_R4P_n)
      if(associated(self%q1_I8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q1_I8P(:,:,:,:,b),q_name=q1_I8P_n)
      if(associated(self%q2_I8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q2_I8P(:,:,:,:,b),q_name=q2_I8P_n)
      if(associated(self%q3_I8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q3_I8P(:,:,:,:,b),q_name=q3_I8P_n)
      if(associated(self%q4_I8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q4_I8P(:,:,:,:,b),q_name=q4_I8P_n)
      if(associated(self%q5_I8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q5_I8P(:,:,:,:,b),q_name=q5_I8P_n)
      if(associated(self%q1_I4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q1_I4P(:,:,:,:,b),q_name=q1_I4P_n)
      if(associated(self%q2_I4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q2_I4P(:,:,:,:,b),q_name=q2_I4P_n)
      if(associated(self%q3_I4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q3_I4P(:,:,:,:,b),q_name=q3_I4P_n)
      if(associated(self%q4_I4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q4_I4P(:,:,:,:,b),q_name=q4_I4P_n)
      if(associated(self%q5_I4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q5_I4P(:,:,:,:,b),q_name=q5_I4P_n)
      if(associated(self%q1_I2P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q1_I2P(:,:,:,:,b),q_name=q1_I2P_n)
      if(associated(self%q2_I2P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q2_I2P(:,:,:,:,b),q_name=q2_I2P_n)
      if(associated(self%q3_I2P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q3_I2P(:,:,:,:,b),q_name=q3_I2P_n)
      if(associated(self%q4_I2P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q4_I2P(:,:,:,:,b),q_name=q4_I2P_n)
      if(associated(self%q5_I2P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q5_I2P(:,:,:,:,b),q_name=q5_I2P_n)
      if(associated(self%q1_I1P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q1_I1P(:,:,:,:,b),q_name=q1_I1P_n)
      if(associated(self%q2_I1P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q2_I1P(:,:,:,:,b),q_name=q2_I1P_n)
      if(associated(self%q3_I1P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q3_I1P(:,:,:,:,b),q_name=q3_I1P_n)
      if(associated(self%q4_I1P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q4_I1P(:,:,:,:,b),q_name=q4_I1P_n)
      if(associated(self%q5_I1P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=q5_I1P(:,:,:,:,b),q_name=q5_I1P_n)
      if(associated(self%s1_R8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s1_R8P(  :,:,:,b),q_name=s1_R8P_n)
      if(associated(self%s2_R8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s2_R8P(  :,:,:,b),q_name=s2_R8P_n)
      if(associated(self%s3_R8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s3_R8P(  :,:,:,b),q_name=s3_R8P_n)
      if(associated(self%s4_R8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s4_R8P(  :,:,:,b),q_name=s4_R8P_n)
      if(associated(self%s5_R8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s5_R8P(  :,:,:,b),q_name=s5_R8P_n)
      if(associated(self%s1_R4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s1_R4P(  :,:,:,b),q_name=s1_R4P_n)
      if(associated(self%s2_R4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s2_R4P(  :,:,:,b),q_name=s2_R4P_n)
      if(associated(self%s3_R4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s3_R4P(  :,:,:,b),q_name=s3_R4P_n)
      if(associated(self%s4_R4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s4_R4P(  :,:,:,b),q_name=s4_R4P_n)
      if(associated(self%s5_R4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s5_R4P(  :,:,:,b),q_name=s5_R4P_n)
      if(associated(self%s1_I8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s1_I8P(  :,:,:,b),q_name=s1_I8P_n)
      if(associated(self%s2_I8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s2_I8P(  :,:,:,b),q_name=s2_I8P_n)
      if(associated(self%s3_I8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s3_I8P(  :,:,:,b),q_name=s3_I8P_n)
      if(associated(self%s4_I8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s4_I8P(  :,:,:,b),q_name=s4_I8P_n)
      if(associated(self%s5_I8P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s5_I8P(  :,:,:,b),q_name=s5_I8P_n)
      if(associated(self%s1_I4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s1_I4P(  :,:,:,b),q_name=s1_I4P_n)
      if(associated(self%s2_I4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s2_I4P(  :,:,:,b),q_name=s2_I4P_n)
      if(associated(self%s3_I4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s3_I4P(  :,:,:,b),q_name=s3_I4P_n)
      if(associated(self%s4_I4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s4_I4P(  :,:,:,b),q_name=s4_I4P_n)
      if(associated(self%s5_I4P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s5_I4P(  :,:,:,b),q_name=s5_I4P_n)
      if(associated(self%s1_I2P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s1_I2P(  :,:,:,b),q_name=s1_I2P_n)
      if(associated(self%s2_I2P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s2_I2P(  :,:,:,b),q_name=s2_I2P_n)
      if(associated(self%s3_I2P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s3_I2P(  :,:,:,b),q_name=s3_I2P_n)
      if(associated(self%s4_I2P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s4_I2P(  :,:,:,b),q_name=s4_I2P_n)
      if(associated(self%s5_I2P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s5_I2P(  :,:,:,b),q_name=s5_I2P_n)
      if(associated(self%s1_I1P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s1_I1P(  :,:,:,b),q_name=s1_I1P_n)
      if(associated(self%s2_I1P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s2_I1P(  :,:,:,b),q_name=s2_I1P_n)
      if(associated(self%s3_I1P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s3_I1P(  :,:,:,b),q_name=s3_I1P_n)
      if(associated(self%s4_I1P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s4_I1P(  :,:,:,b),q_name=s4_I1P_n)
      if(associated(self%s5_I1P))call self%save_field(xh5f=xh5f,block_name=bn,ijk=ijk,nijk=nijk,q=s5_I1P(  :,:,:,b),q_name=s5_I1P_n)
      call xh5f%close_block
   enddo
   endassociate
   call xh5f%close_grid
   call xh5f%close_grid(grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
   call xh5f%close_file
   endsubroutine save_xh5f

   ! private methods
#define KKP R8P
#define VARTYPE real
#define REGISTER_AUX_FIELD_4D_KKP register_aux_field_4D_R8P
#define REGISTER_AUX_FIELD_5D_KKP register_aux_field_5D_R8P
#define SAVE_XH5F_FIELD_4D_KKP save_xh5f_field_4D_R8P
#define SAVE_XH5F_FIELD_5D_KKP save_xh5f_field_5D_R8P
#include "adam_io_agnostic.INC"

#define KKP R4P
#define VARTYPE real
#define REGISTER_AUX_FIELD_4D_KKP register_aux_field_4D_R4P
#define REGISTER_AUX_FIELD_5D_KKP register_aux_field_5D_R4P
#define SAVE_XH5F_FIELD_4D_KKP save_xh5f_field_4D_R4P
#define SAVE_XH5F_FIELD_5D_KKP save_xh5f_field_5D_R4P
#include "adam_io_agnostic.INC"

#define KKP I8P
#define VARTYPE integer
#define REGISTER_AUX_FIELD_4D_KKP register_aux_field_4D_I8P
#define REGISTER_AUX_FIELD_5D_KKP register_aux_field_5D_I8P
#define SAVE_XH5F_FIELD_4D_KKP save_xh5f_field_4D_I8P
#define SAVE_XH5F_FIELD_5D_KKP save_xh5f_field_5D_I8P
#include "adam_io_agnostic.INC"

#define KKP I4P
#define VARTYPE integer
#define REGISTER_AUX_FIELD_4D_KKP register_aux_field_4D_I4P
#define REGISTER_AUX_FIELD_5D_KKP register_aux_field_5D_I4P
#define SAVE_XH5F_FIELD_4D_KKP save_xh5f_field_4D_I4P
#define SAVE_XH5F_FIELD_5D_KKP save_xh5f_field_5D_I4P
#include "adam_io_agnostic.INC"

#define KKP I2P
#define VARTYPE integer
#define REGISTER_AUX_FIELD_4D_KKP register_aux_field_4D_I2P
#define REGISTER_AUX_FIELD_5D_KKP register_aux_field_5D_I2P
#define SAVE_XH5F_FIELD_4D_KKP save_xh5f_field_4D_I2P
#define SAVE_XH5F_FIELD_5D_KKP save_xh5f_field_5D_I2P
#include "adam_io_agnostic.INC"

#define KKP I1P
#define VARTYPE integer
#define REGISTER_AUX_FIELD_4D_KKP register_aux_field_4D_I1P
#define REGISTER_AUX_FIELD_5D_KKP register_aux_field_5D_I1P
#define SAVE_XH5F_FIELD_4D_KKP save_xh5f_field_4D_I1P
#define SAVE_XH5F_FIELD_5D_KKP save_xh5f_field_5D_I1P
#include "adam_io_agnostic.INC"
endmodule adam_io_object
