!< ADAM, IB class FNL (FNL backend of [[ib_object]]).
module adam_fnl_ib_object
!< ADAM, IB class FNL (FNL backend of [[ib_object]]).

! ADAM classes, libraries, parameters
use :: adam_ib_object
! ADAM singleton objects
use :: adam_ib_global,    only : ib
use :: adam_field_global, only : field
use :: adam_grid_global,  only : grid
! ADAM FNL classes, libraries, parameters
use :: adam_fnl_ib_kernels
! ADAM FNL singleton objects
use :: adam_fnl_mpih_global, only : mpih_fnl
! third party modules
use :: fundal
use :: penf

implicit none
private
public :: ib_fnl_object

type :: ib_fnl_object
   !< IB FNL class definition.
   ! device data
   ! `=> null()` is mandatory: dev_assign_to_device / dev_alloc test
   ! `associated(dst)` before allocating, and associated() on a pointer
   ! that was never nullified is undefined behaviour (nvfortran -Mchkptr
   ! traps it as "Null pointer").
   real(R8P), pointer :: q_bcs_vars_gpu(:,:) => null() !< Variables array for immersed boundary on GPU.
   real(R8P), pointer :: phi_gpu(:,:,:,:,:)  => null() !< Distance function on GPU.
   contains
      ! public methods
      procedure, pass(self) :: evolve_eikonal !< Evolve eikonal equation.
      procedure, pass(self) :: initialize     !< Initialize class from global singletons.
      procedure, pass(self) :: invert_eikonal !< Invert momentum eikonal equation.
endtype ib_fnl_object
contains
   ! public methods
   subroutine evolve_eikonal(self, dq_gpu, q_gpu, dxyz_gpu)
   !< Evolve eikonal equation.
   class(ib_fnl_object), intent(in)    :: self       !< IB.
   real(R8P),            intent(inout) :: dq_gpu(1:,         &
                                                  1-grid%ngc:,&
                                                  1-grid%ngc:,&
                                                  1-grid%ngc:,&
                                                  1:) !< State variables variations.
   real(R8P),            intent(inout) :: q_gpu(1:,         &
                                                 1-grid%ngc:,&
                                                 1-grid%ngc:,&
                                                 1-grid%ngc:,&
                                                 1:)  !< Conservative variables.
   real(R8P),            intent(in)    :: dxyz_gpu(:,:) !< Cell deltas on GPU [nb,3].
   integer(I4P)                        :: iib        !< Counter.

   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, ngc=>grid%ngc, nv=>field%nv, &
             dx_gpu=>dxyz_gpu(:,1), dy_gpu=>dxyz_gpu(:,2), dz_gpu=>dxyz_gpu(:,3),                                    &
             solids_number=>ib%solids_number, phi_gpu=>self%phi_gpu)
      do iib=1, solids_number
         call compute_eikonal_dq_phi_dev(ib=iib, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                         dx_gpu=dx_gpu, dy_gpu=dy_gpu, dz_gpu=dz_gpu,                              &
                                         phi_gpu=phi_gpu, dq_gpu=dq_gpu, q_gpu=q_gpu)
         call evolve_eikonal_q_phi_dev(ib=iib, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                       phi_gpu=phi_gpu, dq_gpu=dq_gpu, q_gpu=q_gpu)
      enddo
   endassociate
   endsubroutine evolve_eikonal

   subroutine initialize(self)
   !< Initialize class from program-scope `ib` (adam_ib_global), `field` (adam_field_global)
   !< and `grid` (adam_grid_global) singletons.
   !< Requires `mpih_fnl` (adam_fnl_mpih_global) and the field/grid/ib singletons to be ready.
   class(ib_fnl_object), intent(inout) :: self !< IB FNL object.
   integer(I4P)                        :: ierr !< Error status.

   call mpih_fnl%print_message('ib_fnl_object%initialize start')
   associate(ngc=>grid%ngc, ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, nb=>field%nb, solids_number=>ib%solids_number)
   ! ib%q is allocated only when solids_number > 0; with no immersed solids
   ! it is unallocated, and passing an unallocated allocatable as the
   ! assumed-shape intent(in) src of dev_assign_to_device is illegal Fortran.
   ! Gate on solids_number, matching the phi_gpu allocation below; with no
   ! solids the GPU pointer stays null, the correct "nothing to offload" state.
   if (solids_number>0) then
      call dev_assign_to_device(dst=self%q_bcs_vars_gpu, src=ib%q)
      call dev_alloc(fptr_dev=self%phi_gpu,                              &
                     ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,solids_number+1], &
                     lbounds=[1, 1-ngc, 1-ngc, 1-ngc, 1              ], &
                     ierr=ierr, init_value=-1._R8P)
   endif
   endassociate
   endsubroutine initialize

   subroutine invert_eikonal(self, q_gpu)
   !< Invert momentum eikonal equation.
   class(ib_fnl_object), intent(in)    :: self       !< IB.
   real(R8P),            intent(inout) :: q_gpu(1:,         &
                                                 1-grid%ngc:,&
                                                 1-grid%ngc:,&
                                                 1-grid%ngc:,&
                                                 1:) !< Conservative variables.
   integer(I4P)                        :: iib        !< Counter.

   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, ngc=>grid%ngc, nv=>field%nv, &
             blocks_number=>field%blocks_number,                                   &
             bcs_type=>ib%bc_type, solids_number=>ib%solids_number, phi_gpu=>self%phi_gpu)
      do iib=1, solids_number
         call invert_eikonal_q_phi_dev(BCS_VISCOUS=BCS_VISCOUS, BCS_EULER=BCS_EULER,               &
                                       ib=iib, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv,               &
                                       blocks_number=blocks_number, bcs_type=bcs_type(iib),         &
                                       phi_gpu=phi_gpu, q_gpu=q_gpu)
      enddo
   endassociate
   endsubroutine invert_eikonal
endmodule adam_fnl_ib_object
