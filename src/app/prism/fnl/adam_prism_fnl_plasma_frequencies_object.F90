!< PRISM FNL plasma frequencies: the device maxima behind the PIC time-step limit.
!<
!< Device twin of prism_plasma_frequencies_object%compute_cpu (adam_prism_plasma_frequencies_object.F90, which holds the
!< definitions, the composite frequencies and the dt criterion): max over the particles of Omega_c = |q_p| |B_p| / m_p
!< (electrons q < 0, ions q > 0) and of max_d |v_d|/dx_d, and max over the interior cells of omega_p^2 deposited NGP in the
!< cell of the neighbour list (no ghosts, no fold, no communication). Work array: wp2_gpu [nb,ni,nj,nk,2].
#include "fundal.H"

module adam_prism_fnl_plasma_frequencies_object

use :: adam_common_library
use :: adam_fnl_library
use :: adam_prism_common_library
use :: adam_prism_fnl_pic_object
use :: fundal
use :: mpi
use :: penf

implicit none
private
public :: prism_fnl_plasma_frequencies_object

type :: prism_fnl_plasma_frequencies_object
   real(R8P), pointer :: wp2_gpu(:,:,:,:,:) => null() !< omega_p^2 per interior cell [nb,ni,nj,nk,2] (1 electrons, 2 ions).
   contains
      procedure, pass(self) :: initialize
      procedure, pass(self) :: compute
      procedure, pass(self) :: destroy
endtype prism_fnl_plasma_frequencies_object

contains
   subroutine initialize(self, nb, ni, nj, nk)
   !< Allocate the work array on the block slots.
   class(prism_fnl_plasma_frequencies_object), intent(inout) :: self       !< Device helper.
   integer(I4P),                               intent(in)    :: nb         !< Block slots.
   integer(I4P),                               intent(in)    :: ni, nj, nk !< Block cells.
   integer(I4P)                                              :: ierr       !< Error status.

   call self%destroy
   call dev_alloc(fptr_dev=self%wp2_gpu, ubounds=[nb,ni,nj,nk,2], lbounds=[1,1,1,1,1], ierr=ierr)
   if (ierr /= 0_I4P) call mpih_fnl%error_stop(msg=': failed to allocate wp2_gpu in prism_fnl_plasma_frequencies_object')
   endsubroutine initialize

   subroutine destroy(self)
   !< Free the work array.
   class(prism_fnl_plasma_frequencies_object), intent(inout) :: self !< Device helper.

   if (associated(self%wp2_gpu)) then
      call dev_free(self%wp2_gpu, mydev)
      nullify(self%wp2_gpu)
   endif
   endsubroutine destroy

   subroutine compute(self, pic_fnl, field_fnl, grid, inv_eps, wpe2, wpi2, wce, wci, transit)
   !< Local maxima of this rank (the caller reduces them, prism_plasma_frequencies_object%set_maxima).
   class(prism_fnl_plasma_frequencies_object), intent(inout) :: self      !< Device helper.
   type(prism_fnl_pic_object),                 intent(in)    :: pic_fnl   !< Particles on device.
   type(field_fnl_object),                     intent(in)    :: field_fnl !< Device field helper.
   type(grid_object),                          intent(in)    :: grid      !< Grid.
   real(R8P),                                  intent(in)    :: inv_eps   !< 1/eps0 (1 for the non-dimensional models).
   real(R8P),                                  intent(out)   :: wpe2      !< Max omega_pe^2.
   real(R8P),                                  intent(out)   :: wpi2      !< Max omega_pi^2.
   real(R8P),                                  intent(out)   :: wce       !< Max Omega_ce.
   real(R8P),                                  intent(out)   :: wci       !< Max Omega_ci.
   real(R8P),                                  intent(out)   :: transit   !< Max |v_d|/dx_d.
   real(R8P), pointer                                        :: w_gpu(:,:,:,:,:), q_pic_gpu(:,:), pic_fields_gpu(:,:)
   real(R8P), pointer                                        :: dxyz_gpu(:,:)
   integer(I4P), pointer                                     :: neighbour_list_gpu(:,:)
   integer(I4P)                                              :: n, b, i, j, k, s, nb_, ni, nj, nk, np
   real(R8P)                                                 :: qp, w, val, tr, mce, mci, mpe, mpi

   wpe2 = 0._R8P ; wpi2 = 0._R8P ; wce = 0._R8P ; wci = 0._R8P ; transit = 0._R8P
   if (.not.associated(pic_fnl%q_pic_gpu)) return
   np = pic_fnl%particle_number
   if (np <= 0_I4P) return
   w_gpu => self%wp2_gpu
   nb_ = field_fnl%blocks_number ; ni = grid%ni ; nj = grid%nj ; nk = grid%nk
   q_pic_gpu => pic_fnl%q_pic_gpu ; pic_fields_gpu => pic_fnl%pic_fields_gpu
   neighbour_list_gpu => pic_fnl%neighbour_list_gpu ; dxyz_gpu => field_fnl%dxyz_gpu

   !$acc parallel loop collapse(5) independent DEVICEVAR(w_gpu)
   !$omp OMPLOOP collapse(5) DEVICEPTR(w_gpu)
   do s = 1, 2
      do k = 1, nk
         do j = 1, nj
            do i = 1, ni
               do b = 1, nb_
                  w_gpu(b,i,j,k,s) = 0._R8P
               enddo
            enddo
         enddo
      enddo
   enddo

   tr = 0._R8P ; mce = 0._R8P ; mci = 0._R8P
   !$acc parallel loop independent DEVICEVAR(w_gpu, q_pic_gpu, pic_fields_gpu, neighbour_list_gpu, dxyz_gpu) &
   !$acc& private(qp, w, val, b, i, j, k, s) reduction(max: mce, mci, tr)
   !$omp OMPLOOP DEVICEPTR(w_gpu, q_pic_gpu, pic_fields_gpu, neighbour_list_gpu, dxyz_gpu) &
   !$omp& private(qp, w, val, b, i, j, k, s) reduction(max: mce, mci, tr)
   do n = 1, np
      qp = q_pic_gpu(n,7)
      if (qp == 0._R8P) cycle
      w = abs(qp)*sqrt(pic_fields_gpu(n,4)**2 + pic_fields_gpu(n,5)**2 + pic_fields_gpu(n,6)**2)/q_pic_gpu(n,8)
      if (qp < 0._R8P) then
         mce = max(mce, w)
      else
         mci = max(mci, w)
      endif
      b = neighbour_list_gpu(n,1)
      if (b <= 0_I4P) cycle
      tr = max(tr, abs(q_pic_gpu(n,4))/dxyz_gpu(b,1), abs(q_pic_gpu(n,5))/dxyz_gpu(b,2), abs(q_pic_gpu(n,6))/dxyz_gpu(b,3))
      i = neighbour_list_gpu(n,2) ; j = neighbour_list_gpu(n,3) ; k = neighbour_list_gpu(n,4)
      ! a particle in the ghost halo (wrapped across a periodic face at the end of the step, one rank) is skipped
      if (i < 1 .or. i > ni .or. j < 1 .or. j > nj .or. k < 1 .or. k > nk) cycle
      s = 1_I4P ; if (qp > 0._R8P) s = 2_I4P
      val = qp*qp*inv_eps/(q_pic_gpu(n,8)*dxyz_gpu(b,1)*dxyz_gpu(b,2)*dxyz_gpu(b,3))
      !$acc atomic update
      !$omp atomic update
      w_gpu(b,i,j,k,s) = w_gpu(b,i,j,k,s) + val
   enddo

   mpe = 0._R8P ; mpi = 0._R8P
   !$acc parallel loop collapse(4) independent DEVICEVAR(w_gpu) reduction(max: mpe, mpi)
   !$omp OMPLOOP collapse(4) DEVICEPTR(w_gpu) reduction(max: mpe, mpi)
   do k = 1, nk
      do j = 1, nj
         do i = 1, ni
            do b = 1, nb_
               mpe = max(mpe, w_gpu(b,i,j,k,1))
               mpi = max(mpi, w_gpu(b,i,j,k,2))
            enddo
         enddo
      enddo
   enddo
   wpe2 = mpe ; wpi2 = mpi ; wce = mce ; wci = mci ; transit = tr
   endsubroutine compute
endmodule adam_prism_fnl_plasma_frequencies_object
