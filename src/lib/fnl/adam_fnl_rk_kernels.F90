!< ADAM, RK FNL kernels (FNL backend of [[rk_fnl_object]]).

#include "fundal.H"

module adam_fnl_rk_kernels
!< ADAM, RK FNL kernels (FNL backend of [[rk_fnl_object]]).

! third party modules
use :: fundal
use :: penf, only : I4P, R8P

implicit none
private
public :: rk_assign_stage_dev
public :: rk_compute_stage_dev
public :: rk_compute_stage_ls_dev
public :: rk_initialize_stages_dev
public :: rk_update_q_dev

contains
   ! public procedures
   subroutine rk_assign_stage_dev(ni,nj,nk,ngc,nv,blocks_number,s,phi_gpu,q_gpu,q_rk_gpu)
   !< Assign q to RK stage.
   integer(I4P), intent(in)           :: ni                                      !< Grid cells number in I direction.
   integer(I4P), intent(in)           :: nj                                      !< Grid cells number in J direction.
   integer(I4P), intent(in)           :: nk                                      !< Grid cells number in K direction.
   integer(I4P), intent(in)           :: ngc                                     !< Ghost cells number.
   integer(I4P), intent(in)           :: nv                                      !< Number of conservative varibales.
   integer(I4P), intent(in)           :: blocks_number                           !< Number of blocks.
   integer(I4P), intent(in)           :: s                                       !< Stage index.
   real(R8P),    intent(in), optional ::  phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< IB distance.
   real(R8P),    intent(in)           ::    q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative field.
   real(R8P),    intent(inout)        :: q_rk_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:) !< Conservative field stage.
   integer(I4P)                       :: all_solids                              !< Last phi index, all solids summary.
   integer(I4P)                       :: i, j, k, b, v                           !< Counter.

   if (present(phi_gpu)) then
      all_solids = ubound(phi_gpu, dim=5)
      !$acc parallel loop independent DEVICEVAR(phi_gpu, q_gpu, q_rk_gpu)
      !$omp OMPLOOP DEVICEPTR(phi_gpu, q_gpu, q_rk_gpu)
      do v=1, nv
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do b=1, blocks_number
                     if (phi_gpu(b,i,j,k,all_solids) < 0._R8P) then
                        q_rk_gpu(b,i,j,k,v,s) = q_gpu(b,i,j,k,v)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
   else
      !$acc parallel loop independent DEVICEVAR(q_gpu, q_rk_gpu)
      !$omp OMPLOOP DEVICEPTR(q_gpu, q_rk_gpu)
      do v=1, nv
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do b=1, blocks_number
                     q_rk_gpu(b,i,j,k,v,s) = q_gpu(b,i,j,k,v)
                  enddo
               enddo
            enddo
         enddo
      enddo
   endif
   endsubroutine rk_assign_stage_dev

   subroutine rk_compute_stage_dev(ni,nj,nk,ngc,nv,blocks_number,s,dt,alph,phi_gpu,q_rk_gpu)
   !< Sum RK stages up to s.
   integer(I4P), intent(in)           :: ni                                      !< Grid cells number in I direction.
   integer(I4P), intent(in)           :: nj                                      !< Grid cells number in J direction.
   integer(I4P), intent(in)           :: nk                                      !< Grid cells number in K direction.
   integer(I4P), intent(in)           :: ngc                                     !< Ghost cells number.
   integer(I4P), intent(in)           :: nv                                      !< Number of conservative varibales.
   integer(I4P), intent(in)           :: blocks_number                           !< Number of blocks.
   integer(I4P), intent(in)           :: s                                       !< Current stage.
   real(R8P),    intent(in)           :: dt                                      !< Time step.
   real(R8P),    intent(in)           :: alph(1:,1:)                             !< RK alpha coefficients.
   real(R8P),    intent(in), optional ::  phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< IB distance.
   real(R8P),    intent(inout)        :: q_rk_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:) !< Conservative field stages.
   integer(I4P)                       :: all_solids                              !< Last phi index, all solids summary.
   integer(I4P)                       :: i, j, k, b, v, ss                       !< Counter.
   real(R8P)                          :: dq                                      !< Local increment.

   if (present(phi_gpu)) then
      all_solids = ubound(phi_gpu, dim=5)
      !$acc parallel loop collapse(5) independent DEVICEVAR(phi_gpu, alph, q_rk_gpu)&
      !$acc& private(dq) firstprivate(dt) 
      !$omp OMPLOOP collapse(5) DEVICEPTR(phi_gpu, alph, q_rk_gpu)&
      !$omp& private(dq) firstprivate(dt)
      do v=1, nv
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do b=1, blocks_number
                     if (phi_gpu(b,i,j,k,all_solids) < 0._R8P) then
                        dq = 0._R8P
                        !$acc loop seq
                        do ss=1, s-1
                           dq = dq + dt * alph(s,ss) * q_rk_gpu(b,i,j,k,v,ss)
                        enddo
                        q_rk_gpu(b,i,j,k,v,s) = q_rk_gpu(b,i,j,k,v,s) + dq
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
   else
      !$acc parallel loop collapse(5) independent DEVICEVAR(alph, q_rk_gpu)&
      !$acc& private(dq) firstprivate(dt) 
      !$omp OMPLOOP collapse(5) DEVICEPTR(alph, q_rk_gpu)&
      !$omp& private(dq) firstprivate(dt)
      do v=1, nv
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do b=1, blocks_number
                     dq = 0._R8P
                     !$acc loop seq
                     do ss=1, s-1
                        dq = dq + dt * alph(s,ss) * q_rk_gpu(b,i,j,k,v,ss)
                     enddo
                     q_rk_gpu(b,i,j,k,v,s) = q_rk_gpu(b,i,j,k,v,s) + dq
                  enddo
               enddo
            enddo
         enddo
      enddo
   endif
   endsubroutine rk_compute_stage_dev

   subroutine rk_compute_stage_ls_dev(ni,nj,nk,ngc,nv,blocks_number,dt,ark,brk,crk,phi_gpu,q_n_gpu,dq_gpu,q_rk_gpu)
   !< Compute RK stage, low storage scheme using only q(n) and q.
   integer(I4P), intent(in)           :: ni                                   !< Grid cells number in I direction.
   integer(I4P), intent(in)           :: nj                                   !< Grid cells number in J direction.
   integer(I4P), intent(in)           :: nk                                   !< Grid cells number in K direction.
   integer(I4P), intent(in)           :: ngc                                  !< Ghost cells number.
   integer(I4P), intent(in)           :: nv                                   !< Number of conservative varibales.
   integer(I4P), intent(in)           :: blocks_number                        !< Number of blocks.
   real(R8P),    intent(in)           :: dt                                   !< Time step.
   real(R8P),    intent(in)           :: ark, brk, crk                        !< RK coefficients.
   real(R8P),    intent(in), optional ::  phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< IB distance.
   real(R8P),    intent(in)           ::  q_n_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< RK stage 0, Q^n.
   real(R8P),    intent(in)           ::   dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Residuals.
   real(R8P),    intent(inout)        :: q_rk_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field stage.
   integer(I4P)                       :: all_solids                           !< Last phi index, all solids summary.
   integer(I4P)                       :: i, j, k, b, v                        !< Counter.

   if (present(phi_gpu)) then
      all_solids = ubound(phi_gpu, dim=5)
      !$acc parallel loop independent DEVICEVAR(phi_gpu, q_n_gpu, dq_gpu, q_rk_gpu)
      !$omp OMPLOOP DEVICEPTR(phi_gpu, q_n_gpu, dq_gpu, q_rk_gpu)
      do v=1, nv
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do b=1, blocks_number
                     if (phi_gpu(b,i,j,k,all_solids) < 0._R8P) then
                        q_rk_gpu(b,i,j,k,v) = ark * q_n_gpu(b,i,j,k,v) + brk * q_rk_gpu(b,i,j,k,v) + dt * crk * dq_gpu(b,i,j,k,v)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
   else
      !$acc parallel loop independent DEVICEVAR(q_n_gpu, dq_gpu, q_rk_gpu)
      !$omp OMPLOOP DEVICEPTR(q_n_gpu, dq_gpu, q_rk_gpu)
      do v=1, nv
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do b=1, blocks_number
                     q_rk_gpu(b,i,j,k,v) = ark * q_n_gpu(b,i,j,k,v) + brk * q_rk_gpu(b,i,j,k,v) + dt * crk * dq_gpu(b,i,j,k,v)
                  enddo
               enddo
            enddo
         enddo
      enddo
   endif
   endsubroutine rk_compute_stage_ls_dev

   subroutine rk_initialize_stages_dev(ni,nj,nk,ngc,nv,blocks_number,q_gpu,q_rk_gpu)
   !< Initialize RK stages.
   integer(I4P), intent(in)    :: ni                                      !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                      !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                      !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                     !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                                      !< Number of conservative varibales.
   integer(I4P), intent(in)    :: blocks_number                           !< Number of blocks.
   real(R8P),    intent(in)    ::    q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative field.
   real(R8P),    intent(inout) :: q_rk_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:) !< RK stage.
   integer(I4P)                :: i, j, k, b, v,s                         !< Counter.

   !$acc parallel loop independent DEVICEVAR(q_gpu, q_rk_gpu)
   !$omp OMPLOOP DEVICEPTR(q_gpu, q_rk_gpu)
   do s=lbound(q_rk_gpu,dim=6),ubound(q_rk_gpu,dim=6)
      do v=1, nv
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do b=1, blocks_number
                     q_rk_gpu(b,i,j,k,v,s) = q_gpu(b,i,j,k,v)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine rk_initialize_stages_dev


   subroutine rk_update_q_dev(ni,nj,nk,ngc,nv,blocks_number,nrk,dt,beta,phi_gpu,q_rk_gpu,q_gpu)
   !< Update RK q.
   integer(I4P), intent(in)           :: ni
   integer(I4P), intent(in)           :: nj
   integer(I4P), intent(in)           :: nk
   integer(I4P), intent(in)           :: ngc
   integer(I4P), intent(in)           :: nv
   integer(I4P), intent(in)           :: blocks_number
   integer(I4P), intent(in)           :: nrk
   real(R8P),    intent(in)           :: dt
   real(R8P),    intent(in)           :: beta(1:)
   real(R8P),    intent(in), optional :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(in)           :: q_rk_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:)
   real(R8P),    intent(inout)        :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P)                       :: all_solids
   integer(I4P)                       :: i, j, k, b, v, s
   real(R8P)                          :: dq

   if (present(phi_gpu)) then
      all_solids = ubound(phi_gpu, dim=5)

      !$acc parallel loop collapse(5) independent DEVICEVAR(phi_gpu, beta, q_rk_gpu, q_gpu)&
      !$acc& private(dq) firstprivate(dt) 
      !$omp OMPLOOP collapse(5) DEVICEPTR(phi_gpu, beta, q_rk_gpu, q_gpu)&
      !$omp& private(dq) firstprivate(dt)
      do v=1, nv
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do b=1, blocks_number
                     if (phi_gpu(b,i,j,k,all_solids) < 0._R8P) then
                        dq = 0._R8P
                        !$acc loop seq
                        do s=1, nrk
                           dq = dq + beta(s) * q_rk_gpu(b,i,j,k,v,s)
                        enddo
                        q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + dt * dq
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo

   else

      !$acc parallel loop collapse(5) independent DEVICEVAR(beta, q_rk_gpu, q_gpu)&
      !$acc& private(dq) firstprivate(dt) 
      !$omp OMPLOOP collapse(5) DEVICEPTR(beta, q_rk_gpu, q_gpu)&
      !$omp& private(dq) firstprivate(dt)
      do v=1, nv
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do b=1, blocks_number
                     dq = 0._R8P
                     !$acc loop seq
                     do s=1, nrk
                        dq = dq + beta(s) * q_rk_gpu(b,i,j,k,v,s)
                     enddo
                     q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + dt * dq
                  enddo
               enddo
            enddo
         enddo
      enddo

   endif
   endsubroutine rk_update_q_dev
endmodule adam_fnl_rk_kernels
