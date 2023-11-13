!< ADAM, RK NVF kernels (NVF backend of [[rk_nvf_object]]).
module adam_rk_nvf_kernels
!< ADAM, RK NVF kernels (NVF backend of [[rk_nvf_object]]).

use penf, only : I4P, R8P
use cudafor

implicit none
private
public :: rk_assign_stage_cuf
public :: rk_compute_stage_cuf
public :: rk_compute_stage_ls_cuf
public :: rk_initialize_stages_cuf
public :: rk_update_q_cuf

contains
   ! public procedures
   subroutine rk_assign_stage_cuf(ni,nj,nk,ngc,nv,blocks_number,s,phi_gpu,q_gpu,q_rk_gpu)
   !< Assign q to RK stage.
   integer(I4P), intent(in),    value            :: ni                                      !< Grid cells number in I direction.
   integer(I4P), intent(in),    value            :: nj                                      !< Grid cells number in J direction.
   integer(I4P), intent(in),    value            :: nk                                      !< Grid cells number in K direction.
   integer(I4P), intent(in),    value            :: ngc                                     !< Ghost cells number.
   integer(I4P), intent(in),    value            :: nv                                      !< Number of conservative varibales.
   integer(I4P), intent(in),    value            :: blocks_number                           !< Number of blocks.
   integer(I4P), intent(in),    value            :: s                                       !< Stage index.
   real(R8P),    intent(in),    device, optional ::  phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< IB distance.
   real(R8P),    intent(in),    device           ::    q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative field.
   real(R8P),    intent(inout), device           :: q_rk_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:) !< Conservative field stage.
   integer(I4P)                                  :: all_solids                              !< Last phi index, all solids summary.
   integer(I4P)                                  :: i, j, k, b, v                           !< Counter.
   integer(I4P)                                  :: iercuda                                 !< Error trapping flag for CUDAFortran.

   if (present(phi_gpu)) then
      all_solids = ubound(phi_gpu, dim=5)
      !$cuf kernel do(5) <<<*,*>>>
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
      !@cuf iercuda=cudaDeviceSynchronize()
   else
      !$cuf kernel do(5) <<<*,*>>>
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
      !@cuf iercuda=cudaDeviceSynchronize()
   endif
   endsubroutine rk_assign_stage_cuf

   subroutine rk_compute_stage_cuf(ni,nj,nk,ngc,nv,blocks_number,s,dt,alph,phi_gpu,q_rk_gpu)
   !< Sum RK stages up to s.
   integer(I4P), intent(in),    value            :: ni                                      !< Grid cells number in I direction.
   integer(I4P), intent(in),    value            :: nj                                      !< Grid cells number in J direction.
   integer(I4P), intent(in),    value            :: nk                                      !< Grid cells number in K direction.
   integer(I4P), intent(in),    value            :: ngc                                     !< Ghost cells number.
   integer(I4P), intent(in),    value            :: nv                                      !< Number of conservative varibales.
   integer(I4P), intent(in),    value            :: blocks_number                           !< Number of blocks.
   integer(I4P), intent(in),    value            :: s                                       !< Current stage.
   real(R8P),    intent(in),    value            :: dt                                      !< Time step.
   real(R8P),    intent(in),    device           :: alph(1:,1:)                             !< RK alpha coefficients.
   real(R8P),    intent(in),    device, optional ::  phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< IB distance.
   real(R8P),    intent(inout), device           :: q_rk_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:) !< Conservative field stages.
   integer(I4P)                                  :: all_solids                              !< Last phi index, all solids summary.
   integer(I4P)                                  :: i, j, k, b, v, ss                       !< Counter.
   integer(I4P)                                  :: iercuda                                 !< Error trapping flag for CUDAFortran.

   if (present(phi_gpu)) then
      all_solids = ubound(phi_gpu, dim=5)
      !$cuf kernel do(6) <<<*,*>>>
      do ss=1, s-1
         do v=1, nv
            do k=1, nk
               do j=1, nj
                  do i=1, ni
                     do b=1, blocks_number
                        if (phi_gpu(b,i,j,k,all_solids) < 0._R8P) then
                           q_rk_gpu(b,i,j,k,v,s) = q_rk_gpu(b,i,j,k,v,s) + dt * alph(s,ss) * q_rk_gpu(b,i,j,k,v,ss)
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   else
      !$cuf kernel do(6) <<<*,*>>>
      do ss=1, s-1
         do v=1, nv
            do k=1, nk
               do j=1, nj
                  do i=1, ni
                     do b=1, blocks_number
                        q_rk_gpu(b,i,j,k,v,s) = q_rk_gpu(b,i,j,k,v,s) + dt * alph(s, ss) * q_rk_gpu(b,i,j,k,v,ss)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   endif
   endsubroutine rk_compute_stage_cuf

   subroutine rk_compute_stage_ls_cuf(ni,nj,nk,ngc,nv,blocks_number,dt,ark,brk,crk,phi_gpu,q_n_gpu,dq_gpu,q_rk_gpu)
   !< Compute RK stage, low storage scheme using only q(n) and q.
   integer(I4P), intent(in),    value            :: ni                                   !< Grid cells number in I direction.
   integer(I4P), intent(in),    value            :: nj                                   !< Grid cells number in J direction.
   integer(I4P), intent(in),    value            :: nk                                   !< Grid cells number in K direction.
   integer(I4P), intent(in),    value            :: ngc                                  !< Ghost cells number.
   integer(I4P), intent(in),    value            :: nv                                   !< Number of conservative varibales.
   integer(I4P), intent(in),    value            :: blocks_number                        !< Number of blocks.
   real(R8P),    intent(in),    value            :: dt                                   !< Time step.
   real(R8P),    intent(in),    value            :: ark, brk, crk                        !< RK coefficients.
   real(R8P),    intent(in),    device, optional ::  phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< IB distance.
   real(R8P),    intent(in),    device           ::  q_n_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< RK stage 0, Q^n.
   real(R8P),    intent(in),    device           ::   dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Residuals.
   real(R8P),    intent(inout), device           :: q_rk_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field stage.
   integer(I4P)                                  :: all_solids                           !< Last phi index, all solids summary.
   integer(I4P)                                  :: i, j, k, b, v                        !< Counter.
   integer(I4P)                                  :: iercuda                              !< Error trapping flag for CUDAFortran.

   if (present(phi_gpu)) then
      all_solids = ubound(phi_gpu, dim=5)
      !$cuf kernel do(5) <<<*,*>>>
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
      !@cuf iercuda=cudaDeviceSynchronize()
   else
      !$cuf kernel do(5) <<<*,*>>>
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
      !@cuf iercuda=cudaDeviceSynchronize()
   endif
   endsubroutine rk_compute_stage_ls_cuf

   subroutine rk_initialize_stages_cuf(ni,nj,nk,ngc,nv,blocks_number,q_gpu,q_rk_gpu)
   !< Initialize RK stages.
   integer(I4P), intent(in)            :: ni                                      !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                      !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                      !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                     !< Ghost cells number.
   integer(I4P), intent(in)            :: nv                                      !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number                           !< Number of blocks.
   real(R8P),    intent(in),    device ::    q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative field.
   real(R8P),    intent(inout), device :: q_rk_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:) !< RK stage.
   integer(I4P)                        :: i, j, k, b, v,s                         !< Counter.
   integer(I4P)                        :: iercuda                                 !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(6) <<<*,*>>>
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
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine rk_initialize_stages_cuf

   subroutine rk_update_q_cuf(ni,nj,nk,ngc,nv,blocks_number,nrk,dt,beta,phi_gpu,q_rk_gpu,q_gpu)
   !< Update RK q.
   integer(I4P), intent(in),    value            :: ni                                      !< Grid cells number in I direction.
   integer(I4P), intent(in),    value            :: nj                                      !< Grid cells number in J direction.
   integer(I4P), intent(in),    value            :: nk                                      !< Grid cells number in K direction.
   integer(I4P), intent(in),    value            :: ngc                                     !< Ghost cells number.
   integer(I4P), intent(in),    value            :: nv                                      !< Number of conservative varibales.
   integer(I4P), intent(in),    value            :: blocks_number                           !< Number of blocks.
   integer(I4P), intent(in),    value            :: nrk                                     !< Number of RK stages.
   real(R8P),    intent(in),    value            :: dt                                      !< Time step.
   real(R8P),    intent(in),    device           :: beta(1:)                                !< RK betaa coefficients.
   real(R8P),    intent(in),    device, optional ::  phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< IB distance.
   real(R8P),    intent(in),    device           :: q_rk_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:) !< Conservative field stage.
   real(R8P),    intent(inout), device           ::    q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative field.
   integer(I4P)                                  :: all_solids                              !< Last phi index, all solids summary.
   integer(I4P)                                  :: i, j, k, b, v, s                        !< Counter.
   integer(I4P)                                  :: iercuda                                 !< Error trapping flag for CUDAFortran.

   if (present(phi_gpu)) then
      all_solids = ubound(phi_gpu, dim=5)
      !$cuf kernel do(6) <<<*,*>>>
      do s=1, nrk
         do v=1, nv
            do k=1, nk
               do j=1, nj
                  do i=1, ni
                     do b=1, blocks_number
                        if (phi_gpu(b,i,j,k,all_solids) < 0._R8P) then
                           q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + dt * beta(s) * q_rk_gpu(b,i,j,k,v,s)
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   else
      !$cuf kernel do(6) <<<*,*>>>
      do s=1, nrk
         do v=1, nv
            do k=1, nk
               do j=1, nj
                  do i=1, ni
                     do b=1, blocks_number
                        q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + dt * beta(s) * q_rk_gpu(b,i,j,k,v,s)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   endif
   endsubroutine rk_update_q_cuf
endmodule adam_rk_nvf_kernels
