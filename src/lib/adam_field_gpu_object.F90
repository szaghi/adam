!< ADAM, GPU field class definition.
module adam_field_gpu_object
!< ADAM, GPU field class definition.

use adam_field_object
use PENF
#ifdef _MPI_
use MPI
#endif
use CUDAFOR

implicit none
private
public :: field_gpu_object

type :: field_gpu_object
   !< GPU Field class definition.
   type(field_object), pointer :: field_cpu=>null() !< Pointer to CPU field data.
   integer(I4P)                :: error=0_I4P !< Error traping flag.
   ! GPU data.
   real(R8P), allocatable, device :: u_gpu(:,:,:,:)      !< Field cell centered variables [ni+gc12,nj+gc34,nk+gc56,nv,nb].
   real(R8P), allocatable, device :: u_work_gpu(:,:,:,:) !< Field working buffer.
   real(R8P), allocatable, device :: u_s_gpu(:,:,:,:,:)  !< RK field stages.
   real(R8P), allocatable, device :: alph_gpu(:,:)       !< RK alpha coefficients.
   real(R8P), allocatable, device :: beta_gpu(:)         !< RK beta  coefficients.
   real(R8P), allocatable, device :: gamm_gpu(:)         !< RK gamma coefficients.
   ! MPI data
   integer(I4P) :: mydev=0_I4P !< My GPU rank.
   contains
      ! public methods
      procedure, pass(self) :: allocate_gpu !< Allocate GPU data.
      procedure, pass(self) :: initialize   !< Initialize field.
      procedure, pass(self) :: copy_cpu_gpu !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu !< Copy data from GPU to CPU.
      procedure, pass(self) :: rk_integrate !< Runge Kutta integration of field.
      ! procedure, pass(self) :: residuals    !< Compute residuals of field.
endtype field_gpu_object

contains
      ! public methods
   subroutine initialize(self, field_cpu)
   !< Initialize field.
   class(field_gpu_object), intent(inout)      :: self       !< The field.
   type(field_object),      intent(in), target :: field_cpu  !< CPU field data.
   integer(I4P)                                :: local_comm !< Local communicator.

   self%field_cpu => field_cpu
#ifdef _MPI_
   call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, local_comm, self%error)
   call MPI_COMM_RANK(local_comm, self%mydev,self%error)
#endif
   self%error = cudaSetDevice(self%mydev)
   write(*,*) "MPI rank ", self%field_cpu%myrank, " using GPU ", self%mydev
   call self%allocate_gpu
   endsubroutine initialize

   subroutine allocate_gpu(self)
   !< Allocate GPU data.
   class(field_gpu_object), intent(inout) :: self !< The field.

   associate(grid=>self%field_cpu%grid)
      if (allocated(self%u_gpu   )) deallocate(self%u_gpu   )
      allocate(self%u_gpu(1-grid%gc1:grid%ni+grid%gc2, &
                          1-grid%gc3:grid%nj+grid%gc4, &
                          1-grid%gc5:grid%nk+grid%gc6, 1:self%field_cpu%nb))
      if (allocated(self%u_work_gpu )) deallocate(self%u_work_gpu )
      allocate(self%u_work_gpu(1-grid%gc1:grid%ni+grid%gc2, &
                               1-grid%gc3:grid%nj+grid%gc4, &
                               1-grid%gc5:grid%nk+grid%gc6, 1:self%field_cpu%nb))
      if (allocated(self%u_s_gpu )) deallocate(self%u_s_gpu )
      allocate(self%u_s_gpu(1-grid%gc1:grid%ni+grid%gc2, &
                            1-grid%gc3:grid%nj+grid%gc4, &
                            1-grid%gc5:grid%nk+grid%gc6, 1:self%field_cpu%nb, 1:3))
   endassociate

   if (allocated(self%alph_gpu)) deallocate(self%alph_gpu) ; allocate(self%alph_gpu(3,3))
   if (allocated(self%beta_gpu)) deallocate(self%beta_gpu) ; allocate(self%beta_gpu(3))
   if (allocated(self%gamm_gpu)) deallocate(self%gamm_gpu) ; allocate(self%gamm_gpu(3))
   endsubroutine allocate_gpu

   subroutine copy_cpu_gpu(self)
   !< Copy data from CPU to GPU.
   class(field_gpu_object), intent(inout) :: self !< The field.

   self%u_gpu    = self%field_cpu%u
   self%alph_gpu = self%field_cpu%alph
   self%beta_gpu = self%field_cpu%beta
   self%gamm_gpu = self%field_cpu%gamm
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self)
   !< Copy data from GPU to CPU.
   class(field_gpu_object), intent(inout) :: self !< The field.

   self%field_cpu%u = self%u_gpu
   endsubroutine copy_gpu_cpu

   subroutine rk_integrate(self, t, Dt)
   !< Runge Kutta integration of field.
   class(field_gpu_object), intent(inout) :: self  !< The field.
   real(R8P),               intent(in)    :: t     !< Time.
   real(R8P),               intent(in)    :: Dt    !< Time step.

   call rk_integrate_gpu(u=self%u_gpu,                               &
                         u_work=self%u_work_gpu,                     &
                         u_s=self%u_s_gpu,                           &
                         blocks_number=self%field_cpu%blocks_number, &
                         ni=self%field_cpu%grid%ni,                  &
                         nj=self%field_cpu%grid%nj,                  &
                         nk=self%field_cpu%grid%nk,                  &
                         t=t,                                        &
                         Dt=Dt,                                      &
                         alph=self%alph_gpu,                         &
                         beta=self%beta_gpu,                         &
                         gamm=self%field_cpu%gamm)
   endsubroutine rk_integrate

   subroutine rk_integrate_gpu(u, u_work, u_s, alph, beta, gamm, blocks_number, ni, nj, nk, t, Dt)
   real(R8P),    intent(inout), device :: u(:,:,:,:)        !< Field cell centered variables.
   real(R8P),    intent(inout), device :: u_work(:,:,:,:)   !< Field working buffer.
   real(R8P),    intent(inout), device :: u_s(:,:,:,:,:)    !< RK field stages.
   real(R8P),    intent(in),    device :: alph(:,:)         !< RK alpha coefficients.
   real(R8P),    intent(in),    device :: beta(:)           !< RK beta coefficients.
   real(R8P),    intent(in)            :: gamm(:)           !< RK gamma coefficients.
   integer(I4P), intent(in)            :: blocks_number     !< Number of blocks actually stored.
   integer(I4P), intent(in)            :: ni                !< Number of cell in I direction.
   integer(I4P), intent(in)            :: nj                !< Number of cell in J direction.
   integer(I4P), intent(in)            :: nk                !< Number of cell in K direction.
   real(R8P),    intent(in)            :: t                 !< Time.
   real(R8P),    intent(in)            :: Dt                !< Time step.
   integer(I4P)                        :: i, j, k, b, s, ss !< Counter.
   integer(I4P)                        :: iercuda           !< Error trapping flag for CUDAFortran.

   do s=1, 3
      !$cuf kernel do(4) <<<*,*>>>
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  u_s(i,j,k,b,s) = u(i,j,k,b)
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      do ss=1, s - 1
         !$cuf kernel do(4) <<<*,*>>>
         do b=1, blocks_number
            do k=1, nk
               do j=1, nj
                  do i=1, ni
                     u_s(i,j,k,b,s) = u_s(i,j,k,b,s) + (u_s(i,j,k,b,ss) * (Dt * alph(s, ss)))
                  enddo
               enddo
            enddo
         enddo
         !@cuf iercuda=cudaDeviceSynchronize()
      enddo
      call residuals_gpu(u_work=u_work, u_s=u_s, blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, s=s, t=t + gamm(s) * Dt)
   enddo
   do s=1, 3
      !$cuf kernel do(4) <<<*,*>>>
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  u(i,j,k,b) = u(i,j,k,b) + u_s(i,j,k,b,s) * Dt * beta(s)
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   enddo
   endsubroutine rk_integrate_gpu

   subroutine residuals_gpu(u_work, u_s, blocks_number, ni, nj, nk, s, t)
   real(R8P),    intent(inout), device :: u_work(:,:,:,:)   !< Field working buffer.
   real(R8P),    intent(inout), device :: u_s(:,:,:,:,:)    !< RK field stages.
   integer(I4P), intent(in)            :: blocks_number     !< Number of blocks actually stored.
   integer(I4P), intent(in)            :: ni                !< Number of cell in I direction.
   integer(I4P), intent(in)            :: nj                !< Number of cell in J direction.
   integer(I4P), intent(in)            :: nk                !< Number of cell in K direction.
   integer(I4P), intent(in)            :: s                 !< Working stage.
   real(R8P),    intent(in)            :: t                 !< Time.
   integer(I4P)                        :: b, i, j, k        !< Counter.
   integer(I4P)                        :: im, jm, km        !< Counter for fake bc.
   integer(I4P)                        :: ip, jp, kp        !< Counter for fake bc.
   integer(I4P)                        :: iercuda           !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(4) <<<*,*>>>
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               kp = min(nk, k+1)
               km = max(1, k-1)
               jp = min(nj, j+1)
               jm = max(1, j-1)
               ip = min(ni, i+1)
               im = max(1, i-1)
               u_work(i,j,k,b) = u_s(ip,j, k, b,s) + u_s(im,j, k, b,s) + &
                                 u_s(i, jp,k, b,s) + u_s(i, jm,k, b,s) + &
                                 u_s(i, j, kp,b,s) + u_s(i, j, km,b,s) - &
                        6._R8P * u_s(i, j, k, b,s)

            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   !$cuf kernel do(4) <<<*,*>>>
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               u_s(i,j,k,b,s) = u_work(i,j,k,b)
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine residuals_gpu
endmodule adam_field_gpu_object
