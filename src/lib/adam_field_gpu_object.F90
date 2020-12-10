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
   integer(I4P)                :: error=0_I4P       !< Error traping flag.
   ! GPU data.
   ! RK data, related to field equations
   real(R8P), allocatable, device :: alph_gpu(:,:)       !< RK alpha coefficients.
   real(R8P), allocatable, device :: beta_gpu(:)         !< RK beta  coefficients.
   real(R8P), allocatable, device :: gamm_gpu(:)         !< RK gamma coefficients.
   ! field equations data
   real(R8P), allocatable,    device :: u_gpu(     :,:,:,:  )    !< Field cell centered variables [ni+2gci,nj+2gcj,nk+2gck,nv,nb].
   real(R8P), allocatable,    device :: u_work_gpu(:,:,:,:  )    !< Field working buffer.
   real(R8P), allocatable,    device :: u_s_gpu(   :,:,:,:,:)    !< RK field stages.
   integer(I8P), allocatable, device :: local_map_ghost_gpu(:,:) !< Local map for ghost cells updating.
   ! MPI data
   integer(I4P) :: mydev=0_I4P      !< My GPU rank.
   integer(I4P) :: local_comm=0_I4P !< Local communicator.
   real(R8P),    allocatable, device :: send_buffer_ghost_gpu(:)       !< Send buffer of ghost cells.
   real(R8P),    allocatable, device :: recv_buffer_ghost_gpu(:)       !< Receive buffer of ghost cells.
   integer(I4P), allocatable, device :: comm_map_send_ptr_ghost_gpu(:) !< Communication map, pointers in list to send.
   integer(I4P), allocatable, device :: comm_map_recv_ptr_ghost_gpu(:) !< Communication map, pointers in list to recv.
   integer(I8P), allocatable, device :: comm_map_recv_ghost_gpu(:,:)   !< Communication map, `fec` information.
   integer(I8P), allocatable, device :: comm_map_send_ghost_gpu(:,:)   !< Communication map, `fec` information.
   contains
      ! public methods
      procedure, pass(self) :: allocate_gpu !< Allocate GPU data.
      procedure, pass(self) :: initialize   !< Initialize field.
      procedure, pass(self) :: copy_cpu_gpu !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu !< Copy data from GPU to CPU.
      procedure, pass(self) :: rk_integrate !< Runge Kutta integration of field.
endtype field_gpu_object

contains
      ! public methods
   subroutine initialize(self, field_cpu)
   !< Initialize field.
   class(field_gpu_object), intent(inout)      :: self       !< The field.
   type(field_object),      intent(in), target :: field_cpu  !< CPU field data.

   self%field_cpu => field_cpu
#ifdef _MPI_
   call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, self%local_comm, self%error)
   call MPI_COMM_RANK(self%local_comm, self%mydev,self%error)
#endif
   self%error = CudaSetDevice(self%mydev)
   write(*,*) "MPI rank ", self%field_cpu%myrank, " using GPU ", self%mydev
   call self%allocate_gpu
   endsubroutine initialize

   subroutine allocate_gpu(self)
   !< Allocate GPU data.
   class(field_gpu_object), intent(inout) :: self !< The field.

   associate(grid=>self%field_cpu%grid)
      if (allocated(self%u_gpu)) deallocate(self%u_gpu)
      allocate(self%u_gpu(1-grid%gci:grid%ni+grid%gci, &
                          1-grid%gcj:grid%nj+grid%gcj, &
                          1-grid%gck:grid%nk+grid%gck, 1:self%field_cpu%nb))
      if (allocated(self%u_work_gpu)) deallocate(self%u_work_gpu)
      allocate(self%u_work_gpu(1-grid%gci:grid%ni+grid%gci, &
                               1-grid%gcj:grid%nj+grid%gcj, &
                               1-grid%gck:grid%nk+grid%gck, 1:self%field_cpu%nb))
      if (allocated(self%u_s_gpu)) deallocate(self%u_s_gpu)
      allocate(self%u_s_gpu(1-grid%gci:grid%ni+grid%gci, &
                            1-grid%gcj:grid%nj+grid%gcj, &
                            1-grid%gck:grid%nk+grid%gck, 1:self%field_cpu%nb, 1:3))
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
   if (allocated(self%local_map_ghost_gpu        )) deallocate(self%local_map_ghost_gpu        )
   if (allocated(self%send_buffer_ghost_gpu      )) deallocate(self%send_buffer_ghost_gpu      )
   if (allocated(self%recv_buffer_ghost_gpu      )) deallocate(self%recv_buffer_ghost_gpu      )
   if (allocated(self%comm_map_send_ptr_ghost_gpu)) deallocate(self%comm_map_send_ptr_ghost_gpu)
   if (allocated(self%comm_map_recv_ptr_ghost_gpu)) deallocate(self%comm_map_recv_ptr_ghost_gpu)
   if (allocated(self%comm_map_recv_ghost_gpu    )) deallocate(self%comm_map_recv_ghost_gpu    )
   if (allocated(self%comm_map_send_ghost_gpu    )) deallocate(self%comm_map_send_ghost_gpu    )
   self%local_map_ghost_gpu         = self%field_cpu%local_map_ghost
   self%send_buffer_ghost_gpu       = self%field_cpu%send_buffer_ghost
   self%recv_buffer_ghost_gpu       = self%field_cpu%recv_buffer_ghost
   self%comm_map_send_ptr_ghost_gpu = self%field_cpu%comm_map_send_ptr_ghost
   self%comm_map_recv_ptr_ghost_gpu = self%field_cpu%comm_map_recv_ptr_ghost
   self%comm_map_recv_ghost_gpu     = self%field_cpu%comm_map_recv_ghost
   self%comm_map_send_ghost_gpu     = self%field_cpu%comm_map_send_ghost
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

   call rk_integrate_gpu(u=self%u_gpu,                                                 &
                         u_work=self%u_work_gpu,                                       &
                         u_s=self%u_s_gpu,                                             &
                         blocks_number=self%field_cpu%blocks_number,                   &
                         ni=self%field_cpu%grid%ni,                                    &
                         nj=self%field_cpu%grid%nj,                                    &
                         nk=self%field_cpu%grid%nk,                                    &
                         t=t,                                                          &
                         Dt=Dt,                                                        &
                         alph=self%alph_gpu,                                           &
                         beta=self%beta_gpu,                                           &
                         gamm=self%field_cpu%gamm,                                     &
                         local_map_ghost_gpu=self%local_map_ghost_gpu,                 &
                         comm_map_recv_ptr_ghost_gpu=self%comm_map_recv_ptr_ghost_gpu, &
                         comm_map_send_ptr_ghost_gpu=self%comm_map_send_ptr_ghost_gpu, &
                         comm_map_recv_ghost_gpu=self%comm_map_recv_ghost_gpu,         &
                         comm_map_send_ghost_gpu=self%comm_map_send_ghost_gpu,         &
                         recv_buffer_ghost_gpu=self%recv_buffer_ghost_gpu,             &
                         send_buffer_ghost_gpu=self%send_buffer_ghost_gpu,             &
                         procs_number=self%field_cpu%procs_number)
   endsubroutine rk_integrate

   subroutine rk_integrate_gpu(u, u_work, u_s, alph, beta, gamm, blocks_number, ni, nj, nk, t, Dt, &
                               local_map_ghost_gpu,                                                &
                               comm_map_recv_ptr_ghost_gpu, comm_map_send_ptr_ghost_gpu,           &
                               comm_map_recv_ghost_gpu, comm_map_send_ghost_gpu,                   &
                               recv_buffer_ghost_gpu, send_buffer_ghost_gpu, procs_number)
   real(R8P),    intent(inout), device              :: u(:,:,:,:)                     !< Field cell centered variables.
   real(R8P),    intent(inout), device              :: u_work(:,:,:,:)                !< Field working buffer.
   real(R8P),    intent(inout), device              :: u_s(:,:,:,:,:)                 !< RK field stages.
   real(R8P),    intent(in),    device              :: alph(:,:)                      !< RK alpha coefficients.
   real(R8P),    intent(in),    device              :: beta(:)                        !< RK beta coefficients.
   real(R8P),    intent(in)                         :: gamm(:)                        !< RK gamma coefficients.
   integer(I4P), intent(in)                         :: blocks_number                  !< Number of blocks actually stored.
   integer(I4P), intent(in)                         :: ni                             !< Number of cell in I direction.
   integer(I4P), intent(in)                         :: nj                             !< Number of cell in J direction.
   integer(I4P), intent(in)                         :: nk                             !< Number of cell in K direction.
   real(R8P),    intent(in)                         :: t                              !< Time.
   real(R8P),    intent(in)                         :: Dt                             !< Time step.
   integer(I8P), intent(in),    device, allocatable :: local_map_ghost_gpu(:,:)       !< Local map for ghost cells updating.
   integer(I4P), intent(in),    device, allocatable :: comm_map_send_ptr_ghost_gpu(:) !< Communication map, pointers in list to send.
   integer(I4P), intent(in),    device, allocatable :: comm_map_recv_ptr_ghost_gpu(:) !< Communication map, pointers in list to recv.
   integer(I8P), intent(in),    device, allocatable :: comm_map_recv_ghost_gpu(:,:)   !< Communication map, `fec` information.
   integer(I8P), intent(in),    device, allocatable :: comm_map_send_ghost_gpu(:,:)   !< Communication map, `fec` information.
   real(R8P),    intent(inout), device              :: recv_buffer_ghost_gpu(:)       !< Receive buffer of ghost cells.
   real(R8P),    intent(inout), device              :: send_buffer_ghost_gpu(:)       !< Send buffer of ghost cells.
   integer(I4P), intent(in)                         :: procs_number                   !< Number of MPI processes.
   integer(I4P)                                     :: i, j, k, b, s, ss              !< Counter.
   integer(I4P)                                     :: iercuda                        !< Error trapping flag for CUDAFortran.

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
      call update_ghost_gpu(local_map_ghost_gpu=local_map_ghost_gpu, u_s=u_s, s=s)
      call update_ghost_mpi_gpu(comm_map_recv_ptr_ghost_gpu=comm_map_recv_ptr_ghost_gpu, &
                                comm_map_send_ptr_ghost_gpu=comm_map_send_ptr_ghost_gpu, &
                                comm_map_recv_ghost_gpu=comm_map_recv_ghost_gpu,         &
                                comm_map_send_ghost_gpu=comm_map_send_ghost_gpu,         &
                                recv_buffer_ghost_gpu=recv_buffer_ghost_gpu,             &
                                send_buffer_ghost_gpu=send_buffer_ghost_gpu,             &
                                u_s=u_s, s=s, procs_number=procs_number)
      call compute_residuals_gpu(u_work=u_work, u_s=u_s, block_start=1, block_end=blocks_number, &
                                 ni=ni, nj=nj, nk=nk, s=s, t=t + gamm(s) * Dt)
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

   subroutine compute_residuals_gpu(u_work, u_s, block_start, block_end, ni, nj, nk, s, t)
   real(R8P),    intent(inout), device :: u_work(:,:,:,:) !< Field working buffer.
   real(R8P),    intent(inout), device :: u_s(:,:,:,:,:)  !< RK field stages.
   integer(I4P), intent(in)            :: block_start     !< Index of block to start residuals computation.
   integer(I4P), intent(in)            :: block_end       !< Index of block to end   residuals computation.
   integer(I4P), intent(in)            :: ni              !< Number of cell in I direction.
   integer(I4P), intent(in)            :: nj              !< Number of cell in J direction.
   integer(I4P), intent(in)            :: nk              !< Number of cell in K direction.
   integer(I4P), intent(in)            :: s               !< Working stage.
   real(R8P),    intent(in)            :: t               !< Time.
   integer(I4P)                        :: b, i, j, k      !< Counter.
   integer(I4P)                        :: im, jm, km      !< Counter for fake bc.
   integer(I4P)                        :: ip, jp, kp      !< Counter for fake bc.
   integer(I4P)                        :: iercuda         !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(4) <<<*,*>>>
   do b=block_start, block_end
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
   do b=block_start, block_end
      do k=1, nk
         do j=1, nj
            do i=1, ni
               u_s(i,j,k,b,s) = u_work(i,j,k,b)
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_residuals_gpu

   subroutine update_ghost_gpu(local_map_ghost_gpu, u_s, s)
   !< Update ghost cells.
   integer(I8P), intent(in),    device, allocatable :: local_map_ghost_gpu(:,:) !< Local map for ghost cells updating.
   real(R8P),    intent(inout), device              :: u_s(:,:,:,:,:)           !< RK field stages.
   integer(I4P), intent(in)                         :: s                        !< Stage.
   integer(I4P)                                     :: i, j, k, mf              !< Counter.
   integer(I4P)                                     :: iii, jjj, kkk            !< Counter.
   integer(I4P)                                     :: fec                      !< Direction where ghost cells are updated, fec.
   integer(I4P)                                     :: portion                  !< Portion of fec updated (0=>whole fec).
   integer(I4P)                                     :: b_recv                   !< Index of receiving block.
   integer(I4P)                                     :: b_send                   !< Index of sending block.
   integer(I4P)                                     :: imin                     !< Lower limit of i indexes.
   integer(I4P)                                     :: jmin                     !< Lower limit of j indexes.
   integer(I4P)                                     :: kmin                     !< Lower limit of j indexes.
   integer(I4P)                                     :: imax                     !< Upper limit of i indexes.
   integer(I4P)                                     :: jmax                     !< Upper limit of j indexes.
   integer(I4P)                                     :: kmax                     !< Upper limit of k indexes.
   integer(I4P)                                     :: idelta                   !< Delta offset for ghost-inner cells of i indexes.
   integer(I4P)                                     :: jdelta                   !< Delta offset for ghost-inner cells of j indexes.
   integer(I4P)                                     :: kdelta                   !< Delta offset for ghost-inner cells of k indexes.
   integer(I4P)                                     :: iercuda                  !< Error trapping flag for CUDAFortran.

   if (.not.allocated(local_map_ghost_gpu)) return
   do mf=1, size(local_map_ghost_gpu, dim=1)
      b_recv  = local_map_ghost_gpu(mf, 1 )
      b_send  = local_map_ghost_gpu(mf, 2 )
      fec     = local_map_ghost_gpu(mf, 3 )
      portion = local_map_ghost_gpu(mf, 4 )
      imin    = local_map_ghost_gpu(mf, 5 )
      jmin    = local_map_ghost_gpu(mf, 6 )
      kmin    = local_map_ghost_gpu(mf, 7 )
      imax    = local_map_ghost_gpu(mf, 8 )
      jmax    = local_map_ghost_gpu(mf, 9 )
      kmax    = local_map_ghost_gpu(mf, 10)
      idelta  = local_map_ghost_gpu(mf, 11)
      jdelta  = local_map_ghost_gpu(mf, 12)
      kdelta  = local_map_ghost_gpu(mf, 13)
      if     (portion==0) then
         ! receiving from a block with the same refinement
         !$cuf kernel do(3) <<<*,*>>>
         do k=kmin, kmax
            do j=jmin, jmax
               do i=imin, imax
                  u_s(i,j,k,b_recv,s) = u_s(i+idelta,j+jdelta,k+kdelta,b_send,s)
               enddo
            enddo
         enddo
         !@cuf iercuda=cudaDeviceSynchronize()
      elseif (portion>0) then
         ! receiving from a block finer than me
         !$cuf kernel do(3) <<<*,*>>>
         do k=kmin, kmax
            do j=jmin, jmax
               do i=imin, imax
                  kkk = 2 * k + kdelta
                  jjj = 2 * j + jdelta
                  iii = 2 * i + idelta
                  u_s(i,j,k,b_recv,s) = (u_s(iii,jjj,  kkk,  b_send,s) + u_s(iii+1,jjj,  kkk,  b_send,s) + &
                                         u_s(iii,jjj+1,kkk,  b_send,s) + u_s(iii+1,jjj+1,kkk,  b_send,s) + &
                                         u_s(iii,jjj,  kkk+1,b_send,s) + u_s(iii+1,jjj,  kkk+1,b_send,s) + &
                                         u_s(iii,jjj+1,kkk+1,b_send,s) + u_s(iii+1,jjj+1,kkk+1,b_send,s)) / 8._R8P
               enddo
            enddo
         enddo
         !@cuf iercuda=cudaDeviceSynchronize()
      else
         ! receiving from a block coarser than me
         !$cuf kernel do(3) <<<*,*>>>
         do k=kmin, kmax
            do j=jmin, jmax
               do i=imin, imax
                  kkk = 2 * k + kdelta
                  jjj = 2 * j + jdelta
                  iii = 2 * i + idelta
                  u_s(iii,  jjj,  kkk  ,b_recv,s) = u_s(i,j,k,b_send,s)
                  u_s(iii+1,jjj,  kkk  ,b_recv,s) = u_s(i,j,k,b_send,s)
                  u_s(iii,  jjj+1,kkk  ,b_recv,s) = u_s(i,j,k,b_send,s)
                  u_s(iii+1,jjj+1,kkk  ,b_recv,s) = u_s(i,j,k,b_send,s)
                  u_s(iii,  jjj,  kkk+1,b_recv,s) = u_s(i,j,k,b_send,s)
                  u_s(iii+1,jjj,  kkk+1,b_recv,s) = u_s(i,j,k,b_send,s)
                  u_s(iii,  jjj+1,kkk+1,b_recv,s) = u_s(i,j,k,b_send,s)
                  u_s(iii+1,jjj+1,kkk+1,b_recv,s) = u_s(i,j,k,b_send,s)
               enddo
            enddo
         enddo
         !@cuf iercuda=cudaDeviceSynchronize()
      endif
   enddo
   endsubroutine update_ghost_gpu

   subroutine update_ghost_mpi_gpu(comm_map_recv_ptr_ghost_gpu, comm_map_send_ptr_ghost_gpu, &
                                   comm_map_recv_ghost_gpu, comm_map_send_ghost_gpu,         &
                                   recv_buffer_ghost_gpu, send_buffer_ghost_gpu,             &
                                   u_s, s, procs_number, step)
   !< Update ghost cells within other processes.
   integer(I4P), allocatable, intent(in),    device :: comm_map_send_ptr_ghost_gpu(:) !< Communication map, pointers list to send.
   integer(I4P), allocatable, intent(in),    device :: comm_map_recv_ptr_ghost_gpu(:) !< Communication map, pointers list to recv.
   integer(I8P), allocatable, intent(in),    device :: comm_map_recv_ghost_gpu(:,:)   !< Communication map, `fec` information.
   integer(I8P), allocatable, intent(in),    device :: comm_map_send_ghost_gpu(:,:)   !< Communication map, `fec` information.
   real(R8P),                 intent(inout), device :: recv_buffer_ghost_gpu(:)       !< Receive buffer of ghost cells.
   real(R8P),                 intent(inout), device :: send_buffer_ghost_gpu(:)       !< Send buffer of ghost cells.
   real(R8P),                 intent(inout), device :: u_s(:,:,:,:,:)                 !< RK field stages.
   integer(I4P),              intent(in)            :: s                              !< Stage.
   integer(I4P),              intent(in)            :: procs_number                   !< Number of MPI processes.
   integer(I4P),              intent(in), optional  :: step                           !< Step to be perfordmed.
   logical                                          :: steps(3)                       !< Steps to be performed.
   integer(I4P)                                     :: i, j, k                        !< Counter.
   integer(I4P)                                     :: iii, jjj, kkk                  !< Counter.
   integer(I4P)                                     :: fec, mf, rf, sf, n, p          !< Counter.
   integer(I4P)                                     :: portion                        !< Portion of fec updated (0=>whole fec).
   integer(I4P)                                     :: b_recv                         !< Index of receiving block.
   integer(I4P)                                     :: b_send                         !< Index of sending block.
   integer(I4P)                                     :: imin                           !< Lower limit of i indexes.
   integer(I4P)                                     :: jmin                           !< Lower limit of j indexes.
   integer(I4P)                                     :: kmin                           !< Lower limit of j indexes.
   integer(I4P)                                     :: imax                           !< Upper limit of i indexes.
   integer(I4P)                                     :: jmax                           !< Upper limit of j indexes.
   integer(I4P)                                     :: kmax                           !< Upper limit of k indexes.
   integer(I4P)                                     :: idelta                         !< Delta offset for ghost-inner cells of i.
   integer(I4P)                                     :: jdelta                         !< Delta offset for ghost-inner cells of j.
   integer(I4P)                                     :: kdelta                         !< Delta offset for ghost-inner cells of k.
   integer(I4P)                                     :: ptr_start, ptr_end             !< Counter.
   integer(I4P)                                     :: send_ptr, send_ctr             !< Counter.
   integer(I4P)                                     :: recv_ptr, recv_ctr             !< Counter.
   integer(I4P)                                     :: n_recv, n_send                 !< Counter.
   integer(I4P)                                     :: recv_rank                      !< Rank of receiving block.
   integer(I4P)                                     :: send_rank                      !< Rank of sending block.
   integer(I4P)                                     :: error                          !< Error traping flag.
   integer(I4P), allocatable                        :: req_send_recv(:)               !< MPI request receive flags.
   integer(I4P)                                     :: iercuda                        !< Error trapping flag for CUDAFortran.

   steps = .true.

   if (present(step)) then
      steps = .false.
      steps(step) = .true.
   endif

   if (steps(1)) then
#ifdef _MPI_
      allocate(req_send_recv(0:procs_number*2-1))
      req_send_recv = MPI_REQUEST_NULL
#endif

      if ((.not.allocated(comm_map_recv_ghost_gpu)).and.(.not.allocated(comm_map_send_ghost_gpu))) return

      ! populate send buffer
      !$cuf kernel do(1) <<<*,*>>>
      do sf=1, size(comm_map_send_ghost_gpu, dim=1)
         b_send    = comm_map_send_ghost_gpu(sf, 2 ) ! neighbor-block-index of block
         send_rank = comm_map_send_ghost_gpu(sf, 3 )
         fec       = comm_map_send_ghost_gpu(sf, 4 )
         portion   = comm_map_send_ghost_gpu(sf, 5 )
         imin      = comm_map_send_ghost_gpu(sf, 6 )
         jmin      = comm_map_send_ghost_gpu(sf, 7 )
         kmin      = comm_map_send_ghost_gpu(sf, 8 )
         imax      = comm_map_send_ghost_gpu(sf, 9 )
         jmax      = comm_map_send_ghost_gpu(sf, 10)
         kmax      = comm_map_send_ghost_gpu(sf, 11)
         idelta    = comm_map_send_ghost_gpu(sf, 12)
         jdelta    = comm_map_send_ghost_gpu(sf, 13)
         kdelta    = comm_map_send_ghost_gpu(sf, 14)
         send_ptr  = comm_map_send_ghost_gpu(sf, 15)
         if (portion==0_I4P) then
            ! sending to a block at my level
            send_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     send_buffer_ghost_gpu(send_ptr + send_ctr) = u_s(i+idelta,j+jdelta,k+kdelta,b_send,s)
                     send_ctr = send_ctr + 1
                  enddo
               enddo
            enddo
         elseif (portion<0_I4P) then ! Beware! This is < 0 because the reference is the receiver
            ! sending to a block finer than me
            send_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     do n=1,8
                        send_buffer_ghost_gpu(send_ptr + send_ctr) = u_s(i,j,k,b_send,s)
                        send_ctr = send_ctr + 1
                     enddo
                  enddo
               enddo
            enddo
         else
            ! sending to a block coarser than me, loop is over the coarser grid
            send_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     kkk = 2 * k + kdelta
                     jjj = 2 * j + jdelta
                     iii = 2 * i + idelta
                     send_buffer_ghost_gpu(send_ptr + send_ctr) =         &
                        (u_s(iii,jjj,  kkk,  b_send,s) + u_s(iii+1,jjj,  kkk,  b_send,s) + &
                         u_s(iii,jjj+1,kkk,  b_send,s) + u_s(iii+1,jjj+1,kkk,  b_send,s) + &
                         u_s(iii,jjj,  kkk+1,b_send,s) + u_s(iii+1,jjj,  kkk+1,b_send,s) + &
                         u_s(iii,jjj+1,kkk+1,b_send,s) + u_s(iii+1,jjj+1,kkk+1,b_send,s)) / 8._R8P
                     send_ctr = send_ctr + 1
                  enddo
               enddo
            enddo
         endif
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   endif

   if (steps(2)) then
      ! receive
      do p=0, procs_number - 1_I4P
         ptr_start = comm_map_recv_ptr_ghost_gpu(p) + 1
         ptr_end   = comm_map_recv_ptr_ghost_gpu(p+1)
         n_recv    = ptr_end - ptr_start + 1
         if (n_recv > 0) then
#ifdef _MPI_
            call MPI_IRECV(recv_buffer_ghost_gpu(ptr_start), n_recv, MPI_REAL8, p, 100, MPI_COMM_WORLD, &
                           req_send_recv(p), error)
#endif
         endif
      enddo

      ! send
      do p=0, procs_number - 1_I4P
         ptr_start = comm_map_send_ptr_ghost_gpu(p) + 1
         ptr_end   = comm_map_send_ptr_ghost_gpu(p+1)
         n_send    = ptr_end - ptr_start + 1
         if (n_send > 0) then
#ifdef _MPI_
            call MPI_ISEND(send_buffer_ghost_gpu(ptr_start), n_send, MPI_REAL8, p, 100, MPI_COMM_WORLD, &
                           req_send_recv(p+procs_number), error)
#endif
         endif
      enddo
   endif

   if (steps(3)) then

#ifdef _MPI_
      call MPI_WAITALL(procs_number * 2, req_send_recv, &
                       MPI_STATUSES_IGNORE, error)
#endif

      ! call MPI_BARRIER(MPI_COMM_WORLD, error)

      ! retrive from receive buffer
      !$cuf kernel do(1) <<<*,*>>>
      do rf=1, size(comm_map_recv_ghost_gpu, dim=1)
         b_recv    = comm_map_recv_ghost_gpu(rf, 1 ) ! block-index
         recv_rank = comm_map_recv_ghost_gpu(rf, 3 )
         fec       = comm_map_recv_ghost_gpu(rf, 4 )
         portion   = comm_map_recv_ghost_gpu(rf, 5 )
         imin      = comm_map_recv_ghost_gpu(rf, 6 )
         jmin      = comm_map_recv_ghost_gpu(rf, 7 )
         kmin      = comm_map_recv_ghost_gpu(rf, 8 )
         imax      = comm_map_recv_ghost_gpu(rf, 9 )
         jmax      = comm_map_recv_ghost_gpu(rf, 10)
         kmax      = comm_map_recv_ghost_gpu(rf, 11)
         idelta    = comm_map_recv_ghost_gpu(rf, 12)
         jdelta    = comm_map_recv_ghost_gpu(rf, 13)
         kdelta    = comm_map_recv_ghost_gpu(rf, 14)
         recv_ptr  = comm_map_recv_ghost_gpu(rf, 15)
         if (portion==0_I4P) then
            recv_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     u_s(i,j,k,b_recv,s) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr)
                     recv_ctr = recv_ctr + 1
                  enddo
               enddo
            enddo
         elseif (portion>0_I4P) then
            ! receiving from a block finer than me
            recv_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     u_s(i,j,k,b_recv,s) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr)
                     recv_ctr = recv_ctr + 1
                  enddo
               enddo
            enddo
         else
            ! receiving from a block coarser than me
            recv_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     kkk = 2 * k + kdelta
                     jjj = 2 * j + jdelta
                     iii = 2 * i + idelta
                     u_s(iii,  jjj,  kkk  ,b_recv,s) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr)
                     recv_ctr = recv_ctr + 1
                     u_s(iii+1,jjj,  kkk  ,b_recv,s) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr)
                     recv_ctr = recv_ctr + 1
                     u_s(iii,  jjj+1,kkk  ,b_recv,s) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr)
                     recv_ctr = recv_ctr + 1
                     u_s(iii+1,jjj+1,kkk  ,b_recv,s) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr)
                     recv_ctr = recv_ctr + 1
                     u_s(iii,  jjj,  kkk+1,b_recv,s) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr)
                     recv_ctr = recv_ctr + 1
                     u_s(iii+1,jjj,  kkk+1,b_recv,s) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr)
                     recv_ctr = recv_ctr + 1
                     u_s(iii,  jjj+1,kkk+1,b_recv,s) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr)
                     recv_ctr = recv_ctr + 1
                     u_s(iii+1,jjj+1,kkk+1,b_recv,s) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr)
                     recv_ctr = recv_ctr + 1
                  enddo
               enddo
            enddo
         endif
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   endif
   endsubroutine update_ghost_mpi_gpu
endmodule adam_field_gpu_object
