!< ADAM, 1D convenction equation class definition, GPU backend.
module adam_equation_convect1D_gpu_object
!< ADAM, 1D convenction equation class definition, GPU backend.

use adam_base_gpu_object
use adam_field_object
use adam_parameters
use PENF
use MPI
use CUDAFOR

implicit none
private
public :: equation_convect1D_gpu_object
public :: BC_EXTRAPOLATION
public :: BC_INFLOW

integer(I4P), parameter :: BC_EXTRAPOLATION = 1_I4P
integer(I4P), parameter :: BC_INFLOW        = 2_I4P

type :: equation_convect1D_gpu_object
   !< 1D convenction equation class definition, GPU backend.
   type(field_object), pointer    :: field=>null()         !< The field.
   type(base_gpu_object)          :: base_gpu              !< The base GPU handler.
   integer(I4P)                   :: myrank=0_I4P          !< MPI rank process.
   integer(I4P)                   :: procs_number=1_I4P    !< Number of MPI processes.
   integer(I4P)                   :: error=0_I4P           !< Error traping flag.
   integer(I4P)                   :: ns=3_I4P              !< Runge-Kutta stages number.
   real(R8P), allocatable         :: alph(:,:)             !< RK alpha coefficients.
   real(R8P), allocatable         :: beta(:)               !< RK beta coefficients.
   real(R8P), allocatable         :: gamm(:)               !< RK gamma coefficients.
   ! cuf data
   real(R8P), allocatable, device :: dxyz_gpu(:,:)         !< Space steps.
   real(R8P), allocatable, device :: alph_gpu(:,:)         !< RK alpha coefficients.
   real(R8P), allocatable, device :: beta_gpu(:)           !< RK beta coefficients.
   real(R8P), allocatable, device :: gamm_gpu(:)           !< RK gamma coefficients.
   real(R8P), allocatable, device :: q_gpu(:,:,:,:,:)      !< Field cell centered variables stages.
   real(R8P), allocatable, device :: q_work_gpu(:,:,:,:,:) !< Field cell centered variables stages, working buffer.
   real(R8P), allocatable, device :: q_s_gpu(:,:,:,:,:,:)  !< RK Field cell centered variables stages.
   contains
      ! public methods
      procedure, pass(self) :: copy_cpu_gpu            !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu            !< Copy data from GPU to CPU.
      procedure, pass(self) :: destroy                 !< Destroy the equation.
      procedure, pass(self) :: initialize              !< Initialize the equation.
      procedure, pass(self) :: mark_by_grad_q          !< Mark blocks to be refined/derefined by a `grad(q)` value.
      procedure, pass(self) :: integrate               !< Runge Kutta integration of equation.
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions of equation.
      procedure, pass(self) :: set_initial_conditions  !< Set initial conditions of equation.
      procedure, pass(self) :: update_ghost_gpu        !< Update ghost cells and set boundary conditions.
      ! operators
      generic :: assignment(=) => eq_assign_eq      !< Overload `=`.
      procedure, pass(lhs), private :: eq_assign_eq !< Operator `=`.
endtype equation_convect1D_gpu_object

contains
   ! public methods
   subroutine copy_cpu_gpu(self)
   !< Copy data from CPU to GPU.
   class(equation_convect1D_gpu_object), intent(inout) :: self !< The base backend.

   self%dxyz_gpu = self%field%dxyz
   call self%base_gpu%copy_transpose_cpu_gpu(q_cpu=self%field%q, q_gpu=self%q_gpu)
   call self%base_gpu%copy_cpu_gpu
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self)
   !< Copy data from GPU to CPU.
   class(equation_convect1D_gpu_object), intent(inout) :: self !< The base backend.

   call self%base_gpu%copy_transpose_gpu_cpu(nv=self%field%nv, q_gpu=self%q_gpu, q_cpu=self%field%q)
   endsubroutine copy_gpu_cpu

   subroutine destroy(self)
   !< Destroy the equation.
   class(equation_convect1D_gpu_object), intent(inout) :: self  !< The equation.
   type(equation_convect1D_gpu_object)                 :: fresh !< Fresh equation.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, field, ns)
   !< Initialize the equation.
   class(equation_convect1D_gpu_object), intent(inout)        :: self  !< The equation.
   type(field_object),                   intent(in), target   :: field !< The field.
   integer(I4P),                         intent(in), optional :: ns    !< Runge-Kutta stages number.

   call self%destroy
   self%field => field
   call self%base_gpu%initialize(field=field)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   if (present(ns)) self%ns = ns
   allocate(self%alph(self%ns,self%ns), self%beta(self%ns), self%gamm(self%ns))
   select case(self%ns)
   case(3_I4P)
      self%alph(:,:) = reshape([0._R8P, 1._R8P, 0.25_R8P, &
                                0._R8P, 0._R8P, 0.25_R8P, &
                                0._R8P, 0._R8P,0._R8P], [3,3])
      self%beta(:) = [1._R8P/6._R8P, &
                      1._R8P/6._R8P, &
                      2._R8P/3._R8P]
      self%gamm(:) = [0._R8P, &
                      1._R8P, &
                      0._R8P]
   endselect
   self%alph_gpu = self%alph
   self%beta_gpu = self%beta
   self%gamm_gpu = self%gamm
   allocate(self%q_gpu(1:field%nb,                                    &
                       1-field%grid%gci:field%grid%ni+field%grid%gci, &
                       1-field%grid%gcj:field%grid%nj+field%grid%gcj, &
                       1-field%grid%gck:field%grid%nk+field%grid%gck, 1:field%nv))
   allocate(self%q_work_gpu(1:field%nb,                                    &
                            1-field%grid%gci:field%grid%ni+field%grid%gci, &
                            1-field%grid%gcj:field%grid%nj+field%grid%gcj, &
                            1-field%grid%gck:field%grid%nk+field%grid%gck, 1:field%nv))
   allocate(self%dxyz_gpu(1:3, 1:field%nb))
   allocate(self%q_s_gpu(1:field%nb,                                    &
                         1-field%grid%gci:field%grid%ni+field%grid%gci, &
                         1-field%grid%gcj:field%grid%nj+field%grid%gcj, &
                         1-field%grid%gck:field%grid%nk+field%grid%gck, 1:field%nv, 1:self%ns))
   self%dxyz_gpu = self%field%dxyz
   endsubroutine initialize

   subroutine mark_by_grad_q(self, grad_tol, delta_fine, delta_coarse, threshold)
   !< Mark blocks to be refined/derefined by a `grad(q)` value.
   class(equation_convect1D_gpu_object), intent(inout)        :: self           !< The equation.
   real(R8P),                            intent(in)           :: grad_tol       !< Gradiend tolerance value.
   real(R8P),                            intent(in)           :: delta_fine     !< Maximum cell delta in fine grids.
   real(R8P),                            intent(in)           :: delta_coarse   !< Minimum cell delta in coarse grids.
   real(R8P),                            intent(in), optional :: threshold      !< Threshold for sphere proximity.
   real(R8P)                                                  :: threshold_     !< Threshold for sphere proximity, local var.
   real(R8P)                                                  :: max_cell_delta !< Maximum cell delta.
   real(R8P)                                                  :: grad_q         !< Value (max) of gradient of q.
   integer(I4P)                                               :: b              !< Counter.

   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   self%field%refinements_needed = [(TO_NOT_TOUCH,b=1,self%field%blocks_number)]
   call self%update_ghost_gpu(q_gpu=self%q_gpu)
   associate (gci=>self%field%grid%gci, gcj=>self%field%grid%gcj, gck=>self%field%grid%gck, &
              ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk, dxyz=>self%field%dxyz)
   do b=1, self%field%blocks_number
      grad_q = gradient_cuf(b=b, ni=ni, nj=nj, nk=nk, gci=gci, gcj=gcj, gck=gck, &
                            dx=dxyz(1,b), dy=dxyz(2,b), dz=dxyz(3,b), q_gpu=self%q_gpu)

      max_cell_delta = max_cell_delta_grad(grad=grad_q)

      if (maxval(dxyz(:,b)) > max_cell_delta) then
         self%field%refinements_needed(b) = TO_BE_REFINED
      elseif (maxval(dxyz(:,b)) * threshold_ < max_cell_delta) then
         self%field%refinements_needed(b) = TO_BE_DEREFINED
      else
         self%field%refinements_needed(b) = TO_NOT_TOUCH
      endif
   enddo
   endassociate
   contains
      function gradient_cuf(b, ni, nj, nk, gci, gcj, gck, dx, dy, dz, q_gpu) result(gradient)
      integer(I4P), intent(in)         :: b            !< Counter.
      integer(I4P), intent(in)         :: ni           !< Grid cells number in I direction.
      integer(I4P), intent(in)         :: nj           !< Grid cells number in J direction.
      integer(I4P), intent(in)         :: nk           !< Grid cells number in K direction.
      integer(I4P), intent(in)         :: gci          !< Ghost grid cells number in I direction.
      integer(I4P), intent(in)         :: gcj          !< Ghost grid cells number in J direction.
      integer(I4P), intent(in)         :: gck          !< Ghost grid cells number in K direction.
      real(R8P),    intent(in)         :: dx           !< X space step.
      real(R8P),    intent(in)         :: dy           !< Y space step.
      real(R8P),    intent(in)         :: dz           !< Z space step.
      real(R8P),    intent(in), device :: q_gpu(1:,    &
                                                1-gci:,&
                                                1-gcj:,&
                                                1-gck:,&
                                                1:)    !< Field component to be updated.
      real(R8P)                        :: gradient     !< Maximum gradient of q.
      real(R8P)                        :: grad         !< Current gradient of q.
      integer(I4P)                     :: i, j, k      !< Counter.
      integer(I4P)                     :: iercuda      !< Error trapping flag for CUDAFortran.

      gradient = 0._R8P
      !$cuf kernel do(3) <<<*,*>>>
      do k=1, nk
         do j=1, nj
            do i=1, ni
               grad = sqrt(((q_gpu(b,i+1,j,k,1) - q_gpu(b,i-1,j,k,1))/(2*dx))**2 + &
                           ((q_gpu(b,i,j+1,k,1) - q_gpu(b,i,j-1,k,1))/(2*dy))**2 + &
                           ((q_gpu(b,i,j,k+1,1) - q_gpu(b,i,j,k-1,1))/(2*dz))**2)
               gradient = max(gradient, grad)

            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
      endfunction gradient_cuf

      function max_cell_delta_grad(grad) result(delta)
      !< Return the maximum cell delta given a gradient tollerance.
      real(R8P), intent(in) :: grad  !< Gradient value.
      real(R8P)             :: delta !< Maximum cell delta admissible.

      if (grad > grad_tol) then
         delta = delta_fine
      else
         delta = delta_coarse
      endif
      endfunction max_cell_delta_grad
   endsubroutine mark_by_grad_q

   subroutine integrate(self, t, Dt, do_ghost_syncro, residual)
   !< Runge Kutta integration of field.
   class(equation_convect1D_gpu_object), intent(inout)         :: self             !< The equation.
   real(R8P),                            intent(in)            :: t                !< Time.
   real(R8P),                            intent(in)            :: Dt               !< Time step.
   logical,                              intent(in),  optional :: do_ghost_syncro  !< Flag to do syncrous ghost update.
   real(R8P),                            intent(out), optional :: residual         !< Global residual.
   logical                                                     :: do_ghost_syncro_ !< Flag to do syncrous ghost update, local var.
   integer(I4P)                                                :: s                !< Counter.

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   associate(alph=>self%alph, beta=>self%beta, gamm=>self%gamm,                            &
             ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk,       &
             gci=>self%field%grid%gci, gcj=>self%field%grid%gcj, gck=>self%field%grid%gck, &
             blocks_number=>self%field%blocks_number,                                      &
             inner_blocks_number=>self%field%inner_blocks_number,                          &
             alph_gpu=>self%alph_gpu, beta_gpu=>self%beta_gpu)
   do s=1, self%ns
      call compute_rk_stage_gpu_cuf(ni=ni, nj=nj, nk=nk, gci=gci, gcj=gcj, gck=gck, blocks_number=blocks_number, &
                                    alph_gpu=alph_gpu, dt=dt, s=s, q_gpu=self%q_gpu, q_s_gpu=self%q_s_gpu)
      if (do_ghost_syncro_) then
         call self%update_ghost_gpu(q_gpu=self%q_s_gpu(:,:,:,:,:,s)) ! all ghosts
         call compute_residuals_gpu_cuf(ni=ni, nj=nj, nk=nk, gci=gci, gcj=gcj, gck=gck,  &
                                        blocks_number=blocks_number, t=t + gamm(s) * dt, &
                                        dxyz=self%dxyz_gpu, q_work_gpu=self%q_work_gpu,  &
                                        q_gpu=self%q_s_gpu(:,:,:,:,:,s))
      else
         ! TODO
      endif
      if (present(residual).and.s==self%ns) then
         ! TODO
      endif
   enddo
   call advance_q_gpu_cuf(ni=ni, nj=nj, nk=nk, gci=gci, gcj=gcj, gck=gck, blocks_number=blocks_number, &
                           beta_gpu=beta_gpu, dt=dt, q_s_gpu=self%q_s_gpu, q_gpu=self%q_gpu)
   endassociate
   endsubroutine integrate

   subroutine set_boundary_conditions(self, q_gpu)
   !< Set boundary conditions of equation.
   class(equation_convect1D_gpu_object), intent(in)            :: self                             !< The equation.
   real(R8P),                            intent(inout), device :: q_gpu(1:,                    &
                                                                        1-self%field%grid%gci:,&
                                                                        1-self%field%grid%gcj:,&
                                                                        1-self%field%grid%gck:,1:) !< Field.

   if (allocated(self%base_gpu%local_map_bc_crown_gpu)) call set_bc_fec(nv=self%field%nv,       &
                                                                        gc=self%field%grid%gci, &
                                                                        local_map_bc=self%base_gpu%local_map_bc_crown_gpu)
   contains
      subroutine set_bc_fec(nv, gc, local_map_bc)
      integer(I4P), intent(in)         :: nv                !< Number of variables.
      integer(I4P), intent(in)         :: gc                !< Ghost cells number.
      integer(I8P), intent(in), device :: local_map_bc(:,:,:) !< Local map for BC ghost cells.
      integer(I4P)                     :: b                 !< Counter.
      integer(I4P)                     :: c, i, j, k, v     !< Counter.
      integer(I4P)                     :: idelta            !< IJK delta step for extrapolation.
      integer(I4P)                     :: jdelta            !< IJK delta step for extrapolation.
      integer(I4P)                     :: kdelta            !< IJK delta step for extrapolation.
      integer(I4P)                     :: bc_type           !< Boundary condition type.
      integer(I4P)                     :: crown             !< Crown counter.
      integer(I4P)                     :: iercuda           !< Error trapping flag for CUDAFortran.

      do crown=1, gc
         !$cuf kernel do(1) <<<*,*>>>
         do c=1, size(local_map_bc, dim=1)
            b = local_map_bc(c, 1 ,crown)
            if (b>0) then
               i       = local_map_bc(c, 2 ,crown)
               j       = local_map_bc(c, 3 ,crown)
               k       = local_map_bc(c, 4 ,crown)
               idelta  = local_map_bc(c, 5 ,crown)
               jdelta  = local_map_bc(c, 6 ,crown)
               kdelta  = local_map_bc(c, 7 ,crown)
               bc_type = local_map_bc(c, 8 ,crown)
               if (bc_type == BC_EXTRAPOLATION) then
                  do v=1, nv
                     q_gpu(b,i,j,k,v) = q_gpu(b,i-idelta,j-jdelta,k-kdelta,v)
                  enddo
               else
               endif
            endif
         enddo
         !@cuf iercuda=cudaDeviceSynchronize()
      enddo
      endsubroutine set_bc_fec
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self)
   !< Set initial conditions of field.
   class(equation_convect1D_gpu_object), intent(inout) :: self    !< The equation.
   integer(I4P)                                        :: b       !< Counter.
   integer(I4P)                                        :: i, j, k !< Counter.
   real(R8P)                                           :: a       !< Gaussian amplitude.
   real(R8P)                                           :: sigma_x !< Gaussian x variance.
   real(R8P)                                           :: sigma_y !< Gaussian y variance.
   real(R8P)                                           :: sigma_z !< Gaussian z variance.
   real(R8P)                                           :: x_0     !< Gaussian x center.
   real(R8P)                                           :: y_0     !< Gaussian y center.
   real(R8P)                                           :: z_0     !< Gaussian z center.

   a = 1.0_R8P
   x_0 = (self%field%grid%domain_emax(1) - self%field%grid%domain_emin(1)) / 5.0_R8P
   y_0 = (self%field%grid%domain_emax(2) - self%field%grid%domain_emin(2)) / 2.0_R8P
   z_0 = (self%field%grid%domain_emax(3) - self%field%grid%domain_emin(3)) / 2.0_R8P
   sigma_x = 0.05_R8P
   sigma_y = 0.05_R8P
   sigma_z = 0.05_R8P
   associate(blocks_number=>self%field%blocks_number,                                      &
             q=>self%field%q,                                                              &
             ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk,       &
             gci=>self%field%grid%gci, gcj=>self%field%grid%gcj, gck=>self%field%grid%gck, &
             x_cell=>self%field%x_cell, y_cell=>self%field%y_cell, z_cell=>self%field%z_cell)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               q(1,i,j,k,b) = a * exp(-((x_cell(i,b) - x_0)**2/(2 * sigma_x**2)+&
                                        (y_cell(j,b) - y_0)**2/(2 * sigma_y**2)+&
                                        (z_cell(k,b) - z_0)**2/(2 * sigma_z**2)))
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine set_initial_conditions

   subroutine update_ghost_gpu(self, q_gpu, step)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(equation_convect1D_gpu_object), intent(inout)         :: self            !< The equation.
   real(R8P),                            intent(inout), device :: q_gpu(1:,                    &
                                                                        1-self%field%grid%gci:,&
                                                                        1-self%field%grid%gcj:,&
                                                                        1-self%field%grid%gck:,&
                                                                        1:)       !< Field component to be updated.
   integer(I4P),                         intent(in), optional  :: step            !< Step to be perfordmed in asyncronous comp.
   logical                                                     :: do_local_update !< Flag for triggering local update.
   logical                                                     :: do_set_bc       !< Flag for triggering setting bc.

   ! perform local update if step is not speficied or if first step is selected
   do_local_update = .false.
   do_set_bc       = .false.
   if (.not.present(step)) then
      do_local_update = .true.
      do_set_bc       = .true.
   else
      if (step==1) do_local_update = .true.
      if (step==3) do_set_bc       = .true.
   endif

   if (do_local_update) call self%base_gpu%update_ghost_local_gpu(q_gpu=q_gpu)
                        call self%base_gpu%update_ghost_mpi_gpu(q_gpu=q_gpu, step=step)
   if (do_set_bc)       call self%set_boundary_conditions(q_gpu=q_gpu)
   endsubroutine update_ghost_gpu

   ! operators
   ! =
   subroutine eq_assign_eq(lhs, rhs)
   !< Operator `=`.
   class(equation_convect1D_gpu_object), intent(inout) :: lhs !< Left hand side.
   type(equation_convect1D_gpu_object),  intent(in)    :: rhs !< Right hand side.

   lhs%field => rhs%field
   lhs%base_gpu = rhs%base_gpu
   lhs%myrank = rhs%myrank
   lhs%procs_number = rhs%procs_number
   lhs%error = rhs%error
   lhs%ns = rhs%ns
   call assign_allocatable(lhs=lhs%alph, rhs=rhs%alph)
   call assign_allocatable(lhs=lhs%beta, rhs=rhs%beta)
   call assign_allocatable(lhs=lhs%gamm, rhs=rhs%gamm)
   call assign_allocatable_gpu(lhs=lhs%dxyz_gpu,   rhs=rhs%dxyz_gpu  )
   call assign_allocatable_gpu(lhs=lhs%alph_gpu,   rhs=rhs%alph_gpu  )
   call assign_allocatable_gpu(lhs=lhs%beta_gpu,   rhs=rhs%beta_gpu  )
   call assign_allocatable_gpu(lhs=lhs%gamm_gpu,   rhs=rhs%gamm_gpu  )
   call assign_allocatable_gpu(lhs=lhs%q_gpu,      rhs=rhs%q_gpu     )
   call assign_allocatable_gpu(lhs=lhs%q_work_gpu, rhs=rhs%q_work_gpu)
   call assign_allocatable_gpu(lhs=lhs%q_s_gpu,    rhs=rhs%q_s_gpu   )
   endsubroutine eq_assign_eq

   ! non TBP cuf methods
   subroutine advance_q_gpu_cuf(ni, nj, nk, gci, gcj, gck, blocks_number, beta_gpu, dt, q_s_gpu, q_gpu)
   !< Advance q_gpu by means of RK stages.
   integer(I4P), intent(in)            :: ni             !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj             !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk             !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: gci            !< Ghost grid cells number in I direction.
   integer(I4P), intent(in)            :: gcj            !< Ghost grid cells number in J direction.
   integer(I4P), intent(in)            :: gck            !< Ghost grid cells number in K direction.
   integer(I4P), intent(in)            :: blocks_number  !< Number of blocks.
   real(R8P),    intent(in),    device :: beta_gpu(:)    !< RK betaa coefficients.
   real(R8P),    intent(in)            :: Dt             !< Time step.
   real(R8P),    intent(in),    device :: q_s_gpu(1:,    &
                                                  1-gci:,&
                                                  1-gcj:,&
                                                  1-gck:,&
                                                  1:,1:) !< RK stage.
   real(R8P),    intent(inout), device ::   q_gpu(1:,    &
                                                  1-gci:,&
                                                  1-gcj:,&
                                                  1-gck:,&
                                                  1:)    !< Conservative field.
   integer(I4P)                        :: i, j, k, b, s  !< Counter.
   integer(I4P)                        :: iercuda        !< Error trapping flag for CUDAFortran.

   do s=1, 3
      !$cuf kernel do(4) <<<*,*>>>
      do b=1, blocks_number
         do k=1-gck, nk+gck
            do j=1-gcj, nj+gcj
               do i=1-gci, ni+gci
                  q_gpu(b,i,j,k,1) = q_gpu(b,i,j,k,1) + q_s_gpu(b,i,j,k,1,s) * dt * beta_gpu(s)
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   enddo
   endsubroutine advance_q_gpu_cuf

   subroutine compute_residuals_gpu_cuf(ni, nj, nk, gci, gcj, gck, blocks_number, t, dxyz, q_work_gpu, q_gpu)
   !< Compute residuals of equation.
   integer(I4P), intent(in)            :: ni             !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj             !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk             !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: gci            !< Ghost grid cells number in I direction.
   integer(I4P), intent(in)            :: gcj            !< Ghost grid cells number in J direction.
   integer(I4P), intent(in)            :: gck            !< Ghost grid cells number in K direction.
   integer(I4P), intent(in)            :: blocks_number  !< Number of blocks.
   real(R8P),    intent(in)            :: t              !< Time.
   real(R8P),    intent(in),    device :: dxyz(1:,1:)    !< Space steps.
   real(R8P),    intent(inout), device :: q_work_gpu(1:,    &
                                                     1-gci:,&
                                                     1-gcj:,&
                                                     1-gck:,&
                                                     1:) !< Field component to be updated.
   real(R8P),    intent(inout), device :: q_gpu(1:,    &
                                                1-gci:,&
                                                1-gcj:,&
                                                1-gck:,&
                                                1:)      !< Field component to be updated.
   integer(I4P)                        :: b, i, j, k     !< Counter.
   integer(I4P)                        :: iercuda        !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(4) <<<*,*>>>
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               q_work_gpu(b,i,j,k,1) = (q_gpu(b,i+1,j,k,1) - q_gpu(b,i-1,j,k,1)) / (2 * dxyz(1,b))
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
               q_gpu(b,i,j,k,1) = q_work_gpu(b,i,j,k,1)
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_residuals_gpu_cuf

   subroutine compute_rk_stage_gpu_cuf(ni, nj, nk, gci, gcj, gck, blocks_number, alph_gpu, dt, s, q_gpu, q_s_gpu)
   !< Initialize RK stage with q_gpu.
   integer(I4P), intent(in)            :: ni             !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj             !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk             !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: gci            !< Ghost grid cells number in I direction.
   integer(I4P), intent(in)            :: gcj            !< Ghost grid cells number in J direction.
   integer(I4P), intent(in)            :: gck            !< Ghost grid cells number in K direction.
   integer(I4P), intent(in)            :: blocks_number  !< Number of blocks.
   real(R8P),    intent(in),    device :: alph_gpu(:,:)  !< RK alpha coefficients.
   real(R8P),    intent(in)            :: Dt             !< Time step.
   integer(I4P), intent(in)            :: s              !< Stage to initialize.
   real(R8P),    intent(in),    device ::   q_gpu(1:,    &
                                                  1-gci:,&
                                                  1-gcj:,&
                                                  1-gck:,&
                                                  1:)    !< Conservative field.
   real(R8P),    intent(inout), device :: q_s_gpu(1:,    &
                                                  1-gci:,&
                                                  1-gcj:,&
                                                  1-gck:,&
                                                  1:,1:) !< RK stage.
   integer(I4P)                        :: i, j, k, b, ss !< Counter.
   integer(I4P)                        :: iercuda        !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(4) <<<*,*>>>
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               q_s_gpu(b,i,j,k,1,s) = q_gpu(b,i,j,k,1)
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
                  q_s_gpu(b,i,j,k,1,s) = q_s_gpu(b,i,j,k,1,s) + (q_s_gpu(b,i,j,k,1,ss) * (dt * alph_gpu(s, ss)))
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   enddo
   endsubroutine compute_rk_stage_gpu_cuf

endmodule adam_equation_convect1D_gpu_object
