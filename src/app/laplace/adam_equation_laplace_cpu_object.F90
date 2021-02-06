!< ADAM, Laplace equation class definition, CPU backend.
module adam_equation_laplace_cpu_object
!< ADAM, Laplace equation class definition, CPU backend.

use adam_base_cpu_object
use adam_field_object
use adam_parameters
use PENF
use MPI

implicit none
private
public :: equation_laplace_cpu_object
public :: BC_EXTRAPOLATION
public :: BC_INFLOW

integer(I4P), parameter :: BC_EXTRAPOLATION = 1_I4P
integer(I4P), parameter :: BC_INFLOW        = 2_I4P

type :: equation_laplace_cpu_object
   !< Laplace equation class definition, CPU backend.
   type(field_object), pointer :: field=>null() !< The field.
   type(base_cpu_object)       :: base_cpu      !< The base CPU handler.
   ! MPI data, unrelated to field equations
   integer(I4P) :: error=0_I4P  !< Error traping flag.
   integer(I4P) :: myrank=0_I4P !< MPI rank process.
   ! RK data, related to field equations
   integer(I4P) :: ns=3_I4P                   !< Runge-Kutta stages number.
   real(R8P), allocatable :: alph(:,:)        !< RK alpha coefficients.
   real(R8P), allocatable :: beta(:)          !< RK beta coefficients.
   real(R8P), allocatable :: gamm(:)          !< RK gamma coefficients.
   real(R8P), allocatable :: q_s(:,:,:,:,:,:) !< RK Field cell centered variables stages.
   contains
      ! public methods
      procedure, pass(self) :: destroy                 !< Destroy the equation.
      procedure, pass(self) :: initialize              !< Initialize the equation.
      procedure, pass(self) :: mark_by_grad_q          !< Mark blocks to be refined/derefined by a `grad(q)` value.
      procedure, pass(self) :: integrate               !< Runge Kutta integration of equation.
      procedure, pass(self) :: compute_residuals       !< Compute residuals of equation.
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions of equation.
      procedure, pass(self) :: set_initial_conditions  !< Set initial conditions of equation.
      procedure, pass(self) :: update_ghost            !< Update ghost cells and set boundary conditions.
      ! operators
      generic :: assignment(=) => eq_assign_eq      !< Overload `=`.
      procedure, pass(lhs), private :: eq_assign_eq !< Operator `=`.
endtype equation_laplace_cpu_object

contains
   ! public methods
   subroutine destroy(self)
   !< Destroy the equation.
   class(equation_laplace_cpu_object), intent(inout) :: self  !< The equation.
   type(equation_laplace_cpu_object)                 :: fresh !< Fresh equation.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, field, ns)
   !< Initialize the equation.
   class(equation_laplace_cpu_object), intent(inout)        :: self  !< The equation.
   type(field_object),                 intent(in), target   :: field !< The field.
   integer(I4P),                       intent(in), optional :: ns    !< Runge-Kutta stages number.

   call self%destroy
   self%field => field
   call self%base_cpu%initialize(field=field)
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
   allocate(self%q_s(1-field%grid%ngc:field%grid%ni+field%grid%ngc, &
                     1-field%grid%ngc:field%grid%nj+field%grid%ngc, &
                     1-field%grid%ngc:field%grid%nk+field%grid%ngc, 1:field%nv, 1:field%nb, 1:self%ns))
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   endsubroutine initialize

   subroutine mark_by_grad_q(self, threshold)
   !< Mark blocks to be refined/derefined by a `grad(q)` value.
   class(equation_laplace_cpu_object), intent(inout)        :: self           !< The equation.
   real(R8P),                          intent(in), optional :: threshold      !< Threshold for sphere proximity.
   real(R8P)                                                :: threshold_     !< Threshold for sphere proximity, local var.
   real(R8P)                                                :: max_cell_delta !< Maximum cell delta.
   real(R8P)                                                :: grad_q         !< Value (max) of gradient of q.
   integer(I4P)                                             :: b, i, j, k     !< Counter.

   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   self%field%refinements_needed = [(TO_NOT_TOUCH,b=1,self%field%blocks_number)]
   call self%update_ghost(q=self%field%q)
   associate (ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk, q=>self%field%q, dxyz=>self%field%dxyz)
   do b=1, self%field%blocks_number
      grad_q = 0._R8P
      do k=1, nk
         do j=1, nj
            do i=1, ni
               grad_q = max(grad_q, sqrt(((q(i+1,j,k,1,b) - q(i-1,j,k,1,b))/(2*dxyz(1,b)))**2 + &
                                         ((q(i,j+1,k,1,b) - q(i,j-1,k,1,b))/(2*dxyz(2,b)))**2 + &
                                         ((q(i,j,k+1,1,b) - q(i,j,k-1,1,b))/(2*dxyz(3,b)))**2))

            enddo
         enddo
      enddo

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
      function max_cell_delta_grad(grad) result(delta)
      !< Return the maximum cell delta given a gradient tollerance.
      real(R8P), intent(in) :: grad  !< Gradient value.
      real(R8P)             :: delta !< Maximum cell delta admissible.

      if (grad > 9.2_R8P) then
         delta = 0.004_R8P
      else
         delta = 0.08_R8P
      endif
      endfunction max_cell_delta_grad
   endsubroutine mark_by_grad_q

   subroutine integrate(self, t, Dt, do_ghost_syncro, residual)
   !< Runge Kutta integration of field.
   class(equation_laplace_cpu_object), intent(inout)         :: self             !< The equation.
   real(R8P),                          intent(in)            :: t                !< Time.
   real(R8P),                          intent(in)            :: Dt               !< Time step.
   logical,                            intent(in),  optional :: do_ghost_syncro  !< Flag to do syncrous ghost update.
   real(R8P),                          intent(out), optional :: residual         !< Global residual.
   logical                                                   :: do_ghost_syncro_ !< Flag to do syncrous ghost update, local var.
   integer(I4P)                                              :: b, s, ss         !< Counter.

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   associate(alph=>self%alph, beta=>self%beta, gamm=>self%gamm,                      &
             ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk, &
             blocks_number=>self%field%blocks_number,                                &
             inner_blocks_number=>self%field%inner_blocks_number,                    &
             q=>self%field%q, q_s=>self%q_s)
   do s=1, self%ns
      q_s(1:ni,1:nj,1:nk,1,1:blocks_number,s) = q(1:ni,1:nj,1:nk,1,1:blocks_number)
      do ss=1, s - 1
         q_s(1:ni,1:nj,1:nk,1,1:blocks_number,s) = q_s(1:ni,1:nj,1:nk,1,1:blocks_number,s ) + &
                                                  (q_s(1:ni,1:nj,1:nk,1,1:blocks_number,ss) * (Dt * alph(s, ss)))
      enddo
      if (do_ghost_syncro_) then
         call self%update_ghost(q=q_s(:,:,:,:,:,s)) ! all ghosts
         call self%compute_residuals(q=q_s(:,:,:,:,:,s), t=t + gamm(s) * Dt, block_start=1, block_end=blocks_number)
      else
         call self%update_ghost(q=q_s(:,:,:,:,:,s), step=1) ! local ghosts
         call self%update_ghost(q=q_s(:,:,:,:,:,s), step=2) ! initialize MPI comms
         call self%compute_residuals(q=q_s(:,:,:,:,:,s), t=t + gamm(s) * Dt, block_start=1, block_end=inner_blocks_number)
         call self%update_ghost(q=q_s(:,:,:,:,:,s), step=3) ! complete MPI comms
         call self%compute_residuals(q=q_s(:,:,:,:,:,s),t=t+gamm(s)*Dt,block_start=inner_blocks_number+1,block_end=blocks_number)
      endif
      if (present(residual).and.s==3) then
         residual = 0._R8P
         do b=1, blocks_number
            residual = residual + sum(q_s(1:ni,1:nj,1:nk,1,b,s))/ni/nj/nk
         enddo
         call MPI_ALLREDUCE(MPI_IN_PLACE, residual, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, self%error)
      endif
   enddo
   do s=1, self%ns
      q(1:ni,1:nj,1:nk,1,1:blocks_number) =   q(1:ni,1:nj,1:nk,1,1:blocks_number) + &
                                            q_s(1:ni,1:nj,1:nk,1,1:blocks_number,s) * Dt * beta(s)
   enddo
   endassociate
   endsubroutine integrate

   subroutine compute_residuals(self, q, t, block_start, block_end)
   !< Compute residuals of equation.
   class(equation_laplace_cpu_object), intent(in)    :: self                         !< The equation.
   real(R8P),                          intent(inout) :: q(1-self%field%grid%ngc:,&
                                                          1-self%field%grid%ngc:,&
                                                          1-self%field%grid%ngc:,&
                                                          1:,1:)                     !< Field component to be updated.
   real(R8P),                          intent(in)    :: t                            !< Time.
   integer(I4P),                       intent(in)    :: block_start                  !< Index of block to start residuals comp.
   integer(I4P),                       intent(in)    :: block_end                    !< Index of block to end   residuals comp.
   real(R8P)                                         :: q_work(1:self%field%grid%ni,&
                                                               1:self%field%grid%nj,&
                                                               1:self%field%grid%nk) !< Field component to be updated, working buffer.
   integer(I4P)                                      :: b, i, j, k                   !< Counter.

   associate(ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk, dxyz=>self%field%dxyz)
   do b=block_start, block_end
      do k=1, nk
         do j=1, nj
            do i=1, ni
               q_work(i,j,k) = (q(i+1,j,  k,  1,b) + q(i-1,j,  k,  1,b) - 2 * q(i,j,k,1,b)) / dxyz(1,b)**2 + &
                               (q(i,  j+1,k,  1,b) + q(i,  j-1,k,  1,b) - 2 * q(i,j,k,1,b)) / dxyz(2,b)**2 + &
                               (q(i,  j,  k+1,1,b) + q(i,  j,  k-1,1,b) - 2 * q(i,j,k,1,b)) / dxyz(3,b)**2
            enddo
         enddo
      enddo
      q(1:ni,1:nj,1:nk,1,b) = q_work(1:ni,1:nj,1:nk)
   enddo
   endassociate
   endsubroutine compute_residuals

   subroutine set_boundary_conditions(self, q)
   !< Set boundary conditions of equation.
   class(equation_laplace_cpu_object), intent(in)    :: self                            !< The equation.
   real(R8P),                          intent(inout) :: q(1-self%field%grid%ngc:,&
                                                          1-self%field%grid%ngc:,&
                                                          1-self%field%grid%ngc:,1:,1:) !< Field component to be updated.

   if (allocated(self%field%local_map_bc_face))   call set_bc_fec(local_map_bc=self%field%local_map_bc_face)
   if (allocated(self%field%local_map_bc_edge))   call set_bc_fec(local_map_bc=self%field%local_map_bc_edge)
   if (allocated(self%field%local_map_bc_corner)) call set_bc_fec(local_map_bc=self%field%local_map_bc_corner)
   contains
      subroutine set_bc_fec(local_map_bc)
      integer(I8P), intent(in) :: local_map_bc(:,:) !< Local map for BC ghost cells.
      integer(I4P)             :: b                 !< Counter.
      integer(I4P)             :: f, i, j, k        !< Counter.
      integer(I4P)             :: fec               !< Counter.
      integer(I4P)             :: ijkmin(3)         !< Lower limit of ijk indexes.
      integer(I4P)             :: ijkmax(3)         !< Upper limit of ijk indexes.
      integer(I4P)             :: ijkdelta(3)       !< IJK delta step for extrapolation.
      integer(I4P)             :: bc_type           !< Boundary condition type.

      do f=1, size(local_map_bc, dim=1)
         b        = local_map_bc(f, 1)
         fec      = local_map_bc(f, 2)
         ijkmin   = local_map_bc(f, 3:5)
         ijkmax   = local_map_bc(f, 6:8)
         ijkdelta = local_map_bc(f, 9:11)
         bc_type  = local_map_bc(f, 12)
         if (bc_type == BC_EXTRAPOLATION) then
            do k=ijkmin(3), ijkmax(3), sign(1, ijkmax(3)-ijkmin(3))
               do j=ijkmin(2), ijkmax(2), sign(1, ijkmax(2)-ijkmin(2))
                  do i=ijkmin(1), ijkmax(1), sign(1, ijkmax(1)-ijkmin(1))
                     q(i,j,k,1,b) = q(i-ijkdelta(1), j-ijkdelta(2), k-ijkdelta(3), 1, b)
                  enddo
               enddo
            enddo
         elseif (bc_type == BC_INFLOW) then
            do k=ijkmin(3), ijkmax(3), sign(1, ijkmax(3)-ijkmin(3))
               do j=ijkmin(2), ijkmax(2), sign(1, ijkmax(2)-ijkmin(2))
                  do i=ijkmin(1), ijkmax(1), sign(1, ijkmax(1)-ijkmin(1))
                     ! q(i,j,k,b) = 1._R8P
                     q(i,j,k,1,b) = exp(-((self%field%y_cell(j,b) - 0.5)**2/(2 * 0.2**2)+&
                                          (self%field%z_cell(k,b) - 0.5)**2/(2 * 0.2**2)))
                  enddo
               enddo
            enddo
         endif
      enddo
      endsubroutine set_bc_fec
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self)
   !< Set initial conditions of field.
   class(equation_laplace_cpu_object), intent(inout) :: self    !< The equation.
   integer(I4P)                                      :: b       !< Counter.
   integer(I4P)                                      :: i, j, k !< Counter.
   real(R8P)                                         :: a       !< Gaussian amplitude.
   real(R8P)                                         :: sigma_x !< Gaussian x variance.
   real(R8P)                                         :: sigma_y !< Gaussian y variance.
   real(R8P)                                         :: sigma_z !< Gaussian z variance.
   real(R8P)                                         :: x_0     !< Gaussian x center.
   real(R8P)                                         :: y_0     !< Gaussian y center.
   real(R8P)                                         :: z_0     !< Gaussian z center.

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
             ngc=>self%field%grid%ngc, x_cell=>self%field%x_cell, y_cell=>self%field%y_cell, z_cell=>self%field%z_cell)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               q(i,j,k,1,b) = a * exp(-((x_cell(i,b) - x_0)**2/(2 * sigma_x**2)+&
                                        (y_cell(j,b) - y_0)**2/(2 * sigma_y**2)+&
                                        (z_cell(k,b) - z_0)**2/(2 * sigma_z**2)))
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine set_initial_conditions

   subroutine update_ghost(self, q, step)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(equation_laplace_cpu_object), intent(inout)        :: self            !< The equation.
   real(R8P),                          intent(inout)        :: q(1-self%field%grid%ngc:,&
                                                                 1-self%field%grid%ngc:,&
                                                                 1-self%field%grid%ngc:,&
                                                                 1:,1:)        !< Field component to be updated.
   integer(I4P),                       intent(in), optional :: step            !< Step to be perfordmed in asyncronous comp.
   logical                                                  :: do_local_update !< Flag for triggering local update.
   logical                                                  :: do_set_bc       !< Flag for triggering setting bc.

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

   if (do_local_update) call self%base_cpu%update_ghost_local(q=q)
                        call self%base_cpu%update_ghost_mpi(q=q, step=step)
   if (do_set_bc)       call self%set_boundary_conditions(q=q)
   endsubroutine update_ghost

   ! operators
   ! =
   subroutine eq_assign_eq(lhs, rhs)
   !< Operator `=`.
   class(equation_laplace_cpu_object), intent(inout) :: lhs !< Left hand side.
   type(equation_laplace_cpu_object),  intent(in)    :: rhs !< Right hand side.

   lhs%field => rhs%field
   lhs%base_cpu = rhs%base_cpu
   lhs%error  = rhs%error
   lhs%myrank = rhs%myrank
   lhs%ns     = rhs%ns
   call assign_allocatable(lhs=lhs%alph, rhs=rhs%alph)
   call assign_allocatable(lhs=lhs%beta, rhs=rhs%beta)
   call assign_allocatable(lhs=lhs%gamm, rhs=rhs%gamm)
   call assign_allocatable(lhs=lhs%q_s, rhs=rhs%q_s)
   endsubroutine eq_assign_eq
endmodule adam_equation_laplace_cpu_object
