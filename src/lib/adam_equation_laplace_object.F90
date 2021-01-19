!< ADAM, Laplace equation class definition.
module adam_equation_laplace_object
!< ADAM, Laplace equation class definition.

use adam_field_object
use adam_grid_object
use adam_parameters
use PENF
#ifdef _MPI_
use MPI
#endif

implicit none
private
public :: equation_laplace_object

type :: equation_laplace_object
   !< Laplace equation class definition.
   ! MPI data, unrelated to field equations
   integer(I4P) :: error=0_I4P  !< Error traping flag.
   integer(I4P) :: myrank=0_I4P !< MPI rank process.
   ! RK data, related to field equations
   integer(I4P) :: ns=3_I4P                                            !< Stages number.
   real(R8P)    :: alph(3,3) = reshape([0._R8P, 1._R8P, 0.25_R8P, &
                                        0._R8P, 0._R8P, 0.25_R8P, &
                                        0._R8P, 0._R8P,0._R8P], [3,3]) !< RK alpha coefficients.
   real(R8P)    :: beta(3) = [1._R8P/6._R8P, &
                              1._R8P/6._R8P, &
                              2._R8P/3._R8P]                           !< RK beta coefficients.
   real(R8P)    :: gamm(3) = [0._R8P, &
                              1._R8P, &
                              0._R8P]                                  !< RK gamma coefficients.
   real(R8P), allocatable :: q_s(:,:,:,:,:,:)                          !< RK Field cell centered variables stages.
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
endtype equation_laplace_object

contains
   ! public methods
   subroutine destroy(self)
   !< Destroy the equation.
   class(equation_laplace_object), intent(inout) :: self  !< The equation.
   type(equation_laplace_object)                 :: fresh !< Fresh equation.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, field)
   !< Initialize the equation.
   class(equation_laplace_object), intent(inout) :: self  !< The equation.
   type(field_object),             intent(inout) :: field !< The field.

   call self%destroy
   allocate(self%q_s(1-field%grid%gci:field%grid%ni+field%grid%gci, &
                     1-field%grid%gcj:field%grid%nj+field%grid%gcj, &
                     1-field%grid%gck:field%grid%nk+field%grid%gck, 1:field%nv, 1:field%nb, 1:self%ns))
#ifdef _MPI_
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
#endif
   endsubroutine initialize

   subroutine mark_by_grad_q(self, field, threshold)
   !< Mark blocks to be refined/derefined by a `grad(q)` value.
   class(equation_laplace_object), intent(in)           :: self           !< The equation.
   type(field_object),             intent(inout)        :: field          !< The field.
   real(R8P),                      intent(in), optional :: threshold      !< Threshold for sphere proximity.
   real(R8P)                                            :: threshold_     !< Threshold for sphere proximity, local var.
   real(R8P)                                            :: max_cell_delta !< Maximum cell delta.
   real(R8P)                                            :: grad_q         !< Value (max) of gradient of q.
   integer(I4P)                                         :: b, i, j, k     !< Counter.

   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   if (allocated(field%refinements_needed)) deallocate(field%refinements_needed)
   allocate(field%refinements_needed(field%blocks_number))
   call field%update_ghost(q=field%q)
   call self%set_boundary_conditions(field=field, q=field%q)
   do b=1, field%blocks_number
      grad_q = 0._R8P
      associate (ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk)
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  grad_q = max(grad_q, sqrt(((field%q(i+1,j,k,1,b) - field%q(i-1,j,k,1,b))/(2*field%dxyz(1,b)))**2 + &
                                            ((field%q(i,j+1,k,1,b) - field%q(i,j-1,k,1,b))/(2*field%dxyz(2,b)))**2 + &
                                            ((field%q(i,j,k+1,1,b) - field%q(i,j,k-1,1,b))/(2*field%dxyz(3,b)))**2))

               enddo
            enddo
         enddo
      endassociate

      max_cell_delta = max_cell_delta_grad(grad=grad_q)

      if (maxval(field%dxyz(:,b)) > max_cell_delta) then
         field%refinements_needed(b) = TO_BE_REFINED
      elseif (maxval(field%dxyz(:,b)) * threshold_ < max_cell_delta) then
         field%refinements_needed(b) = TO_BE_DEREFINED
      else
         field%refinements_needed(b) = TO_NOT_TOUCH
      endif
   enddo
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

   subroutine integrate(self, field, t, Dt, do_ghost_syncro, residual)
   !< Runge Kutta integration of field.
   class(equation_laplace_object), intent(inout)         :: self             !< The equation.
   type(field_object),             intent(inout)         :: field            !< The field.
   real(R8P),                      intent(in)            :: t                !< Time.
   real(R8P),                      intent(in)            :: Dt               !< Time step.
   logical,                        intent(in),  optional :: do_ghost_syncro  !< Flag to do syncrous ghost update.
   real(R8P),                      intent(out), optional :: residual         !< Global residual.
   logical                                               :: do_ghost_syncro_ !< Flag to do syncrous ghost update, local var.
   integer(I4P)                                          :: b, s, ss         !< Counter.

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   associate(alph=>self%alph, beta=>self%beta, gamm=>self%gamm,       &
             ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, &
             blocks_number=>field%blocks_number,                      &
             inner_blocks_number=>field%inner_blocks_number,          &
             q=>field%q, q_s=>self%q_s)
   do s=1, 3
      q_s(1:ni,1:nj,1:nk,1,1:blocks_number,s) = q(1:ni,1:nj,1:nk,1,1:blocks_number)
      do ss=1, s - 1
         q_s(1:ni,1:nj,1:nk,1,1:blocks_number,s) = q_s(1:ni,1:nj,1:nk,1,1:blocks_number,s ) + &
                                                  (q_s(1:ni,1:nj,1:nk,1,1:blocks_number,ss) * (Dt * alph(s, ss)))
      enddo
      if (do_ghost_syncro_) then
         call field%update_ghost(q=q_s(:,:,:,:,:,s))
         call self%set_boundary_conditions(field=field, q=q_s(:,:,:,:,:,s))
         call self%compute_residuals(field=field, q=q_s(:,:,:,:,:,s), t=t + gamm(s) * Dt, block_start=1, block_end=blocks_number)
      else
         call field%update_ghost(q=q_s(:,:,:,:,:,s), step=1)
         call field%update_ghost(q=q_s(:,:,:,:,:,s), step=2)
         call self%compute_residuals(field=field,q=q_s(:,:,:,:,:,s), t=t + gamm(s) * Dt, &
                                     block_start=1, block_end=inner_blocks_number)
         call field%update_ghost(q=q_s(:,:,:,:,:,s), step=3)
         call self%set_boundary_conditions(field=field, q=q_s(:,:,:,:,:,s))
         call self%compute_residuals(field=field, q=q_s(:,:,:,:,:,s), t=t + gamm(s) * Dt, &
                                     block_start=inner_blocks_number+1, block_end=blocks_number)
      endif
      if (present(residual).and.s==3) then
         residual = 0._R8P
         do b=1, blocks_number
            residual = residual + sum(q_s(1:ni,1:nj,1:nk,1,b,s))/ni/nj/nk
         enddo
#ifdef _MPI_
         call MPI_ALLREDUCE(MPI_IN_PLACE, residual, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, self%error)
#endif
      endif
   enddo
   do s=1, 3
      q(1:ni,1:nj,1:nk,1,1:blocks_number) =   q(1:ni,1:nj,1:nk,1,1:blocks_number) + &
                                            q_s(1:ni,1:nj,1:nk,1,1:blocks_number,s) * Dt * beta(s)
   enddo
   endassociate
   endsubroutine integrate

   subroutine compute_residuals(self, field, q, t, block_start, block_end)
   !< Compute residuals of equation.
   class(equation_laplace_object), intent(in)    :: self                       !< The equation.
   type(field_object),             intent(in)    :: field                      !< The field.
   real(R8P),                      intent(inout) :: q(1-field%grid%gci:,&
                                                      1-field%grid%gcj:,&
                                                      1-field%grid%gck:,1:,1:) !< Field component to be updated.
   real(R8P),                      intent(in)    :: t                          !< Time.
   integer(I4P),                   intent(in)    :: block_start                !< Index of block to start residuals computation.
   integer(I4P),                   intent(in)    :: block_end                  !< Index of block to end   residuals computation.
   real(R8P)                                     :: q_work(1:field%grid%ni,&
                                                           1:field%grid%nj,&
                                                           1:field%grid%nk)    !< Field component to be updated, working buffer.
   integer(I4P)                                  :: b, i, j, k                 !< Counter.

   do b=block_start, block_end
      do k=1, field%grid%nk
         do j=1, field%grid%nj
            do i=1, field%grid%ni
               ! q_work(i,j,k) = (q(i+1,j,  k,  b) + q(i-1,j,  k,  b) - 2 * q(i,j,k,b)) / field%dxyz(1,b)**2 + &
               !                 (q(i,  j+1,k,  b) + q(i,  j-1,k,  b) - 2 * q(i,j,k,b)) / field%dxyz(2,b)**2 + &
               !                 (q(i,  j,  k+1,b) + q(i,  j,  k-1,b) - 2 * q(i,j,k,b)) / field%dxyz(3,b)**2
               q_work(i,j,k) = (q(i+1,j,k,1,b) - q(i-1,j,k,1,b))/(2*field%dxyz(1,b))
            enddo
         enddo
      enddo
      q(1:field%grid%ni,1:field%grid%nj,1:field%grid%nk,1,b) = q_work(1:field%grid%ni,1:field%grid%nj,1:field%grid%nk)
   enddo
   endsubroutine compute_residuals

   subroutine set_boundary_conditions(self, field, q)
   !< Set boundary conditions of equation.
   class(equation_laplace_object), intent(in)    :: self                       !< The equation.
   type(field_object),             intent(in)    :: field                      !< The field.
   real(R8P),                      intent(inout) :: q(1-field%grid%gci:,&
                                                      1-field%grid%gcj:,&
                                                      1-field%grid%gck:,1:,1:) !< Field component to be updated.

   if (allocated(field%local_map_bc_face))   call set_bc_fec(local_map_bc=field%local_map_bc_face)
   if (allocated(field%local_map_bc_edge))   call set_bc_fec(local_map_bc=field%local_map_bc_edge)
   if (allocated(field%local_map_bc_corner)) call set_bc_fec(local_map_bc=field%local_map_bc_corner)
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
                     q(i,j,k,1,b) = exp(-((field%y_cell(j,b) - 0.5)**2/(2 * 0.2**2)+&
                                          (field%z_cell(k,b) - 0.5)**2/(2 * 0.2**2)))
                  enddo
               enddo
            enddo
         endif
      enddo
      endsubroutine set_bc_fec
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self, field)
   !< Set initial conditions of field.
   class(equation_laplace_object), intent(in)    :: self    !< The equation.
   type(field_object),             intent(inout) :: field   !< The field.
   integer(I4P)                                  :: b       !< Counter.
   integer(I4P)                                  :: i, j, k !< Counter.
   real(R8P)                                     :: a       !< Gaussian amplitude.
   real(R8P)                                     :: sigma_x !< Gaussian x variance.
   real(R8P)                                     :: sigma_y !< Gaussian y variance.
   real(R8P)                                     :: sigma_z !< Gaussian z variance.
   real(R8P)                                     :: x_0     !< Gaussian x center.
   real(R8P)                                     :: y_0     !< Gaussian y center.
   real(R8P)                                     :: z_0     !< Gaussian z center.

   a = 1.0_R8P
   x_0 = (field%grid%domain_emax(1) - field%grid%domain_emin(1)) / 5.0_R8P
   y_0 = (field%grid%domain_emax(2) - field%grid%domain_emin(2)) / 2.0_R8P
   z_0 = (field%grid%domain_emax(3) - field%grid%domain_emin(3)) / 2.0_R8P
   sigma_x = 0.05_R8P
   sigma_y = 0.05_R8P
   sigma_z = 0.05_R8P
   associate(blocks_number=>field%blocks_number,                            &
             q=>field%q,                                                    &
             ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk,       &
             gci=>field%grid%gci, gcj=>field%grid%gcj, gck=>field%grid%gck, &
             x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell)
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

   subroutine update_ghost(self, field, q)
   !< Update ghost cells and set boundary conditions.
   class(equation_laplace_object), intent(inout) :: self                       !< The equation.
   type(field_object),             intent(inout) :: field                      !< The field.
   real(R8P),                      intent(inout) :: q(1-field%grid%gci:,&
                                                      1-field%grid%gcj:,&
                                                      1-field%grid%gck:,1:,1:) !< Field component to be updated.

   if (field%blocks_number > 1) call field%update_ghost(q=q)
   call self%set_boundary_conditions(field=field, q=q)
   endsubroutine update_ghost

   ! operators
   ! =
   subroutine eq_assign_eq(lhs, rhs)
   !< Operator `=`.
   class(equation_laplace_object), intent(inout) :: lhs !< Left hand side.
   type(equation_laplace_object),  intent(in)    :: rhs !< Right hand side.

   lhs%error  = rhs%error
   lhs%myrank = rhs%myrank
   lhs%ns     = rhs%ns
   lhs%alph   = rhs%alph
   lhs%beta   = rhs%beta
   lhs%gamm   = rhs%gamm
   call assign_allocatable(lhs=lhs%q_s, rhs=rhs%q_s)
   endsubroutine eq_assign_eq
endmodule adam_equation_laplace_object
