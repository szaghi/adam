!< ADAM, Euler equations system class definition, CPU backend.
module adam_equation_euler_cpu_object
!< ADAM, Euler equations system class definition, CPU backend.
!<
!< Multifluids is modeled by the standard thermodynamic model.

use adam_base_cpu_object
use adam_field_object
use adam_parameters
use PENF
use MPI
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none
private
public :: equation_euler_cpu_object
public :: BC_EXTRAPOLATION
public :: BC_INFLOW

integer(I4P), parameter :: BC_EXTRAPOLATION = 1_I4P
integer(I4P), parameter :: BC_INFLOW        = 2_I4P

type :: equation_euler_cpu_object
   !< Euler equations system class definition, CPU backend.
   !<
   !< Multifluids is modeled by the standard thermodynamic model.
   !<
   !< The conservative varibales are arranged as follows:
   !<```
   !< q(1):    rho(1)
   !< q(2):    rho(2),
   !< ...
   !< q(ns):   rho(ns), specific density of last specie
   !< q(ns+1): rho * u
   !< q(ns+2): rho * v
   !< q(ns+3): rho * w
   !< q(ns+4): rho * E
   !<```
   !< Where `rho(s)` is the specific density of s-th specie, `rho=sum(rho(s))`, `[u,v,w]` is the velocity
   !< vector and `E` it the total specific internal energy. The auxiliary variables array is arranged as follows:
   !<```
   !< q_aux(1): c(1)
   !< q_aux(2): c(2)
   !< ...
   !< q_aux(ns): c(ns), specie concentration of last specie
   !< q_aux(ns+1): rho=sum(rho(s))
   !< q_aux(ns+2): u
   !< q_aux(ns+3): v
   !< q_aux(ns+4): w
   !< q_aux(ns+5): g
   !< q_aux(ns+6): p
   !<```
   !< Where `c` is the species concentration, i.e. rho/rho(s), `p` is the pressure and
   !< `g` is the specific heat ratio of the mixture, i.e.
   !<```
   !< cp = sum(rho(s)/rho * cp(s))
   !< cv = sum(rho(s)/rho * cv(s))
   !< g = cp / cv
   !<```
   type(field_object), pointer :: field=>null() !< The field.
   type(base_cpu_object)       :: base_cpu      !< The base CPU handler.
   ! equation data
   integer(I4P)           :: ns=1_I4P         !< Number of fluid species.
   real(R8P), allocatable :: cp0(:)           !< Specific heat at constant pressure of initial species.
   real(R8P), allocatable :: cv0(:)           !< Specific heat at constant pressure of initial species.
   real(R8P), allocatable :: q_aux(:,:,:,:,:) !< Auxiliary cell centered variables.
   real(R8P)              :: dt=0._R8P        !< Maximum time step accordingly to CFL criterion.
   real(R8P)              :: CFL=0.3_R8P      !< CFL limit.
   ! Runge-Kutta data
   integer(I4P)           :: nrk=3_I4P        !< Runge-Kutta stages number.
   real(R8P), allocatable :: alph(:,:)        !< RK alpha coefficients.
   real(R8P), allocatable :: beta(:)          !< RK beta coefficients.
   real(R8P), allocatable :: gamm(:)          !< RK gamma coefficients.
   real(R8P), allocatable :: q_s(:,:,:,:,:,:) !< RK Field cell centered variables stages.
   ! MPI data, unrelated to field equations
   integer(I4P) :: error=0_I4P  !< Error traping flag.
   integer(I4P) :: myrank=0_I4P !< MPI rank process.
   contains
      ! public methods
      procedure, pass(self) :: compute_aux             !< Compute auxiliary variables.
      procedure, pass(self) :: compute_dt              !< Compute maximum time step.
      procedure, pass(self) :: destroy                 !< Destroy the equation.
      procedure, pass(self) :: initialize              !< Initialize the equation.
      procedure, pass(self) :: mark_by_grad_rho        !< Mark blocks to be refined/derefined by a `grad(rho)` value.
      procedure, pass(self) :: integrate               !< Runge Kutta integration of equation.
      procedure, pass(self) :: compute_residuals       !< Compute residuals of equation.
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions of equation.
      procedure, pass(self) :: set_initial_conditions  !< Set initial conditions of equation.
      procedure, pass(self) :: update_ghost            !< Update ghost cells and set boundary conditions.
      ! operators
      generic :: assignment(=) => eq_assign_eq      !< Overload `=`.
      procedure, pass(lhs), private :: eq_assign_eq !< Operator `=`.
endtype equation_euler_cpu_object

contains
   ! public methods
   subroutine compute_aux(self, q, q_aux)
   !< Update auxiliary variables.
   class(equation_euler_cpu_object), intent(in)    :: self         !< The equation.
   real(R8P),                        intent(in)    :: q(1-self%field%grid%gci:,&
                                                        1-self%field%grid%gcj:,&
                                                        1-self%field%grid%gck:,&
                                                        1:,1:)     !< Conservative variables.
   real(R8P),                        intent(inout) :: q_aux(1-self%field%grid%gci:,&
                                                            1-self%field%grid%gcj:,&
                                                            1-self%field%grid%gck:,&
                                                            1:,1:) !< Auxiliary variables.
   integer(I4P)                                    :: b, i, j, k   !< Counter.

   associate(blocks_number=>self%field%blocks_number,                                      &
             ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk,       &
             gci=>self%field%grid%gci, gcj=>self%field%grid%gcj, gck=>self%field%grid%gck, &
             ns=>self%ns)
      do b=1, size(q, dim=5)
         do k=1-gck, nk+gck
            do j=1-gcj, nj+gcj
               do i=1-gci, ni+gci
                  q_aux(i,j,k,ns+1,b) = sum(q(i,j,k,1:ns,b))
                  q_aux(i,j,k,1:ns,b) = q(i,j,k,1:ns,b) / q_aux(i,j,k,ns+1,b)
                  q_aux(i,j,k,ns+2,b) = q(i,j,k,ns+1,b) / q_aux(i,j,k,ns+1,b)
                  q_aux(i,j,k,ns+3,b) = q(i,j,k,ns+2,b) / q_aux(i,j,k,ns+1,b)
                  q_aux(i,j,k,ns+4,b) = q(i,j,k,ns+3,b) / q_aux(i,j,k,ns+1,b)
                  q_aux(i,j,k,ns+5,b) = dot_product(q_aux(i,j,k,1:ns,b), self%cp0) / dot_product(q_aux(i,j,k,1:ns,b), self%cv0)
                  q_aux(i,j,k,ns+6,b) = (q(i,j,k,ns+4,b) - 0.5_R8P * q_aux(i,j,k,ns+1,b) * (q_aux(i,j,k,ns+2,b)**2 +   &
                                                                                            q_aux(i,j,k,ns+3,b)**2 +   &
                                                                                            q_aux(i,j,k,ns+4,b)**2)) * &
                                        (q_aux(i,j,k,ns+5,b) - 1._R8P)
               enddo
            enddo
         enddo
      enddo
   endassociate
   endsubroutine compute_aux

   subroutine compute_dt(self)
   !< Compute maximum time step accordingly to CFL stabilty criterion.
   class(equation_euler_cpu_object), intent(inout) :: self       !< The equation.
   real(R8P)                                       :: umax       !< Maximum speed of waves propagation.
   integer(I4P)                                    :: b, i, j, k !< Counter.

   associate(blocks_number=>self%field%blocks_number, dxyz=>self%field%dxyz,         &
             ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk, &
             ns=>self%ns, q=>self%field%q, q_aux=>self%q_aux, dt=>self%dt, CFL=>self%CFL)
   call self%compute_aux(q=self%field%q(:,:,:,:,1:blocks_number), q_aux=self%q_aux(:,:,:,:,1:blocks_number))
   dt = huge(1._R8P)
   do b=1, self%field%blocks_number
      umax = 0._R8P
      do k=1, nk
         do j=1, nj
            do i=1, ni
               umax = max(umax, maxval(abs(q_aux(i,j,k,ns+2:ns+4,b))) + &
                                a(p=q_aux(i,j,k,ns+6,b), r=q_aux(i,j,k,ns+1,b), g=q_aux(i,j,k,ns+5,b)))

            enddo
         enddo
      enddo
      dt = min(dt, minval(dxyz(:,b)) / umax * CFL)
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, self%error)
   endassociate
   endsubroutine compute_dt

   subroutine destroy(self)
   !< Destroy the equation.
   class(equation_euler_cpu_object), intent(inout) :: self  !< The equation.
   type(equation_euler_cpu_object)                 :: fresh !< Fresh equation.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, field, ns, nrk, cp0, cv0, CFL)
   !< Initialize the equation.
   class(equation_euler_cpu_object), intent(inout)        :: self   !< The equation.
   type(field_object),               intent(in), target   :: field  !< The field.
   integer(I4P),                     intent(in), optional :: ns     !< Species number.
   integer(I4P),                     intent(in), optional :: nrk    !< Runge-Kutta stages number.
   real(R8P),                        intent(in), optional :: cp0(:) !< Initial specific heats at constant pressure.
   real(R8P),                        intent(in), optional :: cv0(:) !< Initial specific heats at constant volume.
   real(R8P),                        intent(in), optional :: CFL    !< CFL value.

   call self%destroy
   self%field => field
   call self%base_cpu%initialize(field=field)
   if (present(ns)) self%ns = ns
   if (present(nrk)) self%nrk = nrk
   if (self%field%nv - self%ns /= 4) then
      write(stderr, '(A)') 'ADAM-ERROR: field%nv must be euler%ns+4'
      call MPI_FINALIZE(self%error)
      stop
   endif
   if (present(cp0)) then
      self%cp0 = cp0
   else
      allocate(self%cp0(self%ns))
      self%cp0 = 1040._R8P
   endif
   if (present(cv0)) then
      self%cv0 = cv0
   else
      allocate(self%cv0(self%ns))
      self%cv0 = 742.85_R8P
   endif
   if (present(CFL)) self%CFL = CFL

   allocate(self%q_aux(1-field%grid%gci:field%grid%ni+field%grid%gci, &
                       1-field%grid%gcj:field%grid%nj+field%grid%gcj, &
                       1-field%grid%gck:field%grid%nk+field%grid%gck, 1:self%ns+6, 1:field%nb))
   allocate(self%alph(self%nrk,self%nrk), self%beta(self%nrk), self%gamm(self%nrk))
   select case(self%nrk)
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
   allocate(self%q_s(1-field%grid%gci:field%grid%ni+field%grid%gci, &
                     1-field%grid%gcj:field%grid%nj+field%grid%gcj, &
                     1-field%grid%gck:field%grid%nk+field%grid%gck, 1:field%nv, 1:field%nb, 1:self%nrk))
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   endsubroutine initialize

   subroutine mark_by_grad_rho(self, grad_tol, delta_fine, delta_coarse, threshold)
   !< Mark blocks to be refined/derefined by a `grad(rho)` value.
   class(equation_euler_cpu_object), intent(inout)        :: self           !< The equation.
   real(R8P),                        intent(in)           :: grad_tol       !< Gradiend tolerance value.
   real(R8P),                        intent(in)           :: delta_fine     !< Maximum cell delta in fine grids.
   real(R8P),                        intent(in)           :: delta_coarse   !< Minimum cell delta in coarse grids.
   real(R8P),                        intent(in), optional :: threshold      !< Threshold for sphere proximity.
   real(R8P)                                              :: threshold_     !< Threshold for sphere proximity, local var.
   real(R8P)                                              :: max_cell_delta !< Maximum cell delta.
   real(R8P)                                              :: grad_rho       !< Value (max) of gradient of rho.
   integer(I4P)                                           :: b, i, j, k     !< Counter.

   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   self%field%refinements_needed = [(TO_NOT_TOUCH,b=1,self%field%blocks_number)]
   call self%update_ghost(q=self%field%q)
   associate (ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk, &
              blocks_number=>self%field%blocks_number, ns=>self%ns,                   &
              dxyz=>self%field%dxyz, q_aux=>self%q_aux)
   call self%compute_aux(q=self%field%q(:,:,:,:,1:blocks_number), q_aux=self%q_aux(:,:,:,:,1:blocks_number))
   do b=1, blocks_number
      grad_rho = 0._R8P
      do k=1, nk
         do j=1, nj
            do i=1, ni
               grad_rho = max(grad_rho, sqrt(((q_aux(i+1,j,k,ns+1,b) - q_aux(i-1,j,k,ns+1,b))/(2*dxyz(1,b)))**2 + &
                                             ((q_aux(i,j+1,k,ns+1,b) - q_aux(i,j-1,k,ns+1,b))/(2*dxyz(2,b)))**2 + &
                                             ((q_aux(i,j,k+1,ns+1,b) - q_aux(i,j,k-1,ns+1,b))/(2*dxyz(3,b)))**2))

            enddo
         enddo
      enddo

      max_cell_delta = max_cell_delta_grad(grad=grad_rho)

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

      if (grad > grad_tol) then
         delta = delta_fine
      else
         delta = delta_coarse
      endif
      endfunction max_cell_delta_grad
   endsubroutine mark_by_grad_rho

   subroutine integrate(self, t, do_ghost_syncro, residual)
   !< Runge Kutta integration of field.
   class(equation_euler_cpu_object), intent(inout)         :: self             !< The equation.
   real(R8P),                        intent(in)            :: t                !< Time.
   logical,                          intent(in),  optional :: do_ghost_syncro  !< Flag to do syncrous ghost update.
   real(R8P),                        intent(out), optional :: residual         !< Global residual.
   logical                                                 :: do_ghost_syncro_ !< Flag to do syncrous ghost update, local var.
   integer(I4P)                                            :: b, s, ss         !< Counter.

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   associate(alph=>self%alph, beta=>self%beta, gamm=>self%gamm,                      &
             ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk, &
             blocks_number=>self%field%blocks_number,                                &
             inner_blocks_number=>self%field%inner_blocks_number,                    &
             q=>self%field%q, q_s=>self%q_s, dt=>self%dt)
   do s=1, self%nrk
      q_s(1:ni,1:nj,1:nk,:,1:blocks_number,s) = q(1:ni,1:nj,1:nk,:,1:blocks_number)
      do ss=1, s - 1
         q_s(1:ni,1:nj,1:nk,:,1:blocks_number,s) = q_s(1:ni,1:nj,1:nk,:,1:blocks_number,s ) + &
                                                  (q_s(1:ni,1:nj,1:nk,:,1:blocks_number,ss) * (dt * alph(s, ss)))
      enddo
      if (do_ghost_syncro_) then
         call self%update_ghost(q=q_s(:,:,:,:,:,s)) ! all ghosts
         call self%compute_residuals(q=q_s(:,:,:,:,1:,s), t=t + gamm(s) * dt)
      else
         call self%update_ghost(q=q_s(:,:,:,:,:,s), step=1) ! local ghosts
         call self%update_ghost(q=q_s(:,:,:,:,:,s), step=2) ! initialize MPI comms
         call self%compute_residuals(q=q_s(:,:,:,:,1:inner_blocks_number,s), t=t + gamm(s) * dt)
         call self%update_ghost(q=q_s(:,:,:,:,:,s), step=3) ! complete MPI comms
         call self%compute_residuals(q=q_s(:,:,:,:,inner_blocks_number+1:blocks_number,s),t=t+gamm(s)*dt)
      endif
      if (present(residual).and.s==3) then
         residual = 0._R8P
         do b=1, blocks_number
            residual = residual + sum(q_s(1:ni,1:nj,1:nk,:,b,s))/ni/nj/nk
         enddo
         call MPI_ALLREDUCE(MPI_IN_PLACE, residual, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, self%error)
      endif
   enddo
   do s=1, self%nrk
      q(1:ni,1:nj,1:nk,:,1:blocks_number) =   q(1:ni,1:nj,1:nk,:,1:blocks_number) + &
                                            q_s(1:ni,1:nj,1:nk,:,1:blocks_number,s) * dt * beta(s)
   enddo
   endassociate
   endsubroutine integrate

   subroutine compute_residuals(self, q, t)
   !< Compute residuals of equation.
   class(equation_euler_cpu_object), intent(inout) :: self                         !< The equation.
   real(R8P),                        intent(inout) :: q(1-self%field%grid%gci:,&
                                                        1-self%field%grid%gcj:,&
                                                        1-self%field%grid%gck:,&
                                                        1:,1:)                     !< Conservative variables.
   real(R8P),                        intent(in)    :: t                            !< Time.
   real(R8P)                                       :: fluxes_x(0:self%field%grid%ni,&
                                                               1:self%field%grid%nj,&
                                                               1:self%field%grid%nk,&
                                                               1:self%field%nv)    !< Convective fluxes in x direction.
   real(R8P)                                       :: fluxes_y(1:self%field%grid%ni,&
                                                               0:self%field%grid%nj,&
                                                               1:self%field%grid%nk,&
                                                               1:self%field%nv)    !< Convective fluxes in y direction.
   real(R8P)                                       :: fluxes_z(1:self%field%grid%ni,&
                                                               1:self%field%grid%nj,&
                                                               0:self%field%grid%nk,&
                                                               1:self%field%nv)    !< Convective fluxes in z direction.
   integer(I4P)                                    :: b, i, j, k                   !< Counter.

   associate(ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk,       &
             gci=>self%field%grid%gci, gcj=>self%field%grid%gcj, gck=>self%field%grid%gck, &
             dxyz=>self%field%dxyz, ns=>self%ns, q_aux=>self%q_aux)
   call self%compute_aux(q=q, q_aux=self%q_aux(:,:,:,:,1:size(q,dim=5)))
   do b=1, size(q, dim=5)
      ! convective fluxes along x direction
      do k=1, nk
         do j=1, nj
            call compute_fluxes_convective(gc=gci, n=ni,                          &
                                           c         =    q_aux( :,  j,k,1:ns,b), &
                                           rho       =    q_aux( :,  j,k,ns+1,b), &
                                           un        =    q_aux( :,  j,k,ns+2,b), & ! u
                                           ut1       =    q_aux( :,  j,k,ns+3,b), & ! v
                                           ut2       =    q_aux( :,  j,k,ns+4,b), & ! w
                                           g         =    q_aux( :,  j,k,ns+5,b), &
                                           p         =    q_aux( :,  j,k,ns+6,b), &
                                           f_rho     = fluxes_x(0:ni,j,k,1:ns),   &
                                           f_rho_un  = fluxes_x(0:ni,j,k,ns+1),   & ! u
                                           f_rho_ut1 = fluxes_x(0:ni,j,k,ns+2),   & ! v
                                           f_rho_ut2 = fluxes_x(0:ni,j,k,ns+3),   & ! w
                                           f_rho_e   = fluxes_x(0:ni,j,k,ns+4))
         enddo
      enddo
      ! convective fluxes along y direction
      do k=1, nk
         do i=1, ni
            call compute_fluxes_convective(gc=gcj, n=nj,                          &
                                           c         =    q_aux(i, :,  k,1:ns,b), &
                                           rho       =    q_aux(i, :,  k,ns+1,b), &
                                           un        =    q_aux(i, :,  k,ns+3,b), & ! v
                                           ut1       =    q_aux(i, :,  k,ns+2,b), & ! u
                                           ut2       =    q_aux(i, :,  k,ns+4,b), & ! w
                                           g         =    q_aux(i, :,  k,ns+5,b), &
                                           p         =    q_aux(i, :,  k,ns+6,b), &
                                           f_rho     = fluxes_y(i,0:nj,k,1:ns),   &
                                           f_rho_un  = fluxes_y(i,0:nj,k,ns+2),   & ! v
                                           f_rho_ut1 = fluxes_y(i,0:nj,k,ns+1),   & ! u
                                           f_rho_ut2 = fluxes_y(i,0:nj,k,ns+3),   & ! w
                                           f_rho_e   = fluxes_y(i,0:nj,k,ns+4))
         enddo
      enddo
      ! convective fluxes along z direction
      do j=1, nj
         do i=1, ni
            call compute_fluxes_convective(gc=gck, n=nk,                          &
                                           c         =    q_aux(i,j, :,  1:ns,b), &
                                           rho       =    q_aux(i,j, :,  ns+1,b), &
                                           un        =    q_aux(i,j, :,  ns+4,b), & ! w
                                           ut1       =    q_aux(i,j, :,  ns+2,b), & ! u
                                           ut2       =    q_aux(i,j, :,  ns+3,b), & ! v
                                           g         =    q_aux(i,j, :,  ns+5,b), &
                                           p         =    q_aux(i,j, :,  ns+6,b), &
                                           f_rho     = fluxes_z(i,j,0:nk,1:ns),   &
                                           f_rho_un  = fluxes_z(i,j,0:nk,ns+3),   & ! w
                                           f_rho_ut1 = fluxes_z(i,j,0:nk,ns+1),   & ! u
                                           f_rho_ut2 = fluxes_z(i,j,0:nk,ns+2),   & ! v
                                           f_rho_e   = fluxes_z(i,j,0:nk,ns+4))
         enddo
      enddo
      ! residuals
      do k=1, nk
         do j=1, nj
            do i=1, ni
               q(i,j,k,:,b) = (fluxes_x(i-1,j  ,k  ,:) - fluxes_x(i,j,k,:)) / dxyz(1,b) + &
                              (fluxes_y(i  ,j-1,k  ,:) - fluxes_y(i,j,k,:)) / dxyz(2,b) + &
                              (fluxes_z(i  ,j  ,k-1,:) - fluxes_z(i,j,k,:)) / dxyz(3,b)
               ! q(i,j,k,:,b) = (fluxes_x(i-1,j  ,k  ,:) - fluxes_x(i,j,k,:)) / dxyz(1,b)
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine compute_residuals

   subroutine set_boundary_conditions(self, q)
   !< Set boundary conditions of equation.
   class(equation_euler_cpu_object), intent(in)    :: self                            !< The equation.
   real(R8P),                        intent(inout) :: q(1-self%field%grid%gci:,&
                                                        1-self%field%grid%gcj:,&
                                                        1-self%field%grid%gck:,1:,1:) !< Field component to be updated.

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
                     q(i,j,k,:,b) = q(i-ijkdelta(1), j-ijkdelta(2), k-ijkdelta(3), :, b)
                  enddo
               enddo
            enddo
         elseif (bc_type == BC_INFLOW) then
            do k=ijkmin(3), ijkmax(3), sign(1, ijkmax(3)-ijkmin(3))
               do j=ijkmin(2), ijkmax(2), sign(1, ijkmax(2)-ijkmin(2))
                  do i=ijkmin(1), ijkmax(1), sign(1, ijkmax(1)-ijkmin(1))
                  enddo
               enddo
            enddo
         endif
      enddo
      endsubroutine set_bc_fec
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self)
   !< Set initial conditions of field.
   class(equation_euler_cpu_object), intent(inout) :: self    !< The equation.
   integer(I4P)                                    :: b       !< Counter.
   integer(I4P)                                    :: i, j, k !< Counter.

   associate(blocks_number=>self%field%blocks_number,                                      &
             q=>self%field%q,                                                              &
             ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk,       &
             gci=>self%field%grid%gci, gcj=>self%field%grid%gcj, gck=>self%field%grid%gck, &
             x_cell=>self%field%x_cell, y_cell=>self%field%y_cell, z_cell=>self%field%z_cell)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               if (x_cell(i,b)<0.5_R8P) then
                  q(i,j,k,1,b) = 1._R8P
                  q(i,j,k,2,b) = 0._R8P
                  q(i,j,k,3,b) = 0._R8P
                  q(i,j,k,4,b) = 0._R8P
                  q(i,j,k,5,b) = 1._R8P * E(p=1._R8P, r=1._R8P, u=0._R8P, g=self%cp0(1)/self%cv0(1))
               else
                  q(i,j,k,1,b) = 0.125_R8P
                  q(i,j,k,2,b) = 0._R8P
                  q(i,j,k,3,b) = 0._R8P
                  q(i,j,k,4,b) = 0._R8P
                  q(i,j,k,5,b) = 0.125_R8P * E(p=0.1_R8P, r=0.125_R8P, u=0._R8P, g=self%cp0(1)/self%cv0(1))
               endif
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine set_initial_conditions

   subroutine update_ghost(self, q, step)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(equation_euler_cpu_object), intent(inout)        :: self            !< The equation.
   real(R8P),                        intent(inout)        :: q(1-self%field%grid%gci:,&
                                                               1-self%field%grid%gcj:,&
                                                               1-self%field%grid%gck:,&
                                                               1:,1:)        !< Field component to be updated.
   integer(I4P),                     intent(in), optional :: step            !< Step to be perfordmed in asyncronous comp.
   logical                                                :: do_local_update !< Flag for triggering local update.
   logical                                                :: do_set_bc       !< Flag for triggering setting bc.

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
   class(equation_euler_cpu_object), intent(inout) :: lhs !< Left hand side.
   type(equation_euler_cpu_object),  intent(in)    :: rhs !< Right hand side.

   lhs%field => rhs%field
   lhs%base_cpu = rhs%base_cpu
   lhs%ns = rhs%ns
   call assign_allocatable(lhs=lhs%cp0, rhs=rhs%cp0)
   call assign_allocatable(lhs=lhs%cv0, rhs=rhs%cv0)
   call assign_allocatable(lhs=lhs%q_aux, rhs=rhs%q_aux)
   lhs%nrk = rhs%nrk
   call assign_allocatable(lhs=lhs%alph, rhs=rhs%alph)
   call assign_allocatable(lhs=lhs%beta, rhs=rhs%beta)
   call assign_allocatable(lhs=lhs%gamm, rhs=rhs%gamm)
   call assign_allocatable(lhs=lhs%q_s, rhs=rhs%q_s)
   lhs%error = rhs%error
   lhs%myrank = rhs%myrank
   endsubroutine eq_assign_eq

   ! non type-bound procedures
   subroutine compute_fluxes_convective(gc, n,                      &
                                        c, rho, un, ut1, ut2, g, p, &
                                        f_rho, f_rho_un, f_rho_ut1, f_rho_ut2, f_rho_e)
   !< Compute the conservative fluxes on a slice of cells along a coordinate direction.
   integer(I4P), intent(in)    :: gc                  !< Number of ghost cells used.
   integer(I4P), intent(in)    :: n                   !< Number of cells.
   real(R8P),    intent(in)    ::         c(1-gc:,1:) !< Species concentration        [1-gc:n+gc,1:ns].
   real(R8P),    intent(in)    ::       rho(1-gc:)    !< Density                      [1-gc:n+gc].
   real(R8P),    intent(in)    ::        un(1-gc:)    !< Normal velocity              [1-gc:n+gc].
   real(R8P),    intent(in)    ::       ut1(1-gc:)    !< Tangential velocity 1        [1-gc:n+gc].
   real(R8P),    intent(in)    ::       ut2(1-gc:)    !< Tangential velocity 2        [1-gc:n+gc].
   real(R8P),    intent(in)    ::         g(1-gc:)    !< Specific heats ratio         [1-gc:n+gc].
   real(R8P),    intent(in)    ::         p(1-gc:)    !< Pressure                     [1-gc:n+gc].
   real(R8P),    intent(inout) ::     f_rho(0:,   1:) !< Flux of mass                 [0:n,1:ns].
   real(R8P),    intent(inout) ::  f_rho_un(0:)       !< Flux normal momentums        [0:n].
   real(R8P),    intent(inout) :: f_rho_ut1(0:)       !< Flux of tangential1 momentum [0:n].
   real(R8P),    intent(inout) :: f_rho_ut2(0:)       !< Flux of tangential2 momentum [0:n].
   real(R8P),    intent(inout) ::   f_rho_e(0:)       !< Flux energy                  [0:n].
   real(R8P)                   ::   fluxes(3)         !< 1D fluxes.
   ! real(R8P),    allocatable   :: qr(:,:,:)       !< Reconstructed variables.
   integer(I4P)                :: i               !< Counter.

   ! allocate(pr(1:2,0:N+1))
   ! call reconstruct_interfaces_characteristic
   ! print*, 'cazzo rho ',  rho
   ! print*, 'cazzo un  ',  un
   ! print*, 'cazzo ut1 ',  ut1
   ! print*, 'cazzo ut2 ',  ut2
   ! print*, 'cazzo g   ',  g
   ! print*, 'cazzo p   ',  p
   do i=0, n
      ! computing normal fluxes solving Riemann problem
      call solve_riemann(p1=p(i), r1=rho(i), u1=un(i), g1=g(i), p4=p(i+1), r4=rho(i+1), u4=un(i+1), g4=g(i+1), F=fluxes)
      if (fluxes(1)>0._R8P) then
             f_rho(i,:) = fluxes(1) * c(i,:)
          f_rho_un(i)   = fluxes(2)
         f_rho_ut1(i)   = fluxes(1) * ut1(i)
         f_rho_ut2(i)   = fluxes(1) * ut2(i)
           f_rho_e(i)   = fluxes(3) + 0.5_R8P * fluxes(1) * (ut1(i)**2 + ut2(i)**2)
      else
             f_rho(i,:) = fluxes(1) * c(i+1,:)
          f_rho_un(i)   = fluxes(2)
         f_rho_ut1(i)   = fluxes(1) * ut1(i+1)
         f_rho_ut2(i)   = fluxes(1) * ut2(i+1)
           f_rho_e(i)   = fluxes(3) + 0.5_R8P * fluxes(1) * (ut1(i+1)**2 + ut2(i+1)**2)
      endif
   enddo
   endsubroutine compute_fluxes_convective

   subroutine solve_riemann(p1, r1, u1, g1, p4, r4, u4, g4, F)
   !< Solve the Riemann problem between the state $1$ and $4$ using the (local) Lax Friedrichs (Rusanov) solver.
   real(R8P), intent(in)  :: p1      !< Pressure of state 1.
   real(R8P), intent(in)  :: r1      !< Density of state 1.
   real(R8P), intent(in)  :: u1      !< Velocity of state 1.
   real(R8P), intent(in)  :: g1      !< Specific heats ratio of state 1.
   real(R8P), intent(in)  :: p4      !< Pressure of state 4.
   real(R8P), intent(in)  :: r4      !< Density of state 4.
   real(R8P), intent(in)  :: u4      !< Velocity of state 4.
   real(R8P), intent(in)  :: g4      !< Specific heats ratio of state 4.
   real(R8P), intent(out) :: F(1:3)  !< Resulting fluxes.
   real(R8P)              :: lmax    !< Maximum wave speed estimation.
   real(R8P)              :: F1(1:3) !< State 1 fluxes.
   real(R8P)              :: F4(1:3) !< State 4 fluxes.
   real(R8P)              :: u       !< Velocity of the intermediate states.
   real(R8P)              :: p       !< Pressure of the intermediate states.
   real(R8P)              :: S1      !< Maximum wave speed of state 1 and 4.
   real(R8P)              :: S4      !< Maximum wave speed of state 1 and 4.

   ! evaluating the intermediates states 2 and 3 from the known states U1,U4 using the PVRS approximation
   call compute_inter_states
   ! evalutaing the maximum waves speed
   lmax = max(abs(S1), abs(u), abs(S4))
   ! computing the fluxes of state 1 and 4
   F1 = fluxes(p = p1, r = r1, u = u1, g = g1)
   F4 = fluxes(p = p4, r = r4, u = u4, g = g4)
   ! computing the Lax-Friedrichs fluxes approximation
   F(1) = 0.5_R8P*(F1(1) + F4(1) - lmax*(r4                        - r1                       ))
   F(2) = 0.5_R8P*(F1(2) + F4(2) - lmax*(r4*u4                     - r1*u1                    ))
   F(3) = 0.5_R8P*(F1(3) + F4(3) - lmax*(r4*E(p=p4,r=r4,u=u4,g=g4) - r1*E(p=p1,r=r1,u=u1,g=g1)))
   contains
      pure function fluxes(p, r, u, g) result(Fc)
      !< 1D Euler fluxes from primitive variables.
      real(R8P), intent(IN) :: p       !< Pressure.
      real(R8P), intent(IN) :: r       !< Density.
      real(R8P), intent(IN) :: u       !< Velocity.
      real(R8P), intent(IN) :: g       !< Specific heats ratio.
      real(R8P)             :: Fc(1:3) !< State fluxes.

      Fc(1) = r*u
      Fc(2) = Fc(1)*u + p
      Fc(3) = Fc(1)*H(p=p, r=r, u=u, g=g)
      endfunction fluxes

      subroutine compute_inter_states
      !< Compute inter states (23*-states) from state1 and state4.
      real(R8P)              :: a1             !< Speed of sound of state 1.
      real(R8P)              :: a4             !< Speed of sound of state 4.
      real(R8P)              :: ram            !< Mean value of rho*a.
      real(R8P), parameter   :: toll=1e-10_R_P !< Tollerance.

      ! evaluation of the intermediate states pressure and velocity
      a1  = sqrt(g1 * p1 / r1)                              ! left speed of sound
      a4  = sqrt(g4 * p4 / r4)                              ! right speed of sound
      ram = 0.5_R8P * (r1 + r4) * 0.5_R8P * (a1 + a4)       ! product of mean density for mean speed of sound
      u   = 0.5_R8P * (u1 + u4) - 0.5_R8P * (p4 - p1) / ram ! evaluation of the contact wave speed (velocity of intermediate states)
      p   = 0.5_R8P * (p1 + p4) - 0.5_R8P * (u4 - u1) * ram ! evaluation of the pressure of the intermediate states
      ! evaluation of the left wave speeds
      if (p<=p1*(1._R8P + toll)) then
        ! rarefaction
        S1 = u1 - a1
      else
        ! shock
        S1 = u1 - a1 * sqrt(1._R8P + (g1 + 1._R8P) / (2._R8P * g1) * (p / p1 - 1._R8P))
      endif
      ! evaluation of the right wave speeds
      if (p<=p4 * (1._R8P + toll)) then
        ! rarefaction
        S4 = u4 + a4
      else
        ! shock
        S4 = u4 + a4 * sqrt(1._R8P + (g4 + 1._R8P) / (2._R8P * g4) * ( p / p4 - 1._R8P))
      endif
      endsubroutine compute_inter_states
   endsubroutine solve_riemann

   elemental function p(r, a, g) result(pressure)
   !< Return pressure for an ideal calorically perfect gas.
   real(R8P), intent(in) :: r        !< Density.
   real(R8P), intent(in) :: a        !< Speed of sound.
   real(R8P), intent(in) :: g        !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
   real(R8P)             :: pressure !< Pressure.

   pressure = r*a*a/g
   endfunction p

   elemental function r(p, a, g) result(density)
   !< Return density for an ideal calorically perfect gas.
   real(R8P), intent(in) :: p       !< Pressure.
   real(R8P), intent(in) :: a       !< Speed of sound.
   real(R8P), intent(in) :: g       !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
   real(R8P)             :: density !< Density.

   density = g*p/(a*a)
   endfunction r

   elemental function a(p, r, g) result(ss)
   !< Return speed of sound for an ideal calorically perfect gas.
   real(R8P), intent(in) :: p  !< Pressure.
   real(R8P), intent(in) :: r  !< Density.
   real(R8P), intent(in) :: g  !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
   real(R8P)             :: ss !< Speed of sound.

   ss = sqrt(g*p/r)
   endfunction a

   elemental function E(p, r, u, g) result(energy)
   !< Return total specific energy (per unit of mass).
   !<$$
   !<  E = \frac{p}{{\left( {\g  - 1} \right)\r }} + \frac{{u^2 }}{2}
   !<$$
   real(R8P), intent(in) :: p      !< Pressure.
   real(R8P), intent(in) :: r      !< Density.
   real(R8P), intent(in) :: u      !< Module of velocity vector.
   real(R8P), intent(in) :: g      !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
   real(R8P)             :: energy !< Total specific energy (per unit of mass).

   energy = p/((g - 1._R8P) * r) + 0.5_R8P * u*u
   endfunction E

   elemental function H(p, r, u, g) result(entalpy)
   !< Return total specific entalpy (per unit of mass).
   !<$$
   !<  H = \frac{{\g p}}{{\left( {\g  - 1} \right)\r }} + \frac{{u^2 }}{2}
   !<$$
   real(R_P), intent(in) :: g       !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
   real(R_P), intent(in) :: p       !< Pressure.
   real(R_P), intent(in) :: r       !< Density.
   real(R8P), intent(in) :: u       !< Module of velocity vector.
   real(R_P)             :: entalpy !< Total specific entalpy (per unit of mass).

   entalpy = g * p / ((g - 1._R_P) * r) + 0.5_R_P * u*u
   endfunction H
endmodule adam_equation_euler_cpu_object
