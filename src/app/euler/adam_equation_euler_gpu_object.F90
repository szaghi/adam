!< ADAM, Euler equations system class definition, GPU backend.
module adam_equation_euler_gpu_object
!< ADAM, Euler equations system class definition, GPU backend.

use adam_base_gpu_object
use adam_field_object
use adam_parameters
use PENF
use MPI
use CUDAFOR
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none
private
public :: equation_euler_gpu_object
public :: BC_EXTRAPOLATION
public :: BC_INFLOW

integer(I4P), parameter :: BC_EXTRAPOLATION = 1_I4P
integer(I4P), parameter :: BC_INFLOW        = 2_I4P

type :: equation_euler_gpu_object
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
   type(field_object), pointer :: field=>null()      !< The field.
   type(base_gpu_object)       :: base_gpu           !< The base GPU handler.
   integer(I4P)                :: myrank=0_I4P       !< MPI rank process.
   integer(I4P)                :: procs_number=1_I4P !< Number of MPI processes.
   integer(I4P)                :: error=0_I4P        !< Error traping flag.
   ! equation data
   integer(I4P)           :: ns=1_I4P               !< Number of fluid species.
   real(R8P), allocatable :: cp0(:)                 !< Specific heat at constant pressure of initial species.
   real(R8P), allocatable :: cv0(:)                 !< Specific heat at constant pressure of initial species.
   real(R8P)              :: dt=0._R8P              !< Maximum time step accordingly to CFL criterion.
   real(R8P)              :: CFL=0.3_R8P            !< CFL limit.
   logical                :: null_xyz(3)=[.false.,&
                                          .false.,&
                                          .false.]  !< Flag triggering 1D/2D simulations.
   real(R8P), allocatable :: q_aux(:,:,:,:,:)       !< Auxiliary cell centered variables.
   ! Runge-Kutta data
   integer(I4P)           :: nrk=3_I4P !< Runge-Kutta stages number.
   real(R8P), allocatable :: alph(:,:) !< RK alpha coefficients.
   real(R8P), allocatable :: beta(:)   !< RK beta coefficients.
   real(R8P), allocatable :: gamm(:)   !< RK gamma coefficients.
   ! WENO data
   integer(I4P)           :: weno_s=1_I4P    !< Stencil number.
   real(R8P), allocatable :: weno_c(:,:)     !< Central difference coefficients    [1:2,1:2*S].
   real(R8P), allocatable :: weno_a(:,:)     !< Optimal weights                    [1:2,0:S-1].
   real(R8P), allocatable :: weno_p(:,:,:)   !< Polinomials coefficients           [1:2,0:S-1,0:S-1].
   real(R8P), allocatable :: weno_d(:,:,:)   !< Smoothness indicators coefficients [0:S-1,0:S-1,0:S-1].
   real(R8P)              :: weno_eps=0._R8P !< Parameter for avoiding divided by zero when computing smoothness indicators.
   integer(I4P)           :: weno_odd=1_I4P  !< Constant for distinguishing between odd and even number of stencils (mod(S,2)).
   integer(I4P)           :: weno_exp=0_I4P  !< Exponent for growing the diffusive part of weights.
   ! cuf data
   real(R8P), allocatable, device :: cp0_gpu(:)           !< Specific heat at constant pressure of initial species.
   real(R8P), allocatable, device :: cv0_gpu(:)           !< Specific heat at constant pressure of initial species.
   real(R8P), allocatable, device :: f_i_gpu(:,:,:,:,:)   !< Fluxes of cell centered variables, i direction.
   real(R8P), allocatable, device :: f_j_gpu(:,:,:,:,:)   !< Fluxes of cell centered variables, j direction.
   real(R8P), allocatable, device :: f_k_gpu(:,:,:,:,:)   !< Fluxes of cell centered variables, k direction.
   real(R8P), allocatable, device :: f_rho_gpu(:,:,:,:)   !< Fluxes of cell centered variables, rho, 1D, normal direction.
   real(R8P), allocatable, device :: f_rho_u_gpu(:,:,:,:) !< Fluxes of cell centered variables, rho*u, 1D, normal direction.
   real(R8P), allocatable, device :: f_rho_E_gpu(:,:,:,:) !< Fluxes of cell centered variables, rho*E, 1D, normal direction.
   real(R8P), allocatable, device :: dxyz_gpu(:,:)        !< Space steps.
   real(R8P), allocatable, device :: alph_gpu(:,:)        !< RK alpha coefficients.
   real(R8P), allocatable, device :: beta_gpu(:)          !< RK beta coefficients.
   real(R8P), allocatable, device :: gamm_gpu(:)          !< RK gamma coefficients.
   real(R8P), allocatable, device :: q_aux_gpu(:,:,:,:,:) !< Auxiliary cell centered variables.
   real(R8P), allocatable, device :: q_gpu(:,:,:,:,:)     !< Field cell centered variables stages.
   real(R8P), allocatable, device :: q_s_gpu(:,:,:,:,:,:) !< RK Field cell centered variables stages.
   real(R8P), allocatable, device :: weno_c_gpu(:,:)      !< Central difference coefficients    [1:2,1:2*S].
   real(R8P), allocatable, device :: weno_a_gpu(:,:)      !< Optimal weights                    [1:2,0:S-1].
   real(R8P), allocatable, device :: weno_p_gpu(:,:,:)    !< Polinomials coefficients           [1:2,0:S-1,0:S-1].
   real(R8P), allocatable, device :: weno_d_gpu(:,:,:)    !< Smoothness indicators coefficients [0:S-1,0:S-1,0:S-1].
   contains
      ! public methods
      procedure, pass(self) :: compute_aux             !< Compute auxiliary variables.
      procedure, pass(self) :: compute_dt              !< Compute time step.
      procedure, pass(self) :: copy_cpu_gpu            !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu            !< Copy data from GPU to CPU.
      procedure, pass(self) :: destroy                 !< Destroy the equation.
      procedure, pass(self) :: initialize              !< Initialize the equation.
      procedure, pass(self) :: mark_by_grad_rho        !< Mark blocks to be refined/derefined by a `grad(rho)` value.
      procedure, pass(self) :: integrate               !< Runge Kutta integration of equation.
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions of equation.
      procedure, pass(self) :: set_initial_conditions  !< Set initial conditions of equation.
      procedure, pass(self) :: update_ghost_gpu        !< Update ghost cells and set boundary conditions.
      ! operators
      generic :: assignment(=) => eq_assign_eq      !< Overload `=`.
      procedure, pass(lhs), private :: eq_assign_eq !< Operator `=`.
      ! private methods
      procedure, pass(self) :: weno_initialize !< Initialize WENO data.
endtype equation_euler_gpu_object

contains
   ! public methods
   subroutine compute_aux(self, q_gpu, q_aux_gpu)
   !< Compute auxiliary variables.
   class(equation_euler_gpu_object), intent(in)          :: self          !< The equation.
   real(R8P),                        intent(in),  device :: q_gpu(1:,                    &
                                                                  1-self%field%grid%gci:,&
                                                                  1-self%field%grid%gcj:,&
                                                                  1-self%field%grid%gck:,&
                                                                  1:)     !< Conservative variables.
   real(R8P),                        intent(out), device :: q_aux_gpu(1:,                    &
                                                                      1-self%field%grid%gci:,&
                                                                      1-self%field%grid%gcj:,&
                                                                      1-self%field%grid%gck:,&
                                                                      1:) !< Auxiliary variables.

   associate(blocks_number=>self%field%blocks_number,                                      &
             ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk,       &
             gci=>self%field%grid%gci, gcj=>self%field%grid%gcj, gck=>self%field%grid%gck, &
             ns=>self%ns)
      call compute_aux_cuf(ni=ni, nj=nj, nk=nk, gci=gci, gcj=gcj, gck=gck, ns=ns, &
                           blocks_number=blocks_number,                           &
                           cp0_gpu=self%cp0_gpu, cv0_gpu=self%cv0_gpu,            &
                           q_gpu=q_gpu, q_aux_gpu=q_aux_gpu)
   endassociate
   endsubroutine compute_aux

   subroutine compute_dt(self)
   !< Compute maximum time step accordingly to CFL stabilty criterion.
   class(equation_euler_gpu_object), intent(inout) :: self !< The equation.
   real(R8P)                                       :: umax !< Maximum speed of waves propagation.
   integer(I4P)                                    :: b    !< Counter.

   associate(blocks_number=>self%field%blocks_number, dxyz=>self%field%dxyz,               &
             ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk,       &
             gci=>self%field%grid%gci, gcj=>self%field%grid%gcj, gck=>self%field%grid%gck, &
             ns=>self%ns, q=>self%field%q, dt=>self%dt, CFL=>self%CFL)
      call self%compute_aux(q_gpu=self%q_gpu, q_aux_gpu=self%q_aux_gpu)
      dt = huge(1._R8P)
      do b=1, self%field%blocks_number
         call compute_umax_cuf(b, ni=ni, nj=nj, nk=nk, gci=gci, gcj=gcj, gck=gck, ns=ns, &
                               dx=dxyz(1,b), dy=dxyz(2,b), dz=dxyz(3,b),                 &
                               q_aux_gpu=self%q_aux_gpu, umax=umax)
         dt = min(dt, minval(dxyz(:,b)) / umax * CFL)
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, self%error)
   endassociate
   endsubroutine compute_dt

   subroutine copy_cpu_gpu(self)
   !< Copy data from CPU to GPU.
   class(equation_euler_gpu_object), intent(inout) :: self        !< The base backend.
   real(R8P), allocatable                          :: dxyz_t(:,:) !< Space steps transposed.
   integer(I4P)                                    :: i, b        !< Counter.

   allocate(dxyz_t(1:self%field%nb,3))
   do b=1, self%field%blocks_number
      do i=1, 3
         dxyz_t(b,i) = self%field%dxyz(i,b)
      enddo
   enddo
   self%dxyz_gpu = dxyz_t
   call self%base_gpu%copy_transpose_cpu_gpu(q_cpu=self%field%q, q_gpu=self%q_gpu)
   call self%base_gpu%copy_cpu_gpu
   endsubroutine copy_cpu_gpu

   subroutine copy_gpu_cpu(self, compute_q_aux)
   !< Copy data from GPU to CPU.
   class(equation_euler_gpu_object), intent(inout)        :: self          !< The base backend.
   logical,                          intent(in), optional :: compute_q_aux !< Flag to compute auxiliary variables.

   call self%base_gpu%copy_transpose_gpu_cpu(nv=self%field%nv, q_gpu=self%q_gpu, q_cpu=self%field%q)
   if (present(compute_q_aux)) then
      if (compute_q_aux) then
         call self%compute_aux(q_gpu=self%q_gpu, q_aux_gpu=self%q_aux_gpu)
         call self%base_gpu%copy_transpose_gpu_cpu(nv=self%ns+6, q_gpu=self%q_aux_gpu, q_cpu=self%q_aux)
      endif
   endif
   endsubroutine copy_gpu_cpu

   subroutine destroy(self)
   !< Destroy the equation.
   class(equation_euler_gpu_object), intent(inout) :: self  !< The equation.
   type(equation_euler_gpu_object)                 :: fresh !< Fresh equation.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, field, ns, nrk, cp0, cv0, CFL, null_xyz, weno_s)
   !< Initialize the equation.
   class(equation_euler_gpu_object), intent(inout)        :: self          !< The equation.
   type(field_object),               intent(in), target   :: field         !< The field.
   integer(I4P),                     intent(in), optional :: ns            !< Species number.
   integer(I4P),                     intent(in), optional :: nrk           !< Runge-Kutta stages number.
   real(R8P),                        intent(in), optional :: cp0(:)        !< Initial specific heats at constant pressure.
   real(R8P),                        intent(in), optional :: cv0(:)        !< Initial specific heats at constant volume.
   real(R8P),                        intent(in), optional :: CFL           !< CFL value.
   logical,                          intent(in), optional :: null_xyz(3)   !< Flag triggering 1D/2D simulations.
   integer(I4P),                     intent(in), optional :: weno_s        !< Number of WENO stencils.
   integer(I4P)                                           :: v             !< Counter.

   ! CPU data
   call self%destroy
   self%field => field
   if (present(ns)) self%ns = ns
   if (self%field%nv - self%ns /= 4) then
      write(stderr, '(A)') 'ADAM-ERROR: field%nv must be euler%ns+4'
      call MPI_FINALIZE(self%error)
      stop
   endif
   call self%base_gpu%initialize(field=field, nv_aux=self%ns+6)
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
   if (present(nrk)) self%nrk = nrk
   if (present(CFL)) self%CFL = CFL
   if (present(null_xyz)) self%null_xyz = null_xyz
   allocate(self%q_aux(1:self%ns+6,                                   &
                       1-field%grid%gci:field%grid%ni+field%grid%gci, &
                       1-field%grid%gcj:field%grid%nj+field%grid%gcj, &
                       1-field%grid%gck:field%grid%nk+field%grid%gck, 1:field%nb))
   if (present(weno_s)) self%weno_s = weno_s
   call self%weno_initialize
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
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   ! GPU data
   allocate(self%f_i_gpu(1:field%nb,                                    &
                         0-field%grid%gci:field%grid%ni+field%grid%gci, &
                         1-field%grid%gcj:field%grid%nj+field%grid%gcj, &
                         1-field%grid%gck:field%grid%nk+field%grid%gck, 1:field%nv))
   allocate(self%f_j_gpu(1:field%nb,                                    &
                         1-field%grid%gci:field%grid%ni+field%grid%gci, &
                         0-field%grid%gcj:field%grid%nj+field%grid%gcj, &
                         1-field%grid%gck:field%grid%nk+field%grid%gck, 1:field%nv))
   allocate(self%f_k_gpu(1:field%nb,                                    &
                         1-field%grid%gci:field%grid%ni+field%grid%gci, &
                         1-field%grid%gcj:field%grid%nj+field%grid%gcj, &
                         0-field%grid%gck:field%grid%nk+field%grid%gck, 1:field%nv))
   allocate(self%f_rho_gpu(1:field%nb,                                    &
                           0-field%grid%gci:field%grid%ni+field%grid%gci, &
                           0-field%grid%gcj:field%grid%nj+field%grid%gcj, &
                           0-field%grid%gck:field%grid%nk+field%grid%gck))
   allocate(self%f_rho_u_gpu(1:field%nb,                                    &
                             0-field%grid%gci:field%grid%ni+field%grid%gci, &
                             0-field%grid%gcj:field%grid%nj+field%grid%gcj, &
                             0-field%grid%gck:field%grid%nk+field%grid%gck))
   allocate(self%f_rho_E_gpu(1:field%nb,                                    &
                             0-field%grid%gci:field%grid%ni+field%grid%gci, &
                             0-field%grid%gcj:field%grid%nj+field%grid%gcj, &
                             0-field%grid%gck:field%grid%nk+field%grid%gck))
   allocate(self%dxyz_gpu(1:field%nb, 1:3))
   allocate(self%q_aux_gpu(1:field%nb,                                    &
                           1-field%grid%gci:field%grid%ni+field%grid%gci, &
                           1-field%grid%gcj:field%grid%nj+field%grid%gcj, &
                           1-field%grid%gck:field%grid%nk+field%grid%gck, 1:self%ns+6))
   allocate(self%q_gpu(1:field%nb,                                    &
                       1-field%grid%gci:field%grid%ni+field%grid%gci, &
                       1-field%grid%gcj:field%grid%nj+field%grid%gcj, &
                       1-field%grid%gck:field%grid%nk+field%grid%gck, 1:field%nv))
   allocate(self%q_s_gpu(1:field%nb,                                    &
                         1-field%grid%gci:field%grid%ni+field%grid%gci, &
                         1-field%grid%gcj:field%grid%nj+field%grid%gcj, &
                         1-field%grid%gck:field%grid%nk+field%grid%gck, 1:field%nv, 1:self%nrk))
   ! copy data that is not variable during the simulation
   self%cp0_gpu    = self%cp0
   self%cv0_gpu    = self%cv0
   self%alph_gpu   = self%alph
   self%beta_gpu   = self%beta
   self%gamm_gpu   = self%gamm
   self%weno_c_gpu = self%weno_c
   self%weno_a_gpu = self%weno_a
   self%weno_p_gpu = self%weno_p
   self%weno_d_gpu = self%weno_d
   endsubroutine initialize

   subroutine mark_by_grad_rho(self, grad_tol, delta_fine, delta_coarse, threshold)
   !< Mark blocks to be refined/derefined by a `grad(rho)` value.
   class(equation_euler_gpu_object), intent(inout)        :: self           !< The equation.
   real(R8P),                        intent(in)           :: grad_tol       !< Gradiend tolerance value.
   real(R8P),                        intent(in)           :: delta_fine     !< Maximum cell delta in fine grids.
   real(R8P),                        intent(in)           :: delta_coarse   !< Minimum cell delta in coarse grids.
   real(R8P),                        intent(in), optional :: threshold      !< Threshold for sphere proximity.
   real(R8P)                                              :: threshold_     !< Threshold for sphere proximity, local var.
   real(R8P)                                              :: max_cell_delta !< Maximum cell delta.
   real(R8P)                                              :: grad_rho       !< Value (max) of gradient of rho.
   integer(I4P)                                           :: b              !< Counter.

   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   self%field%refinements_needed = [(TO_NOT_TOUCH,b=1,self%field%blocks_number)]
   call self%update_ghost_gpu(q_gpu=self%q_gpu)
   associate (ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk,       &
              gci=>self%field%grid%gci, gcj=>self%field%grid%gcj, gck=>self%field%grid%gck, &
              blocks_number=>self%field%blocks_number, ns=>self%ns, dxyz=>self%field%dxyz)
      call self%compute_aux(q_gpu=self%q_gpu, q_aux_gpu=self%q_aux_gpu)
      do b=1, blocks_number
         grad_rho = gradient_cuf(b=b, ni=ni, nj=nj, nk=nk, gci=gci, gcj=gcj, gck=gck, ns=ns, &
                                 dx=dxyz(1,b), dy=dxyz(2,b), dz=dxyz(3,b), q_gpu=self%q_aux_gpu)

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
      function gradient_cuf(b, ni, nj, nk, gci, gcj, gck, ns, dx, dy, dz, q_gpu) result(gradient)
      !< Gradient done by CUF threads.
      integer(I4P), intent(in)         :: b            !< Block index.
      integer(I4P), intent(in)         :: ni           !< Grid cells number in I direction.
      integer(I4P), intent(in)         :: nj           !< Grid cells number in J direction.
      integer(I4P), intent(in)         :: nk           !< Grid cells number in K direction.
      integer(I4P), intent(in)         :: gci          !< Ghost grid cells number in I direction.
      integer(I4P), intent(in)         :: gcj          !< Ghost grid cells number in J direction.
      integer(I4P), intent(in)         :: gck          !< Ghost grid cells number in K direction.
      integer(I4P), intent(in)         :: ns           !< Species number.
      real(R8P),    intent(in)         :: dx           !< X space step.
      real(R8P),    intent(in)         :: dy           !< Y space step.
      real(R8P),    intent(in)         :: dz           !< Z space step.
      real(R8P),    intent(in), device :: q_gpu(1:,    &
                                                1-gci:,&
                                                1-gcj:,&
                                                1-gck:,&
                                                1:)    !< Field component to which apply gradient.
      real(R8P)                        :: gradient     !< Maximum gradient of q.
      real(R8P)                        :: grad         !< Current gradient of q.
      integer(I4P)                     :: i, j, k      !< Counter.
      integer(I4P)                     :: iercuda      !< Error trapping flag for CUDAFortran.

      gradient = 0._R8P
      !$cuf kernel do(3) <<<*,*>>>
      do k=1, nk
         do j=1, nj
            do i=1, ni
               grad = sqrt(((q_gpu(b,i+1,j,k,ns+1) - q_gpu(b,i-1,j,k,ns+1))/(2*dx))**2 + &
                           ((q_gpu(b,i,j+1,k,ns+1) - q_gpu(b,i,j-1,k,ns+1))/(2*dy))**2 + &
                           ((q_gpu(b,i,j,k+1,ns+1) - q_gpu(b,i,j,k-1,ns+1))/(2*dz))**2)
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
   endsubroutine mark_by_grad_rho

   subroutine integrate(self, t, do_ghost_syncro, residual)
   !< Runge Kutta integration of field.
   class(equation_euler_gpu_object), intent(inout)         :: self             !< The equation.
   real(R8P),                        intent(in)            :: t                !< Time.
   logical,                          intent(in),  optional :: do_ghost_syncro  !< Flag to do syncrous ghost update.
   real(R8P),                        intent(out), optional :: residual         !< Global residual.
   logical                                                 :: do_ghost_syncro_ !< Flag to do syncrous ghost update, local var.
   integer(I4P)                                            :: s                !< Counter.

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   associate(alph=>self%alph, beta=>self%beta, gamm=>self%gamm,                                      &
             dt=>self%dt,                                                                            &
             ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk,                 &
             gci=>self%field%grid%gci, gcj=>self%field%grid%gcj, gck=>self%field%grid%gck,           &
             nv=>self%field%nv, nrk=>self%nrk, ns=>self%ns, blocks_number=>self%field%blocks_number, &
             inner_blocks_number=>self%field%inner_blocks_number,                                    &
             alph_gpu=>self%alph_gpu, beta_gpu=>self%beta_gpu)
   do s=1, nrk
      call compute_rk_stage_gpu_cuf(ni=ni, nj=nj, nk=nk, gci=gci, gcj=gcj, gck=gck, nv=nv, blocks_number=blocks_number, &
                                    alph_gpu=alph_gpu, dt=dt, s=s, q_gpu=self%q_gpu, q_s_gpu=self%q_s_gpu)
      if (do_ghost_syncro_) then
         call self%update_ghost_gpu(q_gpu=self%q_s_gpu(:,:,:,:,:,s)) ! all ghosts
         call self%compute_aux(q_gpu=self%q_s_gpu(:,:,:,:,:,s), q_aux_gpu=self%q_aux_gpu)
         call compute_residuals_gpu_cuf(ni=ni, nj=nj, nk=nk, gci=gci, gcj=gcj, gck=gck, blocks_number=blocks_number, ns=ns, &
                                        null_x=self%null_xyz(1), null_y=self%null_xyz(2), null_z=self%null_xyz(3),          &
                                        dx_gpu      = self%dxyz_gpu(:,1),                                                   &
                                        dy_gpu      = self%dxyz_gpu(:,2),                                                   &
                                        dz_gpu      = self%dxyz_gpu(:,3),                                                   &
                                        q_aux_gpu   = self%q_aux_gpu,                                                       &
                                        f_i_gpu     = self%f_i_gpu,                                                         &
                                        f_j_gpu     = self%f_j_gpu,                                                         &
                                        f_k_gpu     = self%f_k_gpu,                                                         &
                                        f_rho_gpu   = self%f_rho_gpu,                                                       &
                                        f_rho_u_gpu = self%f_rho_u_gpu,                                                     &
                                        f_rho_E_gpu = self%f_rho_E_gpu,                                                     &
                                        q_gpu       = self%q_s_gpu(:,:,:,:,:,s))
      else
         ! TODO
      endif
      if (present(residual).and.s==self%nrk) then
         ! TODO
      endif
   enddo
   call advance_q_gpu_cuf(ni=ni, nj=nj, nk=nk, gci=gci, gcj=gcj, gck=gck, nrk=nrk, nv=nv, blocks_number=blocks_number, &
                           beta_gpu=beta_gpu, dt=dt, q_s_gpu=self%q_s_gpu, q_gpu=self%q_gpu)
   endassociate
   endsubroutine integrate

   subroutine set_boundary_conditions(self, q_gpu)
   !< Set boundary conditions of equation.
   class(equation_euler_gpu_object), intent(in)            :: self                             !< The equation.
   real(R8P),                        intent(inout), device :: q_gpu(1:,                    &
                                                                    1-self%field%grid%gci:,&
                                                                    1-self%field%grid%gcj:,&
                                                                    1-self%field%grid%gck:,1:) !< Conservative variables.
   integer(I4P)                                            :: nv                               !< Number of cons. varibales.

   nv = self%field%nv
   if (allocated(self%base_gpu%local_map_bc_face_gpu  )) call set_bc_fec(local_map_bc=self%base_gpu%local_map_bc_face_gpu  )
   if (allocated(self%base_gpu%local_map_bc_edge_gpu  )) call set_bc_fec(local_map_bc=self%base_gpu%local_map_bc_edge_gpu  )
   if (allocated(self%base_gpu%local_map_bc_corner_gpu)) call set_bc_fec(local_map_bc=self%base_gpu%local_map_bc_corner_gpu)
   contains
      subroutine set_bc_fec(local_map_bc)
      integer(I8P), intent(in),    device :: local_map_bc(:,:) !< Local map for BC ghost cells.
      integer(I4P)                        :: b                 !< Counter.
      integer(I4P)                        :: f, i, j, k, v     !< Counter.
      integer(I4P)                        :: fec               !< Counter.
      integer(I4P)                        :: imin              !< Lower limit of ijk indexes.
      integer(I4P)                        :: jmin              !< Lower limit of ijk indexes.
      integer(I4P)                        :: kmin              !< Lower limit of ijk indexes.
      integer(I4P)                        :: imax              !< Upper limit of ijk indexes.
      integer(I4P)                        :: jmax              !< Upper limit of ijk indexes.
      integer(I4P)                        :: kmax              !< Upper limit of ijk indexes.
      integer(I4P)                        :: idelta            !< IJK delta step for extrapolation.
      integer(I4P)                        :: jdelta            !< IJK delta step for extrapolation.
      integer(I4P)                        :: kdelta            !< IJK delta step for extrapolation.
      integer(I4P)                        :: bc_type           !< Boundary condition type.
      integer(I4P)                        :: iercuda           !< Error trapping flag for CUDAFortran.

      !$cuf kernel do(1) <<<*,*>>>
      do f=1, size(local_map_bc, dim=1)
         b       = local_map_bc(f, 1 )
         fec     = local_map_bc(f, 2 )
         imin    = local_map_bc(f, 3 )
         jmin    = local_map_bc(f, 4 )
         kmin    = local_map_bc(f, 5 )
         imax    = local_map_bc(f, 6 )
         jmax    = local_map_bc(f, 7 )
         kmax    = local_map_bc(f, 8 )
         idelta  = local_map_bc(f, 9 )
         jdelta  = local_map_bc(f, 10)
         kdelta  = local_map_bc(f, 11)
         bc_type = local_map_bc(f, 12)
         if (bc_type == BC_EXTRAPOLATION) then
            do k=kmin, kmax, sign(1, kmax-kmin)
               do j=jmin, jmax, sign(1, jmax-jmin)
                  do i=imin, imax, sign(1, imax-imin)
                     do v=1, nv
                        q_gpu(b,i,j,k,v) = q_gpu(b, i-idelta, j-jdelta, k-kdelta, v)
                     enddo
                  enddo
               enddo
            enddo
         elseif (bc_type == BC_INFLOW) then
            do k=kmin, kmax, sign(1, kmax-kmin)
               do j=jmin, jmax, sign(1, jmax-jmin)
                  do i=imin, imax, sign(1, imax-imin)
                     do v=1, nv
                     enddo
                  enddo
               enddo
            enddo
         endif
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
      endsubroutine set_bc_fec
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self)
   !< Set initial conditions of field.
   class(equation_euler_gpu_object), intent(inout) :: self    !< The equation.
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
                  q(1,i,j,k,b) = 1._R8P
                  q(2,i,j,k,b) = 0._R8P
                  q(3,i,j,k,b) = 0._R8P
                  q(4,i,j,k,b) = 0._R8P
                  q(5,i,j,k,b) = 1._R8P * E(p=1._R8P, r=1._R8P, u=0._R8P, g=self%cp0(1)/self%cv0(1))
               else
                  q(1,i,j,k,b) = 0.125_R8P
                  q(2,i,j,k,b) = 0._R8P
                  q(3,i,j,k,b) = 0._R8P
                  q(4,i,j,k,b) = 0._R8P
                  q(5,i,j,k,b) = 0.125_R8P * E(p=0.1_R8P, r=0.125_R8P, u=0._R8P, g=self%cp0(1)/self%cv0(1))
               endif
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine set_initial_conditions

   subroutine update_ghost_gpu(self, q_gpu, step)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(equation_euler_gpu_object), intent(inout)         :: self            !< The equation.
   real(R8P),                        intent(inout), device :: q_gpu(1:,                    &
                                                                    1-self%field%grid%gci:,&
                                                                    1-self%field%grid%gcj:,&
                                                                    1-self%field%grid%gck:,&
                                                                    1:)       !< Conservative variables.
   integer(I4P),                     intent(in), optional  :: step            !< Step to be perfordmed in asyncronous comp.
   logical                                                 :: do_local_update !< Flag for triggering local update.
   logical                                                 :: do_set_bc       !< Flag for triggering setting bc.

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
   class(equation_euler_gpu_object), intent(inout) :: lhs !< Left hand side.
   type(equation_euler_gpu_object),  intent(in)    :: rhs !< Right hand side.

   lhs%field => rhs%field
   lhs%base_gpu = rhs%base_gpu
   lhs%myrank = rhs%myrank
   lhs%procs_number = rhs%procs_number
   lhs%error = rhs%error
   lhs%ns = rhs%ns
   lhs%dt = rhs%dt
   lhs%CFL = rhs%CFL
   lhs%null_xyz = rhs%null_xyz
   lhs%nrk = rhs%nrk
   lhs%weno_s = rhs%weno_s
   lhs%weno_eps = rhs%weno_eps
   lhs%weno_odd = rhs%weno_odd
   lhs%weno_exp = rhs%weno_exp
   call assign_allocatable(lhs=lhs%q_aux, rhs=rhs%q_aux )
   call assign_allocatable(lhs=lhs%cp0 ,  rhs=rhs%cp0   )
   call assign_allocatable(lhs=lhs%cv0 ,  rhs=rhs%cv0   )
   call assign_allocatable(lhs=lhs%alph,  rhs=rhs%alph  )
   call assign_allocatable(lhs=lhs%beta,  rhs=rhs%beta  )
   call assign_allocatable(lhs=lhs%gamm,  rhs=rhs%gamm  )
   call assign_allocatable(lhs=lhs%weno_c,rhs=rhs%weno_c)
   call assign_allocatable(lhs=lhs%weno_a,rhs=rhs%weno_a)
   call assign_allocatable(lhs=lhs%weno_p,rhs=rhs%weno_p)
   call assign_allocatable(lhs=lhs%weno_d,rhs=rhs%weno_d)
   call assign_allocatable_gpu(lhs=lhs%cp0_gpu,    rhs=rhs%cp0_gpu    )
   call assign_allocatable_gpu(lhs=lhs%cv0_gpu,    rhs=rhs%cv0_gpu    )
   call assign_allocatable_gpu(lhs=lhs%f_i_gpu,    rhs=rhs%f_i_gpu    )
   call assign_allocatable_gpu(lhs=lhs%f_j_gpu,    rhs=rhs%f_j_gpu    )
   call assign_allocatable_gpu(lhs=lhs%f_k_gpu,    rhs=rhs%f_k_gpu    )
   call assign_allocatable_gpu(lhs=lhs%f_rho_gpu,  rhs=rhs%f_rho_gpu  )
   call assign_allocatable_gpu(lhs=lhs%f_rho_u_gpu,rhs=rhs%f_rho_u_gpu)
   call assign_allocatable_gpu(lhs=lhs%f_rho_E_gpu,rhs=rhs%f_rho_E_gpu)
   call assign_allocatable_gpu(lhs=lhs%dxyz_gpu,   rhs=rhs%dxyz_gpu   )
   call assign_allocatable_gpu(lhs=lhs%alph_gpu,   rhs=rhs%alph_gpu   )
   call assign_allocatable_gpu(lhs=lhs%beta_gpu,   rhs=rhs%beta_gpu   )
   call assign_allocatable_gpu(lhs=lhs%gamm_gpu,   rhs=rhs%gamm_gpu   )
   call assign_allocatable_gpu(lhs=lhs%q_aux_gpu,  rhs=rhs%q_aux_gpu  )
   call assign_allocatable_gpu(lhs=lhs%q_gpu,      rhs=rhs%q_gpu      )
   call assign_allocatable_gpu(lhs=lhs%q_s_gpu,    rhs=rhs%q_s_gpu    )
   call assign_allocatable_gpu(lhs=lhs%weno_c_gpu, rhs=rhs%weno_c_gpu )
   call assign_allocatable_gpu(lhs=lhs%weno_a_gpu, rhs=rhs%weno_a_gpu )
   call assign_allocatable_gpu(lhs=lhs%weno_p_gpu, rhs=rhs%weno_p_gpu )
   call assign_allocatable_gpu(lhs=lhs%weno_d_gpu, rhs=rhs%weno_d_gpu )
   endsubroutine eq_assign_eq

   ! private methods
   subroutine weno_initialize(self)
   !< Initialize WENO data.
   class(equation_euler_gpu_object), intent(inout) :: self !< The equation.

   if (self%weno_s==1) return
   ! initialize weno_exp
   self%weno_exp = self%weno_s
   if (self%weno_s>4) self%weno_exp = self%weno_s - 1
   ! computing weno_odd
   self%weno_odd = mod(self%weno_s,2)
   self%weno_eps = 0.00000000001_R8P
   ! allocating variables
   if (allocated(self%weno_c)) deallocate(self%weno_c) ; allocate(self%weno_c(1:2,1:2*self%weno_s))
   if (allocated(self%weno_a)) deallocate(self%weno_a) ; allocate(self%weno_a(1:2,0:self%weno_s-1))
   if (allocated(self%weno_p)) deallocate(self%weno_p) ; allocate(self%weno_p(1:2,0:self%weno_s-1,0:self%weno_s-1))
   if (allocated(self%weno_d)) deallocate(self%weno_d) ; allocate(self%weno_d(0:self%weno_s-1,0:self%weno_s-1,0:self%weno_s-1))
   associate(s=>self%weno_s, weno_exp=>self%weno_exp, weno_odd=>self%weno_odd, weno_eps=>self%weno_eps, &
             weno_c=>self%weno_c, weno_a=>self%weno_a, weno_p=>self%weno_p, weno_d=>self%weno_d)
   ! inizializing the coefficients
   select case(s)
   case(2) ! 3rd order WENO reconstruction
     ! central difference coefficients
     ! 1 => left interface (i-1/2)
     weno_c(1,1) = -1._R8P/12._R8P ! cell -2
     weno_c(1,2) =  7._R8P/12._R8P ! cell -1
     weno_c(1,3) =  7._R8P/12._R8P ! cell  0
     weno_c(1,4) = -1._R8P/12._R8P ! cell  1
     ! 2 => right interface (i+1/2)
     weno_c(2,1) = -1._R8P/12._R8P ! cell -1
     weno_c(2,2) =  7._R8P/12._R8P ! cell  0
     weno_c(2,3) =  7._R8P/12._R8P ! cell  1
     weno_c(2,4) = -1._R8P/12._R8P ! cell  2

     ! optimal weights
     ! 1 => left interface (i-1/2)
     weno_a(1,0) = 2._R8P/3._R8P ! stencil 0
     weno_a(1,1) = 1._R8P/3._R8P ! stencil 1
     ! 2 => right interface (i+1/2)
     weno_a(2,0) = 1._R8P/3._R8P ! stencil 0
     weno_a(2,1) = 2._R8P/3._R8P ! stencil 1

     ! polinomials coefficients
     ! 1 => left interface (i-1/2)
     !  cell  0               ;    cell  1
     weno_p(1,0,0) =  0.5_R8P ; weno_p(1,1,0) =  0.5_R8P ! stencil 0
     weno_p(1,0,1) = -0.5_R8P ; weno_p(1,1,1) =  1.5_R8P ! stencil 1
     ! 2 => right interface (i+1/2)
     !  cell  0               ;    cell  1
     weno_p(2,0,0) =  1.5_R8P ; weno_p(2,1,0) = -0.5_R8P ! stencil 0
     weno_p(2,0,1) =  0.5_R8P ; weno_p(2,1,1) =  0.5_R8P ! stencil 1

     ! smoothness indicators coefficients
     ! stencil 0
     !      i*i             ;       (i-1)*i
     weno_d(0,0,0) = 1._R8P ; weno_d(1,0,0) =-2._R8P
     !      /               ;       (i-1)*(i-1)
     weno_d(0,1,0) = 0._R8P ; weno_d(1,1,0) = 1._R8P
     ! stencil 1
     !     (i+1)*(i+1)      ;       (i+1)*i
     weno_d(0,0,1) = 1._R8P ; weno_d(1,0,1) =-2._R8P
     !      /               ;        i*i
     weno_d(0,1,1) = 0._R8P ; weno_d(1,1,1) = 1._R8P
   case(3) ! 5th order WENO reconstruction
     ! central difference coefficients
     ! 1 => left interface (i-1/2)
     weno_c(1,1) =  1._R8P/60._R8P ! cell -3
     weno_c(1,2) = -7.5_R8P        ! cell -2
     weno_c(1,3) = 37._R8P/60._R8P ! cell -1
     weno_c(1,4) = 37._R8P/60._R8P ! cell  0
     weno_c(1,5) = -7.5_R8P        ! cell  1
     weno_c(1,6) =  1._R8P/60._R8P ! cell  2
     ! 2 => right interface (i+1/2)
     weno_c(1,1) =  1._R8P/60._R8P ! cell -2
     weno_c(1,2) = -7.5_R8P        ! cell -1
     weno_c(1,3) = 37._R8P/60._R8P ! cell  0
     weno_c(1,4) = 37._R8P/60._R8P ! cell  1
     weno_c(1,5) = -7.5_R8P        ! cell  2
     weno_c(1,6) =  1._R8P/60._R8P ! cell  3

     ! optimal weights
     ! 1 => left interface (i-1/2)
     weno_a(1,0) = 0.3_R8P ! stencil 0
     weno_a(1,1) = 0.6_R8P ! stencil 1
     weno_a(1,2) = 0.1_R8P ! stencil 2
     ! 2 => right interface (i+1/2)
     weno_a(2,0) = 0.1_R8P ! stencil 0
     weno_a(2,1) = 0.6_R8P ! stencil 1
     weno_a(2,2) = 0.3_R8P ! stencil 2

     ! polinomials coefficients
     ! 1 => left interface (i-1/2)
     !  cell  0                     ;    cell  1                     ;    cell  2
     weno_p(1,0,0) =  1._R8P/3._R8P ; weno_p(1,1,0) =  5._R8P/6._R8P ; weno_p(1,2,0) = -1._R8P/6._R8P ! stencil 0
     weno_p(1,0,1) = -1._R8P/6._R8P ; weno_p(1,1,1) =  5._R8P/6._R8P ; weno_p(1,2,1) =  1._R8P/3._R8P ! stencil 1
     weno_p(1,0,2) =  1._R8P/3._R8P ; weno_p(1,1,2) = -7._R8P/6._R8P ; weno_p(1,2,2) = 11._R8P/6._R8P ! stencil 2
     ! 2 => right interface (i+1/2)
     !  cell  0                     ;    cell  1                     ;    cell  2
     weno_p(2,0,0) = 11._R8P/6._R8P ; weno_p(2,1,0) = -7._R8P/6._R8P ; weno_p(2,2,0) =  1._R8P/3._R8P ! stencil 0
     weno_p(2,0,1) =  1._R8P/3._R8P ; weno_p(2,1,1) =  5._R8P/6._R8P ; weno_p(2,2,1) = -1._R8P/6._R8P ! stencil 1
     weno_p(2,0,2) = -1._R8P/6._R8P ; weno_p(2,1,2) =  5._R8P/6._R8P ; weno_p(2,2,2) =  1._R8P/3._R8P ! stencil 2

     ! smoothness indicators coefficients
     ! stencil 0
     !      i*i                      ;       (i-1)*i                   ;       (i-2)*i
     weno_d(0,0,0) =  10._R8P/3._R8P ; weno_d(1,0,0) = -31._R8P/3._R8P ; weno_d(2,0,0) =  11._R8P/3._R8P
     !      /                        ;       (i-1)*(i-1)               ;       (i-2)*(i-1)
     weno_d(0,1,0) =   0._R8P        ; weno_d(1,1,0) =  25._R8P/3._R8P ; weno_d(2,1,0) = -19._R8P/3._R8P
     !      /                        ;        /                        ;       (i-2)*(i-2)
     weno_d(0,2,0) =   0._R8P        ; weno_d(1,2,0) =   0._R8P        ; weno_d(2,2,0) =   4._R8P/3._R8P
     ! stencil 1
     !     (i+1)*(i+1)               ;        i*(i+1)                  ;       (i-1)*(i+1)
     weno_d(0,0,1) =   4._R8P/3._R8P ; weno_d(1,0,1) = -13._R8P/3._R8P ; weno_d(2,0,1) =   5._R8P/3._R8P
     !      /                        ;        i*i                      ;       (i-1)*i
     weno_d(0,1,1) =   0._R8P        ; weno_d(1,1,1) =  13._R8P/3._R8P ; weno_d(2,1,1) = -13._R8P/3._R8P
     !      /                        ;        /                        ;       (i-1)*(i-1)
     weno_d(0,2,1) =   0._R8P        ; weno_d(1,2,1) =   0._R8P        ; weno_d(2,2,1) =   4._R8P/3._R8P
     ! stencil 2
     !     (i+2)*(i+2)               ;       (i+1)*(i+2)               ;        i*(i+2)
     weno_d(0,0,2) =   4._R8P/3._R8P ; weno_d(1,0,2) = -19._R8P/3._R8P ; weno_d(2,0,2) =  11._R8P/3._R8P
     !      /                        ;       (i+1)*(i+1)               ;        i*(i+1)
     weno_d(0,1,2) =   0._R8P        ; weno_d(1,1,2) =  25._R8P/3._R8P ; weno_d(2,1,2) = -31._R8P/3._R8P
     !      /                        ;        /                        ;        i*i
     weno_d(0,2,2) =   0._R8P        ; weno_d(1,2,2) =   0._R8P        ; weno_d(2,2,2) =  10._R8P/3._R8P
   endselect
   endassociate
   endsubroutine weno_initialize

   ! non TBP cuf methods
   subroutine advance_q_gpu_cuf(ni, nj, nk, gci, gcj, gck, nv, nrk, blocks_number, beta_gpu, dt, q_s_gpu, q_gpu)
   !< Advance q_gpu by means of RK stages.
   integer(I4P), intent(in)            :: ni               !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj               !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk               !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: gci              !< Ghost grid cells number in I direction.
   integer(I4P), intent(in)            :: gcj              !< Ghost grid cells number in J direction.
   integer(I4P), intent(in)            :: gck              !< Ghost grid cells number in K direction.
   integer(I4P), intent(in)            :: nv               !< Number of conservative varibales.
   integer(I4P), intent(in)            :: nrk              !< Number of RK stages.
   integer(I4P), intent(in)            :: blocks_number    !< Number of blocks.
   real(R8P),    intent(in),    device :: beta_gpu(:)      !< RK betaa coefficients.
   real(R8P),    intent(in)            :: Dt               !< Time step.
   real(R8P),    intent(in),    device :: q_s_gpu(1:,    &
                                                  1-gci:,&
                                                  1-gcj:,&
                                                  1-gck:,&
                                                  1:,1:)   !< RK stage.
   real(R8P),    intent(inout), device ::   q_gpu(1:,    &
                                                  1-gci:,&
                                                  1-gcj:,&
                                                  1-gck:,&
                                                  1:)      !< Conservative variables.
   integer(I4P)                        :: i, j, k, b, s, v !< Counter.
   integer(I4P)                        :: iercuda          !< Error trapping flag for CUDAFortran.

   do s=1, nrk
      !$cuf kernel do(5) <<<*,*>>>
      do v=1, nv
         do k=1-gck, nk+gck
            do j=1-gcj, nj+gcj
               do i=1-gci, ni+gci
                  do b=1, blocks_number
                     q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + q_s_gpu(b,i,j,k,v,s) * dt * beta_gpu(s)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   enddo
   endsubroutine advance_q_gpu_cuf

   subroutine compute_aux_cuf(ni, nj, nk, gci, gcj, gck, ns, blocks_number, cp0_gpu, cv0_gpu, q_gpu, q_aux_gpu)
   !< Compute auxiliary variables by means of CUF threads.
   integer(I4P), intent(in)          :: ni            !< Grid cells number in I direction.
   integer(I4P), intent(in)          :: nj            !< Grid cells number in J direction.
   integer(I4P), intent(in)          :: nk            !< Grid cells number in K direction.
   integer(I4P), intent(in)          :: gci           !< Ghost grid cells number in I direction.
   integer(I4P), intent(in)          :: gcj           !< Ghost grid cells number in J direction.
   integer(I4P), intent(in)          :: gck           !< Ghost grid cells number in K direction.
   integer(I4P), intent(in)          :: ns            !< Number of fluid species.
   integer(I4P), intent(in)          :: blocks_number !< Number of blocks.
   real(R8P),    intent(in),  device :: cp0_gpu(:)    !< Specific heat at constant pressure of initial species.
   real(R8P),    intent(in),  device :: cv0_gpu(:)    !< Specific heat at constant pressure of initial species.
   real(R8P),    intent(in),  device :: q_gpu(1:,    &
                                              1-gci:,&
                                              1-gcj:,&
                                              1-gck:,&
                                              1:)     !< Conservative variables.
   real(R8P),    intent(out), device :: q_aux_gpu(1:,    &
                                                  1-gci:,&
                                                  1-gcj:,&
                                                  1-gck:,&
                                                  1:) !< Auxiliary variables.
   integer(I4P)                      :: b, i, j, k, s !< Counter.
   real(R8P)                         :: cp, cv        !< Specific heats.
   integer(I4P)                      :: iercuda       !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(4) <<<*,*>>>
   do k=1-gck, nk+gck
      do j=1-gcj, nj+gcj
         do i=1-gci, ni+gci
            do b=1, blocks_number
               q_aux_gpu(b,i,j,k,ns+1) = 0._R8P
               do s=1, ns
                  q_aux_gpu(b,i,j,k,ns+1) = q_aux_gpu(b,i,j,k,ns+1) + q_gpu(b,i,j,k,s)
               enddo
               do s=1, ns
                  q_aux_gpu(b,i,j,k,s) = q_gpu(b,i,j,k,s) / q_aux_gpu(b,i,j,k,ns+1)
               enddo
               q_aux_gpu(b,i,j,k,ns+2) = q_gpu(b,i,j,k,ns+1) / q_aux_gpu(b,i,j,k,ns+1)
               q_aux_gpu(b,i,j,k,ns+3) = q_gpu(b,i,j,k,ns+2) / q_aux_gpu(b,i,j,k,ns+1)
               q_aux_gpu(b,i,j,k,ns+4) = q_gpu(b,i,j,k,ns+3) / q_aux_gpu(b,i,j,k,ns+1)
               cp = 0._R8P
               cv = 0._R8P
               do s=1, ns
                  cp = cp + q_aux_gpu(b,i,j,k,s) * cp0_gpu(s)
                  cv = cv + q_aux_gpu(b,i,j,k,s) * cv0_gpu(s)
               enddo
               q_aux_gpu(b,i,j,k,ns+5) = cp / cv
               q_aux_gpu(b,i,j,k,ns+6) = (q_gpu(b,i,j,k,ns+4) - 0.5_R8P * q_aux_gpu(b,i,j,k,ns+1) *      &
                                                                         (q_aux_gpu(b,i,j,k,ns+2)**2 +   &
                                                                          q_aux_gpu(b,i,j,k,ns+3)**2 +   &
                                                                          q_aux_gpu(b,i,j,k,ns+4)**2)) * &
                                         (q_aux_gpu(b,i,j,k,ns+5) - 1._R8P)
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_aux_cuf

   subroutine compute_residuals_gpu_cuf(ni, nj, nk, gci, gcj, gck, ns, null_x, null_y, null_z, blocks_number, &
                                        dx_gpu, dy_gpu, dz_gpu, q_aux_gpu,                                    &
                                        f_i_gpu, f_j_gpu, f_k_gpu, f_rho_gpu, f_rho_u_gpu, f_rho_E_gpu, q_gpu)
   !< Compute residuals of equation.
   integer(I4P), intent(in)            :: ni                       !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                       !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                       !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: gci                      !< Ghost grid cells number in I direction.
   integer(I4P), intent(in)            :: gcj                      !< Ghost grid cells number in J direction.
   integer(I4P), intent(in)            :: gck                      !< Ghost grid cells number in K direction.
   integer(I4P), intent(in)            :: ns                       !< Number of species.
   logical,      intent(in)            :: null_x                   !< Nullify x direction.
   logical,      intent(in)            :: null_y                   !< Nullify y direction.
   logical,      intent(in)            :: null_z                   !< Nullify z direction.
   integer(I4P), intent(in)            :: blocks_number            !< Number of blocks.
   real(R8P),    intent(in),    device :: dx_gpu(1:)               !< X space steps.
   real(R8P),    intent(in),    device :: dy_gpu(1:)               !< Y space steps.
   real(R8P),    intent(in),    device :: dz_gpu(1:)               !< Z space steps.
   real(R8P),    intent(in),    device :: q_aux_gpu(1:,    &
                                                    1-gci:,&
                                                    1-gcj:,&
                                                    1-gck:,&
                                                    1:)            !< Auxiliary variables.
   real(R8P),    intent(inout), device :: f_i_gpu(1:,0:,1:,1:,1:)  !< Convective fluxes in i direction.
   real(R8P),    intent(inout), device :: f_j_gpu(1:,1:,0:,1:,1:)  !< Convective fluxes in j direction.
   real(R8P),    intent(inout), device :: f_k_gpu(1:,1:,1:,0:,1:)  !< Convective fluxes in k direction.
   real(R8P),    intent(inout), device ::   f_rho_gpu(1:,0:,0:,0:) !< Convective fluxes in normal direction, rho, 1D.
   real(R8P),    intent(inout), device :: f_rho_u_gpu(1:,0:,0:,0:) !< Convective fluxes in normal direction, rho*u, 1D.
   real(R8P),    intent(inout), device :: f_rho_E_gpu(1:,0:,0:,0:) !< Convective fluxes in normal direction, rho*E, 1D.
   real(R8P),    intent(inout), device :: q_gpu(1:,    &
                                                1-gci:,&
                                                1-gcj:,&
                                                1-gck:,&
                                                1:)                !< Conservative variables.
   integer(I4P)                        :: b, i, j, k, s            !< Counter.
   integer(I4P)                        :: iercuda                  !< Error trapping flag for CUDAFortran.

   ! fluxes in i direction
   if (.not.null_x) then
      !$cuf kernel do(4) <<<*,*>>>
      do k=1, nk
         do j=1, nj
            do i=0,ni
               do b=1, blocks_number
                  ! compute 1D normal fluxes
                  call solve_riemann(r1=q_aux_gpu(b,i  ,j,k,ns+1), &
                                     u1=q_aux_gpu(b,i  ,j,k,ns+2), &
                                     g1=q_aux_gpu(b,i  ,j,k,ns+5), &
                                     p1=q_aux_gpu(b,i  ,j,k,ns+6), &
                                     r4=q_aux_gpu(b,i+1,j,k,ns+1), &
                                     u4=q_aux_gpu(b,i+1,j,k,ns+2), &
                                     g4=q_aux_gpu(b,i+1,j,k,ns+5), &
                                     p4=q_aux_gpu(b,i+1,j,k,ns+6), &
                                     f_rho  =  f_rho_gpu(b,i,j,k), &
                                     f_rho_u=f_rho_u_gpu(b,i,j,k), &
                                     f_rho_E=f_rho_E_gpu(b,i,j,k))
                  ! compute 3D fluxes
                  if (f_rho_gpu(b,i,j,k)>0._R8P) then
                     do s=1, ns
                        f_i_gpu(b,i,j,k,s) = f_rho_gpu(b,i,j,k) * q_aux_gpu(b,i,j,k,s)
                     enddo
                     f_i_gpu(b,i,j,k,ns+1) = f_rho_u_gpu(b,i,j,k)
                     f_i_gpu(b,i,j,k,ns+2) =   f_rho_gpu(b,i,j,k) * q_aux_gpu(b,i,j,k,ns+3)
                     f_i_gpu(b,i,j,k,ns+3) =   f_rho_gpu(b,i,j,k) * q_aux_gpu(b,i,j,k,ns+4)
                     f_i_gpu(b,i,j,k,ns+4) = f_rho_E_gpu(b,i,j,k) + 0.5_R8P * f_rho_gpu(b,i,j,k) * &
                                                                  (q_aux_gpu(b,i,j,k,ns+3)**2 + q_aux_gpu(b,i,j,k,ns+4)**2)
                  else
                     do s=1, ns
                        f_i_gpu(b,i,j,k,s) = f_rho_gpu(b,i,j,k) * q_aux_gpu(b,i+1,j,k,s)
                     enddo
                     f_i_gpu(b,i,j,k,ns+1) = f_rho_u_gpu(b,i,j,k)
                     f_i_gpu(b,i,j,k,ns+2) =   f_rho_gpu(b,i,j,k) * q_aux_gpu(b,i+1,j,k,ns+3)
                     f_i_gpu(b,i,j,k,ns+3) =   f_rho_gpu(b,i,j,k) * q_aux_gpu(b,i+1,j,k,ns+4)
                     f_i_gpu(b,i,j,k,ns+4) = f_rho_E_gpu(b,i,j,k) + 0.5_R8P * f_rho_gpu(b,i,j,k) * &
                                                                  (q_aux_gpu(b,i+1,j,k,ns+3)**2 + q_aux_gpu(b,i+1,j,k,ns+4)**2)
                  endif
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   else
      !$cuf kernel do(5) <<<*,*>>>
      do s=1, ns + 4
         do k=1, nk
            do j=1, nj
               do i=0,ni
                  do b=1, blocks_number
                     f_i_gpu(b,i,j,k,s) = 0._R8P
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   endif

   ! fluxes in j direction
   if (.not.null_y) then
      !$cuf kernel do(4) <<<*,*>>>
      do k=1, nk
         do j=0, nj
            do i=1,ni
               do b=1, blocks_number
                  ! compute 1D normal fluxes
                  call solve_riemann(r1=q_aux_gpu(b,i,j  ,k,ns+1), &
                                     u1=q_aux_gpu(b,i,j  ,k,ns+3), &
                                     g1=q_aux_gpu(b,i,j  ,k,ns+5), &
                                     p1=q_aux_gpu(b,i,j  ,k,ns+6), &
                                     r4=q_aux_gpu(b,i,j+1,k,ns+1), &
                                     u4=q_aux_gpu(b,i,j+1,k,ns+3), &
                                     g4=q_aux_gpu(b,i,j+1,k,ns+5), &
                                     p4=q_aux_gpu(b,i,j+1,k,ns+6), &
                                     f_rho  =  f_rho_gpu(b,i,j,k), &
                                     f_rho_u=f_rho_u_gpu(b,i,j,k), &
                                     f_rho_E=f_rho_E_gpu(b,i,j,k))
                  ! compute 3D fluxes
                  if (f_rho_gpu(b,i,j,k)>0._R8P) then
                     do s=1, ns
                        f_j_gpu(b,i,j,k,s) = f_rho_gpu(b,i,j,k) * q_aux_gpu(b,i,j,k,s)
                     enddo
                     f_j_gpu(b,i,j,k,ns+1) =   f_rho_gpu(b,i,j,k) * q_aux_gpu(b,i,j,k,ns+2)
                     f_j_gpu(b,i,j,k,ns+2) = f_rho_u_gpu(b,i,j,k)
                     f_j_gpu(b,i,j,k,ns+3) =   f_rho_gpu(b,i,j,k) * q_aux_gpu(b,i,j,k,ns+4)
                     f_j_gpu(b,i,j,k,ns+4) = f_rho_E_gpu(b,i,j,k) + 0.5_R8P * f_rho_gpu(b,i,j,k) * &
                                                                  (q_aux_gpu(b,i,j,k,ns+2)**2 + q_aux_gpu(b,i,j,k,ns+4)**2)
                  else
                     do s=1, ns
                        f_j_gpu(b,i,j,k,s) = f_rho_gpu(b,i,j,k) * q_aux_gpu(b,i,j+1,k,s)
                     enddo
                     f_j_gpu(b,i,j,k,ns+1) =   f_rho_gpu(b,i,j,k) * q_aux_gpu(b,i,j+1,k,ns+2)
                     f_j_gpu(b,i,j,k,ns+2) = f_rho_u_gpu(b,i,j,k)
                     f_j_gpu(b,i,j,k,ns+3) =   f_rho_gpu(b,i,j,k) * q_aux_gpu(b,i,j+1,k,ns+4)
                     f_j_gpu(b,i,j,k,ns+4) = f_rho_E_gpu(b,i,j,k) + 0.5_R8P * f_rho_gpu(b,i,j,k) * &
                                                                  (q_aux_gpu(b,i,j+1,k,ns+2)**2 + q_aux_gpu(b,i,j+1,k,ns+4)**2)
                  endif
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   else
      !$cuf kernel do(5) <<<*,*>>>
      do s=1, ns + 4
         do k=1, nk
            do j=0, nj
               do i=1,ni
                  do b=1, blocks_number
                     f_j_gpu(b,i,j,k,s) = 0._R8P
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   endif

   ! fluxes in k direction
   if (.not.null_z) then
      !$cuf kernel do(4) <<<*,*>>>
      do k=0, nk
         do j=1, nj
            do i=1,ni
               do b=1, blocks_number
                  ! compute 1D normal fluxes
                  call solve_riemann(r1=q_aux_gpu(b,i,j,k  ,ns+1), &
                                     u1=q_aux_gpu(b,i,j,k  ,ns+4), &
                                     g1=q_aux_gpu(b,i,j,k  ,ns+5), &
                                     p1=q_aux_gpu(b,i,j,k  ,ns+6), &
                                     r4=q_aux_gpu(b,i,j,k+1,ns+1), &
                                     u4=q_aux_gpu(b,i,j,k+1,ns+4), &
                                     g4=q_aux_gpu(b,i,j,k+1,ns+5), &
                                     p4=q_aux_gpu(b,i,j,k+1,ns+6), &
                                     f_rho  =  f_rho_gpu(b,i,j,k), &
                                     f_rho_u=f_rho_u_gpu(b,i,j,k), &
                                     f_rho_E=f_rho_E_gpu(b,i,j,k))
                  ! compute 3D fluxes
                  if (f_rho_gpu(b,i,j,k)>0._R8P) then
                     do s=1, ns
                        f_k_gpu(b,i,j,k,s) = f_rho_gpu(b,i,j,k) * q_aux_gpu(b,i,j,k,s)
                     enddo
                     f_k_gpu(b,i,j,k,ns+1) =   f_rho_gpu(b,i,j,k) * q_aux_gpu(b,i,j,k,ns+2)
                     f_k_gpu(b,i,j,k,ns+2) =   f_rho_gpu(b,i,j,k) * q_aux_gpu(b,i,j,k,ns+3)
                     f_k_gpu(b,i,j,k,ns+3) = f_rho_u_gpu(b,i,j,k)
                     f_k_gpu(b,i,j,k,ns+4) = f_rho_E_gpu(b,i,j,k) + 0.5_R8P * f_rho_gpu(b,i,j,k) * &
                                                                  (q_aux_gpu(b,i,j,k,ns+2)**2 + q_aux_gpu(b,i,j,k,ns+3)**2)
                  else
                     do s=1, ns
                        f_k_gpu(b,i,j,k,s) = f_rho_gpu(b,i,j,k) * q_aux_gpu(b,i,j,k+1,s)
                     enddo
                     f_k_gpu(b,i,j,k,ns+1) =   f_rho_gpu(b,i,j,k) * q_aux_gpu(b,i,j,k+1,ns+2)
                     f_k_gpu(b,i,j,k,ns+2) =   f_rho_gpu(b,i,j,k) * q_aux_gpu(b,i,j,k+1,ns+3)
                     f_k_gpu(b,i,j,k,ns+3) = f_rho_u_gpu(b,i,j,k)
                     f_k_gpu(b,i,j,k,ns+4) = f_rho_E_gpu(b,i,j,k) + 0.5_R8P * f_rho_gpu(b,i,j,k) * &
                                                                  (q_aux_gpu(b,i,j,k+1,ns+2)**2 + q_aux_gpu(b,i,j,k+1,ns+3)**2)
                  endif
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   else
      !$cuf kernel do(5) <<<*,*>>>
      do s=1, ns + 4
         do k=1, nk
            do j=1, nj
               do i=1,ni
                  do b=1, blocks_number
                     f_k_gpu(b,i,j,k,s) = 0._R8P
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   endif

   !$cuf kernel do(5) <<<*,*>>>
   do s=1, ns+4
      do k=1, nk
         do j=1, nj
            do i=1, ni
               do b=1, blocks_number
                  q_gpu(b,i,j,k,s) = (f_i_gpu(b,i-1,j  ,k  ,s) - f_i_gpu(b,i,j,k,s)) / dx_gpu(b) + &
                                     (f_j_gpu(b,i  ,j-1,k  ,s) - f_j_gpu(b,i,j,k,s)) / dy_gpu(b) + &
                                     (f_k_gpu(b,i  ,j  ,k-1,s) - f_k_gpu(b,i,j,k,s)) / dz_gpu(b)
               enddo
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_residuals_gpu_cuf

   subroutine compute_rk_stage_gpu_cuf(ni, nj, nk, gci, gcj, gck, nv, blocks_number, alph_gpu, dt, s, q_gpu, q_s_gpu)
   !< Initialize RK stage with q_gpu.
   integer(I4P), intent(in)            :: ni                !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: gci               !< Ghost grid cells number in I direction.
   integer(I4P), intent(in)            :: gcj               !< Ghost grid cells number in J direction.
   integer(I4P), intent(in)            :: gck               !< Ghost grid cells number in K direction.
   integer(I4P), intent(in)            :: nv                !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number     !< Number of blocks.
   real(R8P),    intent(in),    device :: alph_gpu(:,:)     !< RK alpha coefficients.
   real(R8P),    intent(in)            :: dt                !< Time step.
   integer(I4P), intent(in)            :: s                 !< Stage to initialize.
   real(R8P),    intent(in),    device ::   q_gpu(1:,    &
                                                  1-gci:,&
                                                  1-gcj:,&
                                                  1-gck:,&
                                                  1:)       !< Conservative field.
   real(R8P),    intent(inout), device :: q_s_gpu(1:,    &
                                                  1-gci:,&
                                                  1-gcj:,&
                                                  1-gck:,&
                                                  1:,1:)    !< RK stage.
   integer(I4P)                        :: i, j, k, b, v, ss !< Counter.
   integer(I4P)                        :: iercuda           !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(5) <<<*,*>>>
   do v=1, nv
      do k=1, nk
         do j=1, nj
            do i=1, ni
               do b=1, blocks_number
                  q_s_gpu(b,i,j,k,v,s) = q_gpu(b,i,j,k,v)
               enddo
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   do ss=1, s - 1
      !$cuf kernel do(5) <<<*,*>>>
      do v=1, nv
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do b=1, blocks_number
                     q_s_gpu(b,i,j,k,v,s) = q_s_gpu(b,i,j,k,v,s) + (q_s_gpu(b,i,j,k,v,ss) * (dt * alph_gpu(s, ss)))
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   enddo
   endsubroutine compute_rk_stage_gpu_cuf

   subroutine compute_umax_cuf(b, ni, nj, nk, gci, gcj, gck, ns, dx, dy, dz, q_aux_gpu, umax)
   !< Compute maximum speed by means of CUF threads.
   integer(I4P), intent(in)         :: b             !< Block index.
   integer(I4P), intent(in)         :: ni            !< Grid cells number in I direction.
   integer(I4P), intent(in)         :: nj            !< Grid cells number in J direction.
   integer(I4P), intent(in)         :: nk            !< Grid cells number in K direction.
   integer(I4P), intent(in)         :: gci           !< Ghost grid cells number in I direction.
   integer(I4P), intent(in)         :: gcj           !< Ghost grid cells number in J direction.
   integer(I4P), intent(in)         :: gck           !< Ghost grid cells number in K direction.
   integer(I4P), intent(in)         :: ns            !< Number of species.
   real(R8P),    intent(in)         :: dx            !< X space step.
   real(R8P),    intent(in)         :: dy            !< Y space step.
   real(R8P),    intent(in)         :: dz            !< Z space step.
   real(R8P),    intent(in), device :: q_aux_gpu(1:,    &
                                                 1-gci:,&
                                                 1-gcj:,&
                                                 1-gck:,&
                                                 1:) !< Auxiliary varibales.
   real(R8P),    intent(out)        :: umax          !< Maximum speed.
   real(R8P)                        :: ss            !< Speed of sound.
   integer(I4P)                     :: i, j, k       !< Counter.
   integer(I4P)                     :: iercuda       !< Error trapping flag for CUDAFortran.

   umax = 0._R8P
   !$cuf kernel do(3) <<<*,*>>>
   do k=1, nk
      do j=1, nj
         do i=1, ni
            ss = a(p=q_aux_gpu(b,i,j,k,ns+6), r=q_aux_gpu(b,i,j,k,ns+1), g=q_aux_gpu(b,i,j,k,ns+5))
            umax = max(umax, abs(q_aux_gpu(b,i,j,k,ns+2)) + ss, &
                             abs(q_aux_gpu(b,i,j,k,ns+3)) + ss, &
                             abs(q_aux_gpu(b,i,j,k,ns+4)) + ss)

         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_umax_cuf

   ! non type-bound kernel procedures
   attributes(device) subroutine solve_riemann(r1, u1, p1, g1, r4, u4, p4, g4, f_rho, f_rho_u, f_rho_E)
   !< Solve the Riemann problem between the state $1$ and $4$ using the (local) Lax Friedrichs (Rusanov) solver.
   real(R8P),    intent(in)  :: r1            !< Density of state 1.
   real(R8P),    intent(in)  :: u1            !< Velocity of state 1.
   real(R8P),    intent(in)  :: p1            !< Pressure of state 1.
   real(R8P),    intent(in)  :: g1            !< Specific heats ratio of state 1.
   real(R8P),    intent(in)  :: r4            !< Density of state 4.
   real(R8P),    intent(in)  :: u4            !< Velocity of state 4.
   real(R8P),    intent(in)  :: p4            !< Pressure of state 4.
   real(R8P),    intent(in)  :: g4            !< Specific heats ratio of state 4.
   real(R8P),    intent(out) :: f_rho         !< Flux of rho at interface.
   real(R8P),    intent(out) :: f_rho_u       !< Flux of rho*u at interface.
   real(R8P),    intent(out) :: f_rho_E       !< Flux of rho*E at interface.
   real(R8P)                 :: F1(1:3)       !< State 1 fluxes.
   real(R8P)                 :: F4(1:3)       !< State 4 fluxes.
   real(R8P)                 :: lmax          !< Maximum wave speed estimation.
   real(R8P)                 :: u             !< Velocity of the intermediate states.
   real(R8P)                 :: S1            !< Maximum wave speed of state 1 and 4.
   real(R8P)                 :: S4            !< Maximum wave speed of state 1 and 4.
   integer(I4P)              :: s             !< Species counter.

   ! evaluate the intermediates states 2 and 3 from the known states U1,U4 using the PVRS approximation
   call compute_inter_states(r1=r1, u1=u1, p1=p1, g1=g1, r4=r4, u4=u4, p4=p4, g4=g4, u=u, S1=S1, S4=S4)
   ! evalute the maximum waves speed
   lmax = max(abs(S1), abs(u), abs(S4))
   ! compute the fluxes of state 1 and 4
   F1 = fluxes(p = p1, r = r1, u = u1, g = g1)
   F4 = fluxes(p = p4, r = r4, u = u4, g = g4)
   ! compute the Lax-Friedrichs fluxes approximation
   f_rho   = 0.5_R8P*(F1(1) + F4(1) - lmax*(r4                        - r1                       ))
   f_rho_u = 0.5_R8P*(F1(2) + F4(2) - lmax*(r4*u4                     - r1*u1                    ))
   f_rho_E = 0.5_R8P*(F1(3) + F4(3) - lmax*(r4*E(p=p4,r=r4,u=u4,g=g4) - r1*E(p=p1,r=r1,u=u1,g=g1)))
   endsubroutine solve_riemann

   attributes(device) subroutine compute_inter_states(r1, u1, p1, g1, r4, u4, p4, g4, u, S1, S4)
   !< Compute inter states (23*-states) from state1 and state4.
   real(R8P), intent(in)  :: r1             !< Density of state 1.
   real(R8P), intent(in)  :: u1             !< Velocity of state 1.
   real(R8P), intent(in)  :: p1             !< Pressure of state 1.
   real(R8P), intent(in)  :: g1             !< Specific heats ratio of state 1.
   real(R8P), intent(in)  :: r4             !< Density of state 4.
   real(R8P), intent(in)  :: u4             !< Velocity of state 4.
   real(R8P), intent(in)  :: p4             !< Pressure of state 4.
   real(R8P), intent(in)  :: g4             !< Specific heats ratio of state 4.
   real(R8P), intent(out) :: u              !< Velocity of the intermediate states.
   real(R8P), intent(out) :: S1             !< Maximum wave speed of state 1 and 4.
   real(R8P), intent(out) :: S4             !< Maximum wave speed of state 1 and 4.
   real(R8P)              :: p              !< Pressure of the intermediate states.
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

   attributes(device) function fluxes(p, r, u, g) result(Fc)
   !< 1D Euler fluxes from primitive variables.
   real(R8P), intent(in) :: p       !< Pressure.
   real(R8P), intent(in) :: r       !< Density.
   real(R8P), intent(in) :: u       !< Velocity.
   real(R8P), intent(in) :: g       !< Specific heats ratio.
   real(R8P)             :: Fc(1:3) !< State fluxes.

   Fc(1) = r*u
   Fc(2) = Fc(1)*u + p
   Fc(3) = Fc(1)*H(p=p, r=r, u=u, g=g)
   endfunction fluxes

   attributes(device) function p(r, a, g) result(pressure)
   !< Return pressure for an ideal calorically perfect gas.
   real(R8P), intent(in) :: r        !< Density.
   real(R8P), intent(in) :: a        !< Speed of sound.
   real(R8P), intent(in) :: g        !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
   real(R8P)             :: pressure !< Pressure.

   pressure = r*a*a/g
   endfunction p

   attributes(device) function r(p, a, g) result(density)
   !< Return density for an ideal calorically perfect gas.
   real(R8P), intent(in) :: p       !< Pressure.
   real(R8P), intent(in) :: a       !< Speed of sound.
   real(R8P), intent(in) :: g       !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
   real(R8P)             :: density !< Density.

   density = g*p/(a*a)
   endfunction r

   attributes(device) function a(p, r, g) result(ss)
   !< Return speed of sound for an ideal calorically perfect gas.
   real(R8P), intent(in) :: p  !< Pressure.
   real(R8P), intent(in) :: r  !< Density.
   real(R8P), intent(in) :: g  !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
   real(R8P)             :: ss !< Speed of sound.

   ss = sqrt(g*p/r)
   endfunction a

   attributes(host,device) function E(p, r, u, g) result(energy)
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

   attributes(device) function H(p, r, u, g) result(entalpy)
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

   !attributes(device) function weno_polynomials(s, weno_p, v) result(VP)
   !!< Return WENO polynomials
   !integer(I4P), intent(in) :: s                       !< Number of stencils used.
   !real(R8P),    intent(in) :: weno_p(1:2,0:S-1,0:S-1) !< Polinomials coefficients.
   !real(R8P),    intent(in) ::     v (1:2,1-s:-1+s)    !< Variable to be reconstructed.
   !real(R8P)                ::     vp(1:2,0:s-1   )    !< Polynomial reconstructions.
   !integer(I4P)             :: s1, s2, f               !< Counters.

   !vp = 0._R_P
   !do s1=0,s-1 ! stencil counter
   !   do s2=0,s-1 ! cell counter counter
   !      do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
   !         vp(f,s1) = vp(f,s1) + weno_p(f,s2,s1) * v(f,-s2+s1)
   !      enddo
   !   enddo
   !enddo
   !endfunction weno_polynomials

   !attributes(device) function weno_reconstructed(s, weno_p, v) result(vr)
   !!< Return WENO reconstruction of 2S-1 order.
   !integer(I4P), intent(in) :: s                       !< Number of stencils used.
   !real(R8P),    intent(in) :: weno_p(1:2,0:S-1,0:S-1) !< Polinomials coefficients.
   !real(R8P),    intent(in) ::     v (1:2,1-s:-1+s)    !< Variable to be reconstructed.
   !real(R8P)                ::     vr(1:2         )    !< Left and right (1,2) interface value of reconstructed V.
   !real(R8P)                ::     vp(1:2,0:s-1   )    !< Polynomial reconstructions.
   !real(R8P)                ::     w (1:2,0:s-1   )    !< Weights of the stencils.

   !vp = weno_polynomials(s=s, v=v)        ! compute the polynomials
   !w = weno_weights(s=s, v=v)             ! compute the weights associated to the polynomials
   !vr = weno_convolution(s=s, vp=vp, w=w) ! compute the convultion of reconstructing plynomials
   !endfunction weno_reconstructed

   !attributes(device) function weno_weights(s, weno_eps, weno_a, weno_d, v) result(w)
   !!< Return WENO weights of the polynomial reconstructions.
   !integer(I4P),                     intent(in) :: s                   !< Number of stencils used.
   !real(R8P),                        intent(in) :: v    (1:2,1-s:-1+s) !< Variable to be reconstructed.
   !real(R8P)                                    :: W    (1:2,0:s-1)    !< Weights of the stencils.
   !real(R8P)                                    :: IS   (1:2,0:s-1)    !< Smoothness indicators of the stencils.
   !real(R8P)                                    :: a    (1:2,0:s-1)    !< Alpha coefficients for the weights.
   !real(R8P)                                    :: a_tot(1:2)          !< Summ of the alpha coefficients.
   !integer(I4P)                                 :: s1, s2, s3, f       !< Counters.

   !! compute smoothness indicators
   !do s1=0,S-1 ! stencil counter
   !   do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
   !      IS(f,s1) = 0._R_P
   !      do s2=0,S-1
   !         do s3=0,S-1
   !            IS(f,s1) = IS(f,s1) + weno_d(s3,s2,s1) * v(f,s1-s3) * v(f,s1-s2)
   !         enddo
   !      enddo
   !   enddo
   !enddo
   !! compute alfa coefficients
   !a_tot = 0._R_P
   !do s1=0,S-1
   !   do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
   !      a(f,s1) = weno_a(f,s1) * (1._R_P / (weno_eps + IS(f,s1))**s)
   !      a_tot(f) = a_tot(f) + a(f,s1)
   !   enddo
   !enddo
   !! compute the weights
   !do s1=0,S-1
   !   do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
   !      w(f,s1) = a(f,s1) / a_tot(f)
   !   enddo
   !enddo
   !endfunction weno_weights
endmodule adam_equation_euler_gpu_object
