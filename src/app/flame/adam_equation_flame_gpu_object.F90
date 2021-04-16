!< ADAM, Euler equations system class definition, GPU backend.
module adam_equation_flame_gpu_object
!< ADAM, Euler equations system class definition, GPU backend.

use adam_adam_object
use adam_base_gpu_object
use adam_field_object
use adam_grid_object
use adam_parameters
use adam_weno_library_gpu
use FiNeR
use PENF
use MPI
use CUDAFOR
use cgal_wrappers
use ISO_C_BINDING
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none
private
public :: equation_flame_gpu_object
public :: BC_EXTRAPOLATION
public :: BC_INFLOW
public :: BC_INFLOW_COLD
public :: BC_NRIN_XMIN
public :: BC_NROUT_XMAX
public :: BC_WALL_ISOTERM
public :: BC_WALL_ADIAB
public :: BC_NROUT_XMIN

integer(I4P), parameter :: BC_EXTRAPOLATION = 1_I4P
integer(I4P), parameter :: BC_INFLOW        = 2_I4P
integer(I4P), parameter :: BC_INFLOW_COLD   = 3_I4P
integer(I4P), parameter :: BC_NRIN_XMIN     = 4_I4P
integer(I4P), parameter :: BC_NROUT_XMAX    = 5_I4P
integer(I4P), parameter :: BC_WALL_ISOTERM  = 6_I4P
integer(I4P), parameter :: BC_WALL_ADIAB    = 7_I4P
integer(I4P), parameter :: BC_NROUT_XMIN    = 8_I4P

type :: equation_flame_gpu_object
   !< Flame single-step class definition, GPU backend.
   !<
   !< Arrenhius form for reaction term.
   !<
   !< The conservative variables are arranged as follows:
   !<```
   !< q(1): rho
   !< q(2): rho * u
   !< q(3): rho * v
   !< q(4): rho * w
   !< q(5): rho * E
   !< q(6): rho * Ya, specific density of specie
   !<```
   !< q_aux(1): rho
   !< q_aux(2): u
   !< q_aux(3): v
   !< q_aux(4): w
   !< q_aux(5): ya
   !< q_aux(6): tem
   !< q_aux(7): pressure
   !< q_aux(8): entalpy
   !< q_aux(9): sound speed
   !<```
   type(adam_object),  pointer :: adam=>null()          !< ADAM.
   type(field_object), pointer :: field=>null()         !< The field.
   type(grid_object),  pointer :: grid=>null()          !< The grid.
   integer(I4P),       pointer :: ngc=>null()           !< Number of ghost cells.
   integer(I4P),       pointer :: ni=>null()            !< Number of cells in i direction.
   integer(I4P),       pointer :: nj=>null()            !< Number of cells in j direction.
   integer(I4P),       pointer :: nk=>null()            !< Number of cells in k direction.
   integer(I4P),       pointer :: nb=>null()            !< Total blocks number for MPI.
   integer(I4P),       pointer :: blocks_number=>null() !< Actual blocks number.
   integer(I4P),       pointer :: nv=>null()            !< Number of variables.
   type(base_gpu_object)       :: base_gpu              !< The base GPU handler.
   integer(I4P)                :: myrank=0_I4P          !< MPI rank process.
   integer(I4P)                :: procs_number=1_I4P    !< Number of MPI processes.
   integer(I4P)                :: error=0_I4P           !< Error traping flag.
   ! equation data
   character(32)          :: solids(1)=["sd7003.off"]   !< Number of fluid species.
   type(c_ptr), allocatable :: ptree(:)
   character(32)          :: flow_type="flame1d"     !< Number of fluid species.
   real(R8P),    allocatable :: phi(:,:,:,:,:)       !< Fluxes.
   real(R8P),    allocatable, device :: phi_gpu(:,:,:,:,:)    !< Fluxes.
   integer(I4P)           :: ns=2_I4P                !< Number of fluid species.
   integer(I4P)           :: iweno=2_I4P             !< WENO order.
   integer(I4P)           :: lmax=2_I4P              !< Central convective half stencil.
   integer(I4P)           :: ivis=4_I4P              !< Diffusive terms order.
   integer(I4P)           :: visc_type=0_I4P         !< Diffusivity type (0=constant, 1=power, 2=Sutherland)
   real(R8P), allocatable :: x_cell_t(:,:)           !< First order derivatives coeffs.
   real(R8P), allocatable :: y_cell_t(:,:)           !< First order derivatives coeffs.
   real(R8P), allocatable :: z_cell_t(:,:)           !< First order derivatives coeffs.
   real(R8P), allocatable :: fd_coeff1(:)            !< First order derivatives coeffs.
   real(R8P), allocatable :: fd_coeff2(:)            !< Second order derivatives coeffs.
   real(R8P), allocatable :: fd_conv(:,:)            !< Second order derivatives coeffs.
   real(R8P)              :: Prandtl=0.74_R8P        !< Prandtl number
   real(R8P)              :: q_coeff=0.74_R8P        !< Prandtl number
   real(R8P)              :: flame_center=10._R8P    !< Flame stabilization abscissa.
   real(R8P)              :: flame_x_min=10._R8P     !< Flame stabilization abscissa.
   real(R8P)              :: flame_x_min_old=10._R8P !< Flame stabilization abscissa.
   real(R8P)              :: flame_x_velrel=0._R8P   !< Flame stabilization abscissa.
   real(R8P)              :: flame_x_vel=0._R8P      !< Flame stabilization abscissa.
   real(R8P)              :: Lewis=1._R8P            !< Lewis number.
   real(R8P)              :: Zeldovich=1060._R8P     !< Zeldovich number.
   real(R8P)              :: Damkohler=1800._R8P     !< Damkohler number.
   real(R8P)              :: gamma_fluid=1.32_R8P    !< Gamma.
   real(R8P)              :: dha=10000._R8P          !< Entalpy formation.
   real(R8P)              :: pres_inflow=10000._R8P  !< Inlet and initial pressure.
   real(R8P)              :: u_inflow=0.00524_R8P    !< Inlet u.
   real(R8P)              :: ya_inflow=0.054861_R8P  !< Inlet mass fraction.
   real(R8P)              :: tem_inflow=86.056_R8P   !< Inlet temperature
   real(R8P)              :: tem_outflow=870._R8P    !< Outlet temperature
   real(R8P)              :: tem_stabil=250._R8P     !< Outlet temperature

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
   ! cuf data
   real(R8P),    allocatable, device :: dq_gpu(:,:,:,:,:)    !<
   real(R8P),    allocatable, device :: fl_gpu(:,:,:,:,:)    !< Fluxes.
   real(R8P),    allocatable, device :: fhat_gpu(:,:,:,:,:)  !< Auxiliary fluxes.
   real(R8P),    allocatable, device :: dxyz_gpu(:,:)        !< Space steps.
   real(R8P),    allocatable, device :: alph_gpu(:,:)        !< RK alpha coefficients.
   real(R8P),    allocatable, device :: beta_gpu(:)          !< RK beta coefficients.
   real(R8P),    allocatable, device :: gamm_gpu(:)          !< RK gamma coefficients.
   real(R8P),    allocatable, device :: fd_coeff1_gpu(:)     !< First order derivatives coeffs.
   real(R8P),    allocatable, device :: fd_coeff2_gpu(:)     !< Second order derivatives coeffs.
   real(R8P),    allocatable, device :: fd_conv_gpu(:,:)     !< Second order derivatives coeffs.
   real(R8P),    allocatable, device :: x_cell_gpu(:,:)      !< First order derivatives coeffs.
   real(R8P),    allocatable, device :: y_cell_gpu(:,:)      !< First order derivatives coeffs.
   real(R8P),    allocatable, device :: z_cell_gpu(:,:)      !< First order derivatives coeffs.
   real(R8P),    allocatable, device :: q_aux_gpu(:,:,:,:,:) !< Auxiliary cell centered variables.
   real(R8P),    allocatable, device :: q_gpu(:,:,:,:,:)     !< Field cell centered variables stages.
   real(R8P),    allocatable, device :: q_s_gpu(:,:,:,:,:,:) !< RK Field cell centered variables stages.
   real(R8P),    allocatable, device :: gplus_x(:,:,:,:,:)   !< For weno-x
   real(R8P),    allocatable, device :: gminus_x(:,:,:,:,:)  !< For weno-x
   real(R8P),    allocatable, device :: gplus_y(:,:,:,:,:)   !< For weno-y
   real(R8P),    allocatable, device :: gminus_y(:,:,:,:,:)  !< For weno-y
   real(R8P),    allocatable, device :: gplus_z(:,:,:,:,:)   !< For weno-z
   real(R8P),    allocatable, device :: gminus_z(:,:,:,:,:)  !< For weno-z
   contains
      ! public methods
      procedure, pass(self) :: amr_update              !< Do AMR update.
      procedure, pass(self) :: compute_aux             !< Compute auxiliary variables.
      procedure, pass(self) :: compute_dt              !< Compute time step.
      procedure, pass(self) :: copy_cpu_gpu            !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_gpu_cpu            !< Copy data from GPU to CPU.
      procedure, pass(self) :: destroy                 !< Destroy the equation.
      procedure, pass(self) :: fd_initialize           !< Initialize Finite Difference Coefficients.
      procedure, pass(self) :: initialize              !< Initialize the equation.
      procedure, pass(self) :: mark_by_grad_rho        !< Mark blocks to be refined/derefined by a `grad(rho)` value.
      procedure, pass(self) :: mark_by_geo             !< Mark blocks to be refined/derefined by a `grad(rho)` value.
      procedure, pass(self) :: integrate               !< Runge Kutta integration of equation.
      procedure, pass(self) :: print_progress          !< Print simulation progress.
      procedure, pass(self) :: refine_uniform          !< Refine all blocks uniformly.
      procedure, pass(self) :: update_cell_gpu         !< Refine all blocks uniformly.
      procedure, pass(self) :: update_phi              !< Refine all blocks uniformly.
      procedure, pass(self) :: runge_kutta_initialize  !< Initialize Runge-Kutta data.
      procedure, pass(self) :: save_hdf5               !< Save simulation data in HDF5 format.
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions of equation.
      procedure, pass(self) :: set_bc_rhs              !< Set boundary conditions for rhs.
      procedure, pass(self) :: set_initial_conditions  !< Set initial conditions of equation.
      procedure, pass(self) :: update_ghost_gpu        !< Update ghost cells and set boundary conditions.
      ! operators
      generic :: assignment(=) => eq_assign_eq      !< Overload `=`.
      procedure, pass(lhs), private :: eq_assign_eq !< Operator `=`.
endtype equation_flame_gpu_object

contains
   ! public methods
   subroutine amr_update(self, iterations)
   !< Do AMR update.
   class(equation_flame_gpu_object), intent(inout)        :: self            !< The equation.
   integer(I4P),                     intent(in), optional :: iterations      !< Number of AMR iterations.
   integer(I4P)                                           :: iterations_     !< Number of AMR iterations, local var.
   logical                                                :: is_grid_changed !< Flag to check grid changes.
   integer(I4P)                                           :: i               !< Counter.
   integer(I4P)                                           :: b, j, k         !< Counter.

   iterations_ = 1 ; if (present(iterations)) iterations_ = iterations
   amr: do i=1, iterations_
      if(self%flow_type == "cold") then
         !call self%update_phi()
         call self%mark_by_geo(tol=2000._R8P, delta_fine=0.015_R8P, delta_coarse=0.15_R8P)
      else
          call self%mark_by_grad_rho(grad_tol=2._R8P, delta_fine=0.015_R8P, delta_coarse=0.15_R8P)
          !call self%mark_by_grad_rho(grad_tol=0.05_R8P, delta_fine=0.006_R8P, delta_coarse=0.015_R8P)
      endif
      call self%update_ghost_gpu(q_gpu=self%q_gpu)
      call self%copy_gpu_cpu()
      call self%adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false., is_grid_changed=is_grid_changed)
      if(self%flow_type == "cold") then
         call self%update_phi()
      endif
      call self%copy_cpu_gpu
      if (.not.is_grid_changed) then
          print*,'AMR Grid stabilized after : ',i,' AMR iterations'
          exit amr
      endif
   enddo amr

   print*,'START - prima update_cell_gpu'
   call self%update_cell_gpu()
   print*,'END - prima update_cell_gpu'

   endsubroutine amr_update

   subroutine update_cell_gpu(self)
   !< Update x/y/z_cell_gpu
   class(equation_flame_gpu_object), intent(inout)        :: self       !< The equation.
   integer(I4P)                                           :: b, i, j, k !< Counter.

   associate(blocks_number=>self%blocks_number, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc)
   if(allocated(self%x_cell_t))   deallocate(self%x_cell_t)
   if(allocated(self%y_cell_t))   deallocate(self%y_cell_t)
   if(allocated(self%z_cell_t))   deallocate(self%z_cell_t)
   if(allocated(self%x_cell_gpu)) deallocate(self%x_cell_gpu)
   if(allocated(self%y_cell_gpu)) deallocate(self%y_cell_gpu)
   if(allocated(self%z_cell_gpu)) deallocate(self%z_cell_gpu)
   if(blocks_number > 0) then
      allocate(self%x_cell_t(blocks_number, 1-ngc:ni+ngc),   self%y_cell_t(blocks_number, 1-ngc:nj+ngc),   &
          self%z_cell_t(blocks_number, 1-ngc:nk+ngc))
      allocate(self%x_cell_gpu(blocks_number, 1-ngc:ni+ngc), self%y_cell_gpu(blocks_number, 1-ngc:nj+ngc), &
          self%z_cell_gpu(blocks_number, 1-ngc:nk+ngc))
      do b=1,blocks_number
          do i=1-ngc,ni+ngc
             self%x_cell_t(b,i) = self%field%x_cell(i,b)
          enddo
          do j=1-ngc,nj+ngc
             self%y_cell_t(b,j) = self%field%y_cell(j,b)
          enddo
          do k=1-ngc,nk+ngc
             self%z_cell_t(b,k) = self%field%z_cell(k,b)
          enddo
      enddo
      self%x_cell_gpu = self%x_cell_t
      self%y_cell_gpu = self%y_cell_t
      self%z_cell_gpu = self%z_cell_t
   endif
   endassociate
   !print*,'amr update 1', lbound(self%x_cell_t,1),   ubound(self%x_cell_t,1),   size(self%x_cell_t,1)
   !print*,'amr update 2', lbound(self%x_cell_t,2),   ubound(self%x_cell_t,2),   size(self%x_cell_t,2)
   !print*,'amr update 1', lbound(self%x_cell_gpu,1), ubound(self%x_cell_gpu,1), size(self%x_cell_gpu,1)
   !print*,'amr update 2', lbound(self%x_cell_gpu,2), ubound(self%x_cell_gpu,2), size(self%x_cell_gpu,2)
   endsubroutine update_cell_gpu

   subroutine update_phi(self)
   !< Update x/y/z_cell_gpu
   class(equation_flame_gpu_object), intent(inout) :: self                      !< The equation.
   integer(I4P)                                    :: b, i, j, k, ib            !< Counter.
   real(R8P)                                       :: query_x, query_y, query_z !< Counter.
   real(R8P)                                       :: near_x, near_y, near_z    !< Counter.
   real(R8P)                                       :: distance                  !< Counter.
   logical                                         :: inside

   associate(blocks_number=>self%blocks_number, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
             x_cell=>self%field%x_cell, y_cell=>self%field%y_cell, z_cell=>self%field%z_cell,         &
             ptree => self%ptree, phi=>self%phi, phi_gpu=>self%phi_gpu, n_solids=>size(self%solids,dim=1))
   print*,'updating distances - start'
   do ib=1,n_solids
      do b=1,blocks_number
         print*,'block: ',b,' / ',blocks_number
         do i=1-ngc,ni+ngc
         do j=1-ngc,nj+ngc
         do k=1-ngc,nk+ngc
            query_x = x_cell(i,b)
            query_y = y_cell(j,b)
            query_z = z_cell(k,b)
            call polyhedron_closest(ptree(ib),query_x,query_y,query_z,near_x,near_y,near_z)
            distance = sqrt((near_x-query_x)**2+(near_y-query_y)**2+(near_z-query_z)**2)
            inside   = cgal_polyhedron_inside(ptree(ib),query_x,query_y,query_z)
            if(.not.inside) distance = - distance
            phi(b,i,j,k,ib) = distance
         enddo
         enddo
         enddo
      enddo
   enddo
   print*,'updating distances - end'
   phi_gpu = phi
   endassociate
   endsubroutine update_phi

   subroutine update_cell_gpu_associatefail(self)
   !< Update x/y/z_cell_gpu
   class(equation_flame_gpu_object), intent(inout)        :: self       !< The equation.
   integer(I4P)                                           :: b, i, j, k !< Counter.

   associate(blocks_number=>self%blocks_number, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
             x_cell=>self%field%x_cell, y_cell=>self%field%y_cell, z_cell=>self%field%z_cell, &
             x_cell_t=>self%x_cell_t, y_cell_t=>self%y_cell_t, z_cell_t=>self%z_cell_t, &
             x_cell_gpu=>self%x_cell_gpu, y_cell_gpu=>self%y_cell_gpu, z_cell_gpu=>self%z_cell_gpu)
   if(allocated(x_cell_t)) deallocate(x_cell_t)
   if(allocated(y_cell_t)) deallocate(y_cell_t)
   if(allocated(z_cell_t)) deallocate(z_cell_t)
   if(allocated(x_cell_gpu)) deallocate(x_cell_gpu)
   if(allocated(y_cell_gpu)) deallocate(y_cell_gpu)
   if(allocated(z_cell_gpu)) deallocate(z_cell_gpu)
   allocate(x_cell_t(blocks_number, 1-ngc:ni+ngc), y_cell_t(blocks_number, 1-ngc:nj+ngc), z_cell_t(blocks_number, 1-ngc:nk+ngc))
   allocate(x_cell_gpu(blocks_number, 1-ngc:ni+ngc), y_cell_gpu(blocks_number, 1-ngc:nj+ngc), z_cell_gpu(blocks_number, 1-ngc:nk+ngc))
   do b=1,blocks_number
       do i=1-ngc,ni+ngc
          x_cell_t(b,i) = x_cell(i,b)
       enddo
       do j=1-ngc,nj+ngc
          y_cell_t(b,j) = y_cell(j,b)
       enddo
       do k=1-ngc,nk+ngc
          z_cell_t(b,k) = z_cell(k,b)
       enddo
   enddo
   x_cell_gpu = x_cell_t
   y_cell_gpu = y_cell_t
   z_cell_gpu = z_cell_t
   !print*,'amr update 1', lbound(x_cell_t,1),   ubound(x_cell_t,1),   size(x_cell_t,1)
   !print*,'amr update 2', lbound(x_cell_t,2),   ubound(x_cell_t,2),   size(x_cell_t,2)
   !print*,'amr update 1', lbound(x_cell_gpu,1), ubound(x_cell_gpu,1), size(x_cell_gpu,1)
   !print*,'amr update 2', lbound(x_cell_gpu,2), ubound(x_cell_gpu,2), size(x_cell_gpu,2)
   endassociate
   !print*,'amr update 1', lbound(self%x_cell_t,1),   ubound(self%x_cell_t,1),   size(self%x_cell_t,1)
   !print*,'amr update 2', lbound(self%x_cell_t,2),   ubound(self%x_cell_t,2),   size(self%x_cell_t,2)
   !print*,'amr update 1', lbound(self%x_cell_gpu,1), ubound(self%x_cell_gpu,1), size(self%x_cell_gpu,1)
   !print*,'amr update 2', lbound(self%x_cell_gpu,2), ubound(self%x_cell_gpu,2), size(self%x_cell_gpu,2)
   endsubroutine update_cell_gpu_associatefail

   subroutine compute_aux(self, q_gpu, q_aux_gpu)
   !< Compute auxiliary variables.
   class(equation_flame_gpu_object), intent(in)          :: self          !< The equation.
   real(R8P),                        intent(in),  device :: q_gpu(1:,                    &
                                                                  1-self%ngc:,&
                                                                  1-self%ngc:,&
                                                                  1-self%ngc:,&
                                                                  1:)     !< Conservative variables.
   real(R8P),                        intent(out), device :: q_aux_gpu(1:,                    &
                                                                      1-self%ngc:,&
                                                                      1-self%ngc:,&
                                                                      1-self%ngc:,&
                                                                      1:) !< Auxiliary variables.

   associate(blocks_number=>self%field%blocks_number, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, ns=>self%ns)
      call compute_aux_cuf(ni=ni, nj=nj, nk=nk, ngc=ngc, ns=ns, blocks_number=blocks_number, &
                           gamma_fluid=self%gamma_fluid, dha=self%dha, q_gpu=q_gpu, q_aux_gpu=q_aux_gpu)
   endassociate
   endsubroutine compute_aux

   subroutine compute_dt(self)
   !< Compute maximum time step accordingly to CFL stabilty criterion.
   class(equation_flame_gpu_object), intent(inout) :: self !< The equation.
   real(R8P)                                       :: umax !< Maximum speed of waves propagation.
   integer(I4P)                                    :: b    !< Counter.

   associate(blocks_number=>self%field%blocks_number, dxyz=>self%field%dxyz, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             ngc=>self%ngc, ns=>self%ns, q=>self%field%q, dt=>self%dt, CFL=>self%CFL)
      call self%compute_aux(q_gpu=self%q_gpu, q_aux_gpu=self%q_aux_gpu)
      dt = huge(1._R8P)
      do b=1, self%field%blocks_number
         call compute_umax_cuf(b, ni=ni, nj=nj, nk=nk, ngc=ngc, ns=ns,   &
                               dx=dxyz(1,b), dy=dxyz(2,b), dz=dxyz(3,b), &
                               q_aux_gpu=self%q_aux_gpu, umax=umax)
         dt = min(dt, minval(dxyz(:,b)) / umax * CFL)
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, self%error)
   endassociate
   endsubroutine compute_dt

   subroutine copy_cpu_gpu(self)
   !< Copy data from CPU to GPU.
   class(equation_flame_gpu_object), intent(inout) :: self        !< The base backend.
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
   class(equation_flame_gpu_object), intent(inout)        :: self          !< The base backend.
   logical,                          intent(in), optional :: compute_q_aux !< Flag to compute auxiliary variables.

   call self%base_gpu%copy_transpose_gpu_cpu(nv=self%nv, q_gpu=self%q_gpu, q_cpu=self%field%q)
   if (present(compute_q_aux)) then
      if (compute_q_aux) then
         call self%compute_aux(q_gpu=self%q_gpu, q_aux_gpu=self%q_aux_gpu)
         call self%base_gpu%copy_transpose_gpu_cpu(nv=self%ns+7, q_gpu=self%q_aux_gpu, q_cpu=self%q_aux)
      endif
   endif
   endsubroutine copy_gpu_cpu

   subroutine destroy(self)
   !< Destroy the equation.
   class(equation_flame_gpu_object), intent(inout) :: self  !< The equation.
   type(equation_flame_gpu_object)                 :: fresh !< Fresh equation.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, adam, ns, nrk, CFL, null_xyz, fields_gpu_number, flow_type)
   !< Initialize the equation.
   class(equation_flame_gpu_object), intent(inout)        :: self              !< The equation.
   type(adam_object),                intent(in), target   :: adam              !< ADAM.
   integer(I4P),                     intent(in), optional :: ns                !< Species number.
   integer(I4P),                     intent(in), optional :: nrk               !< Runge-Kutta stages number.
   real(R8P),                        intent(in), optional :: CFL               !< CFL value.
   logical,                          intent(in), optional :: null_xyz(3)       !< Flag triggering 1D/2D simulations.
   integer(I4P),                     intent(in), optional :: fields_gpu_number !< Number of fields allocated on GPU.
   integer(I4P)                                           :: v, i              !< Counter.
   character(32), optional                                :: flow_type
   integer(I4P)                                           :: n_solids          !< Counter.

   ! CPU data
   call self%destroy
   self%adam          => adam
   self%field         => adam%field
   self%grid          => adam%field%grid
   self%ngc           => adam%field%grid%ngc
   self%ni            => adam%field%grid%ni
   self%nj            => adam%field%grid%nj
   self%nk            => adam%field%grid%nk
   self%nb            => adam%field%nb
   self%blocks_number => adam%field%blocks_number
   self%nv            => adam%field%nv
   if (present(flow_type)) self%flow_type = flow_type
   if(trim(self%flow_type) == "cold") then
       print*,"Setting cold case"
       self%ya_inflow   = 0._R8P
       self%dha         = 1._R8P ! otherwise WENO fails because there is a division to dha
       self%Damkohler   = 0._R8P
       self%pres_inflow = 10000._R8P
       self%u_inflow    = 1.5_R8P
       self%tem_inflow  = 300._R8P
       self%tem_outflow = 300._R8P
       !self%q_coeff   = self%Prandtl
       !self%Prandtl   = 1._R8P/100._R8P
       n_solids = size(self%solids,dim=1)
       allocate(self%ptree(n_solids))
       do i=1,n_solids
          print*,'reading solid: ', self%solids(i)
          call cgal_polyhedron_read(self%ptree(i), self%solids(i))
       enddo
   endif
   if(trim(self%flow_type) == "flamechannel") then
       print*,"Setting flame channel case"
       self%u_inflow    = 0._R8P
       self%tem_inflow  = 86._R8P
       self%tem_outflow = 86._R8P
       !self%q_coeff   = self%Prandtl
       !self%Prandtl   = 1._R8P/100._R8P
   endif
   if(trim(self%flow_type) == "sod") then
       print*,"Setting SOD channel case"
       self%gamma_fluid = 1.4_R8P
       self%Prandtl = 0.
       self%iweno = 1
   endif
   if (present(ns)) self%ns = ns
   if (self%nv - self%ns /= 4) then
      write(stderr, '(A)') 'ADAM-ERROR: field%nv must be flame%ns+4'
      call MPI_FINALIZE(self%error)
      stop
   endif
   call self%base_gpu%initialize(field=self%field, nv_aux=self%ns+7, fields_gpu_number=fields_gpu_number)
   if (present(nrk)) self%nrk = nrk
   if (present(CFL)) self%CFL = CFL
   if (present(null_xyz)) self%null_xyz = null_xyz

   call self%runge_kutta_initialize
   call self%fd_initialize
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   ! allocate large array
   associate(nv=>self%nv, ns=>self%ns, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             nb=>self%nb, nrk=>self%nrk, iweno=>self%iweno)
   ! CPU data
   allocate(self%q_aux(1:ns+7,       1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
   ! GPU data
   allocate(self%fl_gpu(1:nb,        1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv))
   allocate(self%fhat_gpu(1:nb,      1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv))
   allocate(self%q_aux_gpu(1:nb,     1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:ns+7))
   allocate(self%q_gpu(1:nb,         1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv))
   allocate(self%dq_gpu(1:nb,        1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv))
   allocate(self%q_s_gpu(1:nb,       1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nv, 1:nrk))
   if(trim(self%flow_type) == "cold") then
       allocate(self%phi(1:nb,        1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:n_solids))
       allocate(self%phi_gpu(1:nb,    1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:n_solids))
   endif
   allocate(self%gplus_x (6, 2*iweno, nj, nk, nb))
   allocate(self%gminus_x(6, 2*iweno, nj, nk, nb))
   allocate(self%gplus_y (6, 2*iweno, ni, nk, nb))
   allocate(self%gminus_y(6, 2*iweno, ni, nk, nb))
   allocate(self%gplus_z (6, 2*iweno, ni, nj, nb))
   allocate(self%gminus_z(6, 2*iweno, ni, nj, nb))
   allocate(self%dxyz_gpu(1:nb, 1:3))
   endassociate
   ! copy data that is not variable during the simulation
   self%fd_coeff1_gpu = self%fd_coeff1
   self%fd_coeff2_gpu = self%fd_coeff2
   self%fd_conv_gpu   = self%fd_conv
   self%alph_gpu = self%alph
   self%beta_gpu = self%beta
   self%gamm_gpu = self%gamm
   endsubroutine initialize

   subroutine mark_by_grad_rho(self, grad_tol, delta_fine, delta_coarse, threshold)
   !< Mark blocks to be refined/derefined by a `grad(rho)` value.
   class(equation_flame_gpu_object), intent(inout)        :: self           !< The equation.
   real(R8P),                        intent(in)           :: grad_tol       !< Gradiend tolerance value.
   real(R8P),                        intent(in)           :: delta_fine     !< Maximum cell delta in fine grids.
   real(R8P),                        intent(in)           :: delta_coarse   !< Minimum cell delta in coarse grids.
   real(R8P),                        intent(in), optional :: threshold      !< Threshold for sphere proximity.
   real(R8P)                                              :: threshold_     !< Threshold for sphere proximity, local var.
   real(R8P)                                              :: max_cell_delta !< Maximum cell delta.
   real(R8P)                                              :: grad_rho       !< Value (max) of gradient of rho.
   integer(I4P)                                           :: b              !< Counter.

   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   self%field%refinements_needed = [(TO_NOT_TOUCH,b=1,self%blocks_number)]
   call self%update_ghost_gpu(q_gpu=self%q_gpu)
   associate (ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
              blocks_number=>self%blocks_number, ns=>self%ns, dxyz=>self%field%dxyz)
      call self%compute_aux(q_gpu=self%q_gpu, q_aux_gpu=self%q_aux_gpu)
      do b=1, blocks_number
         grad_rho = gradient_cuf(b=b, ni=ni, nj=nj, nk=nk, ngc=ngc, ns=ns, &
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
      function gradient_cuf(b, ni, nj, nk, ngc, ns, dx, dy, dz, q_gpu) result(gradient)
      !< Gradient done by CUF threads.
      integer(I4P), intent(in)         :: b                                 !< Block index.
      integer(I4P), intent(in)         :: ni                                !< Grid cells number in I direction.
      integer(I4P), intent(in)         :: nj                                !< Grid cells number in J direction.
      integer(I4P), intent(in)         :: nk                                !< Grid cells number in K direction.
      integer(I4P), intent(in)         :: ngc                               !< Ghost cells number.
      integer(I4P), intent(in)         :: ns                                !< Species number.
      real(R8P),    intent(in)         :: dx                                !< X space step.
      real(R8P),    intent(in)         :: dy                                !< Y space step.
      real(R8P),    intent(in)         :: dz                                !< Z space step.
      real(R8P),    intent(in), device :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field component to which apply gradient.
      real(R8P)                        :: gradient                          !< Maximum gradient of q.
      real(R8P)                        :: grad                              !< Current gradient of q.
      integer(I4P)                     :: i, j, k                           !< Counter.
      integer(I4P)                     :: iercuda                           !< Error trapping flag for CUDAFortran.

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
   endsubroutine mark_by_grad_rho

   subroutine mark_by_geo(self, tol, delta_fine, delta_coarse, threshold)
   !< Mark blocks to be refined/derefined by a `grad(rho)` value.
   class(equation_flame_gpu_object), intent(inout)        :: self           !< The equation.
   real(R8P),                        intent(in)           :: tol            !< Gradiend tolerance value.
   real(R8P),                        intent(in)           :: delta_fine     !< Maximum cell delta in fine grids.
   real(R8P),                        intent(in)           :: delta_coarse   !< Minimum cell delta in coarse grids.
   real(R8P),                        intent(in), optional :: threshold      !< Threshold for sphere proximity.
   real(R8P)                                              :: threshold_     !< Threshold for sphere proximity, local var.
   real(R8P)                                              :: max_cell_delta !< Maximum cell delta.
   real(R8P)                                              :: distance       !< Value (max) of gradient of rho.
   integer(I4P)                                           :: b              !< Counter.

   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   self%field%refinements_needed = [(TO_NOT_TOUCH,b=1,self%blocks_number)]
   associate (ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
              blocks_number=>self%blocks_number, ns=>self%ns, dxyz=>self%field%dxyz, phi=>self%phi)
      do b=1, blocks_number
         distance = 1._R8P
         if(maxval(phi(b,:,:,:,1))*minval(phi(b,:,:,:,1)) < 0._R8P) then
            distance = 0._R8P
         endif
         max_cell_delta = max_cell_delta_dist(distance=distance)

         if (maxval(dxyz(:,b)) > max_cell_delta) then
            self%field%refinements_needed(b) = TO_BE_REFINED
         elseif (maxval(dxyz(:,b)) * threshold_ < max_cell_delta) then
            self%field%refinements_needed(b) = TO_NOT_TOUCH ! TO_BE_DEREFINED
         else
            self%field%refinements_needed(b) = TO_NOT_TOUCH
         endif
      enddo
   endassociate
   contains
      function max_cell_delta_dist(distance) result(delta)
      !< Return the maximum cell delta given a comparison distance.
      real(R8P),          intent(in) :: distance !< Comparison distance.
      real(R8P)                      :: delta    !< Maximum cell delta admissible.

      if (abs(distance) < epsilon(0._R8P)) then
         ! delta = 0.001_R8P
         delta = 0.005_R8P
      else
         delta = huge(0._R8P)
      endif
      endfunction max_cell_delta_dist
   endsubroutine mark_by_geo

   subroutine integrate(self, t, do_ghost_syncro, residual)
   !< Runge Kutta integration of field.
   class(equation_flame_gpu_object), intent(inout)         :: self             !< The equation.
   real(R8P),                        intent(in)            :: t                !< Time.
   logical,                          intent(in),  optional :: do_ghost_syncro  !< Flag to do syncrous ghost update.
   real(R8P),                        intent(out), optional :: residual         !< Global residual.
   logical                                                 :: do_ghost_syncro_ !< Flag to do syncrous ghost update, local var.
   integer(I4P)                                            :: s                !< Counter.

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   associate(alph=>self%alph, beta=>self%beta, gamm=>self%gamm, dt=>self%dt, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             ngc=>self%ngc, nv=>self%nv, nrk=>self%nrk, ns=>self%ns, blocks_number=>self%blocks_number,             &
             inner_blocks_number=>self%field%inner_blocks_number,                                                   &
             alph_gpu=>self%alph_gpu, beta_gpu=>self%beta_gpu)
   do s=1, nrk
      if(trim(self%flow_type) == "cold") then
         call minimal_immersed_bc(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number,                 &
            gamma_fluid=self%gamma_fluid, q_gpu=self%q_gpu(:,:,:,:,:), phi_gpu=self%phi_gpu,                        &
            x_cell_gpu=self%x_cell_gpu, y_cell_gpu=self%y_cell_gpu, z_cell_gpu=self%z_cell_gpu)

         call evolve_eikonal_q_gpu_cuf(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                       phi_gpu=self%phi_gpu,                                             &
                                       dx_gpu=self%dxyz_gpu(:,1),                                        &
                                       dy_gpu=self%dxyz_gpu(:,2),                                        &
                                       dz_gpu=self%dxyz_gpu(:,3),                                        &
                                       dq_gpu=self%dq_gpu,                                               &
                                       q_gpu=self%q_gpu)
      endif
      call compute_rk_stage_gpu_cuf(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                    alph_gpu=alph_gpu, dt=dt, s=s, q_gpu=self%q_gpu, q_s_gpu=self%q_s_gpu)

      if (do_ghost_syncro_) then
         call self%update_ghost_gpu(q_gpu=self%q_s_gpu(:,:,:,:,:,s)) ! all ghosts
         if(trim(self%flow_type) == "cold") then
            call minimal_immersed_bc(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number,              &
               gamma_fluid=self%gamma_fluid, q_gpu=self%q_s_gpu(:,:,:,:,:,s), phi_gpu=self%phi_gpu,                 &
               x_cell_gpu=self%x_cell_gpu, y_cell_gpu=self%y_cell_gpu, z_cell_gpu=self%z_cell_gpu)
         endif
         call self%compute_aux(q_gpu=self%q_s_gpu(:,:,:,:,:,s), q_aux_gpu=self%q_aux_gpu)
         ! in the next q_s_gpu(...,s) is the flux (it was the rk-stage so far)
         call compute_residuals_gpu(ni=ni, nj=nj, nk=nk, ngc=ngc, ns=ns, blocks_number=blocks_number,          &
                                    null_x=self%null_xyz(1), null_y=self%null_xyz(2), null_z=self%null_xyz(3), &
                                    dx_gpu        = self%dxyz_gpu(:,1),                                        &
                                    dy_gpu        = self%dxyz_gpu(:,2),                                        &
                                    dz_gpu        = self%dxyz_gpu(:,3),                                        &
                                    q_aux_gpu     = self%q_aux_gpu,                                            &
                                    fl_gpu        = self%fl_gpu,                                               &
                                    fhat_gpu      = self%fhat_gpu,                                             &
                                    q_gpu         = self%q_s_gpu(:,:,:,:,:,s),                                 &
                                    iweno         = self%iweno,                                                &
                                    lmax          = self%lmax,                                                 &
                                    dha           = self%dha,                                                  &
                                    gamma_fluid   = self%gamma_fluid,                                          &
                                    gplus_x       = self%gplus_x,                                              &
                                    gminus_x      = self%gminus_x,                                             &
                                    Prandtl       = self%Prandtl,                                              &
                                    q_coeff       = self%q_coeff,                                              &
                                    Lewis         = self%Lewis,                                                &
                                    Zeldovich     = self%Zeldovich,                                            &
                                    Damkohler     = self%Damkohler,                                            &
                                    ivis          = self%ivis,                                                 &
                                    visc_type     = self%visc_type,                                            &
                                    fd_conv_gpu   = self%fd_conv_gpu,                                          &
                                    fd_coeff1_gpu = self%fd_coeff1_gpu,                                        &
                                    fd_coeff2_gpu = self%fd_coeff2_gpu)
          if(trim(self%flow_type) /= "cold") call self%set_bc_rhs(q_gpu=self%q_s_gpu(:,:,:,:,:,s), q_aux_gpu=self%q_aux_gpu)
      else
         ! TODO
      endif
      if (present(residual).and.s==self%nrk) then
         ! TODO
      endif
   enddo
   call advance_q_gpu_cuf(ni=ni, nj=nj, nk=nk, ngc=ngc, nrk=nrk, nv=nv, blocks_number=blocks_number, &
                           beta_gpu=beta_gpu, dt=dt, q_s_gpu=self%q_s_gpu, q_gpu=self%q_gpu)

   if(trim(self%flow_type) == "flame1d") then
      call self%compute_aux(q_gpu=self%q_gpu(:,:,:,:,:), q_aux_gpu=self%q_aux_gpu)
      call flame_find_x_v_cuf(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, dt=dt, tem_stabil=self%tem_stabil, &
                              x_cell_gpu=self%x_cell_gpu, y_cell_gpu=self%y_cell_gpu, z_cell_gpu=self%z_cell_gpu, &
                              q_aux_gpu=self%q_aux_gpu, x_min=self%flame_x_min, x_min_old=self%flame_x_min_old, &
                              velrel=self%flame_x_velrel)
   endif
   endassociate
   endsubroutine integrate

   subroutine print_progress(self, t, time, time_max)
   !< Print simulation progress.
   class(equation_flame_gpu_object), intent(in) :: self     !< The equation.
   integer(I4P),                     intent(in) :: t        !< Time iteration.
   real(R8P),                        intent(in) :: time     !< Time.
   real(R8P)                                    :: time_max !< Maximum time of integration.

   print '(A)', 'blocks number: '//trim(str(self%adam%tree%nodes_number, .true.))
   print '(A)', 'time step:     '//trim(str(self%dt, .true.))
   print '(A)', 'time:          '//trim(str(time, .true.))
   print '(A)', 't:             '//trim(str(t,.true.))
   print '(A)', 'progress:      '//trim(str(int(time/time_max * 100), .true.))//'%'
   endsubroutine print_progress

   subroutine refine_uniform(self, refinement_levels)
   !< Refine all blocks uniformly.
   class(equation_flame_gpu_object), intent(inout) :: self              !< The equation.
   integer(I4P),                     intent(in)    :: refinement_levels !< Number of refinement to be performed.
   integer(I4P)                                    :: l                 !< Counter.

   do l=1, refinement_levels
      call self%adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
      call self%adam%amr_update(do_blocks_reorder=.false., do_mpi_redistribute=.true.)
   enddo
   call self%update_cell_gpu() ! should not be needed since it is in amr_update
   call self%update_phi()  ! this is needed so that next equation%amr_updates start correctly
   !print*,'amr update a', lbound(self%x_cell_gpu,1), ubound(self%x_cell_gpu,1), size(self%x_cell_gpu,1)
   !print*,'amr update b', lbound(self%x_cell_gpu,2), ubound(self%x_cell_gpu,2), size(self%x_cell_gpu,2)
   endsubroutine

   subroutine fd_initialize(self)
   !< Initialize Finite-Difference coefficients.
   class(equation_flame_gpu_object), intent(inout) :: self !< The equation.

   allocate(self%fd_conv(4,4), self%fd_coeff1(3), self%fd_coeff2(0:3))

   ! Coefficients for finite difference computations
   associate(fd_conv=>self%fd_conv, fd_coeff1=>self%fd_coeff1,fd_coeff2=>self%fd_coeff2)

   ! Coefficients for computation of convective terms
   fd_conv(1,1) = 1._R8P/2._R8P

   fd_conv(1,2) =  2._R8P/3._R8P
   fd_conv(2,2) = -1._R8P/12._R8P

   fd_conv(1,3) =  3._R8P/4._R8P
   fd_conv(2,3) = -3._R8P/20._R8P
   fd_conv(3,3) =  1._R8P/60._R8P

   fd_conv(1,4) =  4._R8P/5._R8P
   fd_conv(2,4) = -1._R8P/5._R8P
   fd_conv(3,4) =  4._R8P/105._R8P
   fd_conv(4,4) = -1._R8P/280._R8P

   ! Coefficients for computation of viscous terms
   select case (self%ivis/2)
   case (1)
    fd_coeff1(1) = 0.5_R8P
   case (2)
    fd_coeff1(1) = 2._R8P/3._R8P
    fd_coeff1(2) = -1._R8P/12._R8P
   case (3)
    fd_coeff1(1) = 0.75_R8P
    fd_coeff1(2) = -0.15_R8P
    fd_coeff1(3) = 1._R8P/60._R8P
   end select

   select case (self%ivis/2)
   case (1)
    fd_coeff2(0) = -2._R8P
    fd_coeff2(1) =  1._R8P
   case (2)
    fd_coeff2(0) = -2.5_R8P
    fd_coeff2(1) = 4._R8P/3._R8P
    fd_coeff2(2) = -1._R8P/12._R8P
   case (3)
    fd_coeff2(0) = -245._R8P/90._R8P
    fd_coeff2(1) = 1.5_R8P
    fd_coeff2(2) = -0.15_R8P
    fd_coeff2(3) = 1._R8P/90._R8P
   endselect

   endassociate
   endsubroutine fd_initialize

   subroutine runge_kutta_initialize(self)
   !< Initialize Runge-Kutta data.
   class(equation_flame_gpu_object), intent(inout) :: self !< The equation.

   allocate(self%alph(self%nrk,self%nrk), self%beta(self%nrk), self%gamm(self%nrk))
   select case(self%nrk)
   case(1_I4P)
      self%alph(:,:) = reshape([1._R8P], [1,1])
      self%beta(:) = [1._R8P]
      self%gamm(:) = [0._R8P]
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
   endsubroutine runge_kutta_initialize

   subroutine save_hdf5(self, output_basename, t, time)
   !< Save simulation data in HDF5 format.
   class(equation_flame_gpu_object), intent(inout) :: self            !< The equation.
   character(*),                     intent(in)    :: output_basename !< Output base name.
   integer(I4P),                     intent(in)    :: t               !< Time iteration.
   real(R8P),                        intent(in)    :: time            !< Time.

   call self%copy_gpu_cpu(compute_q_aux=.true.)
   call self%adam%save_hdf5(basename=trim(output_basename)//trim(strz(t,9)),  &
                            q=self%field%q,                                   &
                            q_aux=self%q_aux,                                 &
                            q_name=['rho','rhu','rhv','rhw','rhe','rhY'], &
                            q_aux_name=['rhob','u','v','w','ya','tem','pres','ental','c'],  &
                            with_cell_morton=.true.)
   endsubroutine save_hdf5

   subroutine set_boundary_conditions(self, q_gpu)
   !< Set boundary conditions of equation.
   class(equation_flame_gpu_object), intent(in)            :: self                  !< The equation.
   real(R8P),                        intent(inout), device :: q_gpu(1:,         &
                                                                    1-self%ngc:,&
                                                                    1-self%ngc:,&
                                                                    1-self%ngc:,1:) !< Conservative variables.

   if (allocated(self%base_gpu%local_map_bc_crown_gpu)) call set_bc(nv=self%nv, ngc=self%ngc,                          &
                                                                    local_map_bc=self%base_gpu%local_map_bc_crown_gpu, &
                                                                    pres_inflow=self%pres_inflow,                      &
                                                                    tem_inflow=self%tem_inflow,                        &
                                                                    u_inflow=self%u_inflow,                            &
                                                                    gamma_fluid=self%gamma_fluid,                      &
                                                                    ya_inflow=self%ya_inflow,                          &
                                                                    dha=self%dha,                                      &
                                                                    flame_center=self%flame_center,                    &
                                                                    x_min=self%flame_x_min,                            &
                                                                    x_cell_gpu=self%x_cell_gpu,                        &
                                                                    y_cell_gpu=self%y_cell_gpu,                        &
                                                                    z_cell_gpu=self%z_cell_gpu)
   contains
      subroutine set_bc(nv, ngc, local_map_bc, pres_inflow, tem_inflow, &
          u_inflow, gamma_fluid, ya_inflow, dha, flame_center, x_min, &
          x_cell_gpu, y_cell_gpu, z_cell_gpu)
      integer(I4P), intent(in)         :: nv                      !< Number of variables.
      integer(I4P), intent(in)         :: ngc                     !< Ghost cells number.
      integer(I8P), intent(in), device :: local_map_bc(:,:,:)     !< Local map for BC ghost cells.
      real(R8P),    intent(in), device :: x_cell_gpu(1:,1-ngc:)   !< Conservative variables.
      real(R8P),    intent(in), device :: y_cell_gpu(1:,1-ngc:)   !< Conservative variables.
      real(R8P),    intent(in), device :: z_cell_gpu(1:,1-ngc:)   !< Conservative variables.
      integer(I4P)                     :: b                       !< Counter.
      integer(I4P)                     :: c, i, j, k, v           !< Counter.
      integer(I4P)                     :: idelta                  !< IJK delta step for extrapolation.
      integer(I4P)                     :: jdelta                  !< IJK delta step for extrapolation.
      integer(I4P)                     :: kdelta                  !< IJK delta step for extrapolation.
      integer(I4P)                     :: bc_type                 !< Boundary condition type.
      integer(I4P)                     :: crown                   !< Crown counter.
      integer(I4P)                     :: iercuda                 !< Error trapping flag for CUDAFortran.
      real(R8P)                        :: pres_inflow             !< Inflow values.
      real(R8P)                        :: u_inflow                !< Inflow values.
      real(R8P)                        :: tem_inflow              !< Inflow values.
      real(R8P)                        :: gamma_fluid             !< Inflow values.
      real(R8P)                        :: ya_inflow               !< Inflow values.
      real(R8P)                        :: dha                     !< Inflow values.
      real(R8P)                        :: flame_center            !< Inflow values.
      real(R8P)                        :: x_min                   !< Inflow values.
      real(R8P)                        :: u_bc                    !< Inflow values.
      integer :: i_inner, j_inner, k_inner
      real(R8P) :: rho, uuu, vvv, www, rhe, rya, yya, tem, pre, pres

      do crown=1, ngc
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
               if (bc_type == BC_EXTRAPOLATION .or. bc_type == BC_NROUT_XMAX .or. bc_type == BC_NRIN_XMIN .or. bc_type == BC_NROUT_XMIN) then
                  do v=1, nv
                     q_gpu(b,i,j,k,v) = q_gpu(b,i-idelta,j-jdelta,k-kdelta,v)
                  enddo
               else if (bc_type == BC_INFLOW) then
                   !u_bc = q_gpu(b,i,j,k,2)/q_gpu(b,i,j,k,1) - 0.001_R8P*(x_min-flame_center)
                   u_bc = u_inflow
                   q_gpu(b,i,j,k,1) = gamma_fluid/(gamma_fluid-1._R8P)*pres_inflow/tem_inflow
                   q_gpu(b,i,j,k,2) = q_gpu(b,i,j,k,1)*u_bc !inflow
                   q_gpu(b,i,j,k,3) = 0._R8P
                   q_gpu(b,i,j,k,4) = 0._R8P
                   q_gpu(b,i,j,k,5) = q_gpu(b,i,j,k,1)*(1._R8P/gamma_fluid*tem_inflow+0.5*u_bc**2+ya_inflow*dha)
                   q_gpu(b,i,j,k,6) = q_gpu(b,i,j,k,1)*ya_inflow
               else if (bc_type == BC_INFLOW_COLD) then
                   !jetif((y_cell_gpu(b,j)-2.)**2+(z_cell_gpu(b,k)-2.)**2 < 0.5_R8P) then
                   !jet    u_bc = u_inflow
                   !jetelse
                   !jet    u_bc = 0._R8P
                   !jetendif
                   u_bc = u_inflow
                   q_gpu(b,i,j,k,1) = gamma_fluid/(gamma_fluid-1._R8P)*pres_inflow/tem_inflow
                   q_gpu(b,i,j,k,2) = q_gpu(b,i,j,k,1)*u_bc !inflow
                   q_gpu(b,i,j,k,3) = 0._R8P
                   q_gpu(b,i,j,k,4) = 0._R8P
                   q_gpu(b,i,j,k,5) = q_gpu(b,i,j,k,1)*(1._R8P/gamma_fluid*tem_inflow+0.5*u_bc**2)
                   q_gpu(b,i,j,k,6) = 0._R8P
               else if (bc_type == BC_WALL_ISOTERM) then
                   i_inner = i - idelta*crown
                   j_inner = j - jdelta*crown
                   k_inner = k - kdelta*crown
                   rho = q_gpu(b,i_inner,j_inner,k_inner,1)
                   uuu = q_gpu(b,i_inner,j_inner,k_inner,2)/rho
                   vvv = q_gpu(b,i_inner,j_inner,k_inner,3)/rho
                   www = q_gpu(b,i_inner,j_inner,k_inner,4)/rho
                   rhe = q_gpu(b,i_inner,j_inner,k_inner,5)
                   rya = q_gpu(b,i_inner,j_inner,k_inner,6)
                   yya = rya/rho
                   tem = gamma_fluid*((rhe-rya*dha)/rho-0.5*(uuu**2+vvv**2+www**2))
                   pres = (gamma_fluid-1.)/gamma_fluid*rho*tem
                   ! u,v,w,temperature Dirichelet
                   ! ya, pressure Neumann
                   tem = 2._R8P*tem_inflow - tem
                   rho = pres*gamma_fluid/(gamma_fluid-1.)/tem
                   q_gpu(b,i,j,k,1) = rho
                   q_gpu(b,i,j,k,2) = -uuu*rho
                   q_gpu(b,i,j,k,3) = -vvv*rho
                   q_gpu(b,i,j,k,4) = -www*rho
                   q_gpu(b,i,j,k,5) = rho*(1./gamma_fluid*tem+0.5*(uuu**2+vvv**2+www**2)+yya*dha)
                   q_gpu(b,i,j,k,6) = rho*yya
               else if (bc_type == BC_WALL_ADIAB) then
                   i_inner = i - idelta*crown
                   j_inner = j - jdelta*crown
                   k_inner = k - kdelta*crown
                   rho = q_gpu(b,i_inner,j_inner,k_inner,1)
                   uuu = q_gpu(b,i_inner,j_inner,k_inner,2)/rho
                   vvv = q_gpu(b,i_inner,j_inner,k_inner,3)/rho
                   www = q_gpu(b,i_inner,j_inner,k_inner,4)/rho
                   rhe = q_gpu(b,i_inner,j_inner,k_inner,5)
                   rya = q_gpu(b,i_inner,j_inner,k_inner,6)
                   yya = rya/rho
                   tem = gamma_fluid*((rhe-rya*dha)/rho-0.5*(uuu**2+vvv**2+www**2))
                   ! u,v,w Dirichelet
                   ! ya, temperature, pressure Neumann
                   q_gpu(b,i,j,k,1) = rho
                   q_gpu(b,i,j,k,2) = -uuu*rho
                   q_gpu(b,i,j,k,3) = -vvv*rho
                   q_gpu(b,i,j,k,4) = -www*rho
                   q_gpu(b,i,j,k,5) = rho*(1./gamma_fluid*tem+0.5*(uuu**2+vvv**2+www**2)+yya*dha)
                   q_gpu(b,i,j,k,6) = rho*yya
               endif
            endif
         enddo
         !@cuf iercuda=cudaDeviceSynchronize()
      enddo
      endsubroutine set_bc
   endsubroutine set_boundary_conditions

   subroutine set_bc_rhs(self, q_gpu, q_aux_gpu)
   !< Set boundary conditions of equation.
   class(equation_flame_gpu_object), intent(in)            :: self                  !< The equation.
   real(R8P),                        intent(inout), device :: q_gpu(1:,         &
                                                                    1-self%ngc:,&
                                                                    1-self%ngc:,&
                                                                    1-self%ngc:,1:) !< Conservative variables.
   real(R8P),                        intent(inout), device :: q_aux_gpu(1:,     &
                                                                    1-self%ngc:,&
                                                                    1-self%ngc:,&
                                                                    1-self%ngc:,1:) !< Conservative variables.

   if (allocated(self%base_gpu%local_map_bc_crown_gpu)) call set_bc(nv=self%nv, ngc=self%ngc,                          &
                                                                    local_map_bc=self%base_gpu%local_map_bc_crown_gpu, &
                                                                    pres_inflow=self%pres_inflow,                      &
                                                                    tem_inflow=self%tem_inflow,                        &
                                                                    u_inflow=self%u_inflow,                            &
                                                                    gamma_fluid=self%gamma_fluid,                      &
                                                                    ya_inflow=self%ya_inflow,                          &
                                                                    dha=self%dha,                                      &
                                                                    flame_center=self%flame_center,                    &
                                                                    x_min=self%flame_x_min,                            &
                                                                    x_cell_gpu=self%x_cell_gpu,                        &
                                                                    y_cell_gpu=self%y_cell_gpu,                        &
                                                                    z_cell_gpu=self%z_cell_gpu)
   contains
      subroutine set_bc(nv, ngc, local_map_bc, pres_inflow, tem_inflow,   &
          u_inflow, gamma_fluid, ya_inflow, dha, flame_center, x_min,     &
          x_cell_gpu, y_cell_gpu, z_cell_gpu)
      integer(I4P), intent(in)         :: nv                      !< Number of variables.
      integer(I4P), intent(in)         :: ngc                     !< Ghost cells number.
      integer(I8P), intent(in), device :: local_map_bc(:,:,:)     !< Local map for BC ghost cells.
      real(R8P),    intent(in), device :: x_cell_gpu(1:,1-ngc:)   !< Conservative variables.
      real(R8P),    intent(in), device :: y_cell_gpu(1:,1-ngc:)   !< Conservative variables.
      real(R8P),    intent(in), device :: z_cell_gpu(1:,1-ngc:)   !< Conservative variables.
      integer(I4P)                     :: b                       !< Counter.
      integer(I4P)                     :: c, i, j, k, v           !< Counter.
      integer(I4P)                     :: idelta                  !< IJK delta step for extrapolation.
      integer(I4P)                     :: jdelta                  !< IJK delta step for extrapolation.
      integer(I4P)                     :: kdelta                  !< IJK delta step for extrapolation.
      integer(I4P)                     :: bc_type                 !< Boundary condition type.
      integer(I4P)                     :: crown                   !< Crown counter.
      integer(I4P)                     :: iercuda                 !< Error trapping flag for CUDAFortran.
      real(R8P)                        :: pres_inflow             !< Inflow values.
      real(R8P)                        :: u_inflow                !< Inflow values.
      real(R8P)                        :: tem_inflow              !< Inflow values.
      real(R8P)                        :: gamma_fluid             !< Inflow values.
      real(R8P)                        :: ya_inflow               !< Inflow values.
      real(R8P)                        :: dha                     !< Inflow values.
      real(R8P)                        :: flame_center            !< Inflow values.
      real(R8P)                        :: x_min                   !< Inflow values.
      real(R8P)                        :: u_bc                    !< Inflow values.
      real(R8P)                        :: rho, uuu, vvv, www, tem, csp, yya
      real(R8P)                        :: rhu, rhv, rhw, rhe, rya
      real(R8P)                        :: dya, duu, dvv, dww, dee, dth, dr2, dr1, ds, drho
      real(R8P)                        :: rhs_rho, rhs_rhu, rhs_rhv, rhs_rhw, rhs_rhe, rhs_rya

      crown=1

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
            if (bc_type == BC_NROUT_XMAX) then
               ! select only face-like boundaries because corresponding inner point are all selected this way
               ! actually, the next "if" also includes some edges/corners which are globally faces. in that cases,
               ! fluxes are computed here but never used since they are computed in ghost regions
               if(abs(idelta)+abs(jdelta)+abs(kdelta) == 1) then
                  i = i - idelta ; j = j - jdelta ; k = k - kdelta
                  rho      = q_aux_gpu(b,i,j,k,1)
                  uuu      = q_aux_gpu(b,i,j,k,2)
                  vvv      = q_aux_gpu(b,i,j,k,3)
                  www      = q_aux_gpu(b,i,j,k,4)
                  yya      = q_aux_gpu(b,i,j,k,5)
                  tem      = q_aux_gpu(b,i,j,k,6)
                  csp      = q_aux_gpu(b,i,j,k,9)
                  rhu      = rho*uuu
                  rhv      = rho*vvv
                  rhw      = rho*www
                  rya      = rho*yya
                  rhe      = rho*(1./gamma_fluid*tem+yya*dha+0.5*(uuu**2+vvv**2+www**2))
                  rhs_rho  = q_gpu(b,i,j,k,1)
                  rhs_rhu  = q_gpu(b,i,j,k,2)
                  rhs_rhv  = q_gpu(b,i,j,k,3)
                  rhs_rhw  = q_gpu(b,i,j,k,4)
                  rhs_rhe  = q_gpu(b,i,j,k,5)
                  rhs_rya  = q_gpu(b,i,j,k,6)
                  dya      = rhs_rya/rho-rya/rho**2*rhs_rho
                  duu      = rhs_rhu/rho-rhu/rho**2*rhs_rho
                  dvv      = rhs_rhv/rho-rhv/rho**2*rhs_rho
                  dww      = rhs_rhw/rho-rhw/rho**2*rhs_rho
                  dee      = rhs_rhe/rho-rhe/rho**2*rhs_rho
                  dth      = gamma_fluid*(dee-uuu*duu-vvv*dvv-www*dww-dha*dya)

                  dr2      = duu+(gamma_fluid-1.)/(gamma_fluid*csp)*(dth+tem/rho*rhs_rho)
                  ds       = dth/(gamma_fluid*tem) - (gamma_fluid-1.)/(gamma_fluid*rho)*rhs_rho

                  drho     = csp*rho/(2*(gamma_fluid-1.)*tem)*dr2-rho*ds
                  dth      = csp/2.*dr2+tem*ds

                  rhs_rho  = drho
                  rhs_rhu  = 0.5*rho*dr2+uuu*rhs_rho
                  dee      = 1./gamma_fluid*dth+dya*dha+uuu*0.5*dr2+vvv*dvv+www*dww
                  rhs_rhe  = rhe/rho*rhs_rho+rho*dee
                  rhs_rya  = rya/rho*rhs_rho+rho*dya
                  rhs_rhv  = rhv/rho*rhs_rho+rho*dvv
                  rhs_rhw  = rhw/rho*rhs_rho+rho*dww
                  q_gpu(b,i,j,k,1) = rhs_rho
                  q_gpu(b,i,j,k,2) = rhs_rhu
                  q_gpu(b,i,j,k,3) = rhs_rhv
                  q_gpu(b,i,j,k,4) = rhs_rhw
                  q_gpu(b,i,j,k,5) = rhs_rhe
                  q_gpu(b,i,j,k,6) = rhs_rya
               endif
            endif
            if (bc_type == BC_NROUT_XMIN) then
               ! select only face-like boundaries because corresponding inner point are all selected this way
               ! actually, the next "if" also includes some edges/corners which are globally faces. in that cases,
               ! fluxes are computed here but never used since they are computed in ghost regions
               if(abs(idelta)+abs(jdelta)+abs(kdelta) == 1) then
                  i = i - idelta ; j = j - jdelta ; k = k - kdelta
                  rho      = q_aux_gpu(b,i,j,k,1)
                  uuu      = q_aux_gpu(b,i,j,k,2)
                  vvv      = q_aux_gpu(b,i,j,k,3)
                  www      = q_aux_gpu(b,i,j,k,4)
                  yya      = q_aux_gpu(b,i,j,k,5)
                  tem      = q_aux_gpu(b,i,j,k,6)
                  csp      = q_aux_gpu(b,i,j,k,9)
                  rhu      = rho*uuu
                  rhv      = rho*vvv
                  rhw      = rho*www
                  rya      = rho*yya
                  rhe      = rho*(1./gamma_fluid*tem+yya*dha+0.5*(uuu**2+vvv**2+www**2))
                  rhs_rho  = q_gpu(b,i,j,k,1)
                  rhs_rhu  = q_gpu(b,i,j,k,2)
                  rhs_rhv  = q_gpu(b,i,j,k,3)
                  rhs_rhw  = q_gpu(b,i,j,k,4)
                  rhs_rhe  = q_gpu(b,i,j,k,5)
                  rhs_rya  = q_gpu(b,i,j,k,6)
                  dya      = rhs_rya/rho-rya/rho**2*rhs_rho
                  duu      = rhs_rhu/rho-rhu/rho**2*rhs_rho
                  dvv      = rhs_rhv/rho-rhv/rho**2*rhs_rho
                  dww      = rhs_rhw/rho-rhw/rho**2*rhs_rho
                  dee      = rhs_rhe/rho-rhe/rho**2*rhs_rho
                  dth      = gamma_fluid*(dee-uuu*duu-vvv*dvv-www*dww-dha*dya)

                  dr1      = duu-(gamma_fluid-1.)/(gamma_fluid*csp)*(dth+tem/rho*rhs_rho)
                  ds       = dth/(gamma_fluid*tem) - (gamma_fluid-1.)/(gamma_fluid*rho)*rhs_rho

                  drho     = -csp*rho/(2*(gamma_fluid-1.)*tem)*dr1-rho*ds
                  dth      = -csp/2.*dr1+tem*ds

                  rhs_rho  = drho
                  rhs_rhu  = 0.5*rho*dr1+uuu*rhs_rho
                  dee      = 1./gamma_fluid*dth+dya*dha+uuu*0.5*dr1+vvv*dvv+www*dww
                  rhs_rhe  = rhe/rho*rhs_rho+rho*dee
                  rhs_rya  = rya/rho*rhs_rho+rho*dya
                  rhs_rhv  = rhv/rho*rhs_rho+rho*dvv
                  rhs_rhw  = rhw/rho*rhs_rho+rho*dww
                  q_gpu(b,i,j,k,1) = rhs_rho
                  q_gpu(b,i,j,k,2) = rhs_rhu
                  q_gpu(b,i,j,k,3) = rhs_rhv
                  q_gpu(b,i,j,k,4) = rhs_rhw
                  q_gpu(b,i,j,k,5) = rhs_rhe
                  q_gpu(b,i,j,k,6) = rhs_rya
               endif
            endif
            if (bc_type == BC_NRIN_XMIN) then
               ! select only face-like boundaries because corresponding inner point are all selected this way
               if(abs(idelta)+abs(jdelta)+abs(kdelta) == 1) then
                  i = i - idelta ; j = j - jdelta ; k = k - kdelta
                  rho      = q_aux_gpu(b,i,j,k,1)
                  uuu      = q_aux_gpu(b,i,j,k,2)
                  vvv      = q_aux_gpu(b,i,j,k,3)
                  www      = q_aux_gpu(b,i,j,k,4)
                  yya      = q_aux_gpu(b,i,j,k,5)
                  tem      = q_aux_gpu(b,i,j,k,6)
                  csp      = q_aux_gpu(b,i,j,k,9)
                  rhu      = rho*uuu
                  rhv      = rho*vvv
                  rhw      = rho*www
                  rya      = rho*yya
                  rhe      = rho*(1./gamma_fluid*tem+yya*dha+0.5*(uuu**2+vvv**2+www**2))
                  rhs_rho  = q_gpu(b,i,j,k,1)
                  rhs_rhu  = q_gpu(b,i,j,k,2)
                  rhs_rhv  = q_gpu(b,i,j,k,3)
                  rhs_rhw  = q_gpu(b,i,j,k,4)
                  rhs_rhe  = q_gpu(b,i,j,k,5)
                  rhs_rya  = q_gpu(b,i,j,k,6)
                  dya       = rhs_rya/rho-rya/rho**2*rhs_rho
                  duu       = rhs_rhu/rho-rhu/rho**2*rhs_rho
                  dvv       = rhs_rhv/rho-rhv/rho**2*rhs_rho
                  dww       = rhs_rhw/rho-rhw/rho**2*rhs_rho
                  dee       = rhs_rhe/rho-rhe/rho**2*rhs_rho
                  dth       = gamma_fluid*(dee-uuu*duu-vvv*dvv-www*dww-dha*dya)
                  dr1       = duu-(gamma_fluid-1.)/(gamma_fluid*csp)*(dth+tem/rho*rhs_rho)
                  duu       = 0._R8P
                  dvv       = 0._R8P
                  dww       = 0._R8P
                  rhs_rho   = rho*csp/(((gamma_fluid-1.)/gamma_fluid)*tem)*(duu-dr1)
                  dee       = uuu*duu+vvv*dvv+www*dww
                  rhs_rhu   = rho*duu+uuu*rhs_rho
                  rhs_rhe   = rhe/rho*rhs_rho+rho*dee
                  rhs_rya   = yya*rhs_rho
                  rhs_rhv   = vvv*rhs_rho+rho*dvv
                  rhs_rhw   = www*rhs_rho+rho*dww
                  q_gpu(b,i,j,k,1) = rhs_rho
                  q_gpu(b,i,j,k,2) = rhs_rhu
                  q_gpu(b,i,j,k,3) = rhs_rhv
                  q_gpu(b,i,j,k,4) = rhs_rhw
                  q_gpu(b,i,j,k,5) = rhs_rhe
                  q_gpu(b,i,j,k,6) = rhs_rya
               endif
            endif
         endif
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()

      endsubroutine set_bc
   endsubroutine set_bc_rhs

   subroutine set_initial_conditions(self)
   !< Set initial conditions of field.
   class(equation_flame_gpu_object), intent(inout) :: self       !< The equation.
   integer(I4P)                                    :: b, i, j, k !< Counter.
   real(R8P)                                       :: tem        !< Scalar.
   real(R8P)                                       :: yya        !< Scalar.
   real(R8P)                                       :: rho        !< Scalar.
   real(R8P)                                       :: ene        !< Scalar.
   real(R8P)                                       :: gauss      !< Scalar.
   real(R8P)                                       :: pres       !< Scalar.

   associate(blocks_number=>self%blocks_number, q=>self%field%q, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             ngc=>self%ngc, x_cell=>self%field%x_cell, y_cell=>self%field%y_cell, z_cell=>self%field%z_cell, &
             tem_inflow=>self%tem_inflow, ya_inflow=>self%ya_inflow, gamma_fluid=>self%gamma_fluid, &
             pres_inflow=>self%pres_inflow, tem_outflow=>self%tem_outflow, u_inflow=>self%u_inflow, dha=>self%dha, &
             flame_center=>self%flame_center)

      if(self%flow_type == "flamechannel") then
         do b=1, blocks_number
            do k=1, nk
               do j=1, nj
                  do i=1, ni
                     gauss = exp(-0.5*(((x_cell(i,b)-3)/1.)**2+(y_cell(j,b)-2)**2+(z_cell(k,b)-2)**2))
                     tem = tem_inflow + 700*gauss
                     yya = ya_inflow*(1-gauss)
                     rho = gamma_fluid/(gamma_fluid-1._R8P)*pres_inflow/tem
                     ene = 1._R8P/gamma_fluid*tem+0.5*u_inflow**2+yya*dha
                     q(1,i,j,k,b) = rho
                     q(2,i,j,k,b) = 0.
                     q(3,i,j,k,b) = 0._R8P
                     q(4,i,j,k,b) = 0._R8P
                     q(5,i,j,k,b) = rho * ene
                     q(6,i,j,k,b) = rho * yya
                  enddo
               enddo
            enddo
         enddo
      elseif(self%flow_type == "sod") then
         do b=1, blocks_number
            do k=1, nk
               do j=1, nj
                  do i=1, ni
                     if(x_cell(i,b) < 0.5) then
                        rho  = 1.0
                        pres = 1.0
                        tem  = pres*gamma_fluid/(rho*(gamma_fluid-1))
                        q(1,i,j,k,b) = rho
                        q(2,i,j,k,b) = 0.
                        q(3,i,j,k,b) = 0._R8P
                        q(4,i,j,k,b) = 0._R8P
                        q(5,i,j,k,b) = rho/gamma_fluid*tem
                        q(6,i,j,k,b) = 0.
                     else
                        rho  = 0.125
                        pres = 0.1
                        tem  = pres*gamma_fluid/(rho*(gamma_fluid-1))
                        q(1,i,j,k,b) = rho
                        q(2,i,j,k,b) = 0.
                        q(3,i,j,k,b) = 0._R8P
                        q(4,i,j,k,b) = 0._R8P
                        q(5,i,j,k,b) = rho/gamma_fluid*tem
                        q(6,i,j,k,b) = 0.
                     endif
                  enddo
               enddo
            enddo
         enddo
      else
         do b=1, blocks_number
            do k=1, nk
               do j=1, nj
                  do i=1, ni
                     tem = tem_inflow + (tem_outflow-tem_inflow)/acos(-1._R8P)*(-atan(-500._R8P)+atan(50.*(x_cell(i,b)-flame_center)))
                     yya = ya_inflow*(1._R8P-0.5_R8P*(2._R8P/acos(-1._R8P)*atan(x_cell(i,b)-flame_center)+1.))
                     rho = gamma_fluid/(gamma_fluid-1._R8P)*pres_inflow/tem
                     ene = 1._R8P/gamma_fluid*tem+0.5*u_inflow**2+yya*dha
                     q(1,i,j,k,b) = rho
                     if(trim(self%flow_type) == "cold") then
                         q(2,i,j,k,b) = rho * 1.5_R8P
                     else
                         q(2,i,j,k,b) = rho * u_inflow
                     endif
                     q(3,i,j,k,b) = 0._R8P
                     q(4,i,j,k,b) = 0._R8P
                     q(5,i,j,k,b) = rho * ene
                     q(6,i,j,k,b) = rho * yya
                  enddo
               enddo
            enddo
         enddo
      endif

   endassociate
   call self%copy_cpu_gpu
   endsubroutine set_initial_conditions

   subroutine update_ghost_gpu(self, q_gpu, step)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(equation_flame_gpu_object), intent(inout)         :: self            !< The equation.
   real(R8P),                        intent(inout), device :: q_gpu(1:,         &
                                                                    1-self%ngc:,&
                                                                    1-self%ngc:,&
                                                                    1-self%ngc:,&
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
   class(equation_flame_gpu_object), intent(inout) :: lhs !< Left hand side.
   type(equation_flame_gpu_object),  intent(in)    :: rhs !< Right hand side.

   lhs%adam          => rhs%adam
   lhs%field         => rhs%field
   lhs%grid          => rhs%grid
   lhs%ni            => rhs%ni
   lhs%nj            => rhs%nj
   lhs%nk            => rhs%nk
   lhs%ngc           => rhs%ngc
   lhs%nb            => rhs%nb
   lhs%blocks_number => rhs%blocks_number
   lhs%nv            => rhs%nv
   lhs%base_gpu = rhs%base_gpu
   lhs%myrank = rhs%myrank
   lhs%procs_number = rhs%procs_number
   lhs%error = rhs%error
   lhs%ns = rhs%ns
   lhs%dt = rhs%dt
   lhs%CFL = rhs%CFL
   lhs%null_xyz = rhs%null_xyz
   lhs%nrk = rhs%nrk
   call assign_allocatable(lhs=lhs%q_aux, rhs=rhs%q_aux )
   call assign_allocatable(lhs=lhs%alph,  rhs=rhs%alph  )
   call assign_allocatable(lhs=lhs%beta,  rhs=rhs%beta  )
   call assign_allocatable(lhs=lhs%gamm,  rhs=rhs%gamm  )
   call assign_allocatable_gpu(lhs=lhs%fl_gpu,     rhs=rhs%fl_gpu      )
   call assign_allocatable_gpu(lhs=lhs%fhat_gpu,   rhs=rhs%fhat_gpu      )
   call assign_allocatable_gpu(lhs=lhs%dxyz_gpu,   rhs=rhs%dxyz_gpu   )
   call assign_allocatable_gpu(lhs=lhs%alph_gpu,   rhs=rhs%alph_gpu   )
   call assign_allocatable_gpu(lhs=lhs%beta_gpu,   rhs=rhs%beta_gpu   )
   call assign_allocatable_gpu(lhs=lhs%gamm_gpu,   rhs=rhs%gamm_gpu   )
   call assign_allocatable_gpu(lhs=lhs%q_aux_gpu,  rhs=rhs%q_aux_gpu  )
   call assign_allocatable_gpu(lhs=lhs%q_gpu,      rhs=rhs%q_gpu      )
   call assign_allocatable_gpu(lhs=lhs%q_s_gpu,    rhs=rhs%q_s_gpu    )
   endsubroutine eq_assign_eq

   ! non TBP cuf procedures
   subroutine advance_q_gpu_cuf(ni, nj, nk, ngc, nv, nrk, blocks_number, beta_gpu, dt, q_s_gpu, q_gpu)
   !< Advance q_gpu by means of RK stages.
   integer(I4P), intent(in)            :: ni                                     !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                     !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                     !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                    !< Ghost grid number.
   integer(I4P), intent(in)            :: nv                                     !< Number of conservative varibales.
   integer(I4P), intent(in)            :: nrk                                    !< Number of RK stages.
   integer(I4P), intent(in)            :: blocks_number                          !< Number of blocks.
   real(R8P),    intent(in),    device :: beta_gpu(:)                            !< RK betaa coefficients.
   real(R8P),    intent(in)            :: Dt                                     !< Time step.
   real(R8P),    intent(in),    device :: q_s_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:) !< RK stage.
   real(R8P),    intent(inout), device ::   q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative variables.
   integer(I4P)                        :: i, j, k, b, s, v                       !< Counter.
   integer(I4P)                        :: iercuda                                !< Error trapping flag for CUDAFortran.

   do s=1, nrk
      !$cuf kernel do(5) <<<*,*>>>
      do v=1, nv
         do k=1-ngc, nk+ngc
            do j=1-ngc, nj+ngc
               do i=1-ngc, ni+ngc
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

   subroutine flame_find_x_v_cuf(ni, nj, nk, ngc, nv, blocks_number, dt, &
       tem_stabil, x_cell_gpu, y_cell_gpu, z_cell_gpu, q_aux_gpu, x_min, x_min_old, velrel)

   !< Advance q_gpu by means of RK stages.
   integer(I4P), intent(in)            :: ni                                     !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                     !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                     !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                    !< Ghost grid number.
   integer(I4P), intent(in)            :: nv                                     !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number                          !< Number of blocks.
   real(R8P),    intent(in)            :: tem_stabil                             !< Time step.
   real(R8P),    intent(in)            :: Dt                                     !< Time step.
   real(R8P),    intent(in), device    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Conservative variables.
   real(R8P),    intent(in), device    :: x_cell_gpu(1:,1-ngc:)                  !< Conservative variables.
   real(R8P),    intent(in), device    :: y_cell_gpu(1:,1-ngc:)                  !< Conservative variables.
   real(R8P),    intent(in), device    :: z_cell_gpu(1:,1-ngc:)                  !< Conservative variables.
   real(R8P),    intent(inout)         :: x_min                                  !< Conservative variables.
   real(R8P),    intent(inout)         :: x_min_old                              !< Conservative variables.
   real(R8P),    intent(inout)         :: velrel                                 !< Conservative variables.
   integer(I4P)                        :: i, j, k, b, s, v                       !< Counter.
   integer(I4P)                        :: iercuda                                !< Error trapping flag for CUDAFortran.
   integer(I4P)                        :: ierr                                   !< Error trapping flag for CUDAFortran.

   x_min = huge(1.0_R8P)

   !$cuf kernel do(4) <<<*,*>>>
   do k=1-ngc, nk+ngc
      do j=1-ngc, nj+ngc
         do i=1-ngc, ni+ngc
            do b=1, blocks_number
               if(q_aux_gpu(b,i,j,k,6) > tem_stabil) then
                   x_min = min(x_min, x_cell_gpu(b,i))
               endif
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   call MPI_ALLREDUCE(MPI_IN_PLACE, x_min, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
   print*,'flame x_min: ',x_min
   endsubroutine flame_find_x_v_cuf

   subroutine compute_aux_cuf(ni, nj, nk, ngc, ns, blocks_number, gamma_fluid, dha, q_gpu, q_aux_gpu)
   !< Compute auxiliary variables by means of CUF threads.
   integer(I4P), intent(in)          :: ni                                     !< Grid cells number in I direction.
   integer(I4P), intent(in)          :: nj                                     !< Grid cells number in J direction.
   integer(I4P), intent(in)          :: nk                                     !< Grid cells number in K direction.
   integer(I4P), intent(in)          :: ngc                                    !< Ghost cells number.
   integer(I4P), intent(in)          :: ns                                     !< Number of fluid species.
   integer(I4P), intent(in)          :: blocks_number                          !< Number of blocks.
   real(R8P),    intent(in)          :: gamma_fluid                            !< Gamma fluid.
   real(R8P),    intent(in)          :: dha                                    !< Entalpy fluid.
   real(R8P),    intent(in),  device ::     q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Conservative variables.
   real(R8P),    intent(out), device :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Auxiliary variables.
   integer(I4P)                      :: b, i, j, k, s                          !< Counter.
   real(R8P)                         :: rho, uuu, vvv, www, rhe, rya, yya, tem !< Scalars.
   integer(I4P)                      :: iercuda                                !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(4) <<<*,*>>>
   do k=1-ngc, nk+ngc
      do j=1-ngc, nj+ngc
         do i=1-ngc, ni+ngc
            do b=1, blocks_number
               rho = q_gpu(b,i,j,k,1)
               uuu = q_gpu(b,i,j,k,2)/rho
               vvv = q_gpu(b,i,j,k,3)/rho
               www = q_gpu(b,i,j,k,4)/rho
               rhe = q_gpu(b,i,j,k,5)
               rya = q_gpu(b,i,j,k,6)
               yya = rya/rho
               tem = gamma_fluid*((rhe-rya*dha)/rho-0.5*(uuu**2+vvv**2+www**2))

               q_aux_gpu(b,i,j,k,1) = rho                                          ! density
               q_aux_gpu(b,i,j,k,2) = uuu                                          ! velocity x
               q_aux_gpu(b,i,j,k,3) = vvv                                          ! velocity y
               q_aux_gpu(b,i,j,k,4) = www                                          ! velocity z
               q_aux_gpu(b,i,j,k,5) = yya                                          ! mass fraction
               q_aux_gpu(b,i,j,k,6) = tem                                          ! temperature
               q_aux_gpu(b,i,j,k,7) = (gamma_fluid-1._R8P)/gamma_fluid*rho*tem     ! pressure
               q_aux_gpu(b,i,j,k,8) = rhe/rho+(gamma_fluid-1._R8P)/gamma_fluid*tem ! entalpy
               q_aux_gpu(b,i,j,k,9) = sqrt((gamma_fluid-1._R8P)*tem)               ! sound speed
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_aux_cuf

   subroutine evolve_eikonal_q_gpu_cuf(ni, nj, nk, ngc, nv, phi_gpu, dx_gpu, dy_gpu, dz_gpu, blocks_number, dq_gpu, q_gpu)
   !< Initialize RK stage with q_gpu.
   integer(I4P), intent(in)            :: ni                                  !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                  !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                  !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in)            :: nv                                  !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number                       !< Number of blocks.
   real(R8P),    intent(in),    device :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function.
   real(R8P),    intent(in),    device :: dx_gpu(1:)                          !< X space steps.
   real(R8P),    intent(in),    device :: dy_gpu(1:)                          !< Y space steps.
   real(R8P),    intent(in),    device :: dz_gpu(1:)                          !< Z space steps.
   real(R8P),    intent(inout), device :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !<
   real(R8P),    intent(inout), device :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Conservative variables.
   integer(I4P)                        :: i, j, k, b, v                       !< Counter.
   integer(I4P)                        :: iercuda                             !< Error trapping flag for CUDAFortran.
   type(dim3)                          :: grid, tBlock                        !< CUDA grid and block.

   tBlock = dim3(32,8,1) ; grid = dim3(ceiling(real(blocks_number)/tBlock%x),ceiling(real(ni)/tBlock%y),1)
   call compute_eikonal_dq_gpu<<<grid, tBlock>>>(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                                 phi_gpu=phi_gpu, dx_gpu=dx_gpu, dy_gpu=dy_gpu, dz_gpu=dz_gpu,     &
                                                 dq_gpu=dq_gpu, q_gpu=q_gpu)

   !$cuf kernel do(4) <<<*,*>>>
   do k=1, nk
      do j=1, nj
         do i=1,ni
            do b=1, blocks_number
               do v=1, nv
                  if (phi_gpu(b,i,j,k,1) > 0._R8P) q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) - dq_gpu(b,i,j,k,v)
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine evolve_eikonal_q_gpu_cuf

   attributes(global) subroutine compute_eikonal_dq_gpu(ni, nj, nk, ngc, nv, blocks_number, &
                                                        phi_gpu, dx_gpu, dy_gpu, dz_gpu,    &
                                                        q_gpu, dq_gpu)
   !< Initialize RK stage with q_gpu.
   integer(I4P), intent(in), value     :: ni                                  !< Grid cells number in I direction.
   integer(I4P), intent(in), value     :: nj                                  !< Grid cells number in J direction.
   integer(I4P), intent(in), value     :: nk                                  !< Grid cells number in K direction.
   integer(I4P), intent(in), value     :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in), value     :: nv                                  !< Number of conservative varibales.
   integer(I4P), intent(in), value     :: blocks_number                       !< Number of blocks.
   real(R8P),    intent(in),    device :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function.
   real(R8P),    intent(in),    device :: dx_gpu(1:)                          !< X space steps.
   real(R8P),    intent(in),    device :: dy_gpu(1:)                          !< Y space steps.
   real(R8P),    intent(in),    device :: dz_gpu(1:)                          !< Z space steps.
   real(R8P),    intent(in),    device :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Conservative variables.
   real(R8P),    intent(inout), device :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !<
   integer(I4P)                        :: i, j, k, b, v                       !< Counter.
   real(R8P)                           :: n_phi_x, n_phi_y, n_phi_z, n_phi    !<

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   i = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if (b > blocks_number .or. i > ni) return

   do k=1, nk
      do j=1, nj
         if (phi_gpu(b,i,j,k,1) > 0._R8P) then
            n_phi_x = -(phi_gpu(b,i+1,j,k,1) - phi_gpu(b,i-1,j,k,1) ) / (2 * dx_gpu(b))
            n_phi_y = -(phi_gpu(b,i,j+1,k,1) - phi_gpu(b,i,j-1,k,1) ) / (2 * dy_gpu(b))
            n_phi_z = -(phi_gpu(b,i,j,k+1,1) - phi_gpu(b,i,j,k-1,1) ) / (2 * dz_gpu(b))
            n_phi = abs(n_phi_x) + abs(n_phi_y) + abs(n_phi_z) + 10e-12
            n_phi = 0.5_R8P / n_phi
            n_phi_x = n_phi_x * n_phi
            n_phi_y = n_phi_y * n_phi
            n_phi_z = n_phi_z * n_phi
            do v=1, nv
               dq_gpu(b,i,j,k,v) = 0._R8P
            enddo
            if (n_phi_x > 0._R8P) then
               do v=1, nv
                  dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_x) * (q_gpu(b,i,j,k,v) - q_gpu(b,i-1,j,k,v))
               enddo
            else
               do v=1, nv
                  dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_x) * (q_gpu(b,i,j,k,v) - q_gpu(b,i+1,j,k,v))
               enddo
            endif
            if (n_phi_y > 0._R8P) then
               do v=1, nv
                  dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_y) * (q_gpu(b,i,j,k,v) - q_gpu(b,i,j-1,k,v))
               enddo
            else
               do v=1, nv
                  dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_y) * (q_gpu(b,i,j,k,v) - q_gpu(b,i,j+1,k,v))
               enddo
            endif
            if (n_phi_z > 0._R8P) then
               do v=1, nv
                  dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_z) * (q_gpu(b,i,j,k,v) - q_gpu(b,i,j,k-1,v))
               enddo
            else
               do v=1, nv
                  dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_z) * (q_gpu(b,i,j,k,v) - q_gpu(b,i,j,k+1,v))
               enddo
            endif
         endif
      enddo
   enddo
   endsubroutine compute_eikonal_dq_gpu

   subroutine compute_residuals_gpu(ni, nj, nk, ngc, ns, blocks_number, null_x, null_y, null_z, &
                                    dx_gpu, dy_gpu, dz_gpu, q_aux_gpu, fl_gpu, fhat_gpu, q_gpu, &
                                    iweno, lmax, dha, gamma_fluid, gplus_x, gminus_x, &
                                    Prandtl, q_coeff, Lewis, Zeldovich, Damkohler, ivis, visc_type, &
                                    fd_conv_gpu, fd_coeff1_gpu, fd_coeff2_gpu )
   !< Compute residuals of equation.
   integer(I4P), intent(in)            :: ni                                    !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                    !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                    !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                   !< Ghost cells number.
   integer(I4P), intent(in)            :: ns                                    !< Number of species.
   integer(I4P), intent(in)            :: blocks_number                         !< Number of blocks.
   integer(I4P), intent(in)            :: iweno                                 !< Weno order.
   integer(I4P), intent(in)            :: lmax                                  !< Weno order.
   integer(I4P), intent(in)            :: ivis                                  !< Weno order.
   integer(I4P), intent(in)            :: visc_type                             !< Weno order.
   real(R8P),    intent(in)            :: dha                                   !< Formation entalpy.
   real(R8P),    intent(in)            :: gamma_fluid                           !< Gamma.
   real(R8P),    intent(in)            :: Prandtl                               !< Gamma.
   real(R8P),    intent(in)            :: q_coeff                               !< Gamma.
   real(R8P),    intent(in)            :: Lewis                                 !< Gamma.
   real(R8P),    intent(in)            :: Zeldovich                             !< Gamma.
   real(R8P),    intent(in)            :: Damkohler                             !< Gamma.
   logical,      intent(in)            :: null_x                                !< Nullify x direction.
   logical,      intent(in)            :: null_y                                !< Nullify y direction.
   logical,      intent(in)            :: null_z                                !< Nullify z direction.
   real(R8P),    intent(in),    device :: fd_conv_gpu(1:,1:)                    !< Convective coefficients.
   real(R8P),    intent(in),    device :: fd_coeff1_gpu(1:)                     !< Derivatives first coefficients.
   real(R8P),    intent(in),    device :: fd_coeff2_gpu(0:)                     !< Derivatives second coefficients.
   real(R8P),    intent(in),    device :: dx_gpu(1:)                            !< X space steps.
   real(R8P),    intent(in),    device :: dy_gpu(1:)                            !< Y space steps.
   real(R8P),    intent(in),    device :: dz_gpu(1:)                            !< Z space steps.
   real(R8P),    intent(in),    device :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout), device :: gplus_x(1:,1:,1:,1:,1:)               !< Auxiliary variables.
   real(R8P),    intent(inout), device :: gminus_x(1:,1:,1:,1:,1:)              !< Auxiliary variables.
   real(R8P),    intent(inout), device :: fl_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Positive fluxes.
   real(R8P),    intent(inout), device :: fhat_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Negative fluxes.
   real(R8P),    intent(inout), device :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Conservative variables.
   integer(I4P)                        :: b, i, j, k, v                         !< Counter.
   integer(I4P)                        :: iercuda                               !< Error trapping flag for CUDAFortran.
   type(dim3)                          :: grid, tBlock                          !< CUDA grid and block.

   real(R8P) :: dx    , dy     , dz
   real(R8P) :: rho   , uu     , vv    , ww
   real(R8P) :: t     , p      , h
   real(R8P) :: rho_xp, rho_yp , rho_zp
   real(R8P) :: rho_xm, rho_ym , rho_zm
   real(R8P) :: u_xp  , u_yp   , u_zp
   real(R8P) :: u_xm  , u_ym   , u_zm
   real(R8P) :: v_xp  , v_yp   , v_zp
   real(R8P) :: v_xm  , v_ym   , v_zm
   real(R8P) :: w_xp  , w_yp   , w_zp
   real(R8P) :: w_xm  , w_ym   , w_zm
   real(R8P) :: p_xp  , p_yp   , p_zp
   real(R8P) :: p_xm  , p_ym   , p_zm
   real(R8P) :: t_xp  , t_yp   , t_zp
   real(R8P) :: t_xm  , t_ym   , t_zm
   real(R8P) :: h_xp  , h_yp   , h_zp
   real(R8P) :: h_xm  , h_ym   , h_zm
   real(R8P) :: ulap, vlap, wlap, tlap
   real(R8P) :: ux, uy, uz, vx, vy, vz, wx, wy, wz, div3l
   real(R8P) :: sigx, sigy, sigz, sig11, sig12, sig13, sig22, sig23, sig33, sigah, sigqt, sigq

   !fl_gpu = 0.

   !!!!!! minimal navier-stokes - start
   !!!!!associate(q=>q_aux_gpu)
   !!!!!!$cuf kernel do(4) <<<*,*>>>
   !!!!!do k=1, nk
   !!!!!   do j=1, nj
   !!!!!      do i=1,ni
   !!!!!         do b=1, blocks_number
   !!!!!            dx     = dx_gpu(b)      ; dy     = dy_gpu(b)      ; dz     = dz_gpu(b)
   !!!!!            rho    = q(b,i,j,k,1)   ; t      = q(b,i,j,k,6)   ; p      = q(b,i,j,k,7)
   !!!!!            uu     = q(b,i,j,k,2)   ; vv     = q(b,i,j,k,3)   ; ww     = q(b,i,j,k,4)
   !!!!!            rho_xp = q(b,i+1,j,k,1) ; rho_yp = q(b,i,j+1,k,1) ; rho_zp = q(b,i,j,k+1,1)
   !!!!!            rho_xm = q(b,i-1,j,k,1) ; rho_ym = q(b,i,j-1,k,1) ; rho_zm = q(b,i,j,k-1,1)
   !!!!!            u_xp   = q(b,i+1,j,k,2) ; u_yp   = q(b,i,j+1,k,2) ; u_zp   = q(b,i,j,k+1,2)
   !!!!!            u_xm   = q(b,i-1,j,k,2) ; u_ym   = q(b,i,j-1,k,2) ; u_zm   = q(b,i,j,k-1,2)
   !!!!!            v_xp   = q(b,i+1,j,k,3) ; v_yp   = q(b,i,j+1,k,3) ; v_zp   = q(b,i,j,k+1,3)
   !!!!!            v_xm   = q(b,i-1,j,k,3) ; v_ym   = q(b,i,j-1,k,3) ; v_zm   = q(b,i,j,k-1,3)
   !!!!!            w_xp   = q(b,i+1,j,k,4) ; w_yp   = q(b,i,j+1,k,4) ; w_zp   = q(b,i,j,k+1,4)
   !!!!!            w_xm   = q(b,i-1,j,k,4) ; w_ym   = q(b,i,j-1,k,4) ; w_zm   = q(b,i,j,k-1,4)
   !!!!!            p_xp   = q(b,i+1,j,k,7) ; p_yp   = q(b,i,j+1,k,7) ; p_zp   = q(b,i,j,k+1,7)
   !!!!!            p_xm   = q(b,i-1,j,k,7) ; p_ym   = q(b,i,j-1,k,7) ; p_zm   = q(b,i,j,k-1,7)
   !!!!!            t_xp   = q(b,i+1,j,k,6) ; t_yp   = q(b,i,j+1,k,6) ; t_zp   = q(b,i,j,k+1,6)
   !!!!!            t_xm   = q(b,i-1,j,k,6) ; t_ym   = q(b,i,j-1,k,6) ; t_zm   = q(b,i,j,k-1,6)
   !!!!!            h_xp   = q(b,i+1,j,k,8) ; h_yp   = q(b,i,j+1,k,8) ; h_zp   = q(b,i,j,k+1,8)
   !!!!!            h_xm   = q(b,i-1,j,k,8) ; h_ym   = q(b,i,j-1,k,8) ; h_zm   = q(b,i,j,k-1,8)

   !!!!!            ulap = (u_xp-2.*uu+u_xm)/dx**2 + (u_yp-2.*uu+u_ym)/dy**2 + (u_zp-2.*uu+u_zm)/dz**2
   !!!!!            vlap = (v_xp-2.*vv+v_xm)/dx**2 + (v_yp-2.*vv+v_ym)/dy**2 + (v_zp-2.*vv+v_zm)/dz**2
   !!!!!            wlap = (w_xp-2.*ww+w_xm)/dx**2 + (w_yp-2.*ww+w_ym)/dy**2 + (w_zp-2.*ww+w_zm)/dz**2
   !!!!!            tlap = (t_xp-2.*t+t_xm) /dx**2 + (t_yp-2.*t+t_ym) /dy**2 + (t_zp-2.*t+t_zm) /dz**2
   !!!!!            ux = (u_xp-u_xm)/(2*dx) ; uy = (u_yp-u_ym)/(2*dy) ; uz = (u_zp-u_zm)/(2*dz)
   !!!!!            vx = (v_xp-v_xm)/(2*dx) ; vy = (v_yp-v_ym)/(2*dy) ; vz = (v_zp-v_zm)/(2*dz)
   !!!!!            wx = (w_xp-w_xm)/(2*dx) ; wy = (w_yp-w_ym)/(2*dy) ; wz = (w_zp-w_zm)/(2*dz)
   !!!!!            div3l = ux+vy+wz ; div3l   = div3l/3._R8P
   !!!!!            sig11 = 2._R8P*(ux-div3l) ; sig12 = uy+vx ; sig13 = uz+wx
   !!!!!            sig22 = 2._R8P*(vy-div3l) ; sig23 = vz+wy ; sig33 = 2._R8P*(wz-div3l)
   !!!!!            ! Viscosity diffusion
   !!!!!            sigx  = Prandtl*ulap ; sigy  = Prandtl*vlap ; sigz  = Prandtl*wlap
   !!!!!            ! Heat conduction
   !!!!!            sigqt = tlap
   !!!!!            ! Aerodynamic heating
   !!!!!            sigah = (sig11*ux+sig12*uy+sig13*uz+sig12*vx+sig22*vy+sig23*vz+sig13*wx+sig23*wy+sig33*wz)*Prandtl ! Aerodynamic heating
   !!!!!            ! Total energy diffusion
   !!!!!            sigq  = sigx*uu+sigy*vv+sigz*ww+sigah+sigqt

   !!!!!            q_gpu(b,i,j,k,1) = - (rho_xp*u_xp - rho_xm*u_xm)/(2.*dx) - &
   !!!!!                                 (rho_yp*v_yp - rho_ym*v_ym)/(2.*dy) - &
   !!!!!                                 (rho_zp*w_zp - rho_zm*w_zm)/(2.*dz)
   !!!!!            q_gpu(b,i,j,k,2) = - (rho_xp*u_xp*u_xp + p_xp - rho_xm*u_xm*u_xm - p_xm)/(2.*dx) - &
   !!!!!                                 (rho_yp*u_yp*v_yp        - rho_ym*u_ym*v_ym)       /(2.*dy) - &
   !!!!!                                 (rho_zp*u_zp*w_zp        - rho_zm*u_zm*w_zm)       /(2.*dz) + sigx
   !!!!!            q_gpu(b,i,j,k,3) = - (rho_xp*v_xp*u_xp        - rho_xm*v_xm*u_xm       )/(2.*dx) - &
   !!!!!                                 (rho_yp*v_yp*v_yp + p_yp - rho_ym*v_ym*v_ym - p_ym)/(2.*dy) - &
   !!!!!                                 (rho_zp*v_zp*w_zp        - rho_zm*v_zm*w_zm)       /(2.*dz) + sigy
   !!!!!            q_gpu(b,i,j,k,4) = - (rho_xp*w_xp*u_xp        - rho_xm*w_xm*u_xm       )/(2.*dx) - &
   !!!!!                                 (rho_yp*w_yp*v_yp        - rho_ym*w_ym*v_ym)       /(2.*dy) - &
   !!!!!                                 (rho_zp*w_zp*w_zp + p_zp - rho_zm*w_zm*w_zm - p_zm)/(2.*dz) + sigz
   !!!!!            q_gpu(b,i,j,k,5) = - (rho_xp*h_xp*u_xp        - rho_xm*h_xm*u_xm       )/(2.*dx) - &
   !!!!!                                 (rho_yp*h_yp*v_yp        - rho_ym*h_ym*v_ym)       /(2.*dy) - &
   !!!!!                                 (rho_zp*h_zp*w_zp        - rho_zm*h_zm*w_zm       )/(2.*dz) + sigq
   !!!!!
   !!!!!         enddo
   !!!!!      enddo
   !!!!!   enddo
   !!!!!enddo
   !!!!!!@cuf iercuda=cudaDeviceSynchronize()
   !!!!!endassociate

   !!!! minimal navier-stokes - end

   tBlock = dim3(32,8,1) ; grid = dim3(ceiling(real(blocks_number)/tBlock%x),ceiling(real(nj)/tBlock%y),1)
   !call euler_x_central_kernel<<<grid, tBlock>>>(q_gpu, q_aux_gpu, fl_gpu, fhat_gpu, &
   !                       fd_conv_gpu, dx_gpu, blocks_number, ni, nj, nk, ngc, ns+4, lmax)
   call euler_x_kernel<<<grid, tBlock>>>(q_gpu, q_aux_gpu, fl_gpu, fhat_gpu, &
                                         gplus_x, gminus_x, dx_gpu,          &
                                         blocks_number, ni, nj, nk, ngc, ns+4, iweno, dha, gamma_fluid)

   tBlock = dim3(32,8,1) ; grid = dim3(ceiling(real(blocks_number)/tBlock%x),ceiling(real(ni)/tBlock%y),1)
   call euler_y_central_kernel<<<grid, tBlock>>>(q_gpu, q_aux_gpu, fl_gpu, fhat_gpu, &
                          fd_conv_gpu, dy_gpu, blocks_number, ni, nj, nk, ngc, ns+4, lmax)

   tBlock = dim3(32,8,1) ; grid = dim3(ceiling(real(blocks_number)/tBlock%x),ceiling(real(ni)/tBlock%y),1)
   call euler_z_central_kernel<<<grid, tBlock>>>(q_gpu, q_aux_gpu, fl_gpu, fhat_gpu, &
                          fd_conv_gpu, dz_gpu, blocks_number, ni, nj, nk, ngc, ns+4, lmax)

   !@cuf iercuda=cudaDeviceSynchronize()

   if(Prandtl > 0.) call viscous_cuf(ni, nj, nk, ngc, blocks_number, ivis, visc_type, fd_coeff1_gpu, fd_coeff2_gpu, &
      gamma_fluid, Prandtl, q_coeff, Lewis, Zeldovich, Damkohler, dha, q_aux_gpu, dx_gpu, dy_gpu, dz_gpu, fl_gpu)

   !$cuf kernel do(5) <<<*,*>>>
   do v=1, 6 !ns+4
      do k=1, nk
         do j=1, nj
            do i=1,ni
               do b=1, blocks_number
                  q_gpu(b,i,j,k,v) = - fl_gpu(b,i,j,k,v)
                  !q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) - fl_gpu(b,i,j,k,v)
               enddo
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()

   !!!!! Force specie mass fraction to be greater than 0.
   !!!!!$cuf kernel do(4) <<<*,*>>>
   !!!!do k=1, nk
   !!!!   do j=1, nj
   !!!!      do i=1,ni
   !!!!         do b=1, blocks_number
   !!!!            q_gpu(b,i,j,k,6) = max(0._R8P, q_gpu(b,i,j,k,6))
   !!!!         enddo
   !!!!      enddo
   !!!!   enddo
   !!!!enddo
   !!!!!@cuf iercuda=cudaDeviceSynchronize()

   endsubroutine compute_residuals_gpu

   subroutine viscous_cuf(ni, nj, nk, ngc, blocks_number, ivis, visc_type, fd_coeff1_gpu, fd_coeff2_gpu, &
                          gamma_fluid, Prandtl, q_coeff, Lewis, Zeldovich, Damkohler, dha, &
                          q_aux_gpu, dx_gpu, dy_gpu, dz_gpu, fl_gpu)
   ! Evaluation of the viscous fluxes
   integer(I4P), intent(in)          :: ni, nj, nk, ngc, blocks_number, ivis, visc_type
   real(R8P), intent(in), device     :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::    fl_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(in), device     :: dx_gpu(1:), dy_gpu(1:), dz_gpu(1:)
   real(R8P), intent(in), device     :: fd_coeff1_gpu(1:), fd_coeff2_gpu(0:)
   real(R8P), intent(in)             :: gamma_fluid, Prandtl, Lewis, dha, Zeldovich, Damkohler, q_coeff
   integer                           :: i,j,k,l,b
   real(R8P)                         :: ccl,clapl,div3l,drmutdt
   real(R8P)                         :: ri,rmut,rmutx,rmuty,rmutz
   real(R8P)                         :: sig11,sig12,sig13
   real(R8P)                         :: sig21,sig22,sig23
   real(R8P)                         :: sig31,sig32,sig33
   real(R8P)                         :: sigq,sigx,sigy,sigz,tt,sigqt,sigah
   real(R8P)                         :: uu,vv,ww
   real(R8P)                         :: tx,ty,tz
   real(R8P)                         :: ux,uy,uz
   real(R8P)                         :: vx,vy,vz
   real(R8P)                         :: wx,wy,wz
   real(R8P)                         :: ulap,ulapx,ulapy,ulapz
   real(R8P)                         :: vlap,vlapx,vlapy,vlapz
   real(R8P)                         :: wlap,wlapx,wlapy,wlapz
   real(R8P)                         :: tlap,tlapx,tlapy,tlapz
   real(R8P)                         :: sqgmr2,sqgmr2h,tt2,sqrtt,sdivt,sdivt1
   real(R8P)                         :: tem,ya,yalap,yalapx,yalapy,yalapz,sigya,reaction_rate
   real(R8P)                         :: invdx, invdy, invdz, invdx2, invdy2, invdz2
   integer(I4P)                      :: iercuda                 !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(4) <<<*,*>>>
   do k=1, nk
      do j=1, nj
         do i=1,ni
            do b=1, blocks_number
               invdx  = 1._R8P/dx_gpu(b) ; invdy  = 1._R8P/dy_gpu(b) ; invdz  = 1._R8P/dz_gpu(b)
               invdx2 = invdx*invdx      ; invdy2 = invdy*invdy      ; invdz2 = invdz*invdz

               ! Compute first derivatives
               ux = 0._R8P ; vx = 0._R8P ; wx = 0._R8P ; tx = 0._R8P
               uy = 0._R8P ; vy = 0._R8P ; wy = 0._R8P ; ty = 0._R8P
               uz = 0._R8P ; vz = 0._R8P ; wz = 0._R8P ; tz = 0._R8P
               do l=1,ivis/2
                  ccl = fd_coeff1_gpu(l)
                  ux = ux+ccl*(q_aux_gpu(b,i+l,j,k,2)-q_aux_gpu(b,i-l,j,k,2))
                  vx = vx+ccl*(q_aux_gpu(b,i+l,j,k,3)-q_aux_gpu(b,i-l,j,k,3))
                  wx = wx+ccl*(q_aux_gpu(b,i+l,j,k,4)-q_aux_gpu(b,i-l,j,k,4))
                  tx = tx+ccl*(q_aux_gpu(b,i+l,j,k,6)-q_aux_gpu(b,i-l,j,k,6))

                  uy = uy+ccl*(q_aux_gpu(b,i,j+l,k,2)-q_aux_gpu(b,i,j-l,k,2))
                  vy = vy+ccl*(q_aux_gpu(b,i,j+l,k,3)-q_aux_gpu(b,i,j-l,k,3))
                  wy = wy+ccl*(q_aux_gpu(b,i,j+l,k,4)-q_aux_gpu(b,i,j-l,k,4))
                  ty = ty+ccl*(q_aux_gpu(b,i,j+l,k,6)-q_aux_gpu(b,i,j-l,k,6))

                  uz = uz+ccl*(q_aux_gpu(b,i,j,k+l,2)-q_aux_gpu(b,i,j,k-l,2))
                  vz = vz+ccl*(q_aux_gpu(b,i,j,k+l,3)-q_aux_gpu(b,i,j,k-l,3))
                  wz = wz+ccl*(q_aux_gpu(b,i,j,k+l,4)-q_aux_gpu(b,i,j,k-l,4))
                  tz = tz+ccl*(q_aux_gpu(b,i,j,k+l,6)-q_aux_gpu(b,i,j,k-l,6))
               enddo
               ux = ux*invdx ; vx = vx*invdx ; wx = wx*invdx ; tx = tx*invdx
               uy = uy*invdy ; vy = vy*invdy ; wy = wy*invdy ; ty = ty*invdy
               uz = uz*invdz ; vz = vz*invdz ; wz = wz*invdz ; tz = tz*invdz

               ! Compute second derivatives
               uu  = q_aux_gpu(b,i,j,k,2)    ; vv = q_aux_gpu(b,i,j,k,3) ; ww = q_aux_gpu(b,i,j,k,4)
               tem = q_aux_gpu(b,i,j,k,6)    ; ya = q_aux_gpu(b,i,j,k,5)
               ulapx  = fd_coeff2_gpu(0)*uu  ; ulapy  = ulapx  ; ulapz  = ulapx
               vlapx  = fd_coeff2_gpu(0)*vv  ; vlapy  = vlapx  ; vlapz  = vlapx
               wlapx  = fd_coeff2_gpu(0)*ww  ; wlapy  = wlapx  ; wlapz  = wlapx
               tlapx  = fd_coeff2_gpu(0)*tem ; tlapy  = tlapx  ; tlapz  = tlapx
               yalapx = fd_coeff2_gpu(0)*ya  ; yalapy = yalapx ; yalapz = yalapx
               do l=1,ivis/2
                  clapl  = fd_coeff2_gpu(l)
                  ulapx  = ulapx  + clapl*(q_aux_gpu(b,i+l,j,k,2)+q_aux_gpu(b,i-l,j,k,2))
                  ulapy  = ulapy  + clapl*(q_aux_gpu(b,i,j+l,k,2)+q_aux_gpu(b,i,j-l,k,2))
                  ulapz  = ulapz  + clapl*(q_aux_gpu(b,i,j,k+l,2)+q_aux_gpu(b,i,j,k-l,2))
                  vlapx  = vlapx  + clapl*(q_aux_gpu(b,i+l,j,k,3)+q_aux_gpu(b,i-l,j,k,3))
                  vlapy  = vlapy  + clapl*(q_aux_gpu(b,i,j+l,k,3)+q_aux_gpu(b,i,j-l,k,3))
                  vlapz  = vlapz  + clapl*(q_aux_gpu(b,i,j,k+l,3)+q_aux_gpu(b,i,j,k-l,3))
                  wlapx  = wlapx  + clapl*(q_aux_gpu(b,i+l,j,k,4)+q_aux_gpu(b,i-l,j,k,4))
                  wlapy  = wlapy  + clapl*(q_aux_gpu(b,i,j+l,k,4)+q_aux_gpu(b,i,j-l,k,4))
                  wlapz  = wlapz  + clapl*(q_aux_gpu(b,i,j,k+l,4)+q_aux_gpu(b,i,j,k-l,4))
                  tlapx  = tlapx  + clapl*(q_aux_gpu(b,i+l,j,k,6)+q_aux_gpu(b,i-l,j,k,6))
                  tlapy  = tlapy  + clapl*(q_aux_gpu(b,i,j+l,k,6)+q_aux_gpu(b,i,j-l,k,6))
                  tlapz  = tlapz  + clapl*(q_aux_gpu(b,i,j,k+l,6)+q_aux_gpu(b,i,j,k-l,6))
                  yalapx = yalapx + clapl*(q_aux_gpu(b,i+l,j,k,5)+q_aux_gpu(b,i-l,j,k,5))
                  yalapy = yalapy + clapl*(q_aux_gpu(b,i,j+l,k,5)+q_aux_gpu(b,i,j-l,k,5))
                  yalapz = yalapz + clapl*(q_aux_gpu(b,i,j,k+l,5)+q_aux_gpu(b,i,j,k-l,5))
               enddo
               ulapx  = ulapx*invdx2 ; vlapx  = vlapx*invdx2 ; wlapx  = wlapx*invdx2
               tlapx  = tlapx*invdx2 ; yalapx = yalapx*invdx2
               ulapy  = ulapy*invdy2 ; vlapy  = vlapy*invdy2 ; wlapy  = wlapy*invdy2
               tlapy  = tlapy*invdy2 ; yalapy = yalapy*invdy2
               ulapz  = ulapz*invdz2 ; vlapz  = vlapz*invdz2 ; wlapz  = wlapz*invdz2
               tlapz  = tlapz*invdz2 ; yalapz = yalapz*invdz2

               ulap  = ulapx +ulapy +ulapz ; vlap  = vlapx +vlapy +vlapz ; wlap  = wlapx +wlapy +wlapz
               tlap  = tlapx +tlapy +tlapz ; yalap = yalapx+yalapy+yalapz

               div3l = ux+vy+wz ; div3l   = div3l/3._R8P
               sig11 = 2._R8P*(ux-div3l)
               sig12 = uy+vx
               sig13 = uz+wx
               sig22 = 2._R8P*(vy-div3l)
               sig23 = vz+wy
               sig33 = 2._R8P*(wz-div3l)

               if (visc_type==0) then
                  rmut    = Prandtl
                  drmutdt = 0._R8P
               !elseif (visc_type==1) then
               !   tt      = q_aux_gpu(b,i,j,k,6)
               !   rmut    = sqgmr*tt**vtexp
               !   drmutdt = rmut*vtexp/tt
               !else
               !   tt      = q_aux_gpu(b,i,j,k,6)
               !   tt2     = tt*tt
               !   sqrtt   = sqrt(tt)
               !   sdivt   = s2tinf/tt
               !   sdivt1  = 1._R8P+sdivt
               !   rmut    = sqgmr2*sqrtt/sdivt1  ! molecular viscosity
               !   drmutdt = sqgmr2h/sqrtt+rmut*s2tinf/tt2
               !   drmutdt = drmutdt/sdivt1
               endif
               rmutx = drmutdt*tx ; rmuty = drmutdt*ty ; rmutz = drmutdt*tz
               ! Viscosity diffusion
               sigx  = rmutx*sig11 + rmuty*sig12 + rmutz*sig13 + rmut*ulap
               sigy  = rmutx*sig12 + rmuty*sig22 + rmutz*sig23 + rmut*vlap
               sigz  = rmutx*sig13 + rmuty*sig23 + rmutz*sig33 + rmut*wlap
               ! Heat conduction
               sigqt = (rmutx*tx+rmuty*ty+rmutz*tz+rmut*tlap)/Prandtl
               ! Aerodynamic heating
               sigah = (sig11*ux+sig12*uy+sig13*uz+sig12*vx+sig22*vy+sig23*vz+sig13*wx+sig23*wy+sig33*wz)*rmut ! Aerodynamic heating
               ! Specie diffusion
               sigya = yalap/Lewis
               ! Total energy diffusion
               sigq  = sigx*uu+sigy*vv+sigz*ww+sigah+sigqt+dha/Lewis*yalap
               ! Reactive term
               !reaction_rate = max(0._R8P, Damkohler*q_aux_gpu(b,i,j,k,1)*ya*exp(-Zeldovich/tem*gamma_fluid/(gamma_fluid-1._R8P)))
               reaction_rate = Damkohler*q_aux_gpu(b,i,j,k,1)*ya*exp(-Zeldovich/tem*gamma_fluid/(gamma_fluid-1._R8P))
               !reaction_rate = 0._R8P

               fl_gpu(b,i,j,k,2) = fl_gpu(b,i,j,k,2) - sigx
               fl_gpu(b,i,j,k,3) = fl_gpu(b,i,j,k,3) - sigy
               fl_gpu(b,i,j,k,4) = fl_gpu(b,i,j,k,4) - sigz
               fl_gpu(b,i,j,k,5) = fl_gpu(b,i,j,k,5) - sigq
               fl_gpu(b,i,j,k,6) = fl_gpu(b,i,j,k,6) - sigya + reaction_rate

            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()

   endsubroutine viscous_cuf

   attributes(global) subroutine euler_x_central_kernel(q_gpu, q_aux_gpu, fl_gpu, fhat_gpu, fd_conv_gpu, dx_gpu, &
                                                        blocks_number, ni, nj, nk, ngc, nv, lmax)

   real(R8P), intent(in), device     ::     q_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(in), device     :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::    fl_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::  fhat_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(in), device     :: fd_conv_gpu(1:, 1:)
   real(R8P), intent(in), device     :: dx_gpu(1:)
   integer, intent(in), value        :: blocks_number, ni, nj, nk, ngc, nv, lmax

   real(R8P)                         :: rhom, uui, vvi, wwi, ppi, enti, yai
   real(R8P)                         :: uuip, vvip, wwip, ppip, entip, yaip
   real(R8P)                         :: ft1, ft2, ft3, ft4, ft5, ft6, ft7
   real(R8P)                         :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6, uvs7, uv_part
   integer                           :: b, i, j, k, l, v, m

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   j = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if(b > blocks_number .or. j > nj) return

   do k=1,nk

      do i=0,ni ! loop on faces
         ft1 = 0._R8P ; ft2 = 0._R8P ; ft3 = 0._R8P ; ft4 = 0._R8P ; ft5 = 0._R8P ; ft6 = 0._R8P ; ft7 = 0._R8P
         do l=1,lmax
             uvs1 = 0._R8P ; uvs2 = 0._R8P ; uvs3 = 0._R8P ; uvs4 = 0._R8P ; uvs5 = 0._R8P ; uvs6 = 0._R8P ; uvs7 = 0._R8P
             do m=0,l-1
                 rhom    = q_aux_gpu(b,i-m,j,k,1) + q_aux_gpu(b,i-m+l,j,k,1)

                 uui     = q_aux_gpu(b,i-m,j,k,2)
                 vvi     = q_aux_gpu(b,i-m,j,k,3)
                 wwi     = q_aux_gpu(b,i-m,j,k,4)
                 ppi     = q_aux_gpu(b,i-m,j,k,7)
                 enti    = q_aux_gpu(b,i-m,j,k,8)
                 yai     = q_aux_gpu(b,i-m,j,k,5)

                 uuip    = q_aux_gpu(b,i-m+l,j,k,2)
                 vvip    = q_aux_gpu(b,i-m+l,j,k,3)
                 wwip    = q_aux_gpu(b,i-m+l,j,k,4)
                 ppip    = q_aux_gpu(b,i-m+l,j,k,7)
                 entip   = q_aux_gpu(b,i-m+l,j,k,8)
                 yaip    = q_aux_gpu(b,i-m+l,j,k,5)

                 uv_part = (uui+uuip) * rhom
                 uvs1    = uvs1 + uv_part * (2._R8P)
                 uvs2    = uvs2 + uv_part * (uui+uuip)
                 uvs3    = uvs3 + uv_part * (vvi+vvip)
                 uvs4    = uvs4 + uv_part * (wwi+wwip)
                 uvs5    = uvs5 + uv_part * (enti+entip)
                 uvs6    = uvs6 + (2._R8P)*(ppi+ppip)
                 uvs7    = uvs7 + uv_part * (yai+yaip)
             enddo
             ft1 = ft1 + fd_conv_gpu(l,lmax)*uvs1
             ft2 = ft2 + fd_conv_gpu(l,lmax)*uvs2
             ft3 = ft3 + fd_conv_gpu(l,lmax)*uvs3
             ft4 = ft4 + fd_conv_gpu(l,lmax)*uvs4
             ft5 = ft5 + fd_conv_gpu(l,lmax)*uvs5
             ft6 = ft6 + fd_conv_gpu(l,lmax)*uvs6
             ft7 = ft7 + fd_conv_gpu(l,lmax)*uvs7
         enddo
!
         fhat_gpu(b,i,j,k,1) = 0.25_R8P*ft1
         fhat_gpu(b,i,j,k,2) = 0.25_R8P*ft2 + 0.5_R8P*ft6
         fhat_gpu(b,i,j,k,3) = 0.25_R8P*ft3
         fhat_gpu(b,i,j,k,4) = 0.25_R8P*ft4
         fhat_gpu(b,i,j,k,5) = 0.25_R8P*ft5
         fhat_gpu(b,i,j,k,6) = 0.25_R8P*ft7
      enddo

      ! Update net flux
      do i=1,ni ! loop on inner nodes
         do v=1,nv
            fl_gpu(b,i,j,k,v) = (fhat_gpu(b,i,j,k,v)-fhat_gpu(b,i-1,j,k,v))/dx_gpu(b)
         enddo
      enddo

   enddo
   endsubroutine euler_x_central_kernel

   attributes(global) subroutine euler_y_central_kernel(q_gpu, q_aux_gpu, fl_gpu, fhat_gpu, fd_conv_gpu, dy_gpu, &
                                                        blocks_number, ni, nj, nk, ngc, nv, lmax)

   real(R8P), intent(in), device     ::     q_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(in), device     :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::    fl_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::  fhat_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(in), device     :: fd_conv_gpu(1:, 1:)
   real(R8P), intent(in), device     :: dy_gpu(1:)
   integer, intent(in), value        :: blocks_number, ni, nj, nk, ngc, nv, lmax

   real(R8P)                         :: rhom, uui, vvi, wwi, ppi, enti, yai
   real(R8P)                         :: uuip, vvip, wwip, ppip, entip, yaip
   real(R8P)                         :: ft1, ft2, ft3, ft4, ft5, ft6, ft7
   real(R8P)                         :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6, uvs7, uv_part
   integer                           :: b, i, j, k, l, v, m

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   i = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if(b > blocks_number .or. i > ni) return

   do k=1,nk

      do j=0,nj ! loop on faces
         ft1 = 0._R8P ; ft2 = 0._R8P ; ft3 = 0._R8P ; ft4 = 0._R8P ; ft5 = 0._R8P ; ft6 = 0._R8P ; ft7 = 0._R8P
         do l=1,lmax
             uvs1 = 0._R8P ; uvs2 = 0._R8P ; uvs3 = 0._R8P ; uvs4 = 0._R8P ; uvs5 = 0._R8P ; uvs6 = 0._R8P ; uvs7 = 0._R8P
             do m=0,l-1
                 rhom    = q_aux_gpu(b,i,j-m,k,1) + q_aux_gpu(b,i,j-m+l,k,1)

                 uui     = q_aux_gpu(b,i,j-m,k,2)
                 vvi     = q_aux_gpu(b,i,j-m,k,3)
                 wwi     = q_aux_gpu(b,i,j-m,k,4)
                 ppi     = q_aux_gpu(b,i,j-m,k,7)
                 enti    = q_aux_gpu(b,i,j-m,k,8)
                 yai     = q_aux_gpu(b,i,j-m,k,5)

                 uuip    = q_aux_gpu(b,i,j-m+l,k,2)
                 vvip    = q_aux_gpu(b,i,j-m+l,k,3)
                 wwip    = q_aux_gpu(b,i,j-m+l,k,4)
                 ppip    = q_aux_gpu(b,i,j-m+l,k,7)
                 entip   = q_aux_gpu(b,i,j-m+l,k,8)
                 yaip    = q_aux_gpu(b,i,j-m+l,k,5)

                 uv_part = (vvi+vvip) * rhom
                 uvs1    = uvs1 + uv_part * (2._R8P)
                 uvs2    = uvs2 + uv_part * (uui+uuip)
                 uvs3    = uvs3 + uv_part * (vvi+vvip)
                 uvs4    = uvs4 + uv_part * (wwi+wwip)
                 uvs5    = uvs5 + uv_part * (enti+entip)
                 uvs6    = uvs6 + (2._R8P)*(ppi+ppip)
                 uvs7    = uvs7 + uv_part * (yai+yaip)
             enddo
             ft1 = ft1 + fd_conv_gpu(l,lmax)*uvs1
             ft2 = ft2 + fd_conv_gpu(l,lmax)*uvs2
             ft3 = ft3 + fd_conv_gpu(l,lmax)*uvs3
             ft4 = ft4 + fd_conv_gpu(l,lmax)*uvs4
             ft5 = ft5 + fd_conv_gpu(l,lmax)*uvs5
             ft6 = ft6 + fd_conv_gpu(l,lmax)*uvs6
             ft7 = ft7 + fd_conv_gpu(l,lmax)*uvs7
         enddo
!
         fhat_gpu(b,i,j,k,1) = 0.25_R8P*ft1
         fhat_gpu(b,i,j,k,2) = 0.25_R8P*ft2
         fhat_gpu(b,i,j,k,3) = 0.25_R8P*ft3 + 0.5_R8P*ft6
         fhat_gpu(b,i,j,k,4) = 0.25_R8P*ft4
         fhat_gpu(b,i,j,k,5) = 0.25_R8P*ft5
         fhat_gpu(b,i,j,k,6) = 0.25_R8P*ft7
      enddo

      ! Update net flux
      do j=1,nj ! loop on inner nodes
         do v=1,nv
            fl_gpu(b,i,j,k,v) = fl_gpu(b,i,j,k,v) + (fhat_gpu(b,i,j,k,v)-fhat_gpu(b,i,j-1,k,v))/dy_gpu(b)
         enddo
      enddo

   enddo
   endsubroutine euler_y_central_kernel

   attributes(global) subroutine euler_z_central_kernel(q_gpu, q_aux_gpu, fl_gpu, fhat_gpu, fd_conv_gpu, dz_gpu, &
                                                        blocks_number, ni, nj, nk, ngc, nv, lmax)

   real(R8P), intent(in), device     ::     q_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(in), device     :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::    fl_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::  fhat_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(in), device     :: fd_conv_gpu(1:, 1:)
   real(R8P), intent(in), device     :: dz_gpu(1:)
   integer, intent(in), value        :: blocks_number, ni, nj, nk, ngc, nv, lmax

   real(R8P)                         :: rhom, uui, vvi, wwi, ppi, enti, yai
   real(R8P)                         :: uuip, vvip, wwip, ppip, entip, yaip
   real(R8P)                         :: ft1, ft2, ft3, ft4, ft5, ft6, ft7
   real(R8P)                         :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6, uvs7, uv_part
   integer                           :: b, i, j, k, l, v, m

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   i = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if(b > blocks_number .or. i > ni) return

   do j=1,nj

      do k=0,nk ! loop on faces
         ft1 = 0._R8P ; ft2 = 0._R8P ; ft3 = 0._R8P ; ft4 = 0._R8P ; ft5 = 0._R8P ; ft6 = 0._R8P ; ft7 = 0._R8P
         do l=1,lmax
             uvs1 = 0._R8P ; uvs2 = 0._R8P ; uvs3 = 0._R8P ; uvs4 = 0._R8P ; uvs5 = 0._R8P ; uvs6 = 0._R8P ; uvs7 = 0._R8P
             do m=0,l-1
                 rhom    = q_aux_gpu(b,i,j,k-m,1) + q_aux_gpu(b,i,j,k-m+l,1)

                 uui     = q_aux_gpu(b,i,j,k-m,2)
                 vvi     = q_aux_gpu(b,i,j,k-m,3)
                 wwi     = q_aux_gpu(b,i,j,k-m,4)
                 ppi     = q_aux_gpu(b,i,j,k-m,7)
                 enti    = q_aux_gpu(b,i,j,k-m,8)
                 yai     = q_aux_gpu(b,i,j,k-m,5)

                 uuip    = q_aux_gpu(b,i,j,k-m+l,2)
                 vvip    = q_aux_gpu(b,i,j,k-m+l,3)
                 wwip    = q_aux_gpu(b,i,j,k-m+l,4)
                 ppip    = q_aux_gpu(b,i,j,k-m+l,7)
                 entip   = q_aux_gpu(b,i,j,k-m+l,8)
                 yaip    = q_aux_gpu(b,i,j,k-m+l,5)

                 uv_part = (wwi+wwip) * rhom
                 uvs1    = uvs1 + uv_part * (2._R8P)
                 uvs2    = uvs2 + uv_part * (uui+uuip)
                 uvs3    = uvs3 + uv_part * (vvi+vvip)
                 uvs4    = uvs4 + uv_part * (wwi+wwip)
                 uvs5    = uvs5 + uv_part * (enti+entip)
                 uvs6    = uvs6 + (2._R8P)*(ppi+ppip)
                 uvs7    = uvs7 + uv_part * (yai+yaip)
             enddo
             ft1 = ft1 + fd_conv_gpu(l,lmax)*uvs1
             ft2 = ft2 + fd_conv_gpu(l,lmax)*uvs2
             ft3 = ft3 + fd_conv_gpu(l,lmax)*uvs3
             ft4 = ft4 + fd_conv_gpu(l,lmax)*uvs4
             ft5 = ft5 + fd_conv_gpu(l,lmax)*uvs5
             ft6 = ft6 + fd_conv_gpu(l,lmax)*uvs6
             ft7 = ft7 + fd_conv_gpu(l,lmax)*uvs7
         enddo
!
         fhat_gpu(b,i,j,k,1) = 0.25_R8P*ft1
         fhat_gpu(b,i,j,k,2) = 0.25_R8P*ft2
         fhat_gpu(b,i,j,k,3) = 0.25_R8P*ft3
         fhat_gpu(b,i,j,k,4) = 0.25_R8P*ft4 + 0.5_R8P*ft6
         fhat_gpu(b,i,j,k,5) = 0.25_R8P*ft5
         fhat_gpu(b,i,j,k,6) = 0.25_R8P*ft7
      enddo

      ! Update net flux
      do k=1,nk ! loop on inner nodes
         do v=1,nv
            fl_gpu(b,i,j,k,v) = fl_gpu(b,i,j,k,v) + (fhat_gpu(b,i,j,k,v)-fhat_gpu(b,i,j,k-1,v))/dz_gpu(b)
         enddo
      enddo

   enddo
   endsubroutine euler_z_central_kernel

   attributes(global) subroutine euler_x_kernel_fra(q_gpu, q_aux_gpu, fl_gpu, fhat_gpu, gplus, gminus, dx_gpu, &
                                                blocks_number, ni, nj, nk, ngc, nv, iweno, dha, gamma_fluid)

   real(R8P), intent(in), device     ::     q_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(in), device     :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::    fl_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::  fhat_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::     gplus(1:, 1:, 1:, 1:, 1:)
   real(R8P), intent(inout), device  ::    gminus(1:, 1:, 1:, 1:, 1:)
   real(R8P), intent(in), device     ::    dx_gpu(1:)
   integer, intent(in), value        :: blocks_number, ni, nj, nk, ngc, nv, iweno
   real(R8P), intent(in), value      :: dha, gamma_fluid
   integer                           :: b, i, j, k, l, ll, m, mm, v
   ! here 6 is used instead of nv to help the compiler to use registers instead of global memory
   real(R8P)                         :: er(6,6), el(6,6), ev(6), evmax(6), ghat(6), gl(6), gr(6), fi(6)
   real(R8P)                         :: uu, vv, ww, h, ya, qq, c, ci, b1, b2
   real(R8P)                         :: gc, wc

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   j = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if(b > blocks_number .or. j > nj) return

   do k=1,nk

      do i=0,iweno-2 ! loop on faces

         ! Compute Roe average
         call compute_roe_average(q_aux_gpu, dha, gamma_fluid, &
            ngc, b, i, j, k, i+1, j, k, uu, vv, ww, h, ya, qq, c, ci, b1, b2)

         ! Compute right and left eigenvectors matrices (at Roe state)
         ! call compute_eigenvectors()

         ! Find max eigenvalues on the stencil
         do m=1,nv  ! loop on characteristic fields
            evmax(m) = -1._R8P
         enddo
         do l=1,2*ngc ! LLF
            ll = i + l - iweno
            uu = q_aux_gpu(b,ll,j,k,2)
            c  = q_aux_gpu(b,ll,j,k,9)
            ev(1) = abs(uu-c) ; ev(2) = abs(uu) ; ev(3) = abs(uu+c) ; ev(4) = ev(2) ; ev(5) = ev(2) ; ev(6) = ev(2)
            do m=1,6
                evmax(m) = max(ev(m),evmax(m))
            enddo
         enddo

         ! Decompose fluxes as + and -
         do l=1,2*ngc ! loop over the stencil centered at face i
            ll = i + l - iweno
            fi(1)  = q_gpu(b,ll,j,k,2)
            fi(2)  = q_aux_gpu(b,ll,j,k,2) * q_gpu(b,ll,j,k,2) + q_aux_gpu(b,ll,j,k,7)
            fi(3)  = q_aux_gpu(b,ll,j,k,2) * q_gpu(b,ll,j,k,3)
            fi(4)  = q_aux_gpu(b,ll,j,k,2) * q_gpu(b,ll,j,k,4)
            fi(5)  = q_aux_gpu(b,ll,j,k,2) *(q_gpu(b,ll,j,k,5) + q_aux_gpu(b,ll,j,k,7))
            fi(6)  = q_aux_gpu(b,ll,j,k,2) * q_gpu(b,ll,j,k,6)
            do m=1,nv
               wc = 0._R8P
               gc = 0._R8P
               do mm=1,nv
                  wc = wc + el(mm,m) * q_gpu(b,ll,j,k,mm)
                  gc = gc + el(mm,m) * fi(mm)
               enddo
               gplus (m,l,j,k,b) = 0.5_R8P * (gc + evmax(m) * wc)
               gminus(m,l,j,k,b) = gc - gplus(m,l,j,k,b)
            enddo
         enddo

         ! Reconstruction of the + and - fluxes
         call weno_reconstruction(nv, gplus(1,1,j,k,b), gminus(1,1,j,k,b), gl, gr, iweno)
      enddo

      do i=0,ni ! loop on faces
         ! Compute Roe average
         call compute_roe_average(q_aux_gpu, dha, gamma_fluid, &
            ngc, b, i, j, k, i+1, j, k, uu, vv, ww, h, ya, qq, c, ci, b1, b2)

         ! Compute right and left eigenvectors matrices (at Roe state)
         ! call compute_eigenvectors()

         ! Find max eigenvalues on the stencil
         do m=1,nv  ! loop on characteristic fields
            evmax(m) = -1._R8P
         enddo
         do l=1,2*ngc ! LLF
            ll = i + l - iweno
            uu = q_aux_gpu(b,ll,j,k,2)
            c  = q_aux_gpu(b,ll,j,k,9)
            ev(1) = abs(uu-c) ; ev(2) = abs(uu) ; ev(3) = abs(uu+c) ; ev(4) = ev(2) ; ev(5) = ev(2) ; ev(6) = ev(2)
            do m=1,6
                evmax(m) = max(ev(m),evmax(m))
            enddo
         enddo

         ! Decompose fluxes as + and -
         do l=1,2*ngc ! loop over the stencil centered at face i
            ll = i + l - iweno
            fi(1)  = q_gpu(b,ll,j,k,2)
            fi(2)  = q_aux_gpu(b,ll,j,k,2) * q_gpu(b,ll,j,k,2) + q_aux_gpu(b,ll,j,k,7)
            fi(3)  = q_aux_gpu(b,ll,j,k,2) * q_gpu(b,ll,j,k,3)
            fi(4)  = q_aux_gpu(b,ll,j,k,2) * q_gpu(b,ll,j,k,4)
            fi(5)  = q_aux_gpu(b,ll,j,k,2) *(q_gpu(b,ll,j,k,5) + q_aux_gpu(b,ll,j,k,7))
            fi(6)  = q_aux_gpu(b,ll,j,k,2) * q_gpu(b,ll,j,k,6)
            do m=1,nv
               wc = 0._R8P
               gc = 0._R8P
               do mm=1,nv
                  wc = wc + el(mm,m) * q_gpu(b,ll,j,k,mm)
                  gc = gc + el(mm,m) * fi(mm)
               enddo
               gplus (m,l,j,k,b) = 0.5_R8P * (gc + evmax(m) * wc)
               gminus(m,l,j,k,b) = gc - gplus(m,l,j,k,b)
            enddo
         enddo

         ! Reconstruction of the + and - fluxes
         ! call weno_reconstruction(nv, gplus(1,1,j,k,b), gminus(1,1,j,k,b), gl, gr, iweno, beta)

         ! Reassemble + and - characteristic fluxes
         do m=1,6
            ghat(m) = gl(m) + gr(m)
         enddo

         ! Return to conservative fluxes
         do m=1,nv
            fhat_gpu(b,i,j,k,m) = 0._R8P
            do mm=1,nv
               fhat_gpu(b,i,j,k,m) = fhat_gpu(b,i,j,k,m) + er(mm,m) * ghat(mm)
            enddo
         enddo

      enddo

      ! Update net flux
      do i=1,ni ! loop on inner nodes
         do v=1,nv
            fl_gpu(b,i,j,k,v) = (fhat_gpu(b,i,j,k,v)-fhat_gpu(b,i-1,j,k,v))/dx_gpu(b)
         enddo
      enddo

   enddo

   endsubroutine euler_x_kernel_fra

   attributes(global) subroutine euler_x_kernel_stefano(q_gpu, q_aux_gpu, fl_gpu, fhat_gpu, gplus, gminus, dx_gpu, &
                                                blocks_number, ni, nj, nk, ngc, nv, iweno, dha, gamma_fluid)

   real(R8P), intent(in), device     ::     q_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(in), device     :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::    fl_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::  fhat_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::     gplus(1:, 1:, 1:, 1:, 1:)
   real(R8P), intent(inout), device  ::    gminus(1:, 1:, 1:, 1:, 1:)
   real(R8P), intent(in), device     ::    dx_gpu(1:)
   integer, intent(in), value        :: blocks_number, ni, nj, nk, ngc, nv, iweno
   real(R8P), intent(in), value      :: dha, gamma_fluid
   integer                           :: b, i, j, k, l, ll, m, mm, v
   ! here 6 is used instead of nv to help the compiler to use registers instead of global memory
   real(R8P)                         :: er(6,6), el(6,6), ev(6), evmax(6), ghat(6), gl(6), gr(6), fi(6)
   real(R8P)                         :: uu, vv, ww, h, ya, qq, c, ci, b1, b2
   real(R8P)                         :: gc, wc

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   j = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if(b > blocks_number .or. j > nj) return

   do k=1,nk

      do i=0,ni ! loop on faces

         ! Compute Roe average
         call compute_roe_average(q_aux_gpu, dha, gamma_fluid, &
            ngc, b, i, j, k, i+1, j, k, uu, vv, ww, h, ya, qq, c, ci, b1, b2)

         ! Compute right and left eigenvectors matrices (at Roe state)
         ! call compute_eigenvectors()

         ! Find max eigenvalues on the stencil
         do m=1,nv  ! loop on characteristic fields
            evmax(m) = -1._R8P
         enddo
         do l=1,2*ngc ! LLF
            ll = i + l - iweno
            uu = q_aux_gpu(b,ll,j,k,2)
            c  = q_aux_gpu(b,ll,j,k,9)
            ev(1) = abs(uu-c) ; ev(2) = abs(uu) ; ev(3) = abs(uu+c) ; ev(4) = ev(2) ; ev(5) = ev(2) ; ev(6) = ev(2)
            do m=1,6
                evmax(m) = max(ev(m),evmax(m))
            enddo
         enddo

         ! Decompose fluxes as + and -
         do l=1,2*ngc ! loop over the stencil centered at face i
            ll = i + l - iweno
            fi(1)  = q_gpu(b,ll,j,k,2)
            fi(2)  = q_aux_gpu(b,ll,j,k,2) * q_gpu(b,ll,j,k,2) + q_aux_gpu(b,ll,j,k,7)
            fi(3)  = q_aux_gpu(b,ll,j,k,2) * q_gpu(b,ll,j,k,3)
            fi(4)  = q_aux_gpu(b,ll,j,k,2) * q_gpu(b,ll,j,k,4)
            fi(5)  = q_aux_gpu(b,ll,j,k,2) *(q_gpu(b,ll,j,k,5) + q_aux_gpu(b,ll,j,k,7))
            fi(6)  = q_aux_gpu(b,ll,j,k,2) * q_gpu(b,ll,j,k,6)
            do m=1,nv
               wc = 0._R8P
               gc = 0._R8P
               do mm=1,nv
                  wc = wc + el(mm,m) * q_gpu(b,ll,j,k,mm)
                  gc = gc + el(mm,m) * fi(mm)
               enddo
               gplus (m,l,j,k,b) = 0.5_R8P * (gc + evmax(m) * wc)
               gminus(m,l,j,k,b) = gc - gplus(m,l,j,k,b)
            enddo
         enddo

         ! Reconstruction of the + and - fluxes
         ! call weno_beta(nv, gplus(1,1,j,k,b), gminus(1,1,j,k,b), betap(1,1,j,k,b), betam(1,1,j,k,b), gr, iweno)

      enddo

      do i=0,ni
         ! call weno_reconstruction(nv, gplus(1,1,j,k,b), gminus(1,1,j,k,b), gl, gr, iweno, betap, betam)

         ! Reassemble + and - characteristic fluxes
         do m=1,6
            ghat(m) = gl(m) + gr(m)
         enddo

         ! Return to conservative fluxes
         do m=1,nv
            fhat_gpu(b,i,j,k,m) = 0._R8P
            do mm=1,nv
               fhat_gpu(b,i,j,k,m) = fhat_gpu(b,i,j,k,m) + er(mm,m) * ghat(mm)
            enddo
         enddo
      enddo

      ! Update net flux
      do i=1,ni ! loop on inner nodes
         do v=1,nv
            fl_gpu(b,i,j,k,v) = (fhat_gpu(b,i,j,k,v)-fhat_gpu(b,i-1,j,k,v))/dx_gpu(b)
         enddo
      enddo

   enddo

   endsubroutine euler_x_kernel_stefano

   attributes(global) subroutine euler_x_kernel(q_gpu, q_aux_gpu, fl_gpu, fhat_gpu, gplus, gminus, dx_gpu, &
                                                blocks_number, ni, nj, nk, ngc, nv, iweno, dha, gamma_fluid)

   real(R8P), intent(in), device     ::     q_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(in), device     :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::    fl_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::  fhat_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::     gplus(1:, 1:, 1:, 1:, 1:)
   real(R8P), intent(inout), device  ::    gminus(1:, 1:, 1:, 1:, 1:)
   real(R8P), intent(in), device     ::    dx_gpu(1:)
   integer, intent(in), value        :: blocks_number, ni, nj, nk, ngc, nv, iweno
   real(R8P), intent(in), value      :: dha, gamma_fluid
   integer                           :: b, i, j, k, l, ll, m, mm, v
   ! here 6 is used instead of nv to help the compiler to use registers instead of global memory
   real(R8P)                         :: er(6,6), el(6,6), ev(6), evmax(6), ghat(6), gl(6), gr(6), fi(6), vi(6)
   real(R8P)                         :: uu, vv, ww, h, ya, qq, c, ci, b1, b2
   real(R8P)                         :: gc, wc

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   j = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if(b > blocks_number .or. j > nj) return

   do k=1,nk

      do i=0,ni ! loop on faces

         ! Compute Roe average
         call compute_roe_average(q_aux_gpu, dha, gamma_fluid, &
            ngc, b, i, j, k, i+1, j, k, uu, vv, ww, h, ya, qq, c, ci, b1, b2)

         ! Compute right and left eigenvectors matrices (at Roe state)
         er(1,1) = 1._R8P ;  er(1,2) = uu-c    ; er(1,3) = vv     ; er(1,4) = ww     ; er(1,5) = h-uu*c ; er(1,6) = ya
         er(2,1) = 1._R8P ;  er(2,2) = uu      ; er(2,3) = vv     ; er(2,4) = ww     ; er(2,5) = qq     ; er(2,6) = 0._R8P
         er(3,1) = 1._R8P ;  er(3,2) = uu+c    ; er(3,3) = vv     ; er(3,4) = ww     ; er(3,5) = h+uu*c ; er(3,6) = ya
         er(4,1) = 0._R8P ;  er(4,2) = 0._R8P  ; er(4,3) = 1._R8P ; er(4,4) = 0._R8P ; er(4,5) = 0._R8P ; er(4,6) = -vv/dha
         er(5,1) = 0._R8P ;  er(5,2) = 0._R8P  ; er(5,3) = 0._R8P ; er(5,4) = 1._R8P ; er(5,5) = 0._R8P ; er(5,6) = -ww/dha
         er(6,1) = 0._R8P ;  er(6,2) = 0._R8P  ; er(6,3) = 0._R8P ; er(6,4) = 0._R8P ; er(6,5) = 1._R8P ; er(6,6) = 1._R8P/dha

         el(1,1) =  0.5_R8P*(b1+uu*ci) ; el(1,2) = 1._R8P-b1 ; el(1,3) =  0.5_R8P*(b1-uu*ci)
         el(2,1) = -0.5_R8P*(b2*uu+ci) ; el(2,2) = b2*uu     ; el(2,3) = -0.5_R8P*(b2*uu-ci)
         el(3,1) = -0.5_R8P*(b2*vv   ) ; el(3,2) = b2*vv     ; el(3,3) = -0.5_R8P*(b2*vv   )
         el(4,1) = -0.5_R8P*(b2*ww   ) ; el(4,2) = b2*ww     ; el(4,3) = -0.5_R8P*(b2*ww   )
         el(5,1) =  0.5_R8P*b2         ; el(5,2) = -b2       ; el(5,3) =  0.5_R8P*b2
         el(6,1) = -0.5_R8P*b2*dha     ; el(6,2) = b2*dha    ; el(6,3) = -0.5_R8P*b2*dha

         el(1,4) = -vv                 ; el(1,5) = -ww       ; el(1,6) = -ya*dha*b1-vv**2-ww**2
         el(2,4) = 0._R8P              ; el(2,5) = 0._R8P    ; el(2,6) = uu*ya*dha*b2
         el(3,4) = 1._R8P              ; el(3,5) = 0._R8P    ; el(3,6) = vv*(1._R8P+ya*dha*b2)
         el(4,4) = 0._R8P              ; el(4,5) = 1._R8P    ; el(4,6) = ww*(1._R8P+ya*dha*b2)
         el(5,4) = 0._R8P              ; el(5,5) = 0._R8P    ; el(5,6) = -ya*dha*b2
         el(6,4) = 0._R8P              ; el(6,5) = 0._R8P    ; el(6,6) = dha*(1._R8P+ya*dha*b2)

         ! Find max eigenvalues on the stencil
         do m=1,6  ! loop on characteristic fields
            evmax(m) = -1._R8P
         enddo
         do l=1,2*iweno ! LLF
            ll = i + l - iweno
            uu = q_aux_gpu(b,ll,j,k,2)
            c  = q_aux_gpu(b,ll,j,k,9)
            ev(1) = abs(uu-c) ; ev(2) = abs(uu) ; ev(3) = abs(uu+c) ; ev(4) = ev(2) ; ev(5) = ev(2) ; ev(6) = ev(2)
            do m=1,6
                evmax(m) = max(ev(m),evmax(m))
            enddo
         enddo

         ! Decompose fluxes as + and -
         do l=1,2*iweno ! loop over the stencil centered at face i
            ll = i + l - iweno
            vi(1) = q_aux_gpu(b,ll,j,k,1)
            vi(2) = vi(1)*q_aux_gpu(b,ll,j,k,2)
            vi(3) = vi(1)*q_aux_gpu(b,ll,j,k,3)
            vi(4) = vi(1)*q_aux_gpu(b,ll,j,k,4)
            vi(5) = vi(1)*(1._R8P/gamma_fluid*q_aux_gpu(b,ll,j,k,6)+&
                0.5_R8P*(q_aux_gpu(b,ll,j,k,2)**2+q_aux_gpu(b,ll,j,k,3)**2+q_aux_gpu(b,ll,j,k,4)**2)+&
                q_aux_gpu(b,ll,j,k,5)*dha)
            vi(6) = vi(1)*q_aux_gpu(b,ll,j,k,5)
            fi(1) = vi(2)
            fi(2) = fi(1) * q_aux_gpu(b,ll,j,k,2) + q_aux_gpu(b,ll,j,k,7)
            fi(3) = fi(1) * q_aux_gpu(b,ll,j,k,3)
            fi(4) = fi(1) * q_aux_gpu(b,ll,j,k,4)
            fi(5) = fi(1) * vi(5) / vi(1) + q_aux_gpu(b,ll,j,k,7)*q_aux_gpu(b,ll,j,k,2)
            fi(6) = fi(1) * q_aux_gpu(b,ll,j,k,5)
            do m=1,6
               wc = 0._R8P
               gc = 0._R8P
               do mm=1,6
                  wc = wc + el(mm,m) * vi(mm)
                  gc = gc + el(mm,m) * fi(mm)
               enddo
               gplus (m,l,j,k,b) = 0.5_R8P * (gc + evmax(m) * wc)
               gminus(m,l,j,k,b) = gc - gplus(m,l,j,k,b)
            enddo
         enddo

         ! Reconstruction of the + and - fluxes
         call weno_reconstruction(6, gplus(1,1,j,k,b), gminus(1,1,j,k,b), gl, gr, iweno)

         ! Reassemble + and - characteristic fluxes
         do m=1,6
            ghat(m) = gl(m) + gr(m)
         enddo

         ! Return to conservative fluxes
         do m=1,6
            fhat_gpu(b,i,j,k,m) = 0._R8P
            do mm=1,6
               fhat_gpu(b,i,j,k,m) = fhat_gpu(b,i,j,k,m) + er(mm,m) * ghat(mm)
            enddo
         enddo

      enddo

      ! Update net flux
      do i=1,ni ! loop on inner nodes
         do v=1,6
            fl_gpu(b,i,j,k,v) = (fhat_gpu(b,i,j,k,v)-fhat_gpu(b,i-1,j,k,v))/dx_gpu(b)
         enddo
      enddo

   enddo

   endsubroutine euler_x_kernel

   attributes(global) subroutine euler_y_kernel(q_gpu, q_aux_gpu, fl_gpu, fhat_gpu, gplus, gminus, dy_gpu, &
                                                blocks_number, ni, nj, nk, ngc, nv, iweno, dha, gamma_fluid)

   real(R8P), intent(in), device     ::     q_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(in), device     :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::    fl_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::  fhat_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::     gplus(1:, 1:, 1:, 1:, 1:)
   real(R8P), intent(inout), device  ::    gminus(1:, 1:, 1:, 1:, 1:)
   real(R8P), intent(in), device     ::    dy_gpu(1:)
   integer, intent(in), value        :: blocks_number, ni, nj, nk, ngc, nv, iweno
   real(R8P), intent(in), value      :: dha, gamma_fluid
   integer                           :: b, i, j, k, l, ll, m, mm, v
   ! here 6 is used instead of nv to help the compiler to use registers instead of global memory
   real(R8P)                         :: er(6,6), el(6,6), ev(6), evmax(6), ghat(6), gl(6), gr(6), fi(6), vi(6)
   real(R8P)                         :: uu, vv, ww, h, ya, qq, c, ci, b1, b2
   real(R8P)                         :: gc, wc

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   i = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if(b > blocks_number .or. i > ni) return

   do k=1,nk

      do j=0,nj ! loop on faces

         ! Compute Roe average
         call compute_roe_average(q_aux_gpu, dha, gamma_fluid, &
            ngc, b, i, j, k, i, j+1, k, uu, vv, ww, h, ya, qq, c, ci, b1, b2)

         ! Compute right and left eigenvectors matrices (at Roe state)
         er(1,1) = 1._R8P ;  er(1,2) = uu      ; er(1,3) = vv-c   ; er(1,4) = ww     ; er(1,5) = h-vv*c ; er(1,6) = ya
         er(2,1) = 1._R8P ;  er(2,2) = uu      ; er(2,3) = vv     ; er(2,4) = ww     ; er(2,5) = qq     ; er(2,6) = 0._R8P
         er(3,1) = 1._R8P ;  er(3,2) = uu      ; er(3,3) = vv+c   ; er(3,4) = ww     ; er(3,5) = h+vv*c ; er(3,6) = ya
         er(4,1) = 0._R8P ;  er(4,2) = 1._R8P  ; er(4,3) = 0._R8P ; er(4,4) = 0._R8P ; er(4,5) = 0._R8P ; er(4,6) = -uu/dha
         er(5,1) = 0._R8P ;  er(5,2) = 0._R8P  ; er(5,3) = 0._R8P ; er(5,4) = 1._R8P ; er(5,5) = 0._R8P ; er(5,6) = -ww/dha
         er(6,1) = 0._R8P ;  er(6,2) = 0._R8P  ; er(6,3) = 0._R8P ; er(6,4) = 0._R8P ; er(6,5) = 1._R8P ; er(6,6) = 1._R8P/dha

         el(1,1) =  0.5_R8P*(b1+vv*ci) ; el(1,2) = 1._R8P-b1 ; el(1,3) = 0.5_R8P*(b1-vv*ci)
         el(2,1) = -0.5_R8P*(b2*uu)    ; el(2,2) = b2*uu     ; el(2,3) = -0.5_R8P*(b2*uu)
         el(3,1) = -0.5_R8P*(b2*vv+ci) ; el(3,2) = b2*vv     ; el(3,3) = -0.5_R8P*(b2*vv-ci)
         el(4,1) = -0.5_R8P*(b2*ww)    ; el(4,2) = b2*ww     ; el(4,3) = -0.5_R8P*(b2*ww)
         el(5,1) =  0.5_R8P*b2         ; el(5,2) = -b2       ; el(5,3) = 0.5_R8P*b2
         el(6,1) = -0.5_R8P*dha*b2     ; el(6,2) = dha*b2    ; el(6,3) = -0.5_R8P*dha*b2

         el(1,4) = -ww                 ; el(1,5) = uu        ; el(1,6) = -ya*dha*b1-uu**2-ww**2
         el(2,4) = 0._R8P              ; el(2,5) = -1._R8P   ; el(2,6) = uu*(1._R8P+ya*dha*b2)
         el(3,4) = 0._R8P              ; el(3,5) = 0._R8P    ; el(3,6) = ya*dha*vv*b2
         el(4,4) = 1._R8P              ; el(4,5) = 0._R8P    ; el(4,6) = ww*(1._R8P+ya*dha*b2)
         el(5,4) = 0._R8P              ; el(5,5) = 0._R8P    ; el(5,6) = -ya*dha*b2
         el(6,4) = 0._R8P              ; el(6,5) = 0._R8P    ; el(6,6) = dha*(1._R8P+ya*dha*b2)

         ! Find max eigenvalues on the stencil
         do m=1,6  ! loop on characteristic fields
            evmax(m) = -1._R8P
         enddo
         do l=1,2*iweno ! LLF
            ll = j + l - iweno
            uu = q_aux_gpu(b,i,ll,k,3)
            c  = q_aux_gpu(b,i,ll,k,9)
            ev(1) = abs(uu-c) ; ev(2) = abs(uu) ; ev(3) = abs(uu+c) ; ev(4) = ev(2) ; ev(5) = ev(2) ; ev(6) = ev(2)
            do m=1,6
                evmax(m) = max(ev(m),evmax(m))
            enddo
         enddo

         ! Decompose fluxes as + and -
         do l=1,2*iweno ! loop over the stencil centered at face i
            ll = j + l - iweno
            vi(1) = q_aux_gpu(b,i,ll,k,1)
            vi(2) = vi(1)*q_aux_gpu(b,i,ll,k,2)
            vi(3) = vi(1)*q_aux_gpu(b,i,ll,k,3)
            vi(4) = vi(1)*q_aux_gpu(b,i,ll,k,4)
            vi(5) = vi(1)*(1._R8P/gamma_fluid*q_aux_gpu(b,i,ll,k,6)+&
                0.5_R8P*(q_aux_gpu(b,i,ll,k,2)**2+q_aux_gpu(b,i,ll,k,3)**2+q_aux_gpu(b,i,ll,k,4)**2)+&
                q_aux_gpu(b,i,ll,k,5)*dha)
            vi(6) = vi(1)*q_aux_gpu(b,i,ll,k,5)
            fi(1) = vi(3)
            fi(2) = fi(1) * q_aux_gpu(b,i,ll,k,2)
            fi(3) = fi(1) * q_aux_gpu(b,i,ll,k,3) + q_aux_gpu(b,i,ll,k,7)
            fi(4) = fi(1) * q_aux_gpu(b,i,ll,k,4)
            fi(5) = fi(1) * vi(5) / vi(1) + q_aux_gpu(b,i,ll,k,7)*q_aux_gpu(b,i,ll,k,3)
            fi(6) = fi(1) * q_aux_gpu(b,i,ll,k,5)
            do m=1,6
               wc = 0._R8P
               gc = 0._R8P
               do mm=1,6
                  wc = wc + el(mm,m) * vi(mm)
                  gc = gc + el(mm,m) * fi(mm)
               enddo
               gplus (m,l,i,k,b) = 0.5_R8P * (gc + evmax(m) * wc)
               gminus(m,l,i,k,b) = gc - gplus(m,l,i,k,b)
            enddo
         enddo

         ! Reconstruction of the + and - fluxes
         call weno_reconstruction(6, gplus(1,1,i,k,b), gminus(1,1,i,k,b), gl, gr, iweno)

         ! Reassemble + and - characteristic fluxes
         do m=1,6
            ghat(m) = gl(m) + gr(m)
         enddo

         ! Return to conservative fluxes
         do m=1,6
            fhat_gpu(b,i,j,k,m) = 0._R8P
            do mm=1,6
               fhat_gpu(b,i,j,k,m) = fhat_gpu(b,i,j,k,m) + er(mm,m) * ghat(mm)
            enddo
         enddo

      enddo

      ! Update net flux
      do j=1,nj ! loop on inner nodes
         do v=1,6
            fl_gpu(b,i,j,k,v) = fl_gpu(b,i,j,k,v) + (fhat_gpu(b,i,j,k,v)-fhat_gpu(b,i,j-1,k,v))/dy_gpu(b)
         enddo
      enddo

   enddo

   endsubroutine euler_y_kernel

   attributes(global) subroutine euler_z_kernel(q_gpu, q_aux_gpu, fl_gpu, fhat_gpu, gplus, gminus, dz_gpu, &
                                                blocks_number, ni, nj, nk, ngc, nv, iweno, dha, gamma_fluid)

   real(R8P), intent(in), device     ::     q_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(in), device     :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::    fl_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::  fhat_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(inout), device  ::     gplus(1:, 1:, 1:, 1:, 1:)
   real(R8P), intent(inout), device  ::    gminus(1:, 1:, 1:, 1:, 1:)
   real(R8P), intent(in), device     ::    dz_gpu(1:)
   integer, intent(in), value        :: blocks_number, ni, nj, nk, ngc, nv, iweno
   real(R8P), intent(in), value      :: dha, gamma_fluid
   integer                           :: b, i, j, k, l, ll, m, mm, v
   ! here 6 is used instead of nv to help the compiler to use registers instead of global memory
   real(R8P)                         :: er(6,6), el(6,6), ev(6), evmax(6), ghat(6), gl(6), gr(6), fi(6), vi(6)
   real(R8P)                         :: uu, vv, ww, h, ya, qq, c, ci, b1, b2
   real(R8P)                         :: gc, wc

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   i = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if(b > blocks_number .or. i > ni) return

   do j=1,nj

      do k=0,nk ! loop on faces

         ! Compute Roe average
         call compute_roe_average(q_aux_gpu, dha, gamma_fluid, &
            ngc, b, i, j, k, i, j, k+1, uu, vv, ww, h, ya, qq, c, ci, b1, b2)

         ! Compute right and left eigenvectors matrices (at Roe state)
         er(1,1) = 1._R8P ;  er(1,2) = uu      ; er(1,3) = vv     ; er(1,4) = ww-c    ; er(1,5) = h-ww*c ; er(1,6) = ya
         er(2,1) = 1._R8P ;  er(2,2) = uu      ; er(2,3) = vv     ; er(2,4) = ww      ; er(2,5) = qq     ; er(2,6) = 0._R8P
         er(3,1) = 1._R8P ;  er(3,2) = uu      ; er(3,3) = vv     ; er(3,4) = ww+c    ; er(3,5) = h+ww*c ; er(3,6) = ya
         er(4,1) = 0._R8P ;  er(4,2) = 1._R8P  ; er(4,3) = 0._R8P ; er(4,4) = 0._R8P  ; er(4,5) = 0._R8P ; er(4,6) = -uu/dha
         er(5,1) = 0._R8P ;  er(5,2) = 0._R8P  ; er(5,3) = 1._R8P ; er(5,4) = 0._R8P  ; er(5,5) = 0._R8P ; er(5,6) = -vv/dha
         er(6,1) = 0._R8P ;  er(6,2) = 0._R8P  ; er(6,3) = 0._R8P ; er(6,4) = 0._R8P  ; er(6,5) = 1._R8P ; er(6,6) = 1._R8P/dha

         el(1,1) = 0.5_R8P*(b1+ww*ci)          ; el(1,2) = 1._R8P-b1 ; el(1,3) = 0.5_R8P*(b1-ww*ci)
         el(2,1) = -0.5_R8P*(b2*uu)            ; el(2,2) = b2*uu     ; el(2,3) = -0.5_R8P*(b2*uu)
         el(3,1) = -0.5_R8P*(b2*vv)            ; el(3,2) = b2*vv     ; el(3,3) = -0.5_R8P*(b2*vv)
         el(4,1) = -0.5_R8P*(b2*ww+ci)         ; el(4,2) = b2*ww     ; el(4,3) = -0.5_R8P*(b2*ww-ci)
         el(5,1) = 0.5_R8P*b2                  ; el(5,2) = -b2       ; el(5,3) = 0.5_R8P*b2
         el(6,1) = -0.5_R8P*dha*b2             ; el(6,2) = dha*b2    ; el(6,3) = -0.5_R8P*dha*b2

         el(1,4) = -uu                         ; el(1,5) = -vv       ; el(1,6) = -ya*dha*b1-uu**2-vv**2
         el(2,4) = 1._R8P                      ; el(2,5) = 0._R8P    ; el(2,6) = uu*(1+ya*dha*b2)
         el(3,4) = 0._R8P                      ; el(3,5) = 1._R8P    ; el(3,6) = vv*(1+ya*dha*b2)
         el(4,4) = 0._R8P                      ; el(4,5) = 0._R8P    ; el(4,6) = ya*dha*ww*b2
         el(5,4) = 0._R8P                      ; el(5,5) = 0._R8P    ; el(5,6) = -ya*dha*b2
         el(6,4) = 0._R8P                      ; el(6,5) = 0._R8P    ; el(6,6) = dha*(1._R8P+ya*dha*b2)

         ! Find max eigenvalues on the stencil
         do m=1,6  ! loop on characteristic fields
            evmax(m) = -1._R8P
         enddo
         do l=1,2*iweno ! LLF
            ll = k + l - iweno
            uu = q_aux_gpu(b,i,j,ll,4)
            c  = q_aux_gpu(b,i,j,ll,9)
            ev(1) = abs(uu-c) ; ev(2) = abs(uu) ; ev(3) = abs(uu+c) ; ev(4) = ev(2) ; ev(5) = ev(2) ; ev(6) = ev(2)
            do m=1,6
                evmax(m) = max(ev(m),evmax(m))
            enddo
         enddo

         ! Decompose fluxes as + and -
         do l=1,2*iweno ! loop over the stencil centered at face i
            ll = k + l - iweno
            vi(1) = q_aux_gpu(b,i,j,ll,1)
            vi(2) = vi(1)*q_aux_gpu(b,i,j,ll,2)
            vi(3) = vi(1)*q_aux_gpu(b,i,j,ll,3)
            vi(4) = vi(1)*q_aux_gpu(b,i,j,ll,4)
            vi(5) = vi(1)*(1._R8P/gamma_fluid*q_aux_gpu(b,i,j,ll,6)+&
                0.5_R8P*(q_aux_gpu(b,i,j,ll,2)**2+q_aux_gpu(b,i,j,ll,3)**2+q_aux_gpu(b,i,j,ll,4)**2)+&
                q_aux_gpu(b,i,j,ll,5)*dha)
            vi(6) = vi(1)*q_aux_gpu(b,i,j,ll,5)
            fi(1) = vi(4)
            fi(2) = fi(1) * q_aux_gpu(b,i,j,ll,2)
            fi(3) = fi(1) * q_aux_gpu(b,i,j,ll,3)
            fi(4) = fi(1) * q_aux_gpu(b,i,j,ll,4) + q_aux_gpu(b,i,j,ll,7)
            fi(5) = fi(1) * vi(5) / vi(1) + q_aux_gpu(b,i,j,ll,7)*q_aux_gpu(b,i,j,ll,4)
            fi(6) = fi(1) * q_aux_gpu(b,i,j,ll,5)
            do m=1,6
               wc = 0._R8P
               gc = 0._R8P
               do mm=1,6
                  wc = wc + el(mm,m) * vi(mm)
                  gc = gc + el(mm,m) * fi(mm)
               enddo
               gplus (m,l,i,j,b) = 0.5_R8P * (gc + evmax(m) * wc)
               gminus(m,l,i,j,b) = gc - gplus(m,l,i,j,b)
            enddo
         enddo

         ! Reconstruction of the + and - fluxes
         call weno_reconstruction(6, gplus(1,1,i,j,b), gminus(1,1,i,j,b), gl, gr, iweno)

         ! Reassemble + and - characteristic fluxes
         do m=1,6
            ghat(m) = gl(m) + gr(m)
         enddo

         ! Return to conservative fluxes
         do m=1,6
            fhat_gpu(b,i,j,k,m) = 0._R8P
            do mm=1,6
               fhat_gpu(b,i,j,k,m) = fhat_gpu(b,i,j,k,m) + er(mm,m) * ghat(mm)
            enddo
         enddo

      enddo

      ! Update net flux
      do k=1,nk ! loop on inner nodes
         do v=1,6
            fl_gpu(b,i,j,k,v) = fl_gpu(b,i,j,k,v) + (fhat_gpu(b,i,j,k,v)-fhat_gpu(b,i,j,k-1,v))/dz_gpu(b)
         enddo
      enddo

   enddo

   endsubroutine euler_z_kernel

   attributes(device) subroutine compute_roe_average(q_aux_gpu, dha, gamma_fluid, &
      ngc, b, i, j, k, ip, jp, kp, uu, vv, ww, h, ya, qq, c, ci, b1, b2)

   implicit none
   real(R8P), intent(in), device  :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P), intent(in)          :: dha, gamma_fluid
   integer(I4P), intent(in)       :: ngc, b, i, j, k, ip, jp, kp
   real(R8P), intent(out)         :: uu, vv, ww, h, ya, qq, c, ci, b1, b2
   real(R8P)                      :: ri, up, vp, wp, hp, yap, r, rp1, cc
   ! Left state (node i)
   ri        =  1._R8P/q_aux_gpu(b,i,j,k,1)
   uu        =  q_aux_gpu(b,i,j,k,2)
   vv        =  q_aux_gpu(b,i,j,k,3)
   ww        =  q_aux_gpu(b,i,j,k,4)
   h         =  q_aux_gpu(b,i,j,k,8)
   ya        =  q_aux_gpu(b,i,j,k,5)
   ! Right state (node i+1)
   up        =  q_aux_gpu(b,ip,jp,kp,2)
   vp        =  q_aux_gpu(b,ip,jp,kp,3)
   wp        =  q_aux_gpu(b,ip,jp,kp,4)
   hp        =  q_aux_gpu(b,ip,jp,kp,8)
   yap       =  q_aux_gpu(b,ip,jp,kp,5)
   ! Average state
   r         =  sqrt(q_aux_gpu(b,ip,jp,kp,1)*ri)
   rp1       =  1._R8P/(r+1._R8P)
   uu        =  (r*up+uu)*rp1
   vv        =  (r*vp+vv)*rp1
   ww        =  (r*wp+ww)*rp1
   h         =  (r*hp+h)*rp1
   ya        =  (r*yap+ya)*rp1
   qq        =  0.5_R8P * (uu*uu+vv*vv+ww*ww)
   cc        =  (gamma_fluid-1._R8P) * (h - qq - ya*dha)
   c         =  sqrt(cc)
   ci        =  1._R8P/c
   b2        = (gamma_fluid-1._R8P)/cc  ! alias 1/theta
   b1        = b2 * qq                  ! alias q/theta

   endsubroutine compute_roe_average

   attributes(device) subroutine weno_reconstruction(nvar,vp,vm,vminus,vplus,iweno)

   implicit none
   integer, intent(in)                             :: nvar, iweno
   !real(R8P), dimension(nvar,2*iweno), intent(in)  :: vm,vp
   real(R8P), dimension(1:nvar,1:*) :: vm,vp
   real(R8P), dimension(nvar), intent(out)         :: vminus,vplus

   real(R8P), dimension(-1:4) :: dwe               ! linear weights
   real(R8P), dimension(-1:4) :: alfp,alfm         ! alpha_l
   real(R8P), dimension(-1:4) :: alfp_map,alfm_map ! alpha_l
   real(R8P), dimension(-1:4) :: betap,betam       ! beta_l
   real(R8P), dimension(-1:4) :: omp,omm           ! WENO weights
   integer                    :: r,i,j,k,l,m
   real(R8P)                  :: c0,c1,c2,c3,c4,d0,d1,d2,d3,d4,summ,sump
   real(R8P)                  :: x,y,y2

   if (iweno==1) then ! Godunov

       i = iweno ! index of intermediate node to perform reconstruction

       vminus(1:nvar) = vp(1:nvar,i)
       vplus (1:nvar) = vm(1:nvar,i+1)

   elseif (iweno==2) then ! WENO-3

       i = iweno ! index of intermediate node to perform reconstruction

       dwe(1)   = 2._R8P/3._R8P
       dwe(0)   = 1._R8P/3._R8P

       do m=1,nvar

           betap(0)  = (vp(m,i  )-vp(m,i-1))**2
           betap(1)  = (vp(m,i+1)-vp(m,i  ))**2
           betam(0)  = (vm(m,i+2)-vm(m,i+1))**2
           betam(1)  = (vm(m,i+1)-vm(m,i  ))**2

           sump = 0._R8P
           summ = 0._R8P
           do l=0,1
               alfp(l) = dwe(l)/(0.000001_R8P+betap(l))**2
               alfm(l) = dwe(l)/(0.000001_R8P+betam(l))**2
               sump = sump + alfp(l)
               summ = summ + alfm(l)
           enddo
           do l=0,1
               omp(l) = alfp(l)/sump
               omm(l) = alfm(l)/summ
           enddo

           vminus(m) = omp(0) *(-vp(m,i-1)+3*vp(m,i  )) + omp(1) *( vp(m,i  )+ vp(m,i+1))
           vplus(m)  = omm(0) *(-vm(m,i+2)+3*vm(m,i+1)) + omm(1) *( vm(m,i  )+ vm(m,i+1))

       enddo

       do m=1,nvar
           vminus(m) = 0.5_R8P*vminus(m)
           vplus(m)  = 0.5_R8P*vplus(m)
       enddo

     elseif (iweno==3) then ! WENO-5
!
      i = iweno ! index of intermediate node to perform reconstruction
!
      dwe( 0) = 1._R8P/10._R8P
      dwe( 1) = 6._R8P/10._R8P
      dwe( 2) = 3._R8P/10._R8P
!     JS
      d0 = 13._R8P/12._R8P
      d1 = 1._R8P/4._R8P
!     Weights for polynomial reconstructions
      c0 = 1._R8P/3._R8P
      c1 = 5._R8P/6._R8P
      c2 =-1._R8P/6._R8P
      c3 =-7._R8P/6._R8P
      c4 =11._R8P/6._R8P
!
      do m=1,nvar
!
       betap(2) = d0*(     vp(m,i)-2._R8P*vp(m,i+1)+vp(m,i+2))**2+d1*(3._R8P*vp(m,i)-4._R8P*vp(m,i+1)+vp(m,i+2))**2
       betap(1) = d0*(     vp(m,i-1)-2._R8P*vp(m,i)+vp(m,i+1))**2+d1*(     vp(m,i-1)-vp(m,i+1) )**2
       betap(0) = d0*(     vp(m,i)-2._R8P*vp(m,i-1)+vp(m,i-2))**2+d1*(3._R8P*vp(m,i)-4._R8P*vp(m,i-1)+vp(m,i-2))**2
!
       betam(2) = d0*(     vm(m,i+1)-2._R8P*vm(m,i)+vm(m,i-1))**2+d1*(3._R8P*vm(m,i+1)-4._R8P*vm(m,i)+vm(m,i-1))**2
       betam(1) = d0*(     vm(m,i+2)-2._R8P*vm(m,i+1)+vm(m,i))**2+d1*(     vm(m,i+2)-vm(m,i) )**2
       betam(0) = d0*(     vm(m,i+1)-2._R8P*vm(m,i+2)+vm(m,i+3))**2+d1*(3._R8P*vm(m,i+1)-4._R8P*vm(m,i+2)+vm(m,i+3))**2
!
       sump = 0._R8P
       summ = 0._R8P
       do l=0,2
        alfp(l) = dwe(  l)/(0.000001_R8P+betap(l))**2
        alfm(l) = dwe(  l)/(0.000001_R8P+betam(l))**2
        sump = sump + alfp(l)
        summ = summ + alfm(l)
       enddo
       do l=0,2
        omp(l) = alfp(l)/sump
        omm(l) = alfm(l)/summ
       enddo
!
       vminus(m)   = omp(2)*(c0*vp(m,i  )+c1*vp(m,i+1)+c2*vp(m,i+2)) + &
         & omp(1)*(c2*vp(m,i-1)+c1*vp(m,i  )+c0*vp(m,i+1)) + omp(0)*(c0*vp(m,i-2)+c3*vp(m,i-1)+c4*vp(m,i  ))
       vplus(m)   = omm(2)*(c0*vm(m,i+1)+c1*vm(m,i  )+c2*vm(m,i-1)) +  &
         & omm(1)*(c2*vm(m,i+2)+c1*vm(m,i+1)+c0*vm(m,i  )) + omm(0)*(c0*vm(m,i+3)+c3*vm(m,i+2)+c4*vm(m,i+1))
!
      enddo ! end of m-loop
!
   elseif (iweno==4) then ! WENO-7
!
      i = iweno ! index of intermediate node to perform reconstruction
!
      dwe( 0) = 1._R8P/35._R8P
      dwe( 1) = 12._R8P/35._R8P
      dwe( 2) = 18._R8P/35._R8P
      dwe( 3) = 4._R8P/35._R8P
!
!     JS weights
      d1 = 1._R8P/36._R8P
      d2 = 13._R8P/12._R8P
      d3 = 781._R8P/720._R8P
!
      do m=1,nvar
!
       betap(3)= d1*(-11*vp(m,  i)+18*vp(m,i+1)- 9*vp(m,i+2)+ 2*vp(m,i+3))**2+&
       &  d2*(  2*vp(m,  i)- 5*vp(m,i+1)+ 4*vp(m,i+2)-   vp(m,i+3))**2+ &
       & d3*(   -vp(m,  i)+ 3*vp(m,i+1)- 3*vp(m,i+2)+   vp(m,i+3))**2
       betap(2)= d1*(- 2*vp(m,i-1)- 3*vp(m,i  )+ 6*vp(m,i+1)-   vp(m,i+2))**2+&
       &  d2*(    vp(m,i-1)- 2*vp(m,i  )+   vp(m,i+1)             )**2+&
       &  d3*(   -vp(m,i-1)+ 3*vp(m,i  )- 3*vp(m,i+1)+   vp(m,i+2))**2
       betap(1)= d1*(    vp(m,i-2)- 6*vp(m,i-1)+ 3*vp(m,i  )+ 2*vp(m,i+1))**2+&
       &  d2*( vp(m,i-1)- 2*vp(m,i  )+   vp(m,i+1))**2+ &
       &  d3*(   -vp(m,i-2)+ 3*vp(m,i-1)- 3*vp(m,i  )+   vp(m,i+1))**2
       betap(0)= d1*(- 2*vp(m,i-3)+ 9*vp(m,i-2)-18*vp(m,i-1)+11*vp(m,i  ))**2+&
       &  d2*(-   vp(m,i-3)+ 4*vp(m,i-2)- 5*vp(m,i-1)+ 2*vp(m,i  ))**2+&
       &  d3*(   -vp(m,i-3)+ 3*vp(m,i-2)- 3*vp(m,i-1)+   vp(m,i  ))**2
!
       betam(3)= d1*(-11*vm(m,i+1)+18*vm(m,i  )- 9*vm(m,i-1)+ 2*vm(m,i-2))**2+&
       &  d2*(  2*vm(m,i+1)- 5*vm(m,i  )+ 4*vm(m,i-1)-   vm(m,i-2))**2+&
       &  d3*(   -vm(m,i+1)+ 3*vm(m,i  )- 3*vm(m,i-1)+   vm(m,i-2))**2
       betam(2)= d1*(- 2*vm(m,i+2)- 3*vm(m,i+1)+ 6*vm(m,i  )-   vm(m,i-1))**2+&
       &  d2*(    vm(m,i+2)- 2*vm(m,i+1)+   vm(m,i  )             )**2+&
       &  d3*(   -vm(m,i+2)+ 3*vm(m,i+1)- 3*vm(m,i  )+   vm(m,i-1))**2
       betam(1)= d1*(    vm(m,i+3)- 6*vm(m,i+2)+ 3*vm(m,i+1)+ 2*vm(m,i  ))**2+&
       &  d2*(                 vm(m,i+2)- 2*vm(m,i+1)+   vm(m,i  ))**2+&
       &  d3*(   -vm(m,i+3)+ 3*vm(m,i+2)- 3*vm(m,i+1)+   vm(m,i  ))**2
       betam(0)= d1*(- 2*vm(m,i+4)+ 9*vm(m,i+3)-18*vm(m,i+2)+11*vm(m,i+1))**2+&
       &  d2*(-   vm(m,i+4)+ 4*vm(m,i+3)- 5*vm(m,i+2)+ 2*vm(m,i+1))**2+&
       &  d3*(   -vm(m,i+4)+ 3*vm(m,i+3)- 3*vm(m,i+2)+   vm(m,i+1))**2
!
       sump = 0._R8P
       summ = 0._R8P
       do l=0,3
        alfp(l) = dwe(  l)/(0.000001_R8P+betap(l))**2
        alfm(l) = dwe(  l)/(0.000001_R8P+betam(l))**2
        sump = sump + alfp(l)
        summ = summ + alfm(l)
       enddo
       do l=0,3
        omp(l) = alfp(l)/sump
        omm(l) = alfm(l)/summ
       enddo
!
       vminus(m)   = omp(3)*( 6*vp(m,i  )+26*vp(m,i+1)-10*vp(m,i+2)+ 2*vp(m,i+3))+&
        omp(2)*(-2*vp(m,i-1)+14*vp(m,i  )+14*vp(m,i+1)- 2*vp(m,i+2))+&
        omp(1)*( 2*vp(m,i-2)-10*vp(m,i-1)+26*vp(m,i  )+ 6*vp(m,i+1))+&
        omp(0)*(-6*vp(m,i-3)+26*vp(m,i-2)-46*vp(m,i-1)+50*vp(m,i  ))
       vplus(m)   =  omm(3)*( 6*vm(m,i+1)+26*vm(m,i  )-10*vm(m,i-1)+ 2*vm(m,i-2))+&
        omm(2)*(-2*vm(m,i+2)+14*vm(m,i+1)+14*vm(m,i  )- 2*vm(m,i-1))+&
        omm(1)*( 2*vm(m,i+3)-10*vm(m,i+2)+26*vm(m,i+1)+ 6*vm(m,i  ))+&
        omm(0)*(-6*vm(m,i+4)+26*vm(m,i+3)-46*vm(m,i+2)+50*vm(m,i+1))
!
      enddo ! end of m-loop
!
      vminus = vminus/24._R8P
      vplus  = vplus /24._R8P
!
   else
      write(*,*) 'Error! WENO scheme not implemented'
      stop
   endif

   endsubroutine weno_reconstruction

   subroutine compute_rk_stage_gpu_cuf(ni, nj, nk, ngc, nv, blocks_number, alph_gpu, dt, s, q_gpu, q_s_gpu)
   !< Initialize RK stage with q_gpu.
   integer(I4P), intent(in)            :: ni                                     !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                     !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                     !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                    !< Ghost cells number.
   integer(I4P), intent(in)            :: nv                                     !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number                          !< Number of blocks.
   real(R8P),    intent(in),    device :: alph_gpu(:,:)                          !< RK alpha coefficients.
   real(R8P),    intent(in)            :: dt                                     !< Time step.
   integer(I4P), intent(in)            :: s                                      !< Stage to initialize.
   real(R8P),    intent(in),    device ::   q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative field.
   real(R8P),    intent(inout), device :: q_s_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:) !< RK stage.
   integer(I4P)                        :: i, j, k, b, v, ss                      !< Counter.
   integer(I4P)                        :: iercuda                                !< Error trapping flag for CUDAFortran.

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

   subroutine minimal_immersed_bc(ni, nj, nk, ngc, nv, blocks_number, gamma_fluid, &
           q_gpu, phi_gpu, x_cell_gpu, y_cell_gpu, z_cell_gpu)
   !< Initialize RK stage with q_gpu.
   integer(I4P), intent(in)            :: ni                                     !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                     !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                     !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                    !< Ghost cells number.
   integer(I4P), intent(in)            :: nv                                     !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number                          !< Number of blocks.
   real(R8P), intent(in)               :: gamma_fluid
   real(R8P),    intent(inout), device :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< RK stage.
   real(R8P),    intent(inout), device :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< RK stage.
   integer(I4P)                        :: i, j, k, b, v                          !< Counter.
   integer(I4P)                        :: iercuda                                !< Error trapping flag for CUDAFortran.
   real(R8P),    intent(in), device    :: x_cell_gpu(1:,1-ngc:)                  !< Conservative variables.
   real(R8P),    intent(in), device    :: y_cell_gpu(1:,1-ngc:)                  !< Conservative variables.
   real(R8P),    intent(in), device    :: z_cell_gpu(1:,1-ngc:)                  !< Conservative variables.
   real(R8P)                           :: x,y,z,tem,rho,pres

   !$cuf kernel do(4) <<<*,*>>>
   do k=1-ngc, nk+ngc
      do j=1-ngc, nj+ngc
         do i=1-ngc, ni+ngc
            do b=1, blocks_number
               x = x_cell_gpu(b,i)
               y = y_cell_gpu(b,j)
               z = z_cell_gpu(b,k)
               !if(x<2.2_R8P .and. x>1.8_R8P .and. y<2.2_R8P .and. y>1.8_R8P .and. z<2.2_R8P .and. z>1.8_R8P) then
               !if( (x-4._R8P)**2+(y-2._R8P)**2+(z-2._R8P)**2 < 0.1_R8P) then
               if( phi_gpu(b,i,j,k,1) >  0._R8P ) then
                  pres = 10000.
                  tem = 300.
                  rho = pres/tem*gamma_fluid/(gamma_fluid-1._R8P)
                  q_gpu(b,i,j,k,1) = rho
                  q_gpu(b,i,j,k,2) = 0._R8P
                  q_gpu(b,i,j,k,3) = 0._R8P
                  q_gpu(b,i,j,k,4) = 0._R8P
                  q_gpu(b,i,j,k,5) = q_gpu(b,i,j,k,1)*1._R8P/gamma_fluid*tem
               endif
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine minimal_immersed_bc

   subroutine compute_umax_cuf(b, ni, nj, nk, ngc, ns, dx, dy, dz, q_aux_gpu, umax)
   !< Compute maximum speed by means of CUF threads.
   integer(I4P), intent(in)         :: b                                     !< Block index.
   integer(I4P), intent(in)         :: ni                                    !< Grid cells number in I direction.
   integer(I4P), intent(in)         :: nj                                    !< Grid cells number in J direction.
   integer(I4P), intent(in)         :: nk                                    !< Grid cells number in K direction.
   integer(I4P), intent(in)         :: ngc                                   !< Ghost cells number.
   integer(I4P), intent(in)         :: ns                                    !< Number of species.
   real(R8P),    intent(in)         :: dx                                    !< X space step.
   real(R8P),    intent(in)         :: dy                                    !< Y space step.
   real(R8P),    intent(in)         :: dz                                    !< Z space step.
   real(R8P),    intent(in), device :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary varibales.
   real(R8P),    intent(out)        :: umax                                  !< Maximum speed.
   real(R8P)                        :: ss                                    !< Speed of sound.
   integer(I4P)                     :: i, j, k                               !< Counter.
   integer(I4P)                     :: iercuda                               !< Error trapping flag for CUDAFortran.

   umax = 0._R8P
   !$cuf kernel do(3) <<<*,*>>>
   do k=1, nk
      do j=1, nj
         do i=1, ni
            ss = q_aux_gpu(b,i,j,k,9)
            umax = max(umax, abs(q_aux_gpu(b,i,j,k,2)) + ss, &
                             abs(q_aux_gpu(b,i,j,k,3)) + ss, &
                             abs(q_aux_gpu(b,i,j,k,4)) + ss)

         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_umax_cuf

endmodule adam_equation_flame_gpu_object
