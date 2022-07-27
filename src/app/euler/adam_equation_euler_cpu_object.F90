!< ADAM, Euler equations system class definition and solver, base backend.
module adam_equation_euler_cpu_object
!< ADAM, Euler equations system class definition and solver, base backend.

use adam_adam_object,       only : adam_object
use adam_base_cpu_object,   only : base_cpu_object
use adam_eos_ic_cpu_object, only : eos_ic_cpu_object
use adam_field_object,      only : field_object
use adam_grid_object,       only : grid_object
use adam_ib_cpu_object,     only : ib_cpu_object
use adam_mpih_object,       only : mpih_object
use adam_memory_cpu_lib
use adam_parameters
use finer
use penf
use mpi

implicit none
private
public :: equation_euler_cpu_object

! Initial conditions parameters.
integer(I4P), parameter :: IC_UNIFORM         = 1_I4P           !< Uniform IC.
integer(I4P), parameter :: IC_LEFTRIGHT       = 2_I4P           !< Discontinous left/right IC.
integer(I4P), parameter :: IC_VARS_NUMBER(2)  = [6_I4P, 11_I4P] !< Number of IC variables for admitted cases.
integer(I4P), parameter :: IC_VARS_NUMBER_MAX = 11              !< Maximum number of IC variables, maxval(IC_VARS_NUMBER).

! Boundary conditions parameters.
integer(I4P), parameter :: BC_EXTRAPOLATION   = 1_I4P          !< Extrapolation.
integer(I4P), parameter :: BC_INFLOW          = 2_I4P          !< Supersonic inflow.
integer(I4P), parameter :: BC_VARS_NUMBER(2)  = [0_I4P, 5_I4P] !< Number of BC variables for admitted cases.
integer(I4P), parameter :: BC_VARS_NUMBER_MAX = 5_I4P          !< Maximum number of BC variables, maxval(BC_VARS_NUMBER).

! Immersed boundary conditions parameters.
integer(I4P), parameter :: BCS_VISCOUS = 1_I4P !< Visous wall.
integer(I4P), parameter :: BCS_EULER   = 2_I4P !< Inviscid wall.

! Named indexes of q variables.
integer(I4P), parameter :: IRHO  = 1_I4P
integer(I4P), parameter :: IRHOU = 2_I4P
integer(I4P), parameter :: IRHOV = 3_I4P
integer(I4P), parameter :: IRHOW = 4_I4P
integer(I4P), parameter :: IRHOE = 5_I4P

! AMR marker modes parameters.
integer(I4P), parameter :: AMR_GEO  = 1_I4P
integer(I4P), parameter :: AMR_GRAD = 2_I4P

type :: amr_marker_object
   !< AMR marker object.
   integer(I4P) :: mode         !< Marker mode.
   integer(I4P) :: solid        !< Solid number.
   real(R8P)    :: delta_fine   !< Fine cell space step.
   real(R8P)    :: delta_coarse !< Coarse cell space step.
   integer(I4P) :: ivar         !< ivar.
   real(R8P)    :: tol          !< Tolerance.
endtype amr_marker_object

type :: slice_object
   !< Slice object.
   character(99)          :: slice_itype           !< Slice interpolation type.
   integer(I4P)           :: slice_save            !< Iteration interval between subsequent data-slice saves.
   integer(I4P)           :: slice_nijk(3)         !< Slice number of points.
   real(R8P)              :: slice_emin(3)         !< Slice minimum extents.
   real(R8P)              :: slice_emax(3)         !< Slice maximum extents.
   real(R8P), allocatable :: slice_points(:,:,:,:) !< Slice points coordinates [3,ni,nj,nk].
endtype slice_object

type :: equation_euler_cpu_object
   !< Euler equations system class definition and solver, base backend.
   !<
   !< The conservative variables are arranged as follows:
   !<```
   !< q(1): rho
   !< q(2): rho * u
   !< q(3): rho * v
   !< q(4): rho * w
   !< q(5): rho * E
   !<```
   !< The auxiliary variables are arranged as follows:
   !<```
   !< q_aux(1): rho
   !< q_aux(2): u
   !< q_aux(3): v
   !< q_aux(4): w
   !< q_aux(5): E
   !< q_aux(6): p
   !< q_aux(7): T
   !< q_aux(8): H
   !< q_aux(9): a
   !<```
   ! Auxiliary objects.
   type(mpih_object)           :: mpih          !< MPI handler.
   type(file_ini)              :: file_input    !< Nasto input file handler.
   type(adam_object)           :: adam          !< ADAM.
   type(field_object), pointer :: field=>null() !< The field.
   type(grid_object),  pointer :: grid=>null()  !< The grid.
   type(base_cpu_object)       :: base_cpu      !< The base CPU handler.
   ! Pointers to ADAM data for easy handling.
   integer(I4P), pointer :: ngc          =>null() !< Number of ghost cells.
   integer(I4P), pointer :: ni           =>null() !< Number of cells in i direction.
   integer(I4P), pointer :: nj           =>null() !< Number of cells in j direction.
   integer(I4P), pointer :: nk           =>null() !< Number of cells in k direction.
   integer(I4P), pointer :: nb           =>null() !< Total blocks number for MPI.
   integer(I4P), pointer :: blocks_number=>null() !< Actual blocks number.
   ! Equation dimensions.
   integer(I4P)  :: nv=5_I4P     !< Number of variables.
   integer(I4P)  :: nv_aux=9_I4P !< Number of auxiliary variables.
   ! Immersed Boundary.
   type(ib_cpu_object) :: ib !< IBs handler.
   ! AMR.
   integer(I4P)                         :: amr_iters=5_I4P          !< AMR updates iterations number.
   integer(I4P)                         :: amr_frequency=100_I4P    !< AMR update time step frequency.
   integer(I4P)                         :: amr_markers_number=1_I4P !< AMR number of markers.
   type(amr_marker_object), allocatable :: amr_markers(:)           !< AMR array of marker objects.
   ! Initial conditions.
   integer(I4P) :: ic_type                     !< Initial condition type.
   real(R8P)    :: ic_vars(IC_VARS_NUMBER_MAX) !< Variables array for initial conditions.
   ! Boundary conditions.
   integer(I4P), pointer     :: bc_type(:)=>null()             !< Boundary condition type.
   real(R8P)                 :: bc_vars(BC_VARS_NUMBER_MAX, 6) !< Variables' array for boundary conditions.
   integer(I4P), allocatable :: bcs_type(:)                    !< Immersed boundary condition type.
   real(R8P),    allocatable :: bcs_vars(:, :)                 !< Variables' array for immersed boundary conditions.
   ! Slices.
   integer(I4P)                    :: slices_number=0 !< Number of slices to be saved.
   type(slice_object), allocatable :: slice(:)        !< Slices data.
   ! Solver data.
   character(99)  :: space_operator='weno-js-3' !< Space operator algorithm.
   character(99)  :: time_operator='rk-ssp-3'   !< Time operator algorithm.
   logical        :: restart=.false.            !< Restart flag.
   character(999) :: restart_basename           !< Restart file basename.
   integer(I4P)   :: restart_save=100_I4P       !< Restart saves frequency.
   integer(I4P)   :: it_max=-1_I4P              !< Maximum number of integration time steps.
   integer(I4P)   :: it_save=-1_I4P             !< Iterations saves frequncy.
   integer(I4P)   :: it=0_I4P                   !< Time steps counter.
   real(R8P)      :: time_max=1._R8P            !< Maximum integration time.
   real(R8P)      :: time_save=0.1_R8P          !< Time saves frequency.
   real(R8P)      :: time=0._R8P                !< Time.
   real(R8P)      :: CFL=0.3_R8P                !< CFL time limit.
   real(R8P)      :: dt=0.0001_R8P              !< Maximum time step accordingly to CFL criterion.
   character(999) :: output_basename='euler'    !< Output file basename.
   logical        :: save_memory_status=.false. !< Flag to activate memory status saving.
   ! Physics.
   type(eos_ic_cpu_object) :: eos !< Equations of state (cp, cv...).
   ! Large arrays.
   real(R8P), pointer     ::     q(:,:,:,:,:)=>null() !< Conservative cell centerd variables, allocated in field object.
   real(R8P), allocatable :: q_aux(:,:,:,:,:)         !< Auxiliary cell centered variables.
   real(R8P), allocatable ::  q_ib(:,:,:,:,:)         !< Field cell with boundary set on immersed bodies.
   contains
      ! public methods
      procedure, pass(self) :: amr_update              !< Do AMR update.
      procedure, pass(self) :: compute_aux             !< Compute auxiliary variables.
      procedure, pass(self) :: compute_dt              !< Compute time step.
      procedure, pass(self) :: initialize              !< Initialize the equation.
      procedure, pass(self) :: load_restart_files      !< Load restart files.
      procedure, pass(self) :: mark_by_grad_var        !< Mark blocks to be refined/derefined by a `grad(var)` value.
      procedure, pass(self) :: mark_by_geo             !< Mark blocks to be refined/derefined by a `grad(var)` value.
      procedure, pass(self) :: integrate               !< Runge Kutta integration of equation.
      procedure, pass(self) :: print_progress          !< Print simulation progress.
      procedure, pass(self) :: refine_uniform          !< Refine all blocks uniformly.
      ! procedure, pass(self) :: runge_kutta_initialize  !< Initialize Runge-Kutta data.
      procedure, pass(self) :: save_simulation_data    !< Save all simulation data.
      procedure, pass(self) :: save_restart_files      !< Save restart files.
      procedure, pass(self) :: save_hdf5               !< Save simulation data in HDF5 format.
      procedure, pass(self) :: save_slices             !< Save simulation data slices.
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions of equation.
      procedure, pass(self) :: set_initial_conditions  !< Set initial conditions of equation.
      procedure, pass(self) :: solve                   !< Solve Euler system.
      procedure, pass(self) :: update_ghost            !< Update ghost cells.
endtype equation_euler_cpu_object

contains
   ! public methods
   subroutine amr_update(self)
   !< Do AMR update.
   class(equation_euler_cpu_object), intent(inout) :: self                !< The equation.
   logical                                         :: is_grid_changed     !< Flag to check grid changes for each marker.
   logical                                         :: is_grid_changed_all !< Flag to check grid changes for each iter.
   type(amr_marker_object)                         :: amr_marker          !< Current amr marker.
   integer(I4P)                                    :: i, i_marker         !< Counter.

   amr: do i=1, self%amr_iters
      is_grid_changed_all = .false.
      do i_marker=1, self%amr_markers_number
         amr_marker = self%amr_markers(i_marker)
         call self%update_ghost(q=self%q)
         if(amr_marker%mode == AMR_GEO) then
            call self%mark_by_geo(delta_fine=amr_marker%delta_fine, delta_coarse=amr_marker%delta_coarse)
         elseif(amr_marker%mode == AMR_GRAD) then
            call self%mark_by_grad_var(grad_tol=amr_marker%tol, delta_fine=amr_marker%delta_fine, &
                                       delta_coarse=amr_marker%delta_coarse, q=self%q, ivar=amr_marker%ivar)
         endif
         call self%adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false., is_grid_changed=is_grid_changed)
         call self%ib%update_phi
         is_grid_changed_all = is_grid_changed_all.or.is_grid_changed
      enddo
      if (.not.is_grid_changed_all) then
          print '(A)', self%mpih%myrankstr//'AMR Grid stabilized after : '//trim(str(i))//' AMR iterations'
          exit amr
       elseif (i==self%amr_iters) then
          print '(A)', self%mpih%myrankstr//'AMR Grid is NOT stabilized after : '//trim(str(i))//' AMR iterations'
      endif
   enddo amr
   endsubroutine amr_update

   subroutine compute_aux(self, q, q_aux)
   !< Compute auxiliary variables.
   class(equation_euler_cpu_object), intent(in)  :: self             !< The equation.
   real(R8P),                        intent(in)  :: q(1:,         &
                                                      1-self%ngc:,&
                                                      1-self%ngc:,&
                                                      1-self%ngc:,&
                                                      1:)            !< Conservative variables.
   real(R8P),                        intent(out) :: q_aux(1:,         &
                                                          1-self%ngc:,&
                                                          1-self%ngc:,&
                                                          1-self%ngc:,&
                                                          1:)        !< Auxiliary variables.
   real(R8P)                                     :: velocity_sq_norm !< Velocity square root norm.
   integer(I4P)                                  :: b, i, j, k       !< Counter.

   associate(blocks_number=>self%blocks_number, q=>self%q, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, gm1=>self%eos%gm1)
   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               q_aux(1,i,j,k,b) = q(IRHO, i,j,k,b)                                                        ! rho

               q_aux(2,i,j,k,b) = q(IRHOU,i,j,k,b) / q(IRHO, i,j,k,b)                                     ! u
               q_aux(3,i,j,k,b) = q(IRHOV,i,j,k,b) / q(IRHO, i,j,k,b)                                     ! v
               q_aux(4,i,j,k,b) = q(IRHOW,i,j,k,b) / q(IRHO, i,j,k,b)                                     ! w
               velocity_sq_norm = sqrt(q_aux(2,i,j,k,b)**2 + &
                                       q_aux(3,i,j,k,b)**2 + &
                                       q_aux(4,i,j,k,b)**2)                                               ! velocity sq norm

               q_aux(5,i,j,k,b) = q(IRHOE,i,j,k,b) / q(IRHO, i,j,k,b)                                     ! E
               q_aux(6,i,j,k,b) = gm1 * (q(IRHOE,i,j,k,b) - 0.5_R8P * q(IRHO,i,j,k,b) * velocity_sq_norm) ! p
               q_aux(7,i,j,k,b) = self%eos%temperature(energy=q_aux(5,i,j,k,b))                           ! T
               q_aux(8,i,j,k,b) = self%eos%total_entalpy(density= q_aux(1,i,j,k,b), &
                                                         pressure=q_aux(6,i,j,k,b), &
                                                         velocity_sq_norm=velocity_sq_norm)               ! H
               q_aux(9,i,j,k,b) = self%eos%speed_of_sound(density= q_aux(1,i,j,k,b), &
                                                          pressure=q_aux(6,i,j,k,b))                      ! a
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine compute_aux

   subroutine compute_dt(self)
   !< Compute maximum time step accordingly to CFL stabilty criterion.
   class(equation_euler_cpu_object), intent(inout) :: self    !< The equation.
   real(R8P)                                       :: umax    !< Maximum speed of waves propagation.
   integer(I4P)                                    :: b,i,j,k !< Counter.

   associate(blocks_number=>self%blocks_number, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             ngc=>self%ngc, dxyz=>self%field%dxyz, dt=>self%dt, CFL=>self%CFL)
      call self%compute_aux(q=self%q, q_aux=self%q_aux)
      dt = huge(1._R8P)
      do b=1, blocks_number
         umax = 0._R8P
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  umax = max(umax, (abs(self%q_aux(2,i,j,k,b)) + self%q_aux(9,i,j,k,b))/dxyz(1,b), &
                                   (abs(self%q_aux(3,i,j,k,b)) + self%q_aux(9,i,j,k,b))/dxyz(2,b), &
                                   (abs(self%q_aux(4,i,j,k,b)) + self%q_aux(9,i,j,k,b))/dxyz(3,b))
               enddo
            enddo
         enddo
         dt = min(dt, CFL / umax)
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, self%mpih%error)
   endassociate
   if ((self%it_max <= 0).and.(self%time + self%dt > self%time_max)) self%dt = self%time_max - self%time
   endsubroutine compute_dt

   subroutine initialize(self, filename)
   !< Initialize the equation.
   class(equation_euler_cpu_object), intent(inout) :: self         !< The equation.
   character(*),                     intent(in)    :: filename     !< Input file name.
   integer(I8P)                                    :: nodes_number !< Allocated nodes on tree.
   integer(I4P)                                    :: nb           !< Number of allocated blocks.

   call self%mpih%initialize(do_mpi_init=.true.)

   print '(A)', self%mpih%myrankstr//'equation_euler_cpu_object%initialize start'

   call self%base_cpu%initialize_cpu

   call self%file_input%initialize(filename=trim(filename))
   call self%file_input%load

   call self%adam%grid%initialize(file_parameters=self%file_input, verbose=.true.)

   call self%adam%compute_blocks_number(memory_avail=self%base_cpu%memory_avail,  &
                                        fields_number=80,                         &
                                        nb=nb,                                    &
                                        nodes_number=nodes_number)

   call self%adam%initialize(file_parameters=self%file_input, &
                             do_tree_init=.true.,             &
                             do_field_init=.true.,            &
                             nv=self%nv, nb=nb, nodes_number=nodes_number)
   call forward_main_adam_data(grid=self%adam%grid, field=self%adam%field)

   call self%adam%refine_uniform(refinement_levels=self%adam%tree%iu_ref_levels, do_blocks_reorder=.false.)

   call self%adam%prune(ijkl_prune=self%adam%tree%ijkl_prune, do_blocks_reorder=.false.)

   call self%base_cpu%initialize(field=self%adam%field)

   call load_equation_from_ini_file

   call load_schemes_from_ini_file

   call load_physics_from_ini_file

   call load_amr_from_ini_file

   call load_timing_from_ini_file

   call load_ic_from_ini_file

   call load_bc_from_ini_file

   call load_slices_from_ini_file

   call self%ib%initialize(grid=self%adam%grid, field=self%adam%field, file_parameters=self%file_input)

   ! call self%runge_kutta_initialize

   ! allocate large arrays
   associate(nv=>self%nv, nv_aux=>self%nv_aux, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, nb=>self%nb)
   allocate(self%q_aux(1:nv_aux, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
   allocate(self%q_ib( 1:nv,     1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
   endassociate

   print '(A)', self%mpih%myrankstr//'equation_euler_cpu_object%initialize finish'
   contains
      subroutine forward_main_adam_data(grid, field)
      !< Forward main ADAM data to equation for easy handling.
      type(grid_object),  intent(in), target :: grid  !< The grid.
      type(field_object), intent(in), target :: field !< The field.

      self%grid          => grid
      self%field         => field
      self%q             => field%q
      self%blocks_number => field%blocks_number
      self%ni            => field%grid%ni
      self%nj            => field%grid%nj
      self%nk            => field%grid%nk
      self%ngc           => field%grid%ngc
      self%nb            => field%nb
      self%bc_type       => grid%bc_type
      endsubroutine forward_main_adam_data

      subroutine load_equation_from_ini_file
      !< Parse equation setting from input file.
      logical :: buf_LOG !< Logical buffer.

      call self%file_input%get(section_name='equation', option_name='save_memory_status', val=buf_LOG)
      self%save_memory_status = buf_LOG
      print '(A)', self%mpih%myrankstr//'save memory status: '//trim(str(self%save_memory_status))
      endsubroutine load_equation_from_ini_file

      subroutine load_schemes_from_ini_file
      !< Parse schemes setting from input file.
      character(99) :: buf_CHAR !< String buffer.

      call self%file_input%get(section_name='schemes',option_name='space_operator',val=buf_CHAR);self%space_operator=trim(buf_CHAR)
      call self%file_input%get(section_name='schemes',option_name='time_operator', val=buf_CHAR); self%time_operator=trim(buf_CHAR)
      endsubroutine load_schemes_from_ini_file

      subroutine load_physics_from_ini_file
      !< Parse physics setting from input file.
      real(R8P) :: cp, cv !< Constant specific heats.

      call self%file_input%get(section_name='physics', option_name='cp', val=cp)
      call self%file_input%get(section_name='physics', option_name='cv', val=cv)
      call self%eos%initialize(cp=cp, cv=cv)
      endsubroutine load_physics_from_ini_file

      subroutine load_amr_from_ini_file
      !< Parse AMR setting from input file.
      integer(I4P)   :: buf_I4   !< I4 buffer.
      real(R8P)      :: buf_R8   !< R8 buffer.
      character(999) :: sname    !< Section name.
      integer(I4P)   :: i_marker !< Counter.
      integer(I4P)   :: mode     !< AMR mode.

      call self%file_input%get(section_name='amr', option_name='frequency', val=buf_I4) ; self%amr_frequency = buf_I4
      call self%file_input%get(section_name='amr', option_name='iters',     val=buf_I4) ; self%amr_iters     = buf_I4
      call self%file_input%get(section_name='amr', option_name='n_markers', val=buf_I4) ; self%amr_markers_number = buf_I4
      allocate(self%amr_markers(self%amr_markers_number))
      do i_marker=1,self%amr_markers_number
         sname = 'amr_marker_'//trim(str(i_marker,.true.))
         call self%file_input%get(section_name=sname, option_name='mode', val=mode)
         self%amr_markers(i_marker)%mode = mode
         call self%file_input%get(section_name=sname, option_name='delta_fine', val=buf_R8)
         self%amr_markers(i_marker)%delta_fine = buf_R8
         call self%file_input%get(section_name=sname, option_name='delta_coarse', val=buf_R8)
         self%amr_markers(i_marker)%delta_coarse = buf_R8
         if     (mode == AMR_GEO) then
            call self%file_input%get(section_name=sname, option_name='solid', val=buf_I4)
            self%amr_markers(i_marker)%solid = buf_I4
         elseif (mode == AMR_GRAD) then
            call self%file_input%get(section_name=sname, option_name='var', val=buf_I4)
            self%amr_markers(i_marker)%ivar = buf_I4
            call self%file_input%get(section_name=sname, option_name='tol', val=buf_R8)
            self%amr_markers(i_marker)%tol = buf_R8
         endif
      enddo
      endsubroutine load_amr_from_ini_file

      subroutine load_timing_from_ini_file
      !< Parse timing setting from input file.
      logical        :: buf_BOOL !< Logical buffer.
      character(999) :: buf_CHAR !< String buffer.
      integer(I4P)   :: buf_I4   !< I4 buffer.
      real(R8P)      :: buf_R8   !< R8 buffer.

      call self%file_input%get(section_name="time", option_name="restart",         val=buf_BOOL) ; self%restart          = buf_BOOL
      call self%file_input%get(section_name="time", option_name="restart_basename",val=buf_CHAR) ; self%restart_basename = buf_CHAR
      call self%file_input%get(section_name="time", option_name="restart_save",    val=buf_I4)   ; self%restart_save     = buf_I4
      call self%file_input%get(section_name="time", option_name="it_max",          val=buf_I4)   ; self%it_max           = buf_I4
      call self%file_input%get(section_name="time", option_name="it_save",         val=buf_I4)   ; self%it_save          = buf_I4
      call self%file_input%get(section_name="time", option_name="time_max",        val=buf_R8)   ; self%time_max         = buf_R8
      call self%file_input%get(section_name="time", option_name="time_save",       val=buf_R8)   ; self%time_save        = buf_R8
      call self%file_input%get(section_name="time", option_name="output_basename", val=buf_CHAR) ; self%output_basename  = buf_CHAR
      call self%file_input%get(section_name='time', option_name='CFL',             val=buf_R8)   ; self%CFL              = buf_R8
      endsubroutine load_timing_from_ini_file

      subroutine load_ic_from_ini_file
      !< Parse initial conditions setting from input file.
      integer(I4P)   :: buf_I4 !< I4 buffer.
      real(R8P)      :: buf_R8 !< R8 buffer.
      character(999) :: oname  !< Option name.
      integer(I4P)   :: i_var  !< Counter.
      integer(I4P)   :: n_vars !< Number of vars.

      call self%file_input%get(section_name="initial_conditions", option_name='ic_type', val=buf_I4) ; self%ic_type = buf_I4
      n_vars = IC_VARS_NUMBER(self%ic_type)
      do i_var=1,n_vars
         oname = "var"//trim(str(i_var,.true.))
         call self%file_input%get(section_name="initial_conditions", option_name=oname, val=buf_R8) ; self%ic_vars(i_var) = buf_R8
      enddo
      endsubroutine load_ic_from_ini_file

      subroutine load_bc_from_ini_file
      !< Parse boundary conditions setting from input file.
      real(R8P)      :: buf_R8       !< R8 buffer.
      character(999) :: sname        !< Section name.
      character(999) :: snames_bc(6) !< Section names bc.
      integer(I4P)   :: bc_type_item !< Boundary condition type element.
      integer(I4P)   :: i_var, i_bc  !< Counter.
      integer(I4P)   :: n_vars       !< Number of vars.

      snames_bc(:) = ["bc_x_min", "bc_x_max", "bc_y_min", "bc_y_max", "bc_z_min", "bc_z_max"]
      do i_bc=1,6
         sname = snames_bc(i_bc)
         call self%file_input%get(section_name=sname, option_name='type', val=bc_type_item)
         n_vars = BC_VARS_NUMBER(bc_type_item)
         do i_var=1,n_vars
             call self%file_input%get(section_name=sname, option_name="var"//trim(str(i_var,.true.)), val=buf_R8)
             self%bc_vars(i_var, i_bc) = buf_R8
         enddo
      enddo
      endsubroutine load_bc_from_ini_file

      subroutine load_slices_from_ini_file
      !< Parse slices setting from input file.
      character(999) :: buf_CHAR   !< String buffer.
      integer(I4P)   :: buf_I4     !< I4 buffer.
      real(R8P)      :: buf_R8     !< R8 buffer.
      character(999) :: sname      !< Section name.
      real(R8P)      :: dxyz(3)    !< Space steps.
      integer(I4P)   :: i, j, k, s !< Counter.

      call self%file_input%get(section_name="time", option_name="slices_number", val=buf_I4) ; self%slices_number = buf_I4
      if (self%slices_number > 0) then
         allocate(self%slice(self%slices_number))
         do s=1, self%slices_number
            sname = 'slice_'//trim(str(s,.true.))
            call self%file_input%get(section_name=sname, option_name='slice_itype',  val=buf_CHAR)
            self%slice(s)%slice_itype  =buf_CHAR
            call self%file_input%get(section_name=sname, option_name='slice_save',   val=buf_I4)
            self%slice(s)%slice_save   =buf_I4
            call self%file_input%get(section_name=sname, option_name='slice_ni',     val=buf_I4)
            self%slice(s)%slice_nijk(1)=buf_I4
            call self%file_input%get(section_name=sname, option_name='slice_nj',     val=buf_I4)
            self%slice(s)%slice_nijk(2)=buf_I4
            call self%file_input%get(section_name=sname, option_name='slice_nk',     val=buf_I4)
            self%slice(s)%slice_nijk(3)=buf_I4
            call self%file_input%get(section_name=sname, option_name='slice_emin_x', val=buf_R8)
            self%slice(s)%slice_emin(1)=buf_R8
            call self%file_input%get(section_name=sname, option_name='slice_emin_y', val=buf_R8)
            self%slice(s)%slice_emin(2)=buf_R8
            call self%file_input%get(section_name=sname, option_name='slice_emin_z', val=buf_R8)
            self%slice(s)%slice_emin(3)=buf_R8
            call self%file_input%get(section_name=sname, option_name='slice_emax_x', val=buf_R8)
            self%slice(s)%slice_emax(1)=buf_R8
            call self%file_input%get(section_name=sname, option_name='slice_emax_y', val=buf_R8)
            self%slice(s)%slice_emax(2)=buf_R8
            call self%file_input%get(section_name=sname, option_name='slice_emax_z', val=buf_R8)
            self%slice(s)%slice_emax(3)=buf_R8
            allocate(self%slice(s)%slice_points(3,self%slice(s)%slice_nijk(1),&
                                                  self%slice(s)%slice_nijk(2),&
                                                  self%slice(s)%slice_nijk(3)))
            dxyz(1) = (self%slice(s)%slice_emax(1) - self%slice(s)%slice_emin(1)) / self%slice(s)%slice_nijk(1)
            dxyz(2) = (self%slice(s)%slice_emax(2) - self%slice(s)%slice_emin(2)) / self%slice(s)%slice_nijk(2)
            dxyz(3) = (self%slice(s)%slice_emax(3) - self%slice(s)%slice_emin(3)) / self%slice(s)%slice_nijk(3)
            do k=1, self%slice(s)%slice_nijk(3)
               do j=1, self%slice(s)%slice_nijk(2)
                  do i=1, self%slice(s)%slice_nijk(1)
                     self%slice(s)%slice_points(1,i,j,k) = self%slice(s)%slice_emin(1) + (i - 0.5_R8P) * dxyz(1)
                     self%slice(s)%slice_points(2,i,j,k) = self%slice(s)%slice_emin(2) + (j - 0.5_R8P) * dxyz(2)
                     self%slice(s)%slice_points(3,i,j,k) = self%slice(s)%slice_emin(3) + (k - 0.5_R8P) * dxyz(3)
                  enddo
               enddo
            enddo
         enddo
      endif
      endsubroutine load_slices_from_ini_file
   endsubroutine initialize

   subroutine integrate(self, do_ghost_syncro, residual)
   !< Runge Kutta integration of field.
   class(equation_euler_cpu_object), intent(inout)         :: self             !< The equation.
   logical,                          intent(in),  optional :: do_ghost_syncro  !< Flag to do syncrous ghost update.
   real(R8P),                        intent(out), optional :: residual         !< Global residual.
   logical                                                 :: do_ghost_syncro_ !< Flag to do syncrous ghost update, local var.
   integer(I4P)                                            :: s                !< Counter.
   integer(I4P)                                            :: i_eikonal        !< Counter.
   integer(I4P), parameter                                 :: n_eikonal=2      !< Counter.
   real(R8P)                                               :: t_s

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   associate(time=>self%time, ndt=>self%dt, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             ngc=>self%ngc, nv=>self%nv, blocks_number=>self%blocks_number,        &
             inner_blocks_number=>self%field%inner_blocks_number, solids_number=>self%ib%solids_number, bcs_type=>self%bcs_type(1))

   ! do s=1, nrk
      call MPI_Barrier(MPI_COMM_WORLD, self%mpih%error)
      ! t_s = time + dt*(self%ark(s)+self%brk(s))
      ! call compute_rk_prhs_gpu_cuf(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number,    &
                                   ! dt=dt, s=s, q_gpu=self%q_gpu, prhs_gpu=self%prhs_gpu,                &
                                   ! fl_gpu=self%fl_gpu, phi_gpu=self%phi_gpu, qnrk=dt*self%brk(s))
      call self%update_ghost(q=self%q)
      if (solids_number > 0) then
         do i_eikonal=1, n_eikonal
            call MPI_Barrier(MPI_COMM_WORLD, self%mpih%error)
            ! call evolve_eikonal_q(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
            !                       phi=self%phi,                                                     &
            !                       dx=self%field%dxyz(:,1),                                          &
            !                       dy=self%field%dxyz(:,2),                                          &
            !                       dz=self%field%dxyz(:,3),                                          &
            !                       dq=self%dq,                                                       &
            !                       q=self%q)
            call self%update_ghost(q=self%q)
         enddo
         ! call invert_eikonal_field(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
         !                           q=self%q, q_ib=self%q_ib,                                 &
         !                           phi=self%phi, bcs_type=bcs_type)
      else
         ! added for restart debug...
         self%q_ib = self%q
      endif
      call self%compute_aux(q=self%q_ib, q_aux=self%q_aux)

      call MPI_Barrier(MPI_COMM_WORLD, self%mpih%error)
      ! call self%compute_residuals(ni=ni, nj=nj, nk=nk, ngc=ngc, ns=ns, blocks_number=blocks_number, &
      !                             dx = self%field%dxyz(:,1),                                        &
      !                             dy = self%field%dxyz(:,2),                                        &
      !                             dz = self%field%dxyz(:,3),                                        &
      !                             q_aux = self%q_aux,  phi = self%phi,
      ! call compute_rk_linear(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number,      &
                             ! dt=dt, q=self%q,
   ! enddo
   endassociate
   endsubroutine integrate

   subroutine load_restart_files(self, t, time)
   !< Load restart files.
   class(equation_euler_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P),                     intent(out)   :: t    !< Time iteration.
   real(R8P),                        intent(out)   :: time !< Time.

   call self%adam%load_restart_files(basename=self%restart_basename, t=t, time=time)
   call self%adam%make_comm_local_maps_ghost_bc
   endsubroutine load_restart_files

   subroutine mark_by_geo(self, delta_fine, delta_coarse, threshold, do_init)
   !< Mark blocks to be refined/derefined by a geometry proximity.
   class(equation_euler_cpu_object), intent(inout)        :: self           !< The equation.
   real(R8P),                        intent(in)           :: delta_fine     !< Maximum cell delta in fine grids.
   real(R8P),                        intent(in)           :: delta_coarse   !< Minimum cell delta in coarse grids.
   real(R8P),                        intent(in), optional :: threshold      !< Threshold for sphere proximity.
   logical,                          intent(in), optional :: do_init        !< Flag to initialize refinement.
   real(R8P)                                              :: threshold_     !< Threshold for sphere proximity, local var.
   logical                                                :: do_init_       !< Flag to initialize refinement, local var.
   real(R8P)                                              :: max_cell_delta !< Maximum cell delta.
   real(R8P)                                              :: distance       !< Value (max) of gradient of rho.
   integer(I4P)                                           :: b              !< Counter.

   do_init_ = .true.    ; if (present(do_init)) do_init_ = do_init
   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   if(do_init_) self%field%refinements_needed = [(TO_BE_DEREFINED,b=1,self%blocks_number)]
   associate (ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
              blocks_number=>self%blocks_number, dxyz=>self%field%dxyz, phi=>self%ib%phi)
      do b=1, blocks_number
         distance = 1._R8P
         if (maxval(phi(1,:,:,:,b))*minval(phi(1,:,:,:,b)) < 0._R8P) then
            distance = 0._R8P
         endif
         max_cell_delta = max_cell_delta_dist(distance=distance)

         if (maxval(dxyz(:,b)) > max_cell_delta) then
            self%field%refinements_needed(b) = TO_BE_REFINED
         elseif (maxval(dxyz(:,b)) * threshold_ < max_cell_delta) then
            self%field%refinements_needed(b) = max(self%field%refinements_needed(b), TO_BE_DEREFINED)
         else
            self%field%refinements_needed(b) = max(self%field%refinements_needed(b), TO_NOT_TOUCH)
         endif
      enddo
   endassociate
   contains
      function max_cell_delta_dist(distance) result(delta)
      !< Return the maximum cell delta given a comparison distance.
      real(R8P),          intent(in) :: distance !< Comparison distance.
      real(R8P)                      :: delta    !< Maximum cell delta admissible.

      if (abs(distance) < epsilon(0._R8P)) then
         delta = delta_fine
      else
         delta = delta_coarse
      endif
      endfunction max_cell_delta_dist
   endsubroutine mark_by_geo

   subroutine mark_by_grad_var(self, grad_tol, delta_fine, delta_coarse, q, ivar, threshold, do_init)
   !< Mark blocks to be refined/derefined by a `gradient` value.
   !< @Note Field q to which apply gradient musht have ghost cells updated.
   class(equation_euler_cpu_object), intent(inout)        :: self           !< The equation.
   real(R8P),                        intent(in)           :: grad_tol       !< Gradiend tolerance value.
   real(R8P),                        intent(in)           :: delta_fine     !< Maximum cell delta in fine grids.
   real(R8P),                        intent(in)           :: delta_coarse   !< Minimum cell delta in coarse grids.
   real(R8P),                        intent(in)           :: q(1:,          &
                                                               1-self%ngc:, &
                                                               1-self%ngc:, &
                                                               1-self%ngc:, &
                                                               1:)          !< Field component to which apply gradient.
   integer(I4P),                     intent(in)           :: ivar           !< Variable for marking.
   real(R8P),                        intent(in), optional :: threshold      !< Threshold for sphere proximity.
   logical,                          intent(in), optional :: do_init        !< Flag to initialize refinement.
   real(R8P)                                              :: threshold_     !< Threshold for sphere proximity, local var.
   integer(I4P)                                           :: ivar_          !< Variable for marking (local var).
   logical                                                :: do_init_       !< Flag to initialize refinement, local var.
   real(R8P)                                              :: max_cell_delta !< Maximum cell delta.
   real(R8P)                                              :: grad_var       !< Value (max) of gradient of var.
   integer(I4P)                                           :: b              !< Counter.

   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   do_init_   = .true.  ; if (present(do_init))   do_init_   = do_init
   if (do_init_) self%field%refinements_needed = [(TO_BE_DEREFINED,b=1,self%blocks_number)]
   associate (ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, blocks_number=>self%blocks_number, dxyz=>self%field%dxyz)
      do b=1, blocks_number
         grad_var = gradient(b=b, ni=ni, nj=nj, nk=nk, ngc=ngc, dx=dxyz(1,b), dy=dxyz(2,b), dz=dxyz(3,b))

         max_cell_delta = max_cell_delta_grad(grad=grad_var)

         if     (maxval(dxyz(:,b)) > max_cell_delta) then
            self%field%refinements_needed(b) = TO_BE_REFINED
         elseif (maxval(dxyz(:,b)) * threshold_ < max_cell_delta) then
            self%field%refinements_needed(b) = max(self%field%refinements_needed(b), TO_BE_DEREFINED)
         else
            self%field%refinements_needed(b) = max(self%field%refinements_needed(b), TO_NOT_TOUCH)
         endif
      enddo
   endassociate
   contains
      function gradient(b, ni, nj, nk, ngc, dx, dy, dz)
      !< Return gradient of q(ivar).
      integer(I4P), intent(in) :: b          !< Block index.
      integer(I4P), intent(in) :: ni         !< Grid cells number in I direction.
      integer(I4P), intent(in) :: nj         !< Grid cells number in J direction.
      integer(I4P), intent(in) :: nk         !< Grid cells number in K direction.
      integer(I4P), intent(in) :: ngc        !< Ghost cells number.
      real(R8P),    intent(in) :: dx         !< X space step.
      real(R8P),    intent(in) :: dy         !< Y space step.
      real(R8P),    intent(in) :: dz         !< Z space step.
      real(R8P)                :: gradient   !< Maximum gradient of q.
      real(R8P)                :: grad       !< Current gradient of q.
      integer(I4P)             :: i, j, k    !< Counter.
      real(R8P), parameter     :: tol=1.e-12 !< Gradient denominator tolerance.

      gradient = 0._R8P
      do k=1, nk
         do j=1, nj
            do i=1, ni
               grad = sqrt(((q(ivar,i+1,j,k,b) - q(ivar,i-1,j,k,b))/(2*dx))**2 + &
                           ((q(ivar,i,j+1,k,b) - q(ivar,i,j-1,k,b))/(2*dy))**2 + &
                           ((q(ivar,i,j,k+1,b) - q(ivar,i,j,k-1,b))/(2*dz))**2)
               grad = grad/(abs(q(ivar,i,j,k,b))+tol)
               gradient = max(gradient, grad)
            enddo
         enddo
      enddo
      endfunction gradient

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
   endsubroutine mark_by_grad_var

   subroutine print_progress(self)
   !< Print simulation progress.
   class(equation_euler_cpu_object), intent(in) :: self     !< The equation.

   associate(r=>self%mpih%myrankstr)
      print '(A)', r//''
      print '(A)', r//'t:             '//trim(str(self%it,.true.))
      print '(A)', r//'blocks number: '//trim(str(self%adam%tree%nodes_number, .true.))
      print '(A)', r//'time step:     '//trim(str(self%dt, .true.))
      print '(A)', r//'time:          '//trim(str(self%time, .true.))
   if (self%it_max <= 0) then
      print '(A)', r//'progress:      '//trim(str(int(self%time/self%time_max * 100), .true.))//'%'
   else
      print '(A)', r//'progress:      '//trim(str(int((self%it*1._R8P)/self%it_max * 100), .true.))//'%'
   endif
      print '(A)', r//''
   endassociate
   endsubroutine print_progress

   subroutine refine_uniform(self, refinement_levels)
   !< Refine all blocks uniformly.
   class(equation_euler_cpu_object), intent(inout) :: self              !< The equation.
   integer(I4P),                     intent(in)    :: refinement_levels !< Number of refinement to be performed.
   integer(I4P)                                    :: l                 !< Counter.

   do l=1, refinement_levels
      call self%adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
      call self%adam%amr_update(do_blocks_reorder=.false., do_mpi_redistribute=.true.)
   enddo
   endsubroutine

   subroutine solve(self, filename)
   !< Solve Euler system.
   class(equation_euler_cpu_object), intent(inout) :: self             !< The equation.
   character(*),                     intent(in)    :: filename         !< Input file name.
   real(R8P)                                       :: timing(1:2)      !< Tic toc timing.
   real(R8P)                                       :: timing_step(1:2) !< Tic toc timing.

   call self%initialize(filename=filename)
   if (self%restart) then
      print '(A)', self%mpih%myrankstr//'restart simulation from "'//trim(self%restart_basename)//'" files'
      call self%load_restart_files(t=self%it, time=self%time)
      print '(A)', self%mpih%myrankstr//'restart [t, time]: '//trim(str(self%it))//', '//trim(str(self%time))
   else
      call self%set_initial_conditions()
      self%time = 0._R8P
      self%it = 0
   endif
   if (self%ib%solids_number > 0) call self%ib%update_phi

   call self%amr_update

   call self%save_simulation_data

   call self%mpih%barrier(tictoc=.true., timing=timing(1), single=.true.)
   integration: do
      call self%mpih%barrier(tictoc=.true., timing=timing_step(1), single=.true.)
      self%it = self%it + 1

      if (self%save_memory_status) then
         call save_memory_cpu_status(file_name='memory_cpu-'//self%mpih%myrankstr//'.dat', tag=str(self%it,.true.))
      endif

      if (mod(self%it,self%amr_frequency)==0) then
         call self%mpih%barrier(tictoc=.true.)
         call self%amr_update
         call self%mpih%barrier(tictoc=.true.)
         print '(A, F18.10)', self%mpih%myrankstr//'step timing (AMR): ', self%mpih%tictoc_timing()
      endif

      call self%compute_dt

      call self%integrate

      self%time = self%time + self%dt

      call self%print_progress

      call self%save_simulation_data

      if (((self%it_max <= 0).and.(self%time >= self%time_max)).or.((self%it>=self%it_max).and.(self%it_max > 0))) exit integration

      call self%mpih%barrier(tictoc=.true., timing=timing_step(2), single=.true.)
      print '(A, F18.10)', self%mpih%myrankstr//'step timing: ', timing_step(2) - timing_step(1)
   enddo integration
   call self%mpih%barrier(tictoc=.true., timing=timing(2), single=.true.)
   print '(A, F18.10)', self%mpih%myrankstr//'averaged timing: ', (timing(2) - timing(1))/self%it

   call self%save_simulation_data
   endsubroutine solve

   !subroutine runge_kutta_initialize(self)
   !!< Initialize Runge-Kutta data.
   !class(equation_euler_cpu_object), intent(inout) :: self !< The equation.

   !call self%file_input%get(section_name='time', option_name='nrk', val=self%nrk)
   !allocate(self%ark(self%nrk), self%brk(self%nrk))
   !select case(self%nrk)
   !case(3_I4P)
   !   self%ark(:)=[8._R8P  /15._R8P, 5._R8P  /12._R8P, 3._R8P  /4._R8P]
   !   self%brk(:)=[0._R8P, -17._R8P/60._R8P , -5._R8P /12._R8P]
   !case(4_I4P)
   !   self%ark(:) = [8._R8P/17._R8P,17._R8P /60._R8P,5._R8P /12._R8P,3._R8P/4._R8P]
   !   self%brk(:) = [0._R8P,-15._R8P/68._R8P,-17._R8P/60._R8P,-5._R8P/12._R8P]
   !endselect
   !endsubroutine runge_kutta_initialize

   subroutine save_simulation_data(self)
   !< Save all simulation data.
   class(equation_euler_cpu_object), intent(inout) :: self                 !< The equation.
   logical                                         :: is_update_ghost_done !< Flag to minimize ghosts-update-calls for IO.
   integer(I4P)                                    :: s                    !< Slices counter.

   ! update ghost cells if necessary
   is_update_ghost_done = .false.
   if (self%slices_number>0) then
      do s=1, self%slices_number
         if (mod(self%it,self%slice(s)%slice_save)==0.or.self%it==self%it_max.or.&
            (((self%it_max <= 0).and.(self%time >= self%time_max)).or.((self%it>=self%it_max).and.(self%it_max > 0)))) then
            if (.not.is_update_ghost_done) then
               call self%update_ghost(q=self%q)
               call self%update_ghost(q=self%q_aux)
               is_update_ghost_done = .true.
            endif
         endif
      enddo
   endif
   ! save data
   call self%save_hdf5
   call self%save_restart_files
   call self%save_slices
   endsubroutine save_simulation_data

   subroutine save_hdf5(self, output_basename)
   !< Save simulation data in HDF5 format.
   class(equation_euler_cpu_object), intent(inout)        :: self             !< The equation.
   character(*),                     intent(in), optional :: output_basename  !< Output basename.
   character(:), allocatable                              :: output_basename_ !< Output basename, local var.

   if (mod(self%it,self%it_save)==0.or.self%it==self%it_max.or.&
      (((self%it_max <= 0).and.(self%time >= self%time_max)).or.((self%it>=self%it_max).and.(self%it_max > 0)))) then
      call self%mpih%barrier(tictoc=.true.)
      print '(A)', self%mpih%myrankstr//'save HDF5 files t: '//trim(str(self%it,.true.))//', time: '//&
                   trim(str(self%time,.true.))
      output_basename_ = trim(self%output_basename)//'-'//trim(strz(self%it,9))
      if (present(output_basename)) output_basename_ = trim(output_basename)
      call self%adam%save_hdf5(basename=trim(output_basename_),                   &
                               q=self%field%q,                                    &
                               q_aux=self%q_aux,                                  &
                               q_name=['rho','rhu','rhv','rhw','rhe'],            &
                               q_aux_name=['r','u','v','w','y','t','p','h','c'],  &
                               with_cell_morton=.true.)
      call self%mpih%barrier(tictoc=.true.)
      print '(A, F18.10)', self%mpih%myrankstr//'step timing (save HDF5): ', self%mpih%tictoc_timing()
   endif
   endsubroutine save_hdf5

   subroutine save_restart_files(self)
   !< Save restart files.
   class(equation_euler_cpu_object), intent(inout) :: self !< The equation.

   if (mod(self%it,self%restart_save)==0) then
      call self%mpih%barrier(tictoc=.true.)
      print '(A)', self%mpih%myrankstr//'save restart files t: '//trim(str(self%it,.true.))//', time: '//&
                   trim(str(self%time,.true.))
      call self%adam%save_restart_files(basename=self%restart_basename, t=self%it, time=self%time)
      call self%save_hdf5(output_basename=self%restart_basename)
      call self%mpih%barrier(tictoc=.true.)
      print '(A, F18.10)', self%mpih%myrankstr//'step timing (save restart): ', self%mpih%tictoc_timing()
   endif
   endsubroutine save_restart_files

   subroutine save_slices(self)
   !< Save simulation data slices.
   class(equation_euler_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                                    :: s    !< Slices counter.

   if (self%slices_number>0) then
      do s=1, self%slices_number
         if (self%it>0) then
            if (mod(self%it,self%slice(s)%slice_save)==0.or.self%it==self%it_max.or.&
               (((self%it_max <= 0).and.(self%time >= self%time_max)).or.((self%it>=self%it_max).and.(self%it_max > 0)))) then
               call self%mpih%barrier(tictoc=.true.)
               print '(A)', self%mpih%myrankstr//'save slice n: '//trim(str(s,.true.))//&
                            ', t: '//trim(str(self%it,.true.))//', time: '//trim(str(self%time,.true.))
               call self%adam%save_slice(points=self%slice(s)%slice_points,                               &
                                         itype=trim(self%slice(s)%slice_itype),                           &
                                         basename=trim(self%output_basename)//                            &
                                                  '-slice_'//trim(strz(s,2))//'-'//trim(strz(self%it,9)), &
                                         q=self%field%q,                                                  &
                                         q_name=['rho','rhu','rhv','rhw','rhe'],                          &
                                         phi=self%ib%phi(:,:,:,:,1))
               call self%mpih%barrier(tictoc=.true.)
               print '(A, F18.10)', self%mpih%myrankstr//'step timing (save slice-'//trim(str(s,.true.)) //'): ', &
                                    self%mpih%tictoc_timing()
            endif
         endif
      enddo
   endif
   endsubroutine save_slices

   subroutine set_boundary_conditions(self, q)
   !< Set boundary conditions of equation.
   class(equation_euler_cpu_object), intent(in)    :: self              !< The equation.
   real(R8P),                        intent(inout) :: q(1:,         &
                                                        1-self%ngc:,&
                                                        1-self%ngc:,&
                                                        1-self%ngc:,1:) !< Conservative variables.

   !if (allocated(self%base_gpu%local_map_bc_crown_gpu)) call set_bc(nv=self%nv, ngc=self%ngc,                          &
   !                                                                 local_map_bc=self%base_gpu%local_map_bc_crown_gpu, &
   !                                                                 x_cell_gpu=self%base_gpu%x_cell_gpu,               &
   !                                                                 y_cell_gpu=self%base_gpu%y_cell_gpu,               &
   !                                                                 z_cell_gpu=self%base_gpu%z_cell_gpu,               &
   !                                                                 fec_1_6_array_gpu=self%base_gpu%fec_1_6_array_gpu, &
   !                                                                 q_bc_vars_gpu=self%bc_vars_gpu,                    &
   !                                                                 gamma_fluid=self%gamma_fluid,                      &
   !                                                                 dha_star=self%dha_star,                            &
   !                                                                 cv_star=self%cv_star,                              &
   !                                                                 R_star=self%R_star)
   !contains
   !   subroutine set_bc(nv, ngc, local_map_bc,                    &
   !                     x_cell_gpu, y_cell_gpu, z_cell_gpu,       &
   !                     fec_1_6_array_gpu, q_bc_vars_gpu,         &
   !                     gamma_fluid, dha_star,  cv_star, R_star)
   !   integer(I4P), intent(in)         :: nv                      !< Number of variables.
   !   integer(I4P), intent(in)         :: ngc                     !< Ghost cells number.
   !   integer(I8P), intent(in), device :: local_map_bc(:,:,:)     !< Local map for BC ghost cells.
   !   real(R8P),    intent(in), device :: x_cell_gpu(1:,1-ngc:)   !< Conservative variables.
   !   real(R8P),    intent(in), device :: y_cell_gpu(1:,1-ngc:)   !< Conservative variables.
   !   real(R8P),    intent(in), device :: z_cell_gpu(1:,1-ngc:)   !< Conservative variables.
   !   integer(I4P), intent(in), device :: fec_1_6_array_gpu(:)    !< Local map for BC ghost cells.
   !   real(R8P),    intent(in), device :: q_bc_vars_gpu(:,:)      !< Boundary variables.
   !   integer(I4P)                     :: b                       !< Counter.
   !   integer(I4P)                     :: c, i, j, k, v           !< Counter.
   !   integer(I4P)                     :: idelta                  !< IJK delta step for extrapolation.
   !   integer(I4P)                     :: jdelta                  !< IJK delta step for extrapolation.
   !   integer(I4P)                     :: kdelta                  !< IJK delta step for extrapolation.
   !   integer(I4P)                     :: bc_type                 !< Boundary condition type.
   !   integer(I4P)                     :: crown                   !< Crown counter.
   !   integer(I4P)                     :: iercuda                 !< Error trapping flag for CUDAFortran.
   !   integer(I4P)                     :: fec                     !< Boundary fec (1 to 26).
   !   integer(I4P)                     :: fec_1_6                 !< Boundary fec (1 to 6).
   !   real(R8P)                        :: gamma_fluid             !< Gamma.
   !   real(R8P)                        :: dha_star                !< Entalpy formation.
   !   real(R8P)                        :: cv_star                 !< Constant volume specific heat.
   !   real(R8P)                        :: R_star                  !< Gas constant.

   !   do crown=1, ngc
   !      !$cuf kernel do(1) <<<*,*>>>
   !      do c=1, size(local_map_bc, dim=1)
   !         b = local_map_bc(c, 1 ,crown)
   !         if (b>0) then
   !            i       = local_map_bc(c, 2 ,crown)
   !            j       = local_map_bc(c, 3 ,crown)
   !            k       = local_map_bc(c, 4 ,crown)
   !            idelta  = local_map_bc(c, 5 ,crown)
   !            jdelta  = local_map_bc(c, 6 ,crown)
   !            kdelta  = local_map_bc(c, 7 ,crown)
   !            bc_type = local_map_bc(c, 8 ,crown)
   !            fec     = local_map_bc(c, 9 ,crown)
   !            fec_1_6 = fec_1_6_array_gpu(fec)
   !            if (bc_type == BC_EXTRAPOLATION) then
   !               do v=1, nv
   !                  q_gpu(b,i,j,k,v) = q_gpu(b,i-idelta,j-jdelta,k-kdelta,v)
   !               enddo
   !            else if (bc_type == BC_INFLOW) then
   !                q_gpu(b,i,j,k,1) = q_bc_vars_gpu(1, fec_1_6)
   !                q_gpu(b,i,j,k,2) = q_bc_vars_gpu(1, fec_1_6)* q_bc_vars_gpu(2, fec_1_6)
   !                q_gpu(b,i,j,k,3) = q_bc_vars_gpu(1, fec_1_6)* q_bc_vars_gpu(3, fec_1_6)
   !                q_gpu(b,i,j,k,4) = q_bc_vars_gpu(1, fec_1_6)* q_bc_vars_gpu(4, fec_1_6)
   !                q_gpu(b,i,j,k,5) = q_bc_vars_gpu(1, fec_1_6)*                         &
   !                    (cv_star*q_bc_vars_gpu(5, fec_1_6)/(q_bc_vars_gpu(1, fec_1_6)*R_star)+ &
   !                    0.5_R8P*(q_bc_vars_gpu(2, fec_1_6)**2+q_bc_vars_gpu(3, fec_1_6)**2+q_bc_vars_gpu(4, fec_1_6)**2))
   !            endif
   !         endif
   !      enddo
   !      !@cuf iercuda=cudaDeviceSynchronize()
   !   enddo
   !   endsubroutine set_bc
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self)
   !< Set initial conditions of field.
   class(equation_euler_cpu_object), intent(inout) :: self       !< The equation.
   integer(I4P)                                    :: b, i, j, k !< Counter.
   real(R8P)                                       :: x_split    !< Scalar.
   real(R8P)                                       :: uu, vv, ww !< Scalar.
   real(R8P)                                       :: rn         !< Scalar.

   associate(blocks_number=>self%blocks_number, q=>self%q, ni=>self%ni, nj=>self%nj, nk=>self%nk,              &
             ngc=>self%ngc, x_cell=>self%field%x_cell, y_cell=>self%field%y_cell, z_cell=>self%field%z_cell,   &
             g=>self%eos%g, R=>self%eos%R, cv=>self%eos%cv, cp=>self%eos%cp, ic_vars=>self%ic_vars)

      if     (self%ic_type == IC_UNIFORM) then
         do b=1, blocks_number
            do k=1, nk
               do j=1, nj
                  do i=1, ni
                  enddo
               enddo
            enddo
         enddo
      elseif (self%ic_type == IC_LEFTRIGHT) then
         x_split = ic_vars(11)
         do b=1, blocks_number
            do k=1, nk
               do j=1, nj
                  do i=1, ni
                     if(x_cell(i,b) <= x_split) then
                        q(1,i,j,k,b) = ic_vars(1)
                        q(2,i,j,k,b) = ic_vars(1)*ic_vars(2)
                        q(3,i,j,k,b) = ic_vars(1)*ic_vars(3)
                        q(4,i,j,k,b) = ic_vars(1)*ic_vars(4)
                        q(5,i,j,k,b) = ic_vars(1)*(cv*ic_vars(5)/(ic_vars(1)*R)+ &
                            0.5_R8P*(ic_vars(2)**2+ic_vars(3)**2+ic_vars(4)**2))
                     else
                        q(1,i,j,k,b) = ic_vars(6)
                        q(2,i,j,k,b) = ic_vars(6)*ic_vars(7)
                        q(3,i,j,k,b) = ic_vars(6)*ic_vars(8)
                        q(4,i,j,k,b) = ic_vars(6)*ic_vars(9)
                        q(5,i,j,k,b) = ic_vars(6)*(cv*ic_vars(10)/(ic_vars(6)*R)+&
                            0.5_R8P*(ic_vars(7)**2+ic_vars(8)**2+ic_vars(9)**2))
                     endif
                  enddo
               enddo
            enddo
         enddo
      endif
   endassociate
   endsubroutine set_initial_conditions

   subroutine update_ghost(self, q, step)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(equation_euler_cpu_object), intent(inout)        :: self            !< The equation.
   real(R8P),                        intent(inout)        :: q(1:,         &
                                                               1-self%ngc:,&
                                                               1-self%ngc:,&
                                                               1-self%ngc:,&
                                                               1:)           !< Conservative variables.
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

   ! non TBP methods
   !subroutine compute_residuals(self, ni, nj, nk, ngc, ns, blocks_number,   &
   !                             dx , dy , dz ,                  &
   !                             q_aux ,  phi ,
   !!< Compute residuals of equation.
   !class(equation_euler_cpu_object), intent(inout) :: self                      !< The equation.
   !integer(I4P), intent(in)            :: ni                                    !< Grid cells number in I direction.
   !integer(I4P), intent(in)            :: nj                                    !< Grid cells number in J direction.
   !integer(I4P), intent(in)            :: nk                                    !< Grid cells number in K direction.
   !integer(I4P), intent(in)            :: ngc                                   !< Ghost cells number.
   !integer(I4P), intent(in)            :: ns                                    !< Number of species.
   !integer(I4P), intent(in)            :: blocks_number                         !< Number of blocks.
   !real(R8P),    intent(in)            :: dx_gpu(1:)                            !< X space steps.
   !real(R8P),    intent(in)            :: dy_gpu(1:)                            !< Y space steps.
   !real(R8P),    intent(in)            :: dz_gpu(1:)                            !< Z space steps.
   !real(R8P),    intent(in)            :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   !real(R8P),    intent(inout)         :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Conservative variables.
   !real(R8P),    intent(inout)         :: fl_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Positive fluxes.
   !real(R8P),    intent(inout)         :: flx_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Positive fluxes.
   !real(R8P),    intent(inout)         :: fly_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Positive fluxes.
   !real(R8P),    intent(inout)         :: flz_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Positive fluxes.
   !real(R8P),    intent(in)            :: fd_conv_gpu(1:,1:)                    !< Convective coefficients.
   !real(R8P),    intent(in)            :: fd_coeff1_gpu(1:)                     !< Derivatives first coefficients.
   !real(R8P),    intent(in)            :: fd_coeff2_gpu(0:)                     !< Derivatives second coefficients.
   !real(R8P),    intent(inout)         :: gplus_x(1:,1:,1:,1:,1:)               !< Auxiliary variables.
   !real(R8P),    intent(inout)         :: gminus_x(1:,1:,1:,1:,1:)              !< Auxiliary variables.
   !real(R8P),    intent(inout)         :: gplus_y(1:,1:,1:,1:,1:)               !< Auxiliary variables.
   !real(R8P),    intent(inout)         :: gminus_y(1:,1:,1:,1:,1:)              !< Auxiliary variables.
   !real(R8P),    intent(inout)         :: gplus_z(1:,1:,1:,1:,1:)               !< Auxiliary variables.
   !real(R8P),    intent(inout)         :: gminus_z(1:,1:,1:,1:,1:)              !< Auxiliary variables.
   !integer(I4P),  intent(in)           :: euler_scheme                          !< Euler scheme.
   !integer(I4P),  intent(in)           :: visc_scheme                           !< Diffusive terms scheme.
   !integer(I4P), intent(in)            :: lmax                                  !< Conservative stencil size.
   !integer(I4P), intent(in)            :: iweno                                 !< Weno order.
   !integer(I4P), intent(in)            :: visc_order                            !< Diffusive terms' order.
   !integer(I4P), intent(in)            :: visc_law                              !< Viscosity temperature law.
   !real(R8P),    intent(in)            :: cp_star                               !< Constant pressure specific heat.
   !real(R8P),    intent(in)            :: cv_star                               !< Constant volume specific heat.
   !real(R8P),    intent(in)            :: R_star                                !< Gas number.
   !real(R8P),    intent(in)            :: gamma_fluid                           !< Gas gamma..
   !real(R8P),    intent(in)            :: mu_star                               !< Viscosity.
   !real(R8P),    intent(in)            :: k_star                                !< Thermal diffusion.
   !real(R8P),    intent(in)            :: dha_star                              !< Formation entalpy.
   !real(R8P),    intent(in)            :: Lewis                                 !< Lewis number
   !real(R8P),    intent(in)            :: Zeldovich                             !< Zeldovich number.
   !real(R8P),    intent(in)            :: Damkohler                             !< Damkohler number.
   !real(R8P)                           :: ib_eps                                !< Tolerance immersed boundary delta ratio.

   !if (blocks_number > 0) then
   !    if(euler_scheme == 1) then
   !       tBlock = dim3(32,8,1) ; grid = dim3(ceiling(real(blocks_number)/tBlock%x),ceiling(real(nj)/tBlock%y),1)
   !       call euler_x_central_kernel<<<grid, tBlock>>>(q_aux_gpu, flx_gpu, fd_conv_gpu, dx_gpu,     &
   !                                                     blocks_number, ni, nj, nk, ngc, ns+4, lmax)
   !       tBlock = dim3(32,8,1) ; grid = dim3(ceiling(real(blocks_number)/tBlock%x),ceiling(real(ni)/tBlock%y),1)
   !       call euler_y_central_kernel<<<grid, tBlock>>>(q_aux_gpu, fly_gpu, fd_conv_gpu, dy_gpu,     &
   !                                                     blocks_number, ni, nj, nk, ngc, ns+4, lmax)
   !       tBlock = dim3(32,8,1) ; grid = dim3(ceiling(real(blocks_number)/tBlock%x),ceiling(real(ni)/tBlock%y),1)
   !       call euler_z_central_kernel<<<grid, tBlock>>>(q_aux_gpu, flz_gpu, fd_conv_gpu, dz_gpu,     &
   !                                                     blocks_number, ni, nj, nk, ngc, ns+4, lmax)
   !    else
   !       tBlock = dim3(32,8,1) ; grid = dim3(ceiling(real(blocks_number)/tBlock%x),ceiling(real(nj)/tBlock%y),1)
   !       call euler_x_kernel<<<grid, tBlock>>>(q_aux_gpu, flx_gpu, gplus_x, gminus_x,                                       &
   !                                             blocks_number, ni, nj, nk, ngc, ns+4, iweno, dha_star, gamma_fluid, R_star, cv_star)
   !       tBlock = dim3(32,8,1) ; grid = dim3(ceiling(real(blocks_number)/tBlock%x),ceiling(real(ni)/tBlock%y),1)
   !       call euler_y_kernel<<<grid, tBlock>>>(q_aux_gpu, fly_gpu, gplus_y, gminus_y,                                        &
   !                                             blocks_number, ni, nj, nk, ngc, ns+4, iweno, dha_star, gamma_fluid, R_star, cv_star)
   !       tBlock = dim3(32,8,1) ; grid = dim3(ceiling(real(blocks_number)/tBlock%x),ceiling(real(ni)/tBlock%y),1)
   !       call euler_z_kernel<<<grid, tBlock>>>(q_aux_gpu, flz_gpu, gplus_z, gminus_z,                                        &
   !                                             blocks_number, ni, nj, nk, ngc, ns+4, iweno, dha_star, gamma_fluid, R_star, cv_star)
   !    endif
   !endif

   !   self%mpih%error = cudaGetLastError()
   !   if(self%mpih%error /= cudaSuccess) then
   !      print*,'FRA-2 CUDA ERROR ',cudaGetErrorString(self%mpih%error)
   !      call MPI_Abort(MPI_COMM_WORLD, -15,self%mpih%error)
   !      STOP
   !   endif

   !if(visc_scheme == 1) then
   !   !if(mu_star > 0.) call viscous_cuf(ni, nj, nk, ngc, blocks_number, ivis, visc_type, fd_coeff1_gpu, fd_coeff2_gpu, &
   !   !   gamma_fluid, Prandtl, q_coeff, Lewis, Zeldovich, Damkohler, dha, q_aux_gpu, dx_gpu, dy_gpu, dz_gpu, fl_gpu)
   !else
   !   if(mu_star > 0.) call viscous_part(blocks_number, ni, nj, nk, ngc, ns+4, mu_star, k_star, &
   !                                      q_aux_gpu, flx_gpu, fly_gpu, flz_gpu,                           &
   !                                      dx_gpu, dy_gpu, dz_gpu)
   !endif

   !ib_eps = 1.e-12_R8P
   !call compute_flux_diff(blocks_number, ni, nj, nk, ngc, ns+4, &
   !                       fl_gpu, flx_gpu, fly_gpu, flz_gpu, phi_gpu, &
   !                       dx_gpu, dy_gpu, dz_gpu, ib_eps)

   !endsubroutine compute_residuals_gpu

   !subroutine compute_flux_diff(blocks_number, ni, nj, nk, ngc, nv, &
   !                             fl_gpu, flx_gpu, fly_gpu, flz_gpu, phi_gpu, &
   !                             dx_gpu, dy_gpu, dz_gpu, ib_eps)
   !    integer(I4P), intent(in)          :: blocks_number, ni, nj, nk, ngc, nv
   !    real(R8P), intent(inout), device  ::  fl_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   !    real(R8P), intent(in)   , device  :: flx_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   !    real(R8P), intent(in)   , device  :: fly_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   !    real(R8P), intent(in)   , device  :: flz_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   !    real(R8P), intent(in)   , device  :: phi_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   !    real(R8P), intent(in)   , device  :: dx_gpu(1:), dy_gpu(1:), dz_gpu(1:)
   !    real(R8P), intent(in)             :: ib_eps
   !    real(R8P)                         :: delta_x, delta_y, delta_z, dx_locale, dy_locale, dz_locale
   !    integer(I4P)                      :: b, i, j, k, v, iercuda

   !   !$cuf kernel do(4) <<<*,*>>>
   !   do k=1,nk
   !   do j=1,nj
   !   do i=1,ni
   !   do b=1,blocks_number

   !      dx_locale = dx_gpu(b)
   !      ! Update net flux (procedura alternativa all'interpolazione proposta nel paper, utilizza dx_locale).
   !      if(phi_gpu(b,i,j,k,1)<0.) then
   !          if(phi_gpu(b,i+1,j,k,1)*phi_gpu(b,i-1,j,k,1)<0) then
   !              if(phi_gpu(b,i+1,j,k,1)>0.) then
   !                  delta_x = -phi_gpu(b,i,j,k,1)/(phi_gpu(b,i+1,j,k,1)-phi_gpu(b,i,j,k,1)+ib_eps)*dx_gpu(b)
   !                  dx_locale = dx_gpu(b)/2 + delta_x
   !              else !if(phi_gpu(b,i-1,j,k,1)>0) then
   !                  delta_x = -phi_gpu(b,i,j,k,1)/(phi_gpu(b,i-1,j,k,1)-phi_gpu(b,i,j,k,1)+ib_eps)*dx_gpu(b)
   !                  dx_locale = dx_gpu(b)/2 + delta_x
   !              endif
   !          endif
   !      endif

   !      dy_locale = dy_gpu(b)
   !      if(phi_gpu(b,i,j,k,1)<0.) then
   !          if(phi_gpu(b,i,j+1,k,1)*phi_gpu(b,i,j-1,k,1)<0) then
   !              if(phi_gpu(b,i,j+1,k,1)>0.) then
   !                  delta_y = -phi_gpu(b,i,j,k,1)/(phi_gpu(b,i,j+1,k,1)-phi_gpu(b,i,j,k,1)+ib_eps)*dy_gpu(b)
   !                  dy_locale = dy_gpu(b)/2 + delta_y
   !              else !if(phi_gpu(b,i-1,j,k,1)>0) then
   !                  delta_y = -phi_gpu(b,i,j,k,1)/(phi_gpu(b,i,j-1,k,1)-phi_gpu(b,i,j,k,1)+ib_eps)*dy_gpu(b)
   !                  dy_locale = dy_gpu(b)/2 + delta_y
   !              endif
   !          endif
   !      endif

   !      dz_locale = dz_gpu(b)
   !      if(phi_gpu(b,i,j,k,1)<0.) then
   !          if(phi_gpu(b,i,j,k+1,1)*phi_gpu(b,i,j,k-1,1)<0) then
   !              if(phi_gpu(b,i,j,k+1,1)>0.) then
   !                  delta_z = -phi_gpu(b,i,j,k,1)/(phi_gpu(b,i,j,k+1,1)-phi_gpu(b,i,j,k,1)+ib_eps)*dz_gpu(b)
   !                  dz_locale = dz_gpu(b)/2 + delta_z
   !              else !if(phi_gpu(b,i,j,k-1,1)>0) then
   !                  delta_z = -phi_gpu(b,i,j,k,1)/(phi_gpu(b,i,j,k-1,1)-phi_gpu(b,i,j,k,1)+ib_eps)*dz_gpu(b)
   !                  dz_locale = dz_gpu(b)/2 + delta_z
   !              endif
   !          endif
   !      endif

   !      do v=1,nv
   !         fl_gpu(b,i,j,k,v) = - (flx_gpu(b,i,j,k,v)-flx_gpu(b,i-1,j,k,v))/dx_locale &
   !                             - (fly_gpu(b,i,j,k,v)-fly_gpu(b,i,j-1,k,v))/dy_locale &
   !                             - (flz_gpu(b,i,j,k,v)-flz_gpu(b,i,j,k-1,v))/dz_locale
   !      enddo

   !   enddo
   !   enddo
   !   enddo
   !   enddo
   !   !@cuf iercuda=cudaDeviceSynchronize()
   !endsubroutine compute_flux_diff

   subroutine set_bc_ib(ni, nj, nk, ngc, nv, blocks_number, q, q_ib, phi, bcs_type)
   !< Set BC on IB cells.
   integer(I4P), intent(in)    :: ni                               !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                               !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                               !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                              !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                               !< Number of conservative varibales.
   integer(I4P), intent(in)    :: blocks_number                    !< Number of blocks.
   integer(I4P), intent(in)    :: bcs_type                         !< Immersed boundary type.
   real(R8P),    intent(in)    ::  phi(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance field.
   real(R8P),    intent(in)    ::    q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.
   real(R8P),    intent(inout) :: q_ib(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field with BC imposed in IB cells.
   integer(I4P)                :: i, j, k, b, v                    !< Counter.
   real(R8P)                   :: n_phi_x, n_phi_y, n_phi_z        !< Distance function normals.
   real(R8P)                   :: n_phi_mod, un_mod                !< Distance abs normal and normal velocity.

   select case(bcs_type)
   case(BCS_VISCOUS)
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               do b=1, blocks_number
                  if (phi(1,i,j,k,b) < 0) then
                     do v=1, nv
                        q_ib(v,i,j,k,b) = q(v,i,j,k,b)
                     enddo
                  else
                     q_ib(1,i,j,k,b) =   q(1,i,j,k,b)
                     q_ib(2,i,j,k,b) = - q(2,i,j,k,b)
                     q_ib(3,i,j,k,b) = - q(3,i,j,k,b)
                     q_ib(4,i,j,k,b) = - q(4,i,j,k,b)
                     q_ib(5,i,j,k,b) =   q(5,i,j,k,b)
                  endif
               enddo
            enddo
         enddo
      enddo
   case(BCS_EULER)
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               do b=1, blocks_number
                  if (phi(1,i,j,k,b) < 0) then
                     do v=1,nv
                        q_ib(v,i,j,k,b) = q(v,i,j,k,b)
                     enddo
                  else
                     n_phi_x = phi(1,i+1,j,k,b) - phi(1,i-1,j,k,b)
                     n_phi_y = phi(1,i,j+1,k,b) - phi(1,i,j-1,k,b)
                     n_phi_z = phi(1,i,j,k+1,b) - phi(1,i,j,k-1,b)
                     n_phi_mod = sqrt(n_phi_x**2 + n_phi_y**2 + n_phi_z**2)
                     n_phi_x = n_phi_x/n_phi_mod
                     n_phi_y = n_phi_y/n_phi_mod
                     n_phi_z = n_phi_z/n_phi_mod
                     un_mod = q(2,i,j,k,b)*n_phi_x + q(3,i,j,k,b)*n_phi_y + q(4,i,j,k,b)*n_phi_z

                     q_ib(1,i,j,k,b) = q(1,i,j,k,b)
                     q_ib(2,i,j,k,b) = q(2,i,j,k,b) - 2*un_mod*n_phi_x
                     q_ib(3,i,j,k,b) = q(3,i,j,k,b) - 2*un_mod*n_phi_y
                     q_ib(4,i,j,k,b) = q(4,i,j,k,b) - 2*un_mod*n_phi_z
                     q_ib(5,i,j,k,b) = q(5,i,j,k,b)
                  endif
               enddo
            enddo
         enddo
      enddo
   endselect
   endsubroutine set_bc_ib
endmodule adam_equation_euler_cpu_object
