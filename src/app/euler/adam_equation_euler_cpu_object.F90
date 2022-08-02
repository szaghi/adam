!< ADAM, Euler equations system class definition and solver, CPU backend.
module adam_equation_euler_cpu_object
!< ADAM, Euler equations system class definition and solver, CPU backend.

use adam_adam_object,       only : adam_object
use adam_base_cpu_object,   only : base_cpu_object
use adam_eos_ic_cpu_object, only : eos_ic_cpu_object
use adam_field_object,      only : field_object
use adam_grid_object,       only : grid_object
use adam_ib_cpu_object,     only : ib_cpu_object
use adam_mpih_object,       only : mpih_object
use adam_rk_cpu_object,     only : rk_cpu_object
use adam_slices_cpu_object, only : slices_cpu_object
use adam_tree_object,       only : tree_object
use adam_weno_cpu_object,   only : weno_cpu_object
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

! Named indexes of q_aux variables, to be initialized when ns is known.
integer(I4P) :: IRHO = 2_I4P
integer(I4P) :: IU   = 3_I4P
integer(I4P) :: IV   = 4_I4P
integer(I4P) :: IW   = 5_I4P
integer(I4P) :: IG   = 6_I4P
integer(I4P) :: IP   = 7_I4P

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

type :: equation_euler_cpu_object
   !< Euler equations system class definition and solver, base backend.
   !<
   !< The conservative variables are arranged as follows:
   !<```
   !< q(1):    rho(1)   specific density of 1st specie
   !< q(2):    rho(2)   specific density of 2nd specie
   !< ...
   !< q(ns):   rho(ns)  specific density of last specie
   !< q(ns+1): rho * u  momentum x
   !< q(ns+2): rho * v  momentum y
   !< q(ns+3): rho * w  momentum z
   !< q(ns+4): rho * E  total internal energy

   !< q(1): rho
   !< q(2): rho * u
   !< q(3): rho * v
   !< q(4): rho * w
   !< q(5): rho * E
   !<```
   !< The auxiliary variables are arranged as follows:
   !<```
   !< q_aux(1):    c(1)            specie concentration of 1st specie
   !< q_aux(2):    c(2)            specie concentration of 2nd specie
   !< ...
   !< q_aux(ns):   c(ns)           specie concentration of last specie
   !< q_aux(ns+1): rho=sum(rho(s)) density
   !< q_aux(ns+2): u               velocity x
   !< q_aux(ns+3): v               velocity y
   !< q_aux(ns+4): w               velocity z
   !< q_aux(ns+5): g               specific heats ratio
   !< q_aux(ns+6): p               pressure
   !<```
   ! Objects.
   type(mpih_object)           :: mpih          !< MPI handler.
   type(file_ini)              :: file_input    !< Nasto input file handler.
   type(adam_object)           :: adam          !< ADAM.
   type(tree_object),  pointer :: tree=>null()  !< The tree.
   type(field_object), pointer :: field=>null() !< The field.
   type(grid_object),  pointer :: grid=>null()  !< The grid.
   type(base_cpu_object)       :: base_cpu      !< The base CPU handler.
   type(rk_cpu_object)         :: rk            !< Runge Kutta solver.
   type(weno_cpu_object)       :: weno          !< WENO Kutta solver.
   type(ib_cpu_object)         :: ib            !< Immersed Boundary solver.
   type(slices_cpu_object)     :: slices        !< Slices handler.
   ! Pointers to ADAM data for easy handling.
   integer(I4P), pointer :: ngc          =>null() !< Number of ghost cells.
   integer(I4P), pointer :: ni           =>null() !< Number of cells in i direction.
   integer(I4P), pointer :: nj           =>null() !< Number of cells in j direction.
   integer(I4P), pointer :: nk           =>null() !< Number of cells in k direction.
   integer(I4P), pointer :: nb           =>null() !< Total blocks number for MPI.
   integer(I4P), pointer :: blocks_number=>null() !< Actual blocks number.
   ! Physics.
   integer(I4P)                         :: ns=1_I4P     !< Number of species.
   integer(I4P)                         :: nv=5_I4P     !< Number of variables.
   integer(I4P)                         :: nv_aux=7_I4P !< Number of auxiliary variables.
   type(eos_ic_cpu_object), allocatable :: eos(:)       !< Equations of state (cp, cv...) of each specie [1:ns].
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
   ! Solver data.
   logical        :: null_xyz(3)=[.false.,&
                                  .false.,&
                                  .false.]      !< Flag triggering 1D/2D simulations.
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
   ! Large arrays.
   real(R8P), pointer     ::     q(:,:,:,:,:)=>null() !< Conservative cell centerd variables, allocated in field object.
   real(R8P), allocatable :: q_aux(:,:,:,:,:)         !< Auxiliary cell centered variables.
   real(R8P), allocatable ::  q_ib(:,:,:,:,:)         !< Field cell with boundary set on immersed bodies.
   real(R8P), allocatable :: dq_ib(:,:,:,:,:)         !< Difference field cell with boundary set on immersed bodies.
   real(R8P), allocatable ::   q_s(:,:,:,:,:,:)       !< RK stage.
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
      procedure, pass(self) :: save_simulation_data    !< Save all simulation data.
      procedure, pass(self) :: save_restart_files      !< Save restart files.
      procedure, pass(self) :: save_hdf5               !< Save simulation data in HDF5 format.
      procedure, pass(self) :: set_boundary_conditions !< Set boundary conditions of equation.
      procedure, pass(self) :: set_initial_conditions  !< Set initial conditions of equation.
      procedure, pass(self) :: solve                   !< Solve Euler system.
      procedure, pass(self) :: update_ghost            !< Update ghost cells.
      ! private methods
      procedure, pass(self) :: compute_residuals         !< Compute residuals of equation.
      procedure, pass(self) :: compute_fluxes_convective !< Compute the conservative fluxes on a slice along a coordinate direction.
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
   integer(I4P)                                  :: b, i, j, k       !< Counter.

   associate(blocks_number=>self%blocks_number, q=>self%q, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
             ns=>self%ns, eos=>self%eos)
   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               q_aux(ns+1,i,j,k,b) = sum(q(1:ns,i,j,k,b))
               q_aux(1:ns,i,j,k,b) = q(1:ns,i,j,k,b) / q_aux(ns+1,i,j,k,b)
               q_aux(ns+2,i,j,k,b) = q(ns+1,i,j,k,b) / q_aux(ns+1,i,j,k,b)
               q_aux(ns+3,i,j,k,b) = q(ns+2,i,j,k,b) / q_aux(ns+1,i,j,k,b)
               q_aux(ns+4,i,j,k,b) = q(ns+3,i,j,k,b) / q_aux(ns+1,i,j,k,b)
               q_aux(ns+5,i,j,k,b) = dot_product(q_aux(1:ns,i,j,k,b), eos%cp) / dot_product(q_aux(1:ns,i,j,k,b), eos%cv)
               q_aux(ns+6,i,j,k,b) = (q(ns+4,i,j,k,b) - 0.5_R8P * q_aux(ns+1,i,j,k,b) * (q_aux(ns+2,i,j,k,b)**2 +   &
                                                                                         q_aux(ns+3,i,j,k,b)**2 +   &
                                                                                         q_aux(ns+4,i,j,k,b)**2)) * &
                                     (q_aux(ns+5,i,j,k,b) - 1._R8P)
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
   real(R8P)                                       :: ss      !< Speed of sound.
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
                  ss = a(p=self%q_aux(IP,i,j,k,b), r=self%q_aux(IRHO,i,j,k,b), g=self%q_aux(IG,i,j,k,b))
                  umax = max(umax, &
                             (abs(self%q_aux(IU,i,j,k,b)) + ss)/dxyz(1,b), &
                             (abs(self%q_aux(IV,i,j,k,b)) + ss)/dxyz(2,b), &
                             (abs(self%q_aux(IW,i,j,k,b)) + ss)/dxyz(3,b))
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

   call forward_main_adam_data(grid=self%adam%grid, tree=self%adam%tree, field=self%adam%field)

   call self%adam%refine_uniform(refinement_levels=self%adam%tree%iu_ref_levels, do_blocks_reorder=.false.)

   call self%adam%prune(ijkl_prune=self%adam%tree%ijkl_prune, do_blocks_reorder=.false.)

   call self%base_cpu%initialize(tree=self%adam%tree, field=self%adam%field)

   call self%ib%initialize(grid=self%adam%grid, field=self%adam%field, file_parameters=self%file_input)

   call self%rk%initialize(file_parameters=self%file_input)

   call self%weno%initialize(file_parameters=self%file_input)

   call self%slices%initialize(file_parameters=self%file_input)

   call load_equation_from_ini_file

   call load_physics_from_ini_file

   call load_amr_from_ini_file

   call load_timing_from_ini_file

   call load_ic_from_ini_file

   call load_bc_from_ini_file

   ! allocate large arrays
   associate(nv=>self%nv, nv_aux=>self%nv_aux, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, nb=>self%nb, nrk=>self%rk%nrk)
   allocate(self%q_aux(1:nv_aux, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb       ))
   allocate(self%dq_ib(1:nv,     1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb       ))
   allocate(self%q_s(  1:nv,     1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb, 1:nrk))
   endassociate

   print '(A)', self%mpih%myrankstr//'equation_euler_cpu_object%initialize finish'
   contains
      subroutine forward_main_adam_data(grid, tree, field)
      !< Forward main ADAM data to equation for easy handling.
      type(grid_object),  intent(in), target :: grid  !< The grid.
      type(tree_object),  intent(in), target :: tree  !< The tree.
      type(field_object), intent(in), target :: field !< The field.

      self%grid          => grid
      self%tree          => tree
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

      subroutine load_physics_from_ini_file
      !< Parse physics setting from input file.
      integer(I4P) :: ns     !< Number of initial species.
      real(R8P)    :: cp, cv !< Constant specific heats.
      integer(I4P) :: s      !< Counter.

      call self%file_input%get(section_name='physics', option_name='ns', val=ns)
      self%ns = ns
      ! initialize named index of q_aux array
      IRHO = ns + 1
      IU   = ns + 2
      IV   = ns + 3
      IW   = ns + 4
      IG   = ns + 5
      IP   = ns + 6
      allocate(self%eos(1:ns))
      do s=1, self%ns
         call self%file_input%get(section_name='physics_specie_'//trim(str(s,.true.)), option_name='cp', val=cp)
         call self%file_input%get(section_name='physics_specie_'//trim(str(s,.true.)), option_name='cv', val=cv)
         call self%eos(s)%initialize(cp=cp, cv=cv)
      enddo
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

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   associate(time=>self%time, dt=>self%dt, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             ngc=>self%ngc, nv=>self%nv, blocks_number=>self%blocks_number,       &
             ! nrk=>self%rk%nrk, alph=>self%alph, beta=>self%beta,                  &
             inner_blocks_number=>self%field%inner_blocks_number, solids_number=>self%ib%solids_number, bcs_type=>self%bcs_type(1))

   do s=1, self%rk%nrk
      call MPI_Barrier(MPI_COMM_WORLD, self%mpih%error)
      ! call compute_rk_stage(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                            ! alph=alph, dt=dt, s=s, q=self%q, q_s=self%q_s)

      call self%update_ghost(q=self%q_s(:,:,:,:,:,s))
      if (solids_number > 0) then
         do i_eikonal=1, n_eikonal
            call MPI_Barrier(MPI_COMM_WORLD, self%mpih%error)
            call self%ib%evolve_eikonal_q(dq=self%dq_ib, q=self%q_s(:,:,:,:,:,s))
            call self%update_ghost(q=self%q_s(:,:,:,:,:,s))
         enddo
         call set_bc_ib(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                        q=self%q_s(:,:,:,:,:,s), phi=self%ib%phi, bcs_type=bcs_type)
      endif
      call self%compute_residuals(q=self%q_s(:,:,:,:,:,s), time=time)
   enddo
   ! call advance_q(ni=ni, nj=nj, nk=nk, ngc=ngc, nrk=nrk, nv=nv, blocks_number=blocks_number, &
                  ! beta=beta, dt=dt, q_s=self%q_s, q=self%q)
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

   subroutine save_simulation_data(self)
   !< Save all simulation data.
   class(equation_euler_cpu_object), intent(inout) :: self                 !< The equation.
   logical                                         :: is_update_ghost_done !< Flag to minimize ghosts-update-calls for IO.
   integer(I4P)                                    :: s                    !< Slices counter.

   is_update_ghost_done = .false.
   call self%save_hdf5(is_update_ghost_done=is_update_ghost_done)
   call self%save_restart_files(is_update_ghost_done=is_update_ghost_done)
   if (self%slices%is_to_save(it=self%it,it_max=self%it_max,time=self%time,time_max=self%time_max).and.&
       .not.is_update_ghost_done) then
      call self%update_ghost(q=self%q)
      is_update_ghost_done = .true.
    endif
   call self%slices%save_mat(basename=self%output_basename, &
                             it=self%it,                    &
                             it_max=self%it_max,            &
                             time=self%time,                &
                             time_max=self%time_max,        &
                             adam=self%adam,                &
                             q=self%q,                      &
                             q_name=['rho','rhu','rhv','rhw','rhe'])
   endsubroutine save_simulation_data

   subroutine save_hdf5(self, is_update_ghost_done, output_basename)
   !< Save simulation data in HDF5 format.
   class(equation_euler_cpu_object), intent(inout)        :: self                 !< The equation.
   logical,                          intent(inout)        :: is_update_ghost_done !< Flag to minimize ghosts-update-calls for IO.
   character(*),                     intent(in), optional :: output_basename      !< Output basename.
   character(:), allocatable                              :: output_basename_     !< Output basename, local var.

   if (mod(self%it,self%it_save)==0.or.self%it==self%it_max.or.&
      (((self%it_max <= 0).and.(self%time >= self%time_max)).or.((self%it>=self%it_max).and.(self%it_max > 0)))) then
      if (.not.is_update_ghost_done) then
         call self%update_ghost(q=self%q)
         call self%compute_aux(q=self%q, q_aux=self%q_aux)
         is_update_ghost_done = .true.
      endif
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

   subroutine save_restart_files(self, is_update_ghost_done)
   !< Save restart files.
   class(equation_euler_cpu_object), intent(inout) :: self                 !< The equation.
   logical,                          intent(inout) :: is_update_ghost_done !< Flag to minimize ghosts-update-calls for IO.

   if (mod(self%it,self%restart_save)==0) then
      if (.not.is_update_ghost_done) then
         call self%update_ghost(q=self%q)
         is_update_ghost_done = .true.
      endif
      call self%mpih%barrier(tictoc=.true.)
      print '(A)', self%mpih%myrankstr//'save restart files t: '//trim(str(self%it,.true.))//', time: '//&
                   trim(str(self%time,.true.))
      call self%adam%save_restart_files(basename=self%restart_basename, t=self%it, time=self%time)
      call self%save_hdf5(is_update_ghost_done=is_update_ghost_done, output_basename=self%restart_basename)
      call self%mpih%barrier(tictoc=.true.)
      print '(A, F18.10)', self%mpih%myrankstr//'step timing (save restart): ', self%mpih%tictoc_timing()
   endif
   endsubroutine save_restart_files

   subroutine set_boundary_conditions(self, q)
   !< Set boundary conditions of equation.
   class(equation_euler_cpu_object), intent(in)    :: self              !< The equation.
   real(R8P),                        intent(inout) :: q(1:,         &
                                                        1-self%ngc:,&
                                                        1-self%ngc:,&
                                                        1-self%ngc:,1:) !< Conservative variables.
   integer(I4P)                                    :: crown             !< Crown counter.
   integer(I4P)                                    :: fec               !< Boundary fec (1 to 26).
   integer(I4P)                                    :: fec_1_6           !< Boundary fec (1 to 6).
   integer(I4P)                                    :: c, i, j, k, b     !< Counter.
   integer(I4P)                                    :: idelta            !< IJK delta step for extrapolation.
   integer(I4P)                                    :: jdelta            !< IJK delta step for extrapolation.
   integer(I4P)                                    :: kdelta            !< IJK delta step for extrapolation.
   integer(I4P)                                    :: bc_type           !< Boundary condition type.

   associate(ngc=>self%ngc)
   if (allocated(self%tree%local_map_bc_crown)) then
      do crown=1, ngc
         do c=1, size(self%tree%local_map_bc_crown, dim=1)
            b = self%tree%local_map_bc_crown(c, 1 ,crown)
            if (b>0) then
               i       = self%tree%local_map_bc_crown(c, 2 ,crown)
               j       = self%tree%local_map_bc_crown(c, 3 ,crown)
               k       = self%tree%local_map_bc_crown(c, 4 ,crown)
               idelta  = self%tree%local_map_bc_crown(c, 5 ,crown)
               jdelta  = self%tree%local_map_bc_crown(c, 6 ,crown)
               kdelta  = self%tree%local_map_bc_crown(c, 7 ,crown)
               bc_type = self%tree%local_map_bc_crown(c, 8 ,crown)
               fec     = self%tree%local_map_bc_crown(c, 9 ,crown)
               fec_1_6 = FEC_1_6_ARRAY(fec)
               select case(bc_type)
               case(BC_EXTRAPOLATION)
                   q(:,i,j,k,b) = q(:,i-idelta,j-jdelta,k-kdelta,b)
               case(BC_INFLOW)
                   ! q(1,i,j,k,b) = q_bc_vars(1, fec_1_6)
                   ! q(2,i,j,k,b) = q_bc_vars(1, fec_1_6) * q_bc_vars(2, fec_1_6)
                   ! q(3,i,j,k,b) = q_bc_vars(1, fec_1_6) * q_bc_vars(3, fec_1_6)
                   ! q(4,i,j,k,b) = q_bc_vars(1, fec_1_6) * q_bc_vars(4, fec_1_6)
                   ! q(5,i,j,k,b) = q_bc_vars(1, fec_1_6)*                         &
                   !     (cv_star*q_bc_vars(5, fec_1_6)/(q_bc_vars(1, fec_1_6)*R_star)+ &
                   !     0.5_R8P*(q_bc_vars(2, fec_1_6)**2+q_bc_vars(3, fec_1_6)**2+q_bc_vars(4, fec_1_6)**2))
               endselect
            endif
         enddo
      enddo
   endif
   endassociate
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
             R=>self%eos(1)%R, cv=>self%eos(1)%cv, ic_vars=>self%ic_vars)

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

   ! private methods
   subroutine compute_residuals(self, q, time)
   !< Compute residuals of equation.
   class(equation_euler_cpu_object), intent(inout) :: self                          !< The equation.
   real(R8P),                        intent(inout) :: q(1:,                    &
                                                        1-self%field%grid%ngc:,&
                                                        1-self%field%grid%ngc:,&
                                                        1-self%field%grid%ngc:,&
                                                        1:)                         !< Conservative variables.
   real(R8P),                        intent(in)    :: time                          !< Time.
   real(R8P)                                       :: fluxes_x(1:self%field%nv,     &
                                                               0:self%field%grid%ni,&
                                                               1:self%field%grid%nj,&
                                                               1:self%field%grid%nk)!< Convective fluxes in x direction.
   real(R8P)                                       :: fluxes_y(1:self%field%nv,     &
                                                               1:self%field%grid%ni,&
                                                               0:self%field%grid%nj,&
                                                               1:self%field%grid%nk)!< Convective fluxes in y direction.
   real(R8P)                                       :: fluxes_z(1:self%field%nv,     &
                                                               1:self%field%grid%ni,&
                                                               1:self%field%grid%nj,&
                                                               0:self%field%grid%nk)!< Convective fluxes in z direction.
   integer(I4P)                                    :: b, i, j, k                    !< Counter.

   associate(ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk, &
             ngc=>self%field%grid%ngc, blocks_number=>self%field%blocks_number,      &
             dxyz=>self%field%dxyz, ns=>self%ns, q_aux=>self%q_aux)
   call self%compute_aux(q=q, q_aux=q_aux)
   fluxes_x = 0._R8P
   fluxes_y = 0._R8P
   fluxes_z = 0._R8P
   do b=1, blocks_number
      if (.not.self%null_xyz(1)) then ! convective fluxes along x direction
         do k=1, nk
            do j=1, nj
               call self%compute_fluxes_convective(gc=ngc, n=ni,                          &
                                                   c         =    q_aux(1:ns, :,  j,k,b), &
                                                   rho       =    q_aux(ns+1, :,  j,k,b), &
                                                   un        =    q_aux(ns+2, :,  j,k,b), & ! u
                                                   ut1       =    q_aux(ns+3, :,  j,k,b), & ! v
                                                   ut2       =    q_aux(ns+4, :,  j,k,b), & ! w
                                                   g         =    q_aux(ns+5, :,  j,k,b), &
                                                   p         =    q_aux(ns+6, :,  j,k,b), &
                                                   f_rho     = fluxes_x(1:ns,0:ni,j,k),   &
                                                   f_rho_un  = fluxes_x(ns+1,0:ni,j,k),   & ! u
                                                   f_rho_ut1 = fluxes_x(ns+2,0:ni,j,k),   & ! v
                                                   f_rho_ut2 = fluxes_x(ns+3,0:ni,j,k),   & ! w
                                                   f_rho_e   = fluxes_x(ns+4,0:ni,j,k))
            enddo
         enddo
      endif
      if (.not.self%null_xyz(2)) then ! convective fluxes along y direction
         do k=1, nk
            do i=1, ni
               call self%compute_fluxes_convective(gc=ngc, n=nj,                          &
                                                   c         =    q_aux(1:ns,i, :,  k,b), &
                                                   rho       =    q_aux(ns+1,i, :,  k,b), &
                                                   un        =    q_aux(ns+3,i, :,  k,b), & ! v
                                                   ut1       =    q_aux(ns+2,i, :,  k,b), & ! u
                                                   ut2       =    q_aux(ns+4,i, :,  k,b), & ! w
                                                   g         =    q_aux(ns+5,i, :,  k,b), &
                                                   p         =    q_aux(ns+6,i, :,  k,b), &
                                                   f_rho     = fluxes_y(1:ns,i,0:nj,k),   &
                                                   f_rho_un  = fluxes_y(ns+2,i,0:nj,k),   & ! v
                                                   f_rho_ut1 = fluxes_y(ns+1,i,0:nj,k),   & ! u
                                                   f_rho_ut2 = fluxes_y(ns+3,i,0:nj,k),   & ! w
                                                   f_rho_e   = fluxes_y(ns+4,i,0:nj,k))
            enddo
         enddo
      endif
      if (.not.self%null_xyz(3)) then ! convective fluxes along z direction
         do j=1, nj
            do i=1, ni
               call self%compute_fluxes_convective(gc=ngc, n=nk,                          &
                                                   c         =    q_aux(1:ns,i,j, :,  b), &
                                                   rho       =    q_aux(ns+1,i,j, :,  b), &
                                                   un        =    q_aux(ns+4,i,j, :,  b), & ! w
                                                   ut1       =    q_aux(ns+2,i,j, :,  b), & ! u
                                                   ut2       =    q_aux(ns+3,i,j, :,  b), & ! v
                                                   g         =    q_aux(ns+5,i,j, :,  b), &
                                                   p         =    q_aux(ns+6,i,j, :,  b), &
                                                   f_rho     = fluxes_z(1:ns,i,j,0:nk),   &
                                                   f_rho_un  = fluxes_z(ns+3,i,j,0:nk),   & ! w
                                                   f_rho_ut1 = fluxes_z(ns+1,i,j,0:nk),   & ! u
                                                   f_rho_ut2 = fluxes_z(ns+2,i,j,0:nk),   & ! v
                                                   f_rho_e   = fluxes_z(ns+4,i,j,0:nk))
            enddo
         enddo
      endif
      ! residuals
      do k=1, nk
         do j=1, nj
            do i=1, ni
               q(:,i,j,k,b) = (fluxes_x(:,i-1,j  ,k  ) - fluxes_x(:,i,j,k)) / dxyz(1,b) + &
                              (fluxes_y(:,i  ,j-1,k  ) - fluxes_y(:,i,j,k)) / dxyz(2,b) + &
                              (fluxes_z(:,i  ,j  ,k-1) - fluxes_z(:,i,j,k)) / dxyz(3,b)
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine compute_residuals

   subroutine compute_fluxes_convective(self, gc, n,                &
                                        c, rho, un, ut1, ut2, g, p, &
                                        f_rho, f_rho_un, f_rho_ut1, f_rho_ut2, f_rho_e)
   !< Compute the conservative fluxes on a slice of cells along a coordinate direction.
   class(equation_euler_cpu_object), intent(in)    :: self              !< The equation.
   integer(I4P),                     intent(in)    :: gc                !< Number of ghost cells used.
   integer(I4P),                     intent(in)    :: n                 !< Number of cells.
   real(R8P),                        intent(in)    ::      c(1:,1-gc:)  !< Species concentration        [1-gc:n+gc,1:ns].
   real(R8P),                        intent(in)    ::       rho(1-gc:)  !< Density                      [1-gc:n+gc].
   real(R8P),                        intent(in)    ::        un(1-gc:)  !< Normal velocity              [1-gc:n+gc].
   real(R8P),                        intent(in)    ::       ut1(1-gc:)  !< Tangential velocity 1        [1-gc:n+gc].
   real(R8P),                        intent(in)    ::       ut2(1-gc:)  !< Tangential velocity 2        [1-gc:n+gc].
   real(R8P),                        intent(in)    ::         g(1-gc:)  !< Specific heats ratio         [1-gc:n+gc].
   real(R8P),                        intent(in)    ::         p(1-gc:)  !< Pressure                     [1-gc:n+gc].
   real(R8P),                        intent(inout) :: f_rho(1:, 0:)     !< Flux of mass                 [0:n,1:ns].
   real(R8P),                        intent(inout) ::  f_rho_un(0:)     !< Flux normal momentums        [0:n].
   real(R8P),                        intent(inout) :: f_rho_ut1(0:)     !< Flux of tangential1 momentum [0:n].
   real(R8P),                        intent(inout) :: f_rho_ut2(0:)     !< Flux of tangential2 momentum [0:n].
   real(R8P),                        intent(inout) ::   f_rho_e(0:)     !< Flux energy                  [0:n].
   real(R8P)                                       ::   fluxes(3)       !< 1D fluxes.
   real(R8P)                                       :: qr(1:3,1:2,0:n+1) !< Reconstructed variables.
   integer(I4P)                                    :: i                 !< Counter.

   call reconstruct_interfaces
   do i=0, n
      ! computing normal fluxes solving Riemann problem
      call solve_riemann(r1=qr(1,2,i  ), u1=qr(2,2,i  ), p1=qr(3,2,i  ), g1=g(i  ), &
                         r4=qr(1,1,i+1), u4=qr(2,1,i+1), p4=qr(3,1,i+1), g4=g(i+1), F=fluxes)
      if (fluxes(1)>0._R8P) then
           f_rho(:,i) = fluxes(1) * c(:,i)
          f_rho_un(i) = fluxes(2)
         f_rho_ut1(i) = fluxes(1) * ut1(i)
         f_rho_ut2(i) = fluxes(1) * ut2(i)
           f_rho_e(i) = fluxes(3) + 0.5_R8P * fluxes(1) * (ut1(i)**2 + ut2(i)**2)
      else
           f_rho(:,i) = fluxes(1) * c(:,i+1)
          f_rho_un(i) = fluxes(2)
         f_rho_ut1(i) = fluxes(1) * ut1(i+1)
         f_rho_ut2(i) = fluxes(1) * ut2(i+1)
           f_rho_e(i) = fluxes(3) + 0.5_R8P * fluxes(1) * (ut1(i+1)**2 + ut2(i+1)**2)
      endif
   enddo
   contains
      subroutine reconstruct_interfaces
      !< The reconstruction is done in pseudo characteristic variables.
      integer(I4P) :: i, j, f, v            !< Counter.
      real(R8P)    :: qm(1:3,1:2)           !< Mean primitive variables.
      real(R8P)    :: Lqm(1:3,1:3,1:2)      !< Left eigenvalues matrix of mean primitive variables.
      real(R8P)    :: Rqm(1:3,1:3,1:2)      !< Right eigenvalues matrix of mean primitive variables.
      real(R8P)    :: c(1:3,1:2,1-gc:-1+gc) !< Pseudo characteristic variables.
      real(R8P)    :: cr(1:3,1:2)           !< Pseudo characteristic variables reconstructed.

      select case(self%weno%weno_s)
      case(1_I4P)
         do i=0, n+1
            qr(:,1,i) = [rho(i),un(i),p(i)]
            qr(:,2,i) = qr(:,1,i)
         enddo
      case(2_I4P,3_I4P)
      ! compute WENO reconstruction
         do i=0, n+1
            ! compute pseudo characteristic variables
            do f=1, 2
               if (i==0  .and.f==1) cycle
               if (i==n+1.and.f==2) cycle
               qm(:,f) = [0.5_R8P * (rho(i+f-2) + rho(i+f-1)), &
                          0.5_R8P * ( un(i+f-2) +  un(i+f-1)), &
                          0.5_R8P * (  p(i+f-2) +   p(i+f-1))]
            enddo
            do f=1, 2
               if (i==0  .and.f==1) cycle
               if (i==n+1.and.f==2) cycle
               Lqm(:, :, f) =  left_eigenvectors(q=qm(:,f), g=g(i))
               Rqm(:, :, f) = right_eigenvectors(q=qm(:,f), g=g(i))
            enddo
            do j=i+1-gc, i-1+gc
               do f=1, 2
                  if (i==0  .and.f==1) cycle
                  if (i==n+1.and.f==2) cycle
                  do v=1, 3
                     c(v,f,j-i) = dot_product(Lqm(v, :, f), [rho(j), un(j), p(j)])
                  enddo
               enddo
            enddo

            ! compute WENO reconstruction of pseudo charteristic variables
            do v=1, 3
               cr(v,:) = self%weno%reconstructed(s=self%weno%weno_s, v=c(v,:,1-self%weno%weno_s:-1+self%weno%weno_s))
            enddo

            ! trasform back reconstructed pseudo charteristic variables to primitive ones
            do f=1, 2
               if (i==0  .and.f==1) cycle
               if (i==n+1.and.f==2) cycle
               do v=1, 3
                  qr(v,f,i) = dot_product(Rqm(v,:,f), cr(:,f))
               enddo
            enddo
         enddo
      endselect
      endsubroutine reconstruct_interfaces
   endsubroutine compute_fluxes_convective

   ! non TBP methods
   subroutine solve_riemann(r1, u1, p1, g1, r4, u4, p4, g4, F)
   !< Solve the Riemann problem between the state $1$ and $4$ using the (local) Lax Friedrichs (Rusanov) solver.
   real(R8P), intent(in)  :: r1      !< Density of state 1.
   real(R8P), intent(in)  :: u1      !< Velocity of state 1.
   real(R8P), intent(in)  :: p1      !< Pressure of state 1.
   real(R8P), intent(in)  :: g1      !< Specific heats ratio of state 1.
   real(R8P), intent(in)  :: r4      !< Density of state 4.
   real(R8P), intent(in)  :: u4      !< Velocity of state 4.
   real(R8P), intent(in)  :: p4      !< Pressure of state 4.
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
      real(R8P), intent(in) :: p       !< Pressure.
      real(R8P), intent(in) :: r       !< Density.
      real(R8P), intent(in) :: u       !< Velocity.
      real(R8P), intent(in) :: g       !< Specific heats ratio.
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

   pure function left_eigenvectors(q, g) result(eig)
   !< Return the left eigenvectors matrix `L` as `dF/dP = A = R ^ L` `P`` being the primitive variables and `F` the fluxes.
   !<
   !< Primitive variables `q` are: density `q(1)`, normal velocity `q(2)`, pressure `q(3)`.
   real(R8P), intent(in) :: q(3)          !< Primitive variables.
   real(R8P), intent(in) :: g             !< Specific heats ratio.
   real(R8P)             :: eig(1:3, 1:3) !< Eigenvectors.
   real(R8P)             :: gp            !< `g*p`.
   real(R8P)             :: gp_a          !< `g*p/a`.

   gp = g * q(3)
   gp_a = gp / a(p=q(3), r=q(1), g=g)
   eig = 0._R8P
               eig(1, 1) = 0._R8P    ; eig(1, 2) = -gp_a  ; eig(1, 3) =  1._R8P
   if (q(1)>0) eig(2, 1) = gp / q(1) ; eig(2, 2) = 0._R8P ; eig(2, 3) = -1._R8P
               eig(3, 1) = 0._R8P    ; eig(3, 2) =  gp_a  ; eig(3, 3) =  1._R8P
   endfunction left_eigenvectors

   pure function right_eigenvectors(q, g) result(eig)
   !< Return the right eigenvectors matrix `R` as `dF/dP = A = R ^ L` `P`` being the primitive variables and `F` the fluxes.
   !<
   !< Primitive variables `q` are: density `q(1)`, normal velocity `q(2)`, pressure `q(3)`.
   real(R8P), intent(in) :: q(3)          !< Primitive variables.
   real(R8P), intent(in) :: g             !< Specific heats ratio.
   real(R8P)             :: eig(1:3, 1:3) !< Eigenvectors.
   real(R8P)             :: gp            !< `g*p`.
   real(R8P)             :: gp_inv        !< `1/(g*p)`.

   gp = g * q(3)
   gp_inv = 1._R8P / gp
   eig(1, 1) =  0.5_R8P * q(1) * gp_inv                   ; eig(1, 2) = q(1) * gp_inv ; eig(1, 3) =  eig(1, 1)
   eig(2, 1) = -0.5_R8P * a(p=q(3), r=q(1), g=g) * gp_inv ; eig(2, 2) = 0._R8P        ; eig(2, 3) = -eig(2, 1)
   eig(3, 1) =  0.5_R8P                                   ; eig(3, 2) = 0._R8P        ; eig(3, 3) =  eig(3, 1)
   endfunction right_eigenvectors

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
   real(R8P), intent(in) :: g       !< Specific heats ratio \(\frac{{c_p}}{{c_v}}\).
   real(R8P), intent(in) :: p       !< Pressure.
   real(R8P), intent(in) :: r       !< Density.
   real(R8P), intent(in) :: u       !< Module of velocity vector.
   real(R8P)             :: entalpy !< Total specific entalpy (per unit of mass).

   entalpy = g * p / ((g - 1._R_P) * r) + 0.5_R_P * u*u
   endfunction H

   ! RK
   subroutine advance_q(ni, nj, nk, ngc, nv, nrk, blocks_number, beta, dt, q_s, q)
   !< Advance q by means of RK stages.
   integer(I4P), intent(in)    :: ni                                 !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                 !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                 !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                !< Ghost grid number.
   integer(I4P), intent(in)    :: nv                                 !< Number of conservative varibales.
   integer(I4P), intent(in)    :: nrk                                !< Number of RK stages.
   integer(I4P), intent(in)    :: blocks_number                      !< Number of blocks.
   real(R8P),    intent(in)    :: beta(:)                            !< RK betaa coefficients.
   real(R8P),    intent(in)    :: dt                                 !< Time step.
   real(R8P),    intent(in)    :: q_s(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:) !< RK stage.
   real(R8P),    intent(inout) ::   q(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative variables.
   integer(I4P)                :: s                                  !< Counter.

   do s=1, nrk
      q(1:nv,1:ni,1:nj,1:nk,1:blocks_number) =   q(1:nv,1:ni,1:nj,1:nk,1:blocks_number  ) + &
                                               q_s(1:nv,1:ni,1:nj,1:nk,1:blocks_number,s) * dt * beta(s)
   enddo
   endsubroutine advance_q

   subroutine compute_rk_stage(ni, nj, nk, ngc, nv, blocks_number, alph, dt, s, q, q_s)
   !< Initialize RK stage with q.
   integer(I4P), intent(in)    :: ni                                 !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                 !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                 !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                                 !< Number of conservative varibales.
   integer(I4P), intent(in)    :: blocks_number                      !< Number of blocks.
   real(R8P),    intent(in)    :: alph(:,:)                          !< RK alpha coefficients.
   real(R8P),    intent(in)    :: dt                                 !< Time step.
   integer(I4P), intent(in)    :: s                                  !< Stage to initialize.
   real(R8P),    intent(in)    ::   q(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative field.
   real(R8P),    intent(inout) :: q_s(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:) !< RK stage.
   integer(I4P)                :: ss                                 !< Counter.

   q_s(1:nv,1:ni,1:nj,1:nk,1:blocks_number,s) = q(1:nv,1:ni,1:nj,1:nk,1:blocks_number)
   do ss=1, s - 1
      q_s(1:nv,1:ni,1:nj,1:nk,1:blocks_number,s) = q_s(1:nv,1:ni,1:nj,1:nk,1:blocks_number,s ) + &
                                                  (q_s(1:nv,1:ni,1:nj,1:nk,1:blocks_number,ss) * (dt * alph(s, ss)))
   enddo
   endsubroutine compute_rk_stage

   ! IB
   subroutine set_bc_ib(ni, nj, nk, ngc, nv, blocks_number, q, phi, bcs_type)
   !< Set BC on IB cells.
   integer(I4P), intent(in)    :: ni                              !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                              !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                              !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                             !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                              !< Number of conservative varibales.
   integer(I4P), intent(in)    :: blocks_number                   !< Number of blocks.
   integer(I4P), intent(in)    :: bcs_type                        !< Immersed boundary type.
   real(R8P),    intent(in)    :: phi(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance field.
   real(R8P),    intent(inout) ::   q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.
   integer(I4P)                :: i, j, k, b, v                   !< Counter.
   real(R8P)                   :: n_phi_x, n_phi_y, n_phi_z       !< Distance function normals.
   real(R8P)                   :: n_phi_mod, un_mod               !< Distance abs normal and normal velocity.

   select case(bcs_type)
   case(BCS_VISCOUS)
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               do b=1, blocks_number
                  if (phi(1,i,j,k,b) >= 0) then
                     q(2,i,j,k,b) = - q(2,i,j,k,b)
                     q(3,i,j,k,b) = - q(3,i,j,k,b)
                     q(4,i,j,k,b) = - q(4,i,j,k,b)
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
                  if (phi(1,i,j,k,b) >= 0) then
                     n_phi_x = phi(1,i+1,j,k,b) - phi(1,i-1,j,k,b)
                     n_phi_y = phi(1,i,j+1,k,b) - phi(1,i,j-1,k,b)
                     n_phi_z = phi(1,i,j,k+1,b) - phi(1,i,j,k-1,b)
                     n_phi_mod = sqrt(n_phi_x**2 + n_phi_y**2 + n_phi_z**2)
                     n_phi_x = n_phi_x/n_phi_mod
                     n_phi_y = n_phi_y/n_phi_mod
                     n_phi_z = n_phi_z/n_phi_mod
                     un_mod = q(2,i,j,k,b)*n_phi_x + q(3,i,j,k,b)*n_phi_y + q(4,i,j,k,b)*n_phi_z

                     q(2,i,j,k,b) = q(2,i,j,k,b) - 2*un_mod*n_phi_x
                     q(3,i,j,k,b) = q(3,i,j,k,b) - 2*un_mod*n_phi_y
                     q(4,i,j,k,b) = q(4,i,j,k,b) - 2*un_mod*n_phi_z
                  endif
               enddo
            enddo
         enddo
      enddo
   endselect
   endsubroutine set_bc_ib
endmodule adam_equation_euler_cpu_object
