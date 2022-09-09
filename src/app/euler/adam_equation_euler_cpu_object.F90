!< ADAM, Euler equations system class definition and solver, CPU backend.
module adam_equation_euler_cpu_object
!< ADAM, Euler equations system class definition and solver, CPU backend.

use adam_adam_object,              only : adam_object
use adam_amr_cpu_object,           only : amr_cpu_object, amr_marker_cpu_object, AMR_GEO, AMR_GRAD
use adam_base_cpu_object,          only : base_cpu_object
use adam_bc_cpu_object,            only : bc_cpu_object, BC_EXTRAPOLATION, BC_INFLOW
use adam_eos_ic_cpu_object,        only : eos_ic_cpu_object
use adam_field_object,             only : field_object
use adam_euler_physics_cpu_object, only : euler_physics_cpu_object, IR, IU, IV, IW, IG, IP
use adam_grid_object,              only : grid_object
use adam_ib_cpu_object,            only : ib_cpu_object, BCS_VISCOUS, BCS_EULER
use adam_ic_cpu_object,            only : ic_cpu_object
use adam_mpih_object,              only : mpih_object
use adam_rk_cpu_object,            only : rk_cpu_object
use adam_slices_cpu_object,        only : slices_cpu_object
use adam_tree_object,              only : tree_object
use adam_weno_cpu_object,          only : weno_cpu_object
use adam_memory_cpu_lib
use adam_parameters
use finer
use penf
use mpi

implicit none
private
public :: equation_euler_cpu_object

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
   !< q_aux(1):    rho(1)          partial density of 1st specie
   !< q_aux(2):    rho(2)          partial density of 2nd specie
   !< ...
   !< q_aux(ns):   rho(ns)         partial density of last specie
   !< q_aux(ns+1): rho=sum(rho(s)) density
   !< q_aux(ns+2): u               velocity x
   !< q_aux(ns+3): v               velocity y
   !< q_aux(ns+4): w               velocity z
   !< q_aux(ns+5): g               specific heats ratio
   !< q_aux(ns+6): p               pressure
   !<```
   ! Main objects.
   type(mpih_object)              :: mpih            !< MPI handler.
   type(adam_object)              :: adam            !< ADAM.
   type(tree_object),  pointer    :: tree=>null()    !< The tree.
   type(field_object), pointer    :: field=>null()   !< The field.
   type(grid_object),  pointer    :: grid=>null()    !< The grid.
   type(file_ini)                 :: file_parameters !< Input simulation parameters file handler.
   type(base_cpu_object)          :: base_cpu        !< The base CPU handler.
   type(euler_physics_cpu_object) :: physics         !< Boundary conditions handler.
   type(bc_cpu_object)            :: bc              !< Boundary conditions handler.
   type(ic_cpu_object)            :: ic              !< Initial conditions handler.
   type(rk_cpu_object)            :: rk              !< Runge Kutta solver.
   type(weno_cpu_object)          :: weno            !< WENO Kutta solver.
   type(ib_cpu_object)            :: ib              !< Immersed Boundary handler.
   type(amr_cpu_object)           :: amr             !< AMR markers handler.
   type(slices_cpu_object)        :: slices          !< Slices handler.
   ! Other simulation data.
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
   real(R8P), allocatable :: q_aux(:,:,:,:,:) !< Auxiliary cell centered variables.
   real(R8P), allocatable ::  q_ib(:,:,:,:,:) !< Field cell with boundary set on immersed bodies.
   real(R8P), allocatable ::    rq(:,:,:,:,:) !< Field residuals.
   contains
      ! public methods
      procedure, pass(self) :: amr_update              !< Do AMR update.
      procedure, pass(self) :: compute_aux             !< Compute auxiliary variables.
      procedure, pass(self) :: compute_dt              !< Compute time step.
      procedure, pass(self) :: description             !< Return pretty-printed object description.
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
      procedure, pass(self), private :: apply_immersed_boundary   !< Apply immersed boundary to q.
      procedure, pass(self), private :: compute_residuals         !< Compute residuals of equation.
      procedure, pass(self), private :: compute_fluxes_convective !< Compute convective fluxes on a coordinate direction.
endtype equation_euler_cpu_object

contains
   ! public methods
   subroutine amr_update(self)
   !< Do AMR update.
   class(equation_euler_cpu_object), intent(inout) :: self                !< The equation.
   logical                                         :: is_grid_changed     !< Flag to check grid changes for each marker.
   logical                                         :: is_grid_changed_all !< Flag to check grid changes for each iter.
   type(amr_marker_cpu_object)                     :: amr_marker          !< Current amr marker.
   integer(I4P)                                    :: i, i_marker         !< Counter.

   if (self%amr%markers_number>0) then
      amr : do i=1, self%amr%iters
         is_grid_changed_all = .false.
         do i_marker=1, self%amr%markers_number
            amr_marker = self%amr%markers(i_marker)
            call self%update_ghost(q=self%field%q)
            select case(amr_marker%mode)
            case(AMR_GEO)
               call self%mark_by_geo(delta_fine=amr_marker%delta_fine, delta_coarse=amr_marker%delta_coarse)
            case(AMR_GRAD)
               select case(amr_marker%field)
               case(1)
                  call self%mark_by_grad_var(grad_tol=amr_marker%tol, delta_fine=amr_marker%delta_fine, &
                                             delta_coarse=amr_marker%delta_coarse, q=self%field%q, ivar=amr_marker%ivar)
               case(2)
                  call self%compute_aux(q=self%field%q, q_aux=self%q_aux)
                  call self%mark_by_grad_var(grad_tol=amr_marker%tol, delta_fine=amr_marker%delta_fine, &
                                             delta_coarse=amr_marker%delta_coarse, q=self%q_aux, ivar=amr_marker%ivar)
               endselect
            endselect
            call self%adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false., is_grid_changed=is_grid_changed)
            call self%ib%update_phi
            is_grid_changed_all = is_grid_changed_all.or.is_grid_changed
         enddo
         if (.not.is_grid_changed_all) then
             print '(A)', self%mpih%myrankstr//'AMR Grid stabilized after : '//trim(str(i))//' AMR iterations'
             exit amr
          elseif (i==self%amr%iters) then
             print '(A)', self%mpih%myrankstr//'AMR Grid is NOT stabilized after : '//trim(str(i))//' AMR iterations'
         endif
      enddo amr
   endif
   endsubroutine amr_update

   subroutine compute_aux(self, q, q_aux)
   !< Compute auxiliary variables.
   class(equation_euler_cpu_object), intent(in)  :: self                 !< The equation.
   real(R8P),                        intent(in)  :: q(1:,         &
                                                      1-self%grid%ngc:,&
                                                      1-self%grid%ngc:,&
                                                      1-self%grid%ngc:,&
                                                      1:)                !< Conservative variables.
   real(R8P),                        intent(out) :: q_aux(1:,         &
                                                          1-self%grid%ngc:,&
                                                          1-self%grid%ngc:,&
                                                          1-self%grid%ngc:,&
                                                          1:)            !< Auxiliary variables.
   integer(I4P)                                  :: b, i, j, k           !< Counter.
   ! real(R8P)                                     :: c(1:self%physics%ns) !< Species concentration.

   associate(blocks_number=>self%field%blocks_number, q=>self%field%q, &
             ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, ngc=>self%grid%ngc, &
             ns=>self%physics%ns, eos=>self%physics%eos)
   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               q_aux(:,i,j,k,b) = self%physics%conservative2primitive(conservative=q(:,i,j,k,b))
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

   associate(blocks_number=>self%field%blocks_number, ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, &
             ngc=>self%grid%ngc, dxyz=>self%field%dxyz, dt=>self%dt, CFL=>self%CFL)
      call self%compute_aux(q=self%field%q, q_aux=self%q_aux)
      dt = huge(1._R8P)
      do b=1, blocks_number
         umax = 0._R8P
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  ss = a(p=self%q_aux(IP,i,j,k,b), r=self%q_aux(IR,i,j,k,b), g=self%q_aux(IG,i,j,k,b))
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

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(equation_euler_cpu_object), intent(in) :: self             !< IC.
   character(len=:), allocatable                :: desc             !< Description.
   character(len=1), parameter                  :: NL=new_line('a') !< New line character.
   integer(I4P)                                 :: r, v             !< Counter.

   desc =       self%mpih%myrankstr//'equation main data'//NL
   desc = desc//self%mpih%myrankstr//'  restart:            '//trim(str(self%restart           ))//NL
   desc = desc//self%mpih%myrankstr//'  restart basename:   '//trim(    self%restart_basename   )//NL
   desc = desc//self%mpih%myrankstr//'  restart save:       '//trim(str(self%restart_save      ))//NL
   desc = desc//self%mpih%myrankstr//'  it max:             '//trim(str(self%it_max            ))//NL
   desc = desc//self%mpih%myrankstr//'  it save:            '//trim(str(self%it_save           ))//NL
   desc = desc//self%mpih%myrankstr//'  time max:           '//trim(str(self%time_max          ))//NL
   desc = desc//self%mpih%myrankstr//'  time save:          '//trim(str(self%time_save         ))//NL
   desc = desc//self%mpih%myrankstr//'  output basename:    '//trim(    self%output_basename    )//NL
   desc = desc//self%mpih%myrankstr//'  CFL:                '//trim(str(self%CFL               ))//NL
   desc = desc//self%mpih%myrankstr//'  save memory status: '//trim(str(self%save_memory_status))//NL
   desc = desc//self%mpih%myrankstr//'  null X direction:   '//trim(str(self%null_xyz(1)       ))//NL
   desc = desc//self%mpih%myrankstr//'  null Y direction:   '//trim(str(self%null_xyz(2)       ))//NL
   desc = desc//self%mpih%myrankstr//'  null Z direction:   '//trim(str(self%null_xyz(3)       ))
   endfunction description

   subroutine initialize(self, filename)
   !< Initialize the equation.
   class(equation_euler_cpu_object), intent(inout) :: self         !< The equation.
   character(*),                     intent(in)    :: filename     !< File name of input simulation parameters..
   integer(I8P)                                    :: nodes_number !< Allocated nodes on tree.
   integer(I4P)                                    :: nb           !< Number of allocated blocks.

   call self%mpih%initialize(do_mpi_init=.true.)
   print '(A)', self%mpih%myrankstr//'equation_euler_cpu%initialize start'

   call self%base_cpu%initialize_cpu

   call self%file_parameters%initialize(filename=trim(filename))
   call self%file_parameters%load

   call self%physics%initialize(file_parameters=self%file_parameters)

   call self%adam%grid%initialize(file_parameters=self%file_parameters, verbose=.true.)

   call self%adam%compute_blocks_number(memory_avail=self%base_cpu%memory_avail,  &
                                        fields_number=80,                         &
                                        nb=nb,                                    &
                                        nodes_number=nodes_number)

   call self%adam%initialize(nb=nb, file_parameters=self%file_parameters, &
                             do_tree_init=.true.,                         &
                             do_field_init=.true.,                        &
                             nv=self%physics%nv, nodes_number=nodes_number)

   call self%base_cpu%initialize(tree=self%adam%tree, field=self%adam%field)

   call self%bc%initialize(grid=self%adam%grid, file_parameters=self%file_parameters)

   call self%ic%initialize(file_parameters=self%file_parameters)

   call self%ib%initialize(grid=self%adam%grid, field=self%adam%field, file_parameters=self%file_parameters)

   call self%rk%initialize(file_parameters=self%file_parameters, grid=self%adam%grid, field=self%adam%field)

   call self%weno%initialize(file_parameters=self%file_parameters)

   call self%slices%initialize(file_parameters=self%file_parameters)

   call self%amr%initialize(file_parameters=self%file_parameters)

   call self%adam%refine_uniform(refinement_levels=self%adam%tree%iu_ref_levels, do_blocks_reorder=.false.)

   call self%adam%prune(ijkl_prune=self%adam%tree%ijkl_prune, do_blocks_reorder=.false.)

   call load_simulation_from_ini_file

   call link_objects_data(grid=self%adam%grid, tree=self%adam%tree, field=self%adam%field)

   ! allocate large arrays
   associate(nv=>self%physics%nv, nv_aux=>self%physics%nv_aux, &
             ngc=>self%grid%ngc, ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, &
             nb=>self%field%nb, nrk=>self%rk%nrk)
   allocate(self%q_aux(1:nv_aux, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
   allocate(self%q_ib( 1:nv,     1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
   allocate(self%rq(   1:nv,     1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
   endassociate
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'equation_euler_cpu%initialize finish'
   contains
      subroutine link_objects_data(grid, tree, field)
      !< Forward main ADAM data to equation for easy handling.
      type(grid_object),  intent(in), target :: grid  !< The grid.
      type(tree_object),  intent(in), target :: tree  !< The tree.
      type(field_object), intent(in), target :: field !< The field.

      self%grid  => grid
      self%tree  => tree
      self%field => field
      endsubroutine link_objects_data

      subroutine load_simulation_from_ini_file
      !< Load (other) simulation configs from input file.

      call self%file_parameters%get(section_name="simulation", option_name="restart",            val=self%restart           )
      call self%file_parameters%get(section_name="simulation", option_name="restart_basename",   val=self%restart_basename  )
      call self%file_parameters%get(section_name="simulation", option_name="restart_save",       val=self%restart_save      )
      call self%file_parameters%get(section_name="simulation", option_name="it_max",             val=self%it_max            )
      call self%file_parameters%get(section_name="simulation", option_name="it_save",            val=self%it_save           )
      call self%file_parameters%get(section_name="simulation", option_name="time_max",           val=self%time_max          )
      call self%file_parameters%get(section_name="simulation", option_name="time_save",          val=self%time_save         )
      call self%file_parameters%get(section_name="simulation", option_name="output_basename",    val=self%output_basename   )
      call self%file_parameters%get(section_name="simulation", option_name='CFL',                val=self%CFL               )
      call self%file_parameters%get(section_name="simulation", option_name='save_memory_status', val=self%save_memory_status)
      call self%file_parameters%get(section_name="simulation", option_name='null_x',             val=self%null_xyz(1)       )
      call self%file_parameters%get(section_name="simulation", option_name='null_y',             val=self%null_xyz(2)       )
      call self%file_parameters%get(section_name="simulation", option_name='null_z',             val=self%null_xyz(3)       )
      endsubroutine load_simulation_from_ini_file
   endsubroutine initialize

   subroutine integrate(self, do_ghost_syncro, residual)
   !< Runge Kutta integration of field.
   class(equation_euler_cpu_object), intent(inout)         :: self             !< The equation.
   logical,                          intent(in),  optional :: do_ghost_syncro  !< Flag to do syncrous ghost update.
   real(R8P),                        intent(out), optional :: residual         !< Global residual.
   logical                                                 :: do_ghost_syncro_ !< Flag to do syncrous ghost update, local var.
   integer(I4P)                                            :: s                !< Counter.

   do_ghost_syncro_ = .true. ; if (present(do_ghost_syncro)) do_ghost_syncro_ = do_ghost_syncro
   call self%update_ghost(q=self%field%q)
   associate(ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, ngc=>self%grid%ngc, blocks_number=>self%field%blocks_number, &
             nv=>self%physics%nv, dt=>self%dt)
   select case(self%rk%scheme(1:5))
   case('rk-ls') ! low storage schemes
      call self%rk%initialize_stage(ngc=self%grid%ngc, q=self%field%q)
      do s=1, self%rk%nrk
         call self%apply_immersed_boundary(q=self%field%q, q_ib=self%q_ib)
         call self%mpih%barrier
         call self%update_ghost(q=self%q_ib)
         call self%compute_aux(q=self%q_ib, q_aux=self%q_aux)
         call self%compute_residuals(q_aux=self%q_aux, rq=self%rq)
         call self%rk%compute_stage(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, dt=dt, s=s, &
                                    phi=self%ib%phi, rq=self%rq, q=self%field%q)
      enddo
   ! stop
   case('rk-ns') ! normal storage schemes
      call self%rk%initialize_stage(ngc=self%grid%ngc, q=self%field%q)
      do s=1, self%rk%nrk
         ! call self%apply_immersed_boundary(q=self%rk%q_s(:,:,:,:,:,s), q_ib=self%q_ib)
         call self%mpih%barrier
         call self%rk%compute_stage(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, dt=dt, s=s, &
                                    phi=self%ib%phi, rq=self%rq, q=self%rk%q_s(:,:,:,:,:,s))
         call self%compute_aux(q=self%rk%q_s(:,:,:,:,:,s), q_aux=self%q_aux)
         call self%compute_residuals(q_aux=self%q_aux, rq=self%rk%q_s(:,:,:,:,:,s))
      enddo
      call self%rk%sum_stages(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, dt=dt, q=self%field%q)
   endselect
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
   if(do_init_) self%field%refinements_needed = [(TO_BE_DEREFINED,b=1,self%field%blocks_number)]
   associate (ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, ngc=>self%grid%ngc, &
              blocks_number=>self%field%blocks_number, dxyz=>self%field%dxyz, phi=>self%ib%phi)
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
   !< @note Field q to which apply gradient musht have ghost cells updated.
   class(equation_euler_cpu_object), intent(inout)        :: self           !< The equation.
   real(R8P),                        intent(in)           :: grad_tol       !< Gradiend tolerance value.
   real(R8P),                        intent(in)           :: delta_fine     !< Maximum cell delta in fine grids.
   real(R8P),                        intent(in)           :: delta_coarse   !< Minimum cell delta in coarse grids.
   real(R8P),                        intent(in)           :: q(1:,          &
                                                               1-self%grid%ngc:, &
                                                               1-self%grid%ngc:, &
                                                               1-self%grid%ngc:, &
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
   if (do_init_) self%field%refinements_needed = [(TO_BE_DEREFINED,b=1,self%field%blocks_number)]
   associate (ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, ngc=>self%grid%ngc, blocks_number=>self%field%blocks_number, &
              dxyz=>self%field%dxyz)
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
      print '(A)', r//'t:                     '//trim(str(self%it,.true.))
      print '(A)', r//'blocks number:         '//trim(str(self%adam%tree%nodes_number, .true.))
      print '(A)', r//'process blocks number: '//trim(str(self%adam%field%blocks_number, .true.))
      print '(A)', r//'time step:             '//trim(str(self%dt, .true.))
      print '(A)', r//'time:                  '//trim(str(self%time, .true.))
   if (self%it_max <= 0) then
      print '(A)', r//'progress:              '//trim(str(int(self%time/self%time_max * 100), .true.))//'%'
   else
      print '(A)', r//'progress:              '//trim(str(int((self%it*1._R8P)/self%it_max * 100), .true.))//'%'
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
      call self%update_ghost(q=self%field%q)
      is_update_ghost_done = .true.
    endif
   call self%slices%save_mat(basename=self%output_basename, &
                             it=self%it,                    &
                             it_max=self%it_max,            &
                             time=self%time,                &
                             time_max=self%time_max,        &
                             adam=self%adam,                &
                             q=self%field%q,                &
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
         call self%update_ghost(q=self%field%q)
         call self%compute_aux(q=self%field%q, q_aux=self%q_aux)
         is_update_ghost_done = .true.
      endif
      call self%mpih%barrier(tictoc=.true.)
      print '(A)', self%mpih%myrankstr//'save HDF5 files t: '//trim(str(self%it,.true.))//', time: '//&
                   trim(str(self%time,.true.))
      output_basename_ = trim(self%output_basename)//'-'//trim(strz(self%it,9))
      if (present(output_basename)) output_basename_ = trim(output_basename)
      call self%adam%save_hdf5(basename=trim(output_basename_),          &
                               q=self%field%q,                           &
                               q_aux=self%q_aux,                         &
                               q_name=['rho','rhu','rhv','rhw','rhe'],   &
                               q_aux_name=['c','r','u','v','w','g','p'], &
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
         call self%update_ghost(q=self%field%q)
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
   class(equation_euler_cpu_object), intent(in)    :: self                      !< The equation.
   real(R8P),                        intent(inout) :: q(1:,         &
                                                        1-self%grid%ngc:,&
                                                        1-self%grid%ngc:,&
                                                        1-self%grid%ngc:,1:)    !< Conservative variables.
   integer(I4P)                                    :: crown                     !< Crown counter.
   integer(I4P)                                    :: fec                       !< Boundary fec (1 to 26).
   integer(I4P)                                    :: fec_1_6                   !< Boundary fec (1 to 6).
   integer(I4P)                                    :: c, i, j, k, b             !< Counter.
   integer(I4P)                                    :: idelta                    !< IJK delta step for extrapolation.
   integer(I4P)                                    :: jdelta                    !< IJK delta step for extrapolation.
   integer(I4P)                                    :: kdelta                    !< IJK delta step for extrapolation.
   integer(I4P)                                    :: bc_type                   !< Boundary condition type.
   integer(I4P)                                    :: s_bc                      !< Specie index of BC.
   real(R8P)                                       :: p_bc(self%physics%nv_aux) !< Primitive variables of BC.

   associate(ngc=>self%grid%ngc, q_bc=>self%bc%q)
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
                  ! initial specie index
                  s_bc = int(q_bc(6, fec_1_6))
                  ! build primitive variables array
                  p_bc = 0._R8P
                  p_bc(s_bc) = q_bc(1, fec_1_6)
                  p_bc(IR)   = q_bc(1, fec_1_6)
                  p_bc(IU)   = q_bc(2, fec_1_6)
                  p_bc(IV)   = q_bc(3, fec_1_6)
                  p_bc(IW)   = q_bc(4, fec_1_6)
                  p_bc(IP)   = q_bc(5, fec_1_6)
                  p_bc(IG)   = self%physics%eos(s_bc)%g
                  ! set BC in conservative variables
                  q(:,i,j,k,b) = self%physics%primitive2conservative(primitive=p_bc)
               endselect
            endif
         enddo
      enddo
   endif
   endassociate
   endsubroutine set_boundary_conditions

   subroutine set_initial_conditions(self)
   !< Set initial conditions of field.
   class(equation_euler_cpu_object), intent(inout) :: self                      !< The equation.
   integer(I4P)                                    :: b, i, j, k, r             !< Counter.
   integer(I4P)                                    :: s_ic                      !< Specie index of IC.
   real(R8P)                                       :: p_ic(self%physics%nv_aux) !< Primitive variables of IC.

   associate(blocks_number=>self%field%blocks_number, ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, &
             q=>self%field%q, q_ic=>self%ic%q, &
             x_cell=>self%field%x_cell, y_cell=>self%field%y_cell, z_cell=>self%field%z_cell)

   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               do r=1, self%ic%regions_number
                  if ((x_cell(i,b) > self%ic%emin(1,r).and.x_cell(i,b) <= self%ic%emax(1,r)).and. &
                      (y_cell(j,b) > self%ic%emin(2,r).and.y_cell(j,b) <= self%ic%emax(2,r)).and. &
                      (z_cell(k,b) > self%ic%emin(3,r).and.z_cell(k,b) <= self%ic%emax(3,r))) then
                     ! initial specie index
                     s_ic = int(q_ic(6, r))
                     ! build primitive variables array
                     p_ic = 0._R8P
                     p_ic(s_ic) = q_ic(1, r)
                     p_ic(IR)   = q_ic(1, r)
                     p_ic(IU)   = q_ic(2, r)
                     p_ic(IV)   = q_ic(3, r)
                     p_ic(IW)   = q_ic(4, r)
                     p_ic(IP)   = q_ic(5, r)
                     p_ic(IG)   = self%physics%eos(s_ic)%g
                     ! set IC in conservative variables
                     q(:,i,j,k,b) = self%physics%primitive2conservative(primitive=p_ic)
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine set_initial_conditions

   subroutine solve(self, filename)
   !< Solve Euler system.
   class(equation_euler_cpu_object), intent(inout) :: self             !< The equation.
   character(*),                     intent(in)    :: filename         !< File name of input simulation parameters..
   real(R8P)                                       :: timing(1:2)      !< Tic toc timing.
   real(R8P)                                       :: timing_step(1:2) !< Tic toc timing.

   call self%initialize(filename=filename)
   if (self%restart) then
      print '(A)', self%mpih%myrankstr//'restart simulation from "'//trim(self%restart_basename)//'" files'
      call self%load_restart_files(t=self%it, time=self%time)
      print '(A)', self%mpih%myrankstr//'restart [t, time]: '//trim(str(self%it))//', '//trim(str(self%time))
   else
      call self%set_initial_conditions
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

      if (mod(self%it,self%amr%frequency)==0) then
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
   call self%mpih%finalize
   endsubroutine solve

   subroutine update_ghost(self, q, step)
   !< Update ghost cells.
   !< If not specified all steps are perfermod, syncronous computation
   class(equation_euler_cpu_object), intent(inout)        :: self            !< The equation.
   real(R8P),                        intent(inout)        :: q(1:,         &
                                                               1-self%grid%ngc:,&
                                                               1-self%grid%ngc:,&
                                                               1-self%grid%ngc:,&
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
   subroutine apply_immersed_boundary(self, q, q_ib)
   !< Apply immersed boundary to q.
   class(equation_euler_cpu_object), intent(inout) :: self        !< The equation.
   real(R8P),                        intent(in   ) :: q(1:,              &
                                                        1-self%grid%ngc:,&
                                                        1-self%grid%ngc:,&
                                                        1-self%grid%ngc:,&
                                                        1:)       !< Conservative variables.
   real(R8P),                        intent(inout) :: q_ib(1:,              &
                                                           1-self%grid%ngc:,&
                                                           1-self%grid%ngc:,&
                                                           1-self%grid%ngc:,&
                                                           1:)    !< Conservative variables with applied IB.
   integer(I4P)                                    :: i_eikonal   !< Counter.
   integer(I4P), parameter                         :: n_eikonal=2 !< Counter.

   associate(ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, ngc=>self%grid%ngc, blocks_number=>self%field%blocks_number, &
             ns=>self%physics%ns)
   q_ib = q
   if (self%ib%solids_number > 0) then
      do i_eikonal=1, n_eikonal
         call MPI_BARRIER(MPI_COMM_WORLD, self%mpih%error)
         call self%ib%evolve_eikonal_q(q=q_ib)
         call self%update_ghost(q=q_ib)
      enddo
      call set_bc_ib(ni=ni, nj=nj, nk=nk, ngc=ngc, ns=ns, blocks_number=blocks_number, &
                     bcs_type=self%ib%bcs_type, phi=self%ib%phi, q=q_ib)
   endif
   endassociate
   endsubroutine apply_immersed_boundary

   subroutine compute_fluxes_convective(self, gc, n, ns, np,           &
                                        rhos, rho, un, ut1, ut2, g, p, &
                                        f_rho, f_rho_un, f_rho_ut1, f_rho_ut2, f_rho_e)
   !< Compute convective fluxes on coordinate direction.
   class(equation_euler_cpu_object), intent(in)    :: self                !< The equation.
   integer(I4P),                     intent(in)    :: gc                  !< Number of ghost cells used.
   integer(I4P),                     intent(in)    :: n                   !< Number of cells.
   integer(I4P),                     intent(in)    :: ns                  !< Number of species.
   integer(I4P),                     intent(in)    :: np                  !< Number of 1D primitive varibales.
   real(R8P),                        intent(in)    :: rhos(1:,1-gc:)      !< Partial densities       [1:ns,1-gc:n+gc].
   real(R8P),                        intent(in)    ::     rho(1-gc:)      !< Density                      [1-gc:n+gc].
   real(R8P),                        intent(in)    ::      un(1-gc:)      !< Normal velocity              [1-gc:n+gc].
   real(R8P),                        intent(in)    ::     ut1(1-gc:)      !< Tangential velocity 1        [1-gc:n+gc].
   real(R8P),                        intent(in)    ::     ut2(1-gc:)      !< Tangential velocity 2        [1-gc:n+gc].
   real(R8P),                        intent(in)    ::       g(1-gc:)      !< Specific heats ratio         [1-gc:n+gc].
   real(R8P),                        intent(in)    ::       p(1-gc:)      !< Pressure                     [1-gc:n+gc].
   real(R8P),                        intent(inout) :: f_rho(1:, 0:)       !< Flux of mass                 [0:n,1:ns].
   real(R8P),                        intent(inout) ::  f_rho_un(0:)       !< Flux normal momentums        [0:n].
   real(R8P),                        intent(inout) :: f_rho_ut1(0:)       !< Flux of tangential1 momentum [0:n].
   real(R8P),                        intent(inout) :: f_rho_ut2(0:)       !< Flux of tangential2 momentum [0:n].
   real(R8P),                        intent(inout) ::   f_rho_e(0:)       !< Flux energy                  [0:n].
   real(R8P)                                       ::   fluxes(3)         !< 1D fluxes.
   real(R8P)                                       :: q (1:np, 1-gc:n+gc) !< 1D primitive variables.
   real(R8P)                                       :: qr(1:np,1:2,0:n+1)  !< Reconstructed 1D primitive variables.
   integer(I4P)                                    :: i                   !< Counter.

   ! assembly 1D primitive variables array
   do i=1-gc, n+gc
      q(:, i) = [rhos(:,i), un(i), p(i), rho(i), g(i)]
   enddo
   call reconstruct_interfaces
   do i=0, n
      ! computing normal fluxes solving Riemann problem
      call solve_riemann(r1=qr(ns+3,2,i  ), u1=qr(ns+1,2,i  ), p1=qr(ns+2,2,i  ), g1=qr(ns+4,2,i  ), &
                         r4=qr(ns+3,1,i+1), u4=qr(ns+1,1,i+1), p4=qr(ns+2,1,i+1), g4=qr(ns+4,1,i+1), F=fluxes)
      if (fluxes(1)>0._R8P) then
           f_rho(:,i) = fluxes(1) * qr(1:ns,2,i) / qr(ns+3,2,i)
          f_rho_un(i) = fluxes(2)
         f_rho_ut1(i) = fluxes(1) * ut1(i)
         f_rho_ut2(i) = fluxes(1) * ut2(i)
           f_rho_e(i) = fluxes(3) + 0.5_R8P * fluxes(1) * (ut1(i)**2 + ut2(i)**2)
      else
           f_rho(:,i) = fluxes(1) * qr(1:ns,1,i+1) / qr(ns+3,1,i+1)
          f_rho_un(i) = fluxes(2)
         f_rho_ut1(i) = fluxes(1) * ut1(i+1)
         f_rho_ut2(i) = fluxes(1) * ut2(i+1)
           f_rho_e(i) = fluxes(3) + 0.5_R8P * fluxes(1) * (ut1(i+1)**2 + ut2(i+1)**2)
      endif
   enddo
   contains
      subroutine reconstruct_interfaces
      !< The reconstruction is done in pseudo characteristic variables.
      real(R8P)    :: qm(1:np,1:2)             !< Mean primitive variables.
      real(R8P)    :: Lqm(1:ns+2,1:ns+2,1:2)   !< Left eigenvalues matrix of mean primitive variables.
      real(R8P)    :: Rqm(1:ns+2,1:ns+2,1:2)   !< Right eigenvalues matrix of mean primitive variables.
      real(R8P)    :: c(1:ns+2,1:2,1-gc:-1+gc) !< Pseudo characteristic variables.
      real(R8P)    :: cr(1:ns+2,1:2)           !< Pseudo characteristic variables reconstructed.
      integer(I4P) :: i, j, f, v               !< Counter.

      select case(self%weno%weno_s)
      case(1_I4P)
         do i=0, n+1
            qr(:,1,i) = q(:,i)
            qr(:,2,i) = qr(:,1,i)
         enddo
      case(2_I4P,3_I4P,4_I4P)
      ! compute WENO reconstruction
         do i=0, n+1
            ! compute pseudo characteristic variables
            do f=1, 2
               if (i==0  .and.f==1) cycle
               if (i==n+1.and.f==2) cycle
               qm(:,f) = 0.5_R8P * (q(:,i+f-2) + q(:,i+f-1))
            enddo
            do f=1, 2
               if (i==0  .and.f==1) cycle
               if (i==n+1.and.f==2) cycle
               Lqm(:, :, f) =  left_eigenvectors(ns=ns, np=np, q=qm(:,f))
               Rqm(:, :, f) = right_eigenvectors(ns=ns, np=np, q=qm(:,f))
            enddo
            do j=i+1-gc, i-1+gc
               do f=1, 2
                  if (i==0  .and.f==1) cycle
                  if (i==n+1.and.f==2) cycle
                  do v=1, ns+2
                     c(v,f,j-i) = dot_product(Lqm(v,1:ns+2,f), q(1:ns+2,j))
                  enddo
               enddo
            enddo

            ! compute WENO reconstruction of pseudo charteristic variables
            do v=1, ns+2
               cr(v,:) = self%weno%reconstructed(s=self%weno%weno_s, v=c(v,:,1-self%weno%weno_s:-1+self%weno%weno_s))
            enddo

            ! trasform back reconstructed pseudo charteristic variables to primitive ones
            do f=1, 2
               if (i==0  .and.f==1) cycle
               if (i==n+1.and.f==2) cycle
               do v=1, ns+2
                  qr(v,f,i) = dot_product(Rqm(v,1:ns+2,f), cr(1:ns+2,f))
               enddo
               qr(ns+3,f,i) = sum(qr(1:ns, f, i))
               qr(ns+4,f,i) = dot_product(qr(1:ns,f,i) / qr(ns+3,f,i), self%physics%eos%cp) / &
                              dot_product(qr(1:ns,f,i) / qr(ns+3,f,i), self%physics%eos%cv)
            enddo
         enddo
      endselect
      endsubroutine reconstruct_interfaces
   endsubroutine compute_fluxes_convective

   subroutine compute_residuals(self, q_aux, rq)
   !< Compute residuals of equation.
   class(equation_euler_cpu_object), intent(inout) :: self                    !< The equation.
   real(R8P),                        intent(in)    :: q_aux(1:,         &
                                                            1-self%grid%ngc:,&
                                                            1-self%grid%ngc:,&
                                                            1-self%grid%ngc:,&
                                                            1:)               !< Auxiliary variables.
   real(R8P),                        intent(inout) ::    rq(1:,         &
                                                            1-self%grid%ngc:,&
                                                            1-self%grid%ngc:,&
                                                            1-self%grid%ngc:,&
                                                            1:)               !< Residuals.
   real(R8P)                                       :: fluxes_x(1:self%physics%nv,&
                                                               0:self%grid%ni,&
                                                               1:self%grid%nj,&
                                                               1:self%grid%nk)!< Convective fluxes in x direction.
   real(R8P)                                       :: fluxes_y(1:self%physics%nv,&
                                                               1:self%grid%ni,&
                                                               0:self%grid%nj,&
                                                               1:self%grid%nk)!< Convective fluxes in y direction.
   real(R8P)                                       :: fluxes_z(1:self%physics%nv,&
                                                               1:self%grid%ni,&
                                                               1:self%grid%nj,&
                                                               0:self%grid%nk)!< Convective fluxes in z direction.
   integer(I4P)                                    :: b, i, j, k              !< Counter.
   integer(I4P)                                    :: ii(2), jj(2), kk(2)     !< Loops bounds.

   associate(ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk, &
             ngc=>self%field%grid%ngc, blocks_number=>self%field%blocks_number,      &
             dxyz=>self%field%dxyz, ns=>self%physics%ns, np=>self%physics%np)
   ii = [1,ni] ; if (self%null_xyz(1)) ii = [1,1]
   jj = [1,nj] ; if (self%null_xyz(2)) jj = [1,1]
   kk = [1,nk] ; if (self%null_xyz(3)) kk = [1,1]
   rq = 0._R8P
   do b=1, blocks_number
      fluxes_x = 0._R8P
      fluxes_y = 0._R8P
      fluxes_z = 0._R8P
      if (.not.self%null_xyz(1)) then ! convective fluxes along x direction
         do k=kk(1), kk(2)
            do j=jj(1), jj(2)
               call self%compute_fluxes_convective(gc=ngc, n=ni, ns=ns, np=np,            &
                                                   rhos      =    q_aux(1:ns, :,  j,k,b), &
                                                   rho       =    q_aux(IR  , :,  j,k,b), &
                                                   un        =    q_aux(IU  , :,  j,k,b), & ! u
                                                   ut1       =    q_aux(IV  , :,  j,k,b), & ! v
                                                   ut2       =    q_aux(IW  , :,  j,k,b), & ! w
                                                   g         =    q_aux(IG  , :,  j,k,b), &
                                                   p         =    q_aux(IP  , :,  j,k,b), &
                                                   f_rho     = fluxes_x(1:ns,0:ni,j,k),   &
                                                   f_rho_un  = fluxes_x(ns+1,0:ni,j,k),   & ! u
                                                   f_rho_ut1 = fluxes_x(ns+2,0:ni,j,k),   & ! v
                                                   f_rho_ut2 = fluxes_x(ns+3,0:ni,j,k),   & ! w
                                                   f_rho_e   = fluxes_x(ns+4,0:ni,j,k))
            enddo
         enddo
      endif
      if (.not.self%null_xyz(2)) then ! convective fluxes along y direction
         do k=kk(1), kk(2)
            do i=ii(1), ii(2)
               call self%compute_fluxes_convective(gc=ngc, n=nj, ns=ns, np=np,            &
                                                   rhos      =    q_aux(1:ns,i, :,  k,b), &
                                                   rho       =    q_aux(IR  ,i, :,  k,b), &
                                                   un        =    q_aux(IV  ,i, :,  k,b), & ! v
                                                   ut1       =    q_aux(IU  ,i, :,  k,b), & ! u
                                                   ut2       =    q_aux(IW  ,i, :,  k,b), & ! w
                                                   g         =    q_aux(IG  ,i, :,  k,b), &
                                                   p         =    q_aux(IP  ,i, :,  k,b), &
                                                   f_rho     = fluxes_y(1:ns,i,0:nj,k),   &
                                                   f_rho_un  = fluxes_y(ns+2,i,0:nj,k),   & ! v
                                                   f_rho_ut1 = fluxes_y(ns+1,i,0:nj,k),   & ! u
                                                   f_rho_ut2 = fluxes_y(ns+3,i,0:nj,k),   & ! w
                                                   f_rho_e   = fluxes_y(ns+4,i,0:nj,k))
            enddo
         enddo
      endif
      if (.not.self%null_xyz(3)) then ! convective fluxes along z direction
         do j=jj(1), jj(2)
            do i=ii(1), ii(2)
               call self%compute_fluxes_convective(gc=ngc, n=nk, ns=ns, np=np,            &
                                                   rhos      =    q_aux(1:ns,i,j, :,  b), &
                                                   rho       =    q_aux(IR  ,i,j, :,  b), &
                                                   un        =    q_aux(IW  ,i,j, :,  b), & ! w
                                                   ut1       =    q_aux(IU  ,i,j, :,  b), & ! u
                                                   ut2       =    q_aux(IV  ,i,j, :,  b), & ! v
                                                   g         =    q_aux(IG  ,i,j, :,  b), &
                                                   p         =    q_aux(IP  ,i,j, :,  b), &
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
               rq(:,i,j,k,b) = (fluxes_x(:,i-1,j  ,k  ) - fluxes_x(:,i,j,k)) / dxyz(1,b) + &
                               (fluxes_y(:,i  ,j-1,k  ) - fluxes_y(:,i,j,k)) / dxyz(2,b) + &
                               (fluxes_z(:,i  ,j  ,k-1) - fluxes_z(:,i,j,k)) / dxyz(3,b)
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine compute_residuals

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
      real(R8P), parameter   :: toll=1e-10_R8P !< Tollerance.

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

   pure function left_eigenvectors(ns, np, q) result(eig)
   !< Return the left eigenvectors matrix `L` as `dF/dP = A = R ^ L` `P`` being the primitive variables and `F` the fluxes.
   integer(I4P), intent(in) :: ns                  !< Number of species.
   integer(I4P), intent(in) :: np                  !< Number of 1D primitive varibales.
   real(R8P),    intent(in) :: q(1:np)             !< Primitive variables.
   real(R8P)                :: eig(1:ns+2, 1:ns+2) !< Eigenvectors.
   real(R8P)                :: gp                  !< `g*p`.
   real(R8P)                :: gp_a                !< `g*p/a`.
   integer(I4P)             :: s                   !< Counter.

   gp = q(ns+4) * q(ns+2)
   gp_a = gp / a(p=q(ns+2), r=q(ns+3), g=q(ns+4))
   eig = 0._R8P

      eig(1,    ns+1) = -gp_a      ; eig(1,    ns+2) =  1._R8P
   do s=2, ns+1
      eig(s,    s-1 ) =  gp/q(s-1) ; eig(s,    ns+2) = -1._R8P
   enddo
      eig(ns+2, ns+1) =  gp_a      ; eig(ns+2, ns+2) =  1._R8P
   endfunction left_eigenvectors

   pure function right_eigenvectors(ns, np, q) result(eig)
   !< Return the right eigenvectors matrix `R` as `dF/dP = A = R ^ L` `P`` being the primitive variables and `F` the fluxes.
   integer(I4P), intent(in) :: ns                  !< Number of species.
   integer(I4P), intent(in) :: np                  !< Number of 1D primitive varibales.
   real(R8P),    intent(in) :: q(1:np)             !< Primitive variables.
   real(R8P)                :: eig(1:ns+2, 1:ns+2) !< Eigenvectors.
   real(R8P)                :: gp                  !< `g*p`.
   real(R8P)                :: gp_inv              !< `1/(g*p)`.
   integer(I4P)             :: s                   !< Counter.

   gp = q(ns+4) * q(ns+2)
   gp_inv = 1._R8P / gp
   eig = 0._R8P

   do s=1, ns
     eig(s,   1) =  0.5_R8P * q(s) * gp_inv                             ; eig(s,s+1) = q(s) * gp_inv ; eig(s,   ns+2) =  eig(s,   1)
   enddo
     eig(ns+1,1) = -0.5_R8P * a(p=q(ns+2),r=q(ns+3),g=q(ns+4)) * gp_inv ;                              eig(ns+1,ns+2) = -eig(ns+1,1)
     eig(ns+2,1) =  0.5_R8P                                             ;                              eig(ns+2,ns+2) =  eig(ns+2,1)
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

   entalpy = g * p / ((g - 1._R8P) * r) + 0.5_R8P * u*u
   endfunction H

   ! IB
   subroutine set_bc_ib(ni, nj, nk, ngc, ns, blocks_number, bcs_type, phi, q)
   !< Set BC on IB cells.
   integer(I4P), intent(in)    :: ni                              !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                              !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                              !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                             !< Ghost cells number.
   integer(I4P), intent(in)    :: ns                              !< Number of species.
   integer(I4P), intent(in)    :: blocks_number                   !< Number of blocks.
   integer(I4P), intent(in)    :: bcs_type(1:)                    !< Immersed boundary type.
   real(R8P),    intent(in)    :: phi(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance field.
   real(R8P),    intent(inout) ::   q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.
   integer(I4P)                :: i, j, k, b, v, s                !< Counter.
   real(R8P)                   :: n_phi_x, n_phi_y, n_phi_z       !< Distance function normals.
   real(R8P)                   :: n_phi_mod, un_mod               !< Distance abs normal and normal velocity.

   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               solids_loop : do s=1, size(bcs_type, dim=1)
                  if (phi(s,i,j,k,b) >= 0) then
                     select case(bcs_type(s))
                     case(BCS_VISCOUS)
                        q(ns+1,i,j,k,b) = - q(ns+1,i,j,k,b)
                        q(ns+2,i,j,k,b) = - q(ns+2,i,j,k,b)
                        q(ns+3,i,j,k,b) = - q(ns+3,i,j,k,b)
                     case(BCS_EULER)
                        n_phi_x = phi(s,i+1,j,k,b) - phi(s,i-1,j,k,b)
                        n_phi_y = phi(s,i,j+1,k,b) - phi(s,i,j-1,k,b)
                        n_phi_z = phi(s,i,j,k+1,b) - phi(s,i,j,k-1,b)
                        n_phi_mod = sqrt(n_phi_x**2 + n_phi_y**2 + n_phi_z**2)
                        n_phi_x = n_phi_x/n_phi_mod
                        n_phi_y = n_phi_y/n_phi_mod
                        n_phi_z = n_phi_z/n_phi_mod
                        un_mod = q(ns+1,i,j,k,b)*n_phi_x + q(ns+2,i,j,k,b)*n_phi_y + q(ns+3,i,j,k,b)*n_phi_z

                        q(ns+1,i,j,k,b) = q(ns+1,i,j,k,b) - 2 * un_mod*n_phi_x
                        q(ns+2,i,j,k,b) = q(ns+2,i,j,k,b) - 2 * un_mod*n_phi_y
                        q(ns+3,i,j,k,b) = q(ns+3,i,j,k,b) - 2 * un_mod*n_phi_z
                     endselect
                     exit solids_loop
                  endif
               enddo solids_loop
            enddo
         enddo
      enddo
   enddo
   endsubroutine set_bc_ib
endmodule adam_equation_euler_cpu_object
