!< ADAM, PRISM (Plasma Research usIng Simulation Methods) equations system class definition, CPU backend.
module adam_prism_cpu_object
!< ADAM, PRISM (Plasma Research usIng Simulation Methods) equations system class definition, CPU backend.

! ADAM classes, libraries, parameters
use :: adam_common_library
! PRISM modules
use :: adam_prism_common_library
! third party modules
use :: penf
use :: mpi

implicit none
private
public :: prism_cpu_object

! pointer (abstract) procedures
procedure(compute_convective_fluxes_interface), pointer :: compute_fluxes_maxwell=>null() !< Compute convective fluxes.

integer(I4P), parameter :: PML_FACE_X_M = 1_I4P
integer(I4P), parameter :: PML_FACE_X_P = 2_I4P
integer(I4P), parameter :: PML_FACE_Y_M = 3_I4P
integer(I4P), parameter :: PML_FACE_Y_P = 4_I4P
integer(I4P), parameter :: PML_FACE_Z_M = 5_I4P
integer(I4P), parameter :: PML_FACE_Z_P = 6_I4P
real(R8P),    parameter :: GRMS_3DB_RATIO = 10.0_R8P**(-3.0_R8P/20.0_R8P)

type, extends(prism_common_object) :: prism_cpu_object !commentate procedure AMR e IB
   !< Maxwell equations system class definition, CPU backend.
   real(R8P), allocatable :: flxyz_c(:,:,:,:,:,:,:) !< Fluxes at cell center with +/- decomposition for all directions.
   real(R8P), allocatable ::   flx_f(:,:,:,:,:    ) !< Fluxes along x at cell face.
   real(R8P), allocatable ::   fly_f(:,:,:,:,:    ) !< Fluxes along y at cell face.
   real(R8P), allocatable ::   flz_f(:,:,:,:,:    ) !< Fluxes along z at cell face.
   logical               :: fv_add_phi_damping = .false. !< Apply phi damping in FV residuals.
   logical               :: fv_add_psi_damping = .false. !< Apply psi damping in FV residuals.
   integer(I4P)          :: fv_ivar_phi        = 0_I4P   !< Phi slot index in FV residuals.
   integer(I4P)          :: fv_ivar_psi        = 0_I4P   !< Psi slot index in FV residuals.
   !< Pointer (abstract) TBP.
   procedure(compute_residuals_interface), pass(self),pointer :: compute_residuals=>null()!< Compute residuals, space operator.
   procedure(integrate_interface),         pass(self),pointer :: integrate        =>null()!< Integrate, time operator.
   contains
      ! AMR methods (amr_update + mark_by_geometry promoted to prism_common_object, issue #22 F0)
      procedure, pass(self) :: mark_by_j_vec_total_variation !< Override: the working host TV marker (common default stops).
      ! auxiliary methods
      procedure, pass(self) :: allocate_cpu     !< Allocate CPU data.
      procedure, pass(self) :: initialize_prism !< Initialize PRSIM equation.
      ! IO methods
      procedure, pass(self) :: save_residuals       !< Save residuals history.
      procedure, pass(self) :: save_simulation_data !< Save all simulation data.
      ! IC/BC/sources
      procedure, pass(self) :: apply_fWL_correction       !< Apply fWLayer correction (if present)
      procedure, pass(self) :: compute_coils_current      !< Compute current coils sources.
      procedure, pass(self) :: set_boundary_conditions    !< Set boundary conditions of equation.
      !procedure, pass(self) :: compute_residuals_BC
      !procedure, pass(self) :: update_q_BC
      procedure, pass(self) :: set_initial_conditions     !< Set initial conditions (and coils) of equation.
      procedure, pass(self) :: update_ghost               !< Update ghost cells and set boundary conditions.
      procedure, pass(self) :: compute_field_mean_value   !< Compute field mean value.
      ! forest orchestrator contract methods overridings
      procedure, pass(self) :: initialize_forest            !< Invoked by forest%initialize per realm at startup.
      procedure, pass(self) :: compute_local_dt_forest      !< Invoked by forest%compute_global_dt during the min reduction.
      procedure, pass(self) :: advance_one_step_forest      !< Invoked by forest%evolve_one_step per realm per timestep.
      procedure, pass(self) :: stages_per_step_forest       !< Number of integrator stages this realm exposes per step.
      procedure, pass(self) :: open_step_forest             !< Per-step prologue (multi-realm path).
      procedure, pass(self) :: begin_stage_forest           !< Begin one integrator stage on this realm (multi-realm path).
      procedure, pass(self) :: end_stage_forest             !< End the stage: residuals + assignment (multi-realm path).
      procedure, pass(self) :: close_step_forest            !< Per-step epilogue (multi-realm path).
      procedure, pass(self) :: post_step_forest             !< Invoked by forest%post_step per realm per timestep.
      procedure, pass(self) :: is_done_forest               !< Invoked by forest%is_done during the termination reduction.
      procedure, pass(self) :: finalize_forest              !< Invoked by forest%finalize per realm at shutdown.
      procedure, pass(self) :: fill_seam_from_peer_forest   !< Copy peer's interior into self's ghosts for peer slot p_idx.
      procedure, pass(self) :: apply_reflux_to_stage_forest !< Apply end-of-step Berger-Colella reflux to self's committed q.
      ! numerical methods
      procedure, pass(self) :: compute_dt                   !< Compute time step.
      procedure, pass(self) :: compute_energy               !< Compute energy.
      procedure, pass(self) :: compute_energy_error         !< Compute energy error.
      procedure, pass(self) :: compute_grms                 !< Compute Grms of the rotating magnetic-field amplitude.
      procedure, pass(self) :: impose_div_free              !< Impose divergence-free property.
      procedure, pass(self) :: simulate                     !< Perform the simulation.
endtype prism_cpu_object

interface
   subroutine compute_residuals_interface(self, q, dq, s, flux_register)
   !< Compute residuals of equation, space operator.
   !<
   !< Inter-realm seam ghost cells are filled by the forest BEFORE this
   !< routine fires (Phase 2 of the forest's substage loop, for β seams;
   !< Phase 5 of the previous step, for α seams), so the stencil reads
   !< valid halo data without any peer-realm access here.
   import :: prism_cpu_object, R8P, I4P, flux_register_object
   class(prism_cpu_object),     intent(inout)           :: self          !< The equation.
   real(R8P),                   intent(inout)           :: q(1:,         &
                                                             1-self%ngc:,&
                                                             1-self%ngc:,&
                                                             1-self%ngc:,&
                                                             1:)         !< Conservative variables.
   real(R8P),                   intent(inout)           :: dq(1:,         &
                                                              1-self%ngc:,&
                                                              1-self%ngc:,&
                                                              1-self%ngc:,&
                                                              1:)        !< Residuals.
   integer(I4P),                intent(in),    optional :: s             !< Stage counter.
   class(flux_register_object), intent(inout), optional :: flux_register !< Forest's flux register for FV reflux.
   endsubroutine compute_residuals_interface

   subroutine integrate_interface(self)
   !< Integrate equation, time operator.
   import :: prism_cpu_object, R8P
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   endsubroutine integrate_interface
endinterface

contains
   ! (restrict_fine_face_to_quadrant moved to adam_flux_register_object, issue #23 R3)

   ! AMR methods (amr_update + mark_by_geometry live on prism_common_object since issue #22 F0)
   subroutine mark_by_j_vec_total_variation(self, tv_tol, delta_type, delta_fine, delta_coarse, threshold, do_init)
   !< Mark blocks to be refined/derefined by the value of total variation of j_vec.
   class(prism_cpu_object), intent(inout)        :: self                     !< The equation.
   real(R8P),               intent(in)           :: tv_tol                   !< Total variation tolerance value.
   character(*),            intent(in)           :: delta_type               !< Delta criterion type.
   real(R8P),               intent(in)           :: delta_fine               !< Maximum cell delta in fine grids.
   real(R8P),               intent(in)           :: delta_coarse             !< Minimum cell delta in coarse grids.
   real(R8P),               intent(in), optional :: threshold                !< Threshold for sphere proximity.
   logical,                 intent(in), optional :: do_init                  !< Re-initialize refinements queries.
   logical                                       :: do_init_                 !< Re-initialize refinements queries, local var.
   real(R8P)                                     :: threshold_               !< Threshold for sphere proximity, local var.
   real(R8P)                                     :: max_cell_delta           !< Maximum cell delta.
   real(R8P)                                     :: max_total_variation      !< Total variation of j_vec, max on all coils.
   real(R8P)                                     :: total_variation          !< Total variation of j_vec.
   integer(I4P)                                  :: b,c                      !< Counter.
   real(R8P)                                     :: dc(1:self%blocks_number) !< Delta criterion.

   do_init_   = .true.  ; if (present(do_init  )) do_init_   = do_init
   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   if (do_init_) self%adam%field%refinements_needed = [(TO_BE_DEREFINED,b=1,self%blocks_number)]
   associate (ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
              blocks_number=>self%blocks_number, hs=>self%fdv_half_stencils(1), dxyz=>self%adam%field%dxyz)
      select case(delta_type)
      case(AMR_DELTA_T_X) ; dc(1:blocks_number) = dxyz(1,1:blocks_number)
      case(AMR_DELTA_T_Y) ; dc(1:blocks_number) = dxyz(2,1:blocks_number)
      case(AMR_DELTA_T_Z) ; dc(1:blocks_number) = dxyz(3,1:blocks_number)
      case(AMR_DELTA_T_MAX)
         do b=1, blocks_number
            dc(b) = maxval(dxyz(:,b))
         enddo
      endselect
      call self%update_ghost(q=self%q)
      call self%initialize_coils
      do b=1, blocks_number
         max_total_variation = -huge(1._R8P)
         do c=1, 1!self%coil%total_coils_number
            call self%compute_block_total_variation(hs=hs, dxyz=dxyz(:,b), ivar=1,            &
                                                    q=self%coil%j_vec(:,:,:,:,b,c),           &
                                                    tot_var_field=self%divergence(4,:,:,:,b), &
                                                    total_variation=total_variation)
            max_total_variation = max(max_total_variation,total_variation)
         enddo
         call mpih%print_message('Block '//trim(str(b))//': max_total_variation = '//trim(str(max_total_variation)))
         max_cell_delta = max_cell_delta_tv(tv=max_total_variation)
         if ((dc(b)) > max_cell_delta) then
            self%adam%field%refinements_needed(b) = TO_BE_REFINED
         elseif ((dc(b)) * threshold_ < max_cell_delta) then
            self%adam%field%refinements_needed(b) = max(self%adam%field%refinements_needed(b), TO_BE_DEREFINED)
         else
            self%adam%field%refinements_needed(b) = max(self%adam%field%refinements_needed(b), TO_NOT_TOUCH)
         endif
      enddo
   endassociate
   contains
      function max_cell_delta_tv(tv) result(delta)
      !< Return the maximum cell delta given a total variation tollerance.
      real(R8P), intent(in) :: tv    !< Total variation value.
      real(R8P)             :: delta !< Maximum cell delta admissible.

      if (tv > tv_tol) then
         delta = delta_fine
      else
         delta = delta_coarse
      endif
      endfunction max_cell_delta_tv
   endsubroutine mark_by_j_vec_total_variation

   ! auxiliary methods
   subroutine allocate_cpu(self)
   !< Allocate CPU data.
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   character(:), allocatable              :: msg_ !< Allocating message base.
   character(:), allocatable              :: msg  !< Allocating message.

   call mpih%print_message('prism_cpu_object%allocate_cpu start')
   msg_ = mpih%myrankstr//'prism_cpu_object%allocate_cpu '
   associate(nv=>self%nv, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, nb=>self%nb)
   msg = msg_//' flxyz_c '
   call allocate_variable(var=self%flxyz_c,ulb=reshape([1,nv, &
                                                        1,3,  & ! f at center, f positive at center, f negative at center
                                                        1,3,  & ! direction x, y, z
                                                        1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,nb],[2,7]),msg=msg)
   self%flxyz_c = 0._R8P
   msg = msg_//' flx_f '
   call allocate_variable(var=self%flx_f,ulb=reshape([1,nv,0,ni,1,nj,1,nk,1,nb],[2,5]),msg=msg)
   self%flx_f = 0._R8P
   msg = msg_//' fly_f '
   call allocate_variable(var=self%fly_f,ulb=reshape([1,nv,1,ni,0,nj,1,nk,1,nb],[2,5]),msg=msg)
   self%fly_f = 0._R8P
   msg = msg_//' flz_f '
   call allocate_variable(var=self%flz_f,ulb=reshape([1,nv,1,ni,1,nj,0,nk,1,nb],[2,5]),msg=msg)
   self%flz_f = 0._R8P
   endassociate
   call mpih%print_message('prism_cpu_object%allocate_cpu finish')
   endsubroutine allocate_cpu

   subroutine initialize_prism(self, filename, realms_number)
   !< Initialize PRSIM equation.
   class(prism_cpu_object), intent(inout)        :: self          !< The equation.
   character(*),            intent(in)           :: filename      !< Input file name.
   integer(I4P),            intent(in), optional :: realms_number !< Forest realm count; divides the per-process budget (default 1).
   real(R8P)                                     :: memory_avail_ !< Per-realm memory budget after the forest split.
   integer(I4P)                                  :: realms_number_!< Local realm count (>=1).

   realms_number_ = 1_I4P ; if (present(realms_number)) realms_number_ = max(realms_number, 1_I4P)
   call mpih%initialize(verbose=.true.)
   call mpih%print_message('prism_cpu_object%initialize start')
   memory_avail_ = mpih%memory_avail / real(realms_number_, R8P)
   call self%prism_common_object%initialize(filename=filename,memory_avail=memory_avail_,verbose=.true.)
   call self%allocate_cpu

   ! set pointer (abstract) TBP
   if (self%physics%physical_model == EM_PHYSICAL_MODEL .or. &
         self%physics%physical_model == ADIM_EM_PHYSICAL_MODEL) then
      select case(self%numerics%scheme_time)
      case(NUM_SCHEME_TIME_BLANES_MOAN)        ; self%integrate => integrate_blanesmoan
      case(NUM_SCHEME_TIME_CFM)                ; self%integrate => integrate_cfm
      case(NUM_SCHEME_TIME_LEAPFROG)           ; self%integrate => integrate_leapfrog
      case(NUM_SCHEME_TIME_RUNGE_KUTTA)
         select case(self%rk%scheme)
         case(RK_1, RK_2, RK_3)                ; self%integrate => integrate_rk_ls
         case(RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54) ; self%integrate => integrate_rk_ssp
         case(RK_YOSHIDA)                      ; self%integrate => integrate_rk_yoshida
         endselect
      endselect
   elseif (self%physics%physical_model == PIC_PHYSICAL_MODEL) then !Metterei qualche error stop sulle combinazioni non valide
      select case(self%numerics%scheme_time)
      case(NUM_SCHEME_TIME_LEAPFROG)
         select case(self%pic%scheme_time)
         case(NUM_SCHEME_TIME_PIC_LEAPFROG)
            self%integrate => integrate_leapfrog_pic
         case(NUM_SCHEME_TIME_PIC_RUNGE_KUTTA)
            !self%integrate =>
         endselect
      case(NUM_SCHEME_TIME_RUNGE_KUTTA)
         select case(self%pic%scheme_time)
         case(NUM_SCHEME_TIME_PIC_LEAPFROG)
            self%integrate => integrate_leapfrog_pic
         case(NUM_SCHEME_TIME_PIC_RUNGE_KUTTA)
         select case(self%rk_pic%scheme)
         case(RK_1, RK_2, RK_3)
            !self%integrate => integrate_rk_ls_pic
         case(RK_SSP_22, RK_SSP_33, RK_SSP_54)
            self%integrate => integrate_rk_ssp_pic
         case(RK_YOSHIDA)
            !self%integrate => integrate_rk_yoshida_pic
         endselect
         endselect
      endselect
   endif

   select case(self%numerics%scheme_space)
   case(NUM_SCHEME_SPACE_WENO)        ; self%compute_residuals => compute_residuals_weno
   case(NUM_SCHEME_SPACE_FD_CENTERED) ; self%compute_residuals => compute_residuals_fd_centered
   case(NUM_SCHEME_SPACE_FV_CENTERED) ; self%compute_residuals => compute_residuals_fv_centered
   endselect

   self%fv_add_phi_damping = .false.
   self%fv_add_psi_damping = .false.
   self%fv_ivar_phi        = 0_I4P
   self%fv_ivar_psi        = 0_I4P
   select case(self%physics%physical_model)
   case(ADIM_EM_PHYSICAL_MODEL)
      select case(self%numerics%div_corr_var)
      case(DIV_CORR_VAR_HYPER)
         if (self%numerics%constrained_transport_D .and. .not.self%numerics%constrained_transport_B) then
            compute_fluxes_maxwell  => compute_convective_fluxes_maxwell_adim_div_d
            self%fv_add_phi_damping = .true.
            self%fv_ivar_phi        = self%nv_c
         elseif (.not.self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
            compute_fluxes_maxwell  => compute_convective_fluxes_maxwell_adim_div_b
            self%fv_add_psi_damping = .true.
            self%fv_ivar_psi        = self%nv_c
         elseif (self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
            compute_fluxes_maxwell  => compute_convective_fluxes_maxwell_adim_div_d_b
            self%fv_add_phi_damping = .true.
            self%fv_add_psi_damping = .true.
            self%fv_ivar_phi        = self%nv_c - 1_I4P
            self%fv_ivar_psi        = self%nv_c
         else
            compute_fluxes_maxwell => compute_convective_fluxes_maxwell_adim
         endif
      case default
         compute_fluxes_maxwell => compute_convective_fluxes_maxwell_adim
      endselect
   case default
      select case(self%numerics%div_corr_var)
      case(DIV_CORR_VAR_HYPER)
         if (self%numerics%constrained_transport_D .and. .not.self%numerics%constrained_transport_B) then
            compute_fluxes_maxwell  => compute_convective_fluxes_maxwell_div_d
            self%fv_add_phi_damping = .true.
            self%fv_ivar_phi        = self%nv_c
         elseif (.not.self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
            compute_fluxes_maxwell  => compute_convective_fluxes_maxwell_div_b
            self%fv_add_psi_damping = .true.
            self%fv_ivar_psi        = self%nv_c
         elseif (self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
            compute_fluxes_maxwell  => compute_convective_fluxes_maxwell_div_d_b
            self%fv_add_phi_damping = .true.
            self%fv_add_psi_damping = .true.
            self%fv_ivar_phi        = self%nv_c - 1_I4P
            self%fv_ivar_psi        = self%nv_c
         else
            compute_fluxes_maxwell => compute_convective_fluxes_maxwell
         endif
      case default
         compute_fluxes_maxwell => compute_convective_fluxes_maxwell
      endselect
   endselect

   call check_pml_configuration()

   print '(A)', mpih%description()
   call mpih%print_message('prism_cpu_object%initialize finish')
   contains
      subroutine check_pml_configuration()
      !< The first CPU implementation supports only classic/Bermudez/CFS ADE-PML on the pure-Maxwell SSP FD path.
      logical :: is_supported_ssp

      if (.not. self%pml%enabled) return

      is_supported_ssp = trim(self%rk%scheme) == RK_SSP_11 .or. trim(self%rk%scheme) == RK_SSP_22 .or. &
                         trim(self%rk%scheme) == RK_SSP_33 .or. trim(self%rk%scheme) == RK_SSP_54

      if (trim(self%pml%pml_type) /= 'CLASSIC' .and. trim(self%pml%pml_type) /= 'BERMUDEZ' .and. &
          trim(self%pml%pml_type) /= 'CFS') then
         call mpih%error_stop(msg=': CPU PML time integration is currently implemented only for PML_type = CLASSIC, BERMUDEZ or CFS')
      endif
      if (trim(self%numerics%scheme_space) /= NUM_SCHEME_SPACE_FD_CENTERED) then
         call mpih%error_stop(msg=': CPU PML is currently implemented only for scheme_space = fd_centered')
      endif
      if (trim(adjustl(self%numerics%div_corr_var)) /= 'No') then
         call mpih%error_stop(msg=': CPU PML is currently implemented only without divergence correction')
      endif
      if (self%numerics%constrained_transport_D .or. self%numerics%constrained_transport_B) then
         call mpih%error_stop(msg=': CPU PML is currently implemented only for pure Maxwell states without CT variables')
      endif
      if (trim(self%numerics%scheme_time) /= NUM_SCHEME_TIME_RUNGE_KUTTA .or. .not. is_supported_ssp) then
         call mpih%error_stop(msg=': CPU PML is currently implemented only for SSP Runge-Kutta time integration')
      endif
   endsubroutine check_pml_configuration
   endsubroutine initialize_prism

   ! IO methods
   subroutine save_residuals(self)
   !< Save residuals history.
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: v    !< Counter.

   if (self%time%is_to_save(it_save=self%io%residuals_save)) then
      call self%adam%field%compute_normL2_residuals(grid=self%adam%grid, dq=self%dq, norm=self%adam%field%residuals)
      do v=1, self%nv
         call MPI_ALLREDUCE(MPI_IN_PLACE, self%adam%field%residuals(v), 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
         self%adam%field%residuals(v) = sqrt(self%adam%field%residuals(v))/sqrt(real(self%ni*self%nj*self%nk, R8P))
      enddo
      if (mpih%myrank==0) call self%io%save_residuals(it=self%time%it, time=self%time%time, &
                                                      blocks_number=self%blocks_number, residuals=self%adam%field%residuals)
   endif
   endsubroutine save_residuals

   subroutine save_simulation_data(self)
   !< Save all simulation data.
   class(prism_cpu_object), intent(inout) :: self    !< The equation.

   if ((self%time%is_to_save(it_save=self%io%it_save)).or.      &
       (self%time%is_to_save(it_save=self%io%restart_save)).or. &
       (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max))) then
      call self%update_ghost(q=self%q)
      call self%compute_auxiliary_fields

      if (self%time%is_to_save(it_save=self%io%it_save)) call self%save_xh5f(with_ghost=.true.)
      if (mod(self%time%it,self%io%restart_save)==0) call self%save_restart_files
      ! if (self%slices%is_to_save(it=self%time%it,it_max=self%time%it_max,time=self%time%time,time_max=self%time%time_max)) then
      !    call self%slices%save_mat(basename=self%io%output_basename, &
      !                              it=self%time%it,                  &
      !                              it_max=self%time%it_max,          &
      !                              time=self%time%time,              &
      !                              time_max=self%time%time_max,      &
      !                              adam=self%adam,                   &
      !                              q=self%q,                         &
      !                              q_name=self%q_name)

      ! endif
   endif
   if (self%pic%problem_type == SINGLE_PARTICLE_TYPE_PROBLEM) then
      call write_single_particle_output(filename='single_particle_output.dat', time=self%time%time, q_pic=self%q_pic)
   endif
   endsubroutine save_simulation_data

   ! IC/BC/sources
   subroutine compute_coils_current(self, q, gamma)
   !< Compute current coils sources (DC/AC with smooth envelope).
   class(prism_cpu_object), intent(inout)        :: self
   real(R8P),               intent(inout)        :: q(1:,          &
                                                      1-self%ngc:, &
                                                      1-self%ngc:, &
                                                      1-self%ngc:, &
                                                      1:)
   real(R8P),               intent(in), optional :: gamma
   character(len=128)                            :: fname
   real(R8P)                                     :: current_density
   real(R8P)                                     :: current_density_o
   real(R8P)                                     :: time_s
   real(R8P)                                     :: s, g
   real(R8P)                                     :: phi_rad, omega, theta
   real(R8P)                                     :: f_abs
   integer(I4P)                                  :: w_ac  ! =1 AC, =0 DC (branchless)
   integer(I4P)                                  :: coil_id
   integer(I4P)                                  :: i,j,k,b,n
   real(R8P), parameter                          :: f_tol = 1.0e-30_R8P

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, blocks_number=>self%blocks_number,                     &
             time=>self%time%time, dt=>self%time%dt, td=>self%coil%td,                                     &
             A=>self%coil%coil_amplitude, f=>self%coil%f, phase=>self%coil%phase, J_vec=>self%coil%J_vec,  &
             var_Jx=>self%physics%var_Jx, var_Jy=>self%physics%var_Jy, var_Jz=>self%physics%var_Jz)

      if (present(gamma)) then
         time_s = time + dt*gamma
      else
         time_s = time
      end if
      if (self%coil%total_coils_number >= 1_I4P) then

         ! Azzeramento limitato al supporto coil: preserva il contributo PIC fuori dalle bobine.
         do b=1, blocks_number
            do k=1-self%ngc, nk+self%ngc
               do j=1-self%ngc, nj+self%ngc
                  do i=1-self%ngc, ni+self%ngc
                     if (any(J_vec(:,i,j,k,b,:) /= 0._R8P)) then
                        q(var_Jx,i,j,k,b) = 0._R8P
                        q(var_Jy,i,j,k,b) = 0._R8P
                        q(var_Jz,i,j,k,b) = 0._R8P
                     endif
                  enddo
               enddo
            enddo
         enddo

         ! Envelope C^2: clamp(s) in [0,1], g(0)=0, g(1)=1, g'(0)=g'(1)=0, g''(0)=g''(1)=0
         if (td > 0._R8P) then
            s = time_s / td
         else
            s = 1._R8P
         end if
         s = max(0._R8P, min(1._R8P, s))
         g = 10._R8P*s**3 - 15._R8P*s**4 + 6._R8P*s**5
         do n=1, self%coil%total_coils_number
            coil_id = n

            phi_rad = phase(coil_id) * PI / 180._R8P
            omega   = 2._R8P * PI * f(coil_id)

            ! Se f ~ 0 -> DC (w_ac=0). Se f != 0 -> AC (w_ac=1). Branchless.
            f_abs = abs(f(coil_id))
            w_ac  = nint( (sign(1._R8P, f_abs - f_tol) + 1._R8P) * 0.5_R8P )

            ! Theta: DC -> theta = phi ; AC -> theta = omega*(t-td) + phi
            theta = w_ac * omega * (time_s - td) + phi_rad

            ! Unica formula: DC e AC
            current_density = A(coil_id) * g * cos(theta)

            do b=1, blocks_number
               do k=1-self%ngc, nk+self%ngc
                  do j=1-self%ngc, nj+self%ngc
                     do i=1-self%ngc, ni+self%ngc
                        q(var_Jx,i,j,k,b) = q(var_Jx,i,j,k,b) + current_density * J_vec(1,i,j,k,b,n)
                        q(var_Jy,i,j,k,b) = q(var_Jy,i,j,k,b) + current_density * J_vec(2,i,j,k,b,n)
                        q(var_Jz,i,j,k,b) = q(var_Jz,i,j,k,b) + current_density * J_vec(3,i,j,k,b,n)
                     enddo
                  enddo
               enddo
            enddo
         enddo
         ! Diagnostiche
         do n = 1, self%coil%total_coils_number
            theta = nint( (sign(1._R8P, abs(f(n)) - f_tol) + 1._R8P) * 0.5_R8P ) * (2._R8P*PI*f(n)) * (time_s - td) + &
                    phase(n)*PI/180._R8P
            current_density_o = A(n) * g * cos(theta)
            write(fname,'(A,SS,I0,A)') 'current_density_coil_', n, '.dat'
            call write_current_behavior_tab(trim(fname), time=time_s, current_density=current_density_o)
         enddo
      endif
   endassociate
   endsubroutine compute_coils_current

   subroutine apply_fWL_correction(self, q)
   !< Apply correction if a fWL is present
   class(prism_cpu_object), intent(inout) :: self            !< The equation.
   real(R8P),               intent(inout) :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1:)           !< Conservative variables.
   integer(I4P)                           :: face            !< Counter

   associate(ngc=>self%ngc, blocks_number=>self%blocks_number, dxyz=>self%adam%field%dxyz, layer=>self%fWLayer%layer, &
            ni_fWL=>self%fWLayer%ni_fWL, emin=>self%adam%field%emin, emax=>self%adam%field%emax,              &
            profile_extent=>self%fWLayer%profile_extent, profile_cells=>self%fWLayer%profile_cells,           &
            nj_fWL=>self%fWLayer%nj_fWL, nk_fWL=>self%fWLayer%nk_fWL, n=>self%fWLayer%n, s2=>self%fWLayer%s2, &
            alfa_D=>self%fWLayer%alfa_D, alfa_B=>self%fWLayer%alfa_B, beta_D=>self%fWLayer%beta_D,            &
            beta_B=>self%fWLayer%beta_B)
   if (allocated(self%fWLayer%C)) then
      do face=1, 6
         if (layer(face)) call apply_fWL_correction_fun(blocks_number = blocks_number,        &
                                                        ngc           = ngc,                  &
                                                        ni            = self%ni,             &
                                                        nj            = self%nj,             &
                                                        nk            = self%nk,             &
                                                        face          = face,                 &
                                                        profile_extent= profile_extent(face), &
                                                        profile_cells = profile_cells(face),  &
                                                        ni_fWL        = ni_fWL(:,:,face),    &
                                                        nj_fWL        = nj_fWL(:,:,face),    &
                                                        nk_fWL        = nk_fWL(:,:,face),    &
                                                        n             = n(face),             &
                                                        s2            = s2(face),            &
                                                        alfa_D        = alfa_D(face),        &
                                                        beta_D        = beta_D(face),        &
                                                        alfa_B        = alfa_B(face),        &
                                                        beta_B        = beta_B(face),        &
                                                        domain_emin   = self%adam%grid%domain_emin, &
                                                        domain_emax   = self%adam%grid%domain_emax, &
                                                        emin          = emin,                &
                                                        emax          = emax,                &
                                                        dxyz          = dxyz,                &
                                                        q             = q)
      enddo
   endif
   endassociate
   endsubroutine apply_fWL_correction

   subroutine set_boundary_conditions(self, q, s)
   !< Set boundary conditions of equation.
   class(prism_cpu_object), intent(inout) :: self                         !< The equation.
   real(R8P),               intent(inout) :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,1:)            !< Conservative variables.
   integer(I4P),  optional, intent(in)    :: s                            !< Stage counter.
   integer(I4P)                           :: b, c, i, j, k, v             !< Counter.
   integer(I4P)                           :: idelta,jdelta,kdelta         !< IJK delta step for extrapolation.
   integer(I4P)                           :: bc_type                      !< Boundary condition type.
   integer(I4P)                           :: crown                        !< Crown counter.
   integer(I4P)                           :: fec                          !< Boundary fec (1 to 26).
   integer(I4P)                           :: fec_1_6                      !< Boundary fec (1 to 6).
   integer(I4P)                           :: iref, jref, kref            !< Interior reference indexes for face BCs.
   integer(I4P)                           :: alfa_D, beta_D, gamma_D      !< Indici alfa beta gamma come in Barbas.
   integer(I4P)                           :: alfa_B, beta_B, gamma_B      !< Indici alfa beta gamma come in Barbas.
   real(R8P)                              :: s1                           !< Coefficiente pari a +-1.
   real(R8P)                              :: ngc_r, crown_r               !< Numero di gc totale, reale
   real(R8P)                              :: ref(1:self%nv)               !< Vettore di stato di riferimento per assegnazione gc.

   associate(local_map_bc_crown=>self%adam%maps%local_map_bc_crown,                                                              &
             nv=>self%nv, ngc=>self%ngc, q_bc_vars=>self%bc%q, dx=>self%adam%field%dxyz(1,:), dy=>self%adam%field%dxyz(2,:),     &
             dz=>self%adam%field%dxyz(3,:), ni=>self%ni, nj=>self%nj, nk=>self%nk, dt=>self%time%dt, chi=>self%physics%chi,      &
             nv_c=>self%physics%nv_c, nv_cl=>self%physics%nv_cl, constrained_transport_B=>self%numerics%constrained_transport_B, &
             constrained_transport_D=>self%numerics%constrained_transport_D)

   if (allocated(self%adam%maps%local_map_bc_crown)) then
      do crown=1, ngc
         do c=1, size(local_map_bc_crown, dim=1)
            b = local_map_bc_crown(c, 1 ,crown)
            if (b>0) then
               i       = local_map_bc_crown(c, 2 ,crown)
               j       = local_map_bc_crown(c, 3 ,crown)
               k       = local_map_bc_crown(c, 4 ,crown)
               idelta  = local_map_bc_crown(c, 5 ,crown)
               jdelta  = local_map_bc_crown(c, 6 ,crown)
               kdelta  = local_map_bc_crown(c, 7 ,crown)
               bc_type = local_map_bc_crown(c, 8 ,crown)
               fec     = local_map_bc_crown(c, 9 ,crown) !da qua la faccia e quindi la normale
               if (fec > 6_I4P) cycle
               fec_1_6 = fec_1_6_array(fec)
               if (bc_type == BC_EXTRAPOLATION) then
                  do v=1, nv!(nv_c-nv_cl)
                     q(v,i,j,k,b) = q(v,i-idelta,j-jdelta,k-kdelta,b) !ni,j,k coordinate della cella da cui prendo i valori
                  enddo
               elseif (bc_type == BC_NEUMANN) then
                  call compute_face_mirror_indexes(face=fec_1_6, ni=ni, nj=nj, nk=nk, i_gc=i, j_gc=j, k_gc=k, &
                                                   idelta=idelta, jdelta=jdelta, kdelta=kdelta,               &
                                                   i_d=iref, j_d=jref, k_d=kref)
                  do v=1, nv!(nv_c-nv_cl)
                     q(v,i,j,k,b) = q(v,iref,jref,kref,b)
                  enddo
               elseif (bc_type == BC_SILVER_MULLER) then
                  ! With outward normal n, the Silver-Muller conditions are
                  ! B_t = (n x E_d) / c, B_n = B_n,d, E_t = c (B_d x n), E_n = E_n,d.
                  call compute_face_mirror_indexes(face=fec_1_6, ni=ni, nj=nj, nk=nk, i_gc=i, j_gc=j, k_gc=k, &
                                                   idelta=idelta, jdelta=jdelta, kdelta=kdelta,               &
                                                   i_d=iref, j_d=jref, k_d=kref)
                  ref = q(:,iref,jref,kref,b)
                  select case(fec_1_6)
                  case(1)
                     s1 = -1.0_R8P
                     alfa_D = 2_I4P
                     beta_D = 3_I4P
                     gamma_D = 1_I4P
                     alfa_B = 5_I4P
                     beta_B = 6_I4P
                     gamma_B = 4_I4P
                  case(2)
                     s1 = 1.0_R8P
                     alfa_D = 2_I4P
                     beta_D = 3_I4P
                     gamma_D = 1_I4P
                     alfa_B = 5_I4P
                     beta_B = 6_I4P
                     gamma_B = 4_I4P
                  case(3)
                     s1 = -1.0_R8P
                     alfa_D = 3_I4P
                     beta_D = 1_I4P
                     gamma_D = 2_I4P
                     alfa_B = 6_I4P
                     beta_B = 4_I4P
                     gamma_B = 5_I4P
                  case(4)
                     s1 = 1.0_R8P
                     alfa_D = 3_I4P
                     beta_D = 1_I4P
                     gamma_D = 2_I4P
                     alfa_B = 6_I4P
                     beta_B = 4_I4P
                     gamma_B = 5_I4P
                  case(5)
                     s1 = -1.0_R8P
                     alfa_D = 1_I4P
                     beta_D = 2_I4P
                     gamma_D = 3_I4P
                     alfa_B = 4_I4P
                     beta_B = 5_I4P
                     gamma_B = 6_I4P
                  case(6)
                     s1 = 1.0_R8P
                     alfa_D = 1_I4P
                     beta_D = 2_I4P
                     gamma_D = 3_I4P
                     alfa_B = 4_I4P
                     beta_B = 5_I4P
                     gamma_B = 6_I4P
                  endselect
                  q(alfa_D, i,j,k,b) =  s1*C0*ref(beta_B)*EPS0
                  q(beta_D, i,j,k,b) = -s1*C0*ref(alfa_B)*EPS0
                  q(gamma_D,i,j,k,b) = ref(gamma_D)
                  q(alfa_B, i,j,k,b) = -s1/C0*ref(beta_D)/EPS0
                  q(beta_B, i,j,k,b) =  s1/C0*ref(alfa_D)/EPS0
                  q(gamma_B,i,j,k,b) = ref(gamma_B)

                  do v=(nv_c-nv_cl+1), nv
                     q(v,i,j,k,b) = q(v,iref,jref,kref,b)
                  enddo
               elseif (bc_type == BC_PEC) then
                  call compute_face_mirror_indexes(face=fec_1_6, ni=ni, nj=nj, nk=nk, i_gc=i, j_gc=j, k_gc=k, &
                                                   idelta=idelta, jdelta=jdelta, kdelta=kdelta,               &
                                                   i_d=iref, j_d=jref, k_d=kref)
                  ref = q(:,iref,jref,kref,b)
                  q(:,i,j,k,b) = ref
                  select case(fec_1_6)
                  case(1, 2)
                     q(VAR_DX,i,j,k,b) =  ref(VAR_DX)
                     q(VAR_DY,i,j,k,b) = -ref(VAR_DY)
                     q(VAR_DZ,i,j,k,b) = -ref(VAR_DZ)
                     q(VAR_BX,i,j,k,b) = -ref(VAR_BX)
                     q(VAR_BY,i,j,k,b) =  ref(VAR_BY)
                     q(VAR_BZ,i,j,k,b) =  ref(VAR_BZ)
                  case(3, 4)
                     q(VAR_DX,i,j,k,b) = -ref(VAR_DX)
                     q(VAR_DY,i,j,k,b) =  ref(VAR_DY)
                     q(VAR_DZ,i,j,k,b) = -ref(VAR_DZ)
                     q(VAR_BX,i,j,k,b) =  ref(VAR_BX)
                     q(VAR_BY,i,j,k,b) = -ref(VAR_BY)
                     q(VAR_BZ,i,j,k,b) =  ref(VAR_BZ)
                  case(5, 6)
                     q(VAR_DX,i,j,k,b) = -ref(VAR_DX)
                     q(VAR_DY,i,j,k,b) = -ref(VAR_DY)
                     q(VAR_DZ,i,j,k,b) =  ref(VAR_DZ)
                     q(VAR_BX,i,j,k,b) =  ref(VAR_BX)
                     q(VAR_BY,i,j,k,b) =  ref(VAR_BY)
                     q(VAR_BZ,i,j,k,b) = -ref(VAR_BZ)
                  endselect
               elseif (bc_type == BC_DIRICHLET) then
                     do v=1, nv
                        q(v,i,j,k,b) = 0._R8P
                     enddo
               elseif (bc_type == BC_PERIOD) then
                  q(:,i,j,k,b) = 0._R8P
                  select case(fec_1_6)
                  case(1)
                     q(1:nv_c,i,j,k,b) = q(1:nv_c,ni+i,j   ,k   ,b)
                  case(2)
                     q(1:nv_c,i,j,k,b) = q(1:nv_c,i-ni,j   ,k   ,b)
                  case(3)
                     q(1:nv_c,i,j,k,b) = q(1:nv_c,i   ,nj+j,k   ,b)
                  case(4)
                     q(1:nv_c,i,j,k,b) = q(1:nv_c,i   ,j-nj,k   ,b)
                  case(5)
                     q(1:nv_c,i,j,k,b) = q(1:nv_c,i   ,j   ,nk+k,b)
                  case(6)
                     q(1:nv_c,i,j,k,b) = q(1:nv_c,i   ,j   ,k-nk,b)
                  endselect
               endif
            endif
         enddo
      enddo
   endif

   call enforce_silver_muller_normal_bc_cpu(self=self, q=q)

   if (self%bc%bc_type(1) == BC_radiative .or. self%bc%bc_type(2) == BC_radiative &
       .or. self%bc%bc_type(3) == BC_radiative .or. self%bc%bc_type(4) == BC_radiative &
       .or. self%bc%bc_type(5) == BC_radiative .or. self%bc%bc_type(6) == BC_radiative) then
                                                                                             !Al momento scritta per funzionare solo
                                                                                             !con un secondo ordine
      call mpih%error_stop(msg='radiative BC not already implemented')

      !if (present(s)) then
      !   if (s==1_I4P) call self%rk_bc%initialize_stages(field=self%adam%field, q=q)
      !   if (self%ib%solids_number>0) then !calcolo stadio per le BC
      !      call self%rk_bc%compute_stage(field=self%adam%field, s=s, dt=self%time%dt, phi=self%ib%phi)
      !   else
      !      call self%rk_bc%compute_stage(field=self%adam%field, s=s, dt=self%time%dt)
      !   endif
      !   !Calcolo i residui per l'integrazione temporale delle BC (in un futuro da allineare con operatore spaziale qualsiasi)
      !   call self%compute_residuals_BC(s=s)
      !   !Imponi effettivamente la BC su q: unico punto del ciclo in cui si "uniscono"
      !   !Quindi basta cambiare gli indici di quel do per imporlo su una sola variabile, eventualmente
      !   !(O cambiare i cicli da 1:nv_c a nv_c-nv_cl+1:nv_c)
      !   if (allocated(self%adam%maps%local_map_bc_crown)) then
      !      do crown=1, ngc
      !         do c=1, size(local_map_bc_crown, dim=1)
      !            b = local_map_bc_crown(c, 1 ,crown)
      !            if (b>0) then
      !               i       = local_map_bc_crown(c, 2 ,crown)
      !               j       = local_map_bc_crown(c, 3 ,crown)
      !               k       = local_map_bc_crown(c, 4 ,crown)
      !               idelta  = local_map_bc_crown(c, 5 ,crown)
      !               jdelta  = local_map_bc_crown(c, 6 ,crown)
      !               kdelta  = local_map_bc_crown(c, 7 ,crown)
      !               bc_type = local_map_bc_crown(c, 8 ,crown)
      !               fec     = local_map_bc_crown(c, 9 ,crown) !da qua la faccia e quindi la normale
      !               fec_1_6 = fec_1_6_array(fec)
      !               if (bc_type == BC_radiative) then
      !                  do v=1, nv_c
      !                     q(v,i,j,k,b) = 2*self%rk_bc%q_bc_rk(v,i,j,k,b,s)-q(v,i-idelta,j-jdelta,k-kdelta,b)
      !                  enddo
      !               endif
      !            endif
      !         enddo
      !      enddo
      !   endif
      !   !Concludi assegnando lo stadio
      !   if (self%ib%solids_number>0) then
      !      call self%rk_bc%assign_stage(field=self%adam%field, s=s, phi=self%ib%phi)
      !   else
      !      call self%rk_bc%assign_stage(field=self%adam%field, s=s)
      !   endif
      !else !Mi serve solo per il t0, tanto ic è il vuoto praticamente sempre
      !   !q(v,i,j,k,b) = 0.0_R8P this is bugged, which are v,i,j,k,b? below the fix
      !   if (allocated(self%adam%maps%local_map_bc_crown)) then
      !      do crown=1, ngc
      !         do c=1, size(local_map_bc_crown, dim=1)
      !            b = local_map_bc_crown(c, 1 ,crown)
      !            if (b>0) then
      !               i = local_map_bc_crown(c, 2, crown)
      !               j = local_map_bc_crown(c, 3, crown)
      !               k = local_map_bc_crown(c, 4, crown)
      !               do v=1, nv_c
      !                  q(v,i,j,k,b) = 0.0_R8P
      !               enddo
      !            endif
      !         enddo
      !      enddo
      !   endif
      !endif
   endif
   endassociate
   endsubroutine set_boundary_conditions

   subroutine compute_face_mirror_indexes(face, ni, nj, nk, i_gc, j_gc, k_gc, idelta, jdelta, kdelta, i_d, j_d, k_d)
   !< Return the donor indexes mirrored across the selected boundary face.
   integer(I4P), intent(in)  :: face                     !< Face index in 1:6 numbering.
   integer(I4P), intent(in)  :: ni, nj, nk               !< Grid dimensions.
   integer(I4P), intent(in)  :: i_gc, j_gc, k_gc         !< Ghost-cell indexes.
   integer(I4P), intent(in)  :: idelta, jdelta, kdelta   !< One-step inward deltas.
   integer(I4P), intent(out) :: i_d, j_d, k_d            !< Mirrored donor indexes.

   i_d = i_gc - idelta
   j_d = j_gc - jdelta
   k_d = k_gc - kdelta

   select case(face)
   case(1)
      i_d = 1_I4P - i_gc
   case(2)
      i_d = 2_I4P * ni + 1_I4P - i_gc
   case(3)
      j_d = 1_I4P - j_gc
   case(4)
      j_d = 2_I4P * nj + 1_I4P - j_gc
   case(5)
      k_d = 1_I4P - k_gc
   case(6)
      k_d = 2_I4P * nk + 1_I4P - k_gc
   case default
      call mpih%error_stop(msg='compute_face_mirror_indexes: invalid face index '//trim(str(face)))
   endselect
   endsubroutine compute_face_mirror_indexes

   subroutine enforce_silver_muller_normal_bc_cpu(self, q)
   class(prism_cpu_object), intent(in)    :: self
   real(R8P),               intent(inout) :: q(1:,         &
                                               1-self%ngc:,&
                                               1-self%ngc:,&
                                               1-self%ngc:,1:)
   integer(I4P), parameter                :: sm_face_sweeps = 2_I4P
   integer(I4P)                           :: b, c, face, face_stage, hs, sweep, var_rho
   integer(I4P)                           :: i, j, k
   logical                                :: has_rho

   if (.not. allocated(self%adam%maps%local_map_bc_crown)) return
   if (.not. any(self%bc%bc_type == BC_SILVER_MULLER)) return

   hs = self%fdv_half_stencils(1)
   if (hs <= 0) return
   if (self%ngc < hs) call mpih%error_stop(msg='Silver_Muller requires ngc >= fdv_half_stencils(1)')

   has_rho = self%physics%physical_model == PIC_PHYSICAL_MODEL
   if (has_rho) then
      var_rho = self%nv
   else
      var_rho = 0_I4P
   endif

   do sweep=1, sm_face_sweeps
      do face_stage=1, 6
         do c=1, size(self%adam%maps%local_map_bc_crown, dim=1)
            b = self%adam%maps%local_map_bc_crown(c, 1, 1)
            if (b <= 0) cycle
            if (self%adam%maps%local_map_bc_crown(c, 8, 1) /= BC_SILVER_MULLER) cycle
            if (self%adam%maps%local_map_bc_crown(c, 9, 1) > 6_I4P) cycle

            i = self%adam%maps%local_map_bc_crown(c, 2, 1)
            j = self%adam%maps%local_map_bc_crown(c, 3, 1)
            k = self%adam%maps%local_map_bc_crown(c, 4, 1)
            face = fec_1_6_array(self%adam%maps%local_map_bc_crown(c, 9, 1))
            if (face /= face_stage) cycle
            if (.not. is_face_line_seed_cpu(face=face, i=i, j=j, k=k, ni=self%ni, nj=self%nj, nk=self%nk)) cycle

            call solve_silver_muller_normal_line_cpu(q=q, ngc=self%ngc, ni=self%ni, nj=self%nj, nk=self%nk,      &
                                                     b=b, face=face, i_seed=i, j_seed=j, k_seed=k, hs=hs,         &
                                                     dxyz=self%adam%field%dxyz(:,b), has_rho=has_rho, var_rho=var_rho)
         enddo
      enddo
      do face_stage=6, 1, -1
         do c=1, size(self%adam%maps%local_map_bc_crown, dim=1)
            b = self%adam%maps%local_map_bc_crown(c, 1, 1)
            if (b <= 0) cycle
            if (self%adam%maps%local_map_bc_crown(c, 8, 1) /= BC_SILVER_MULLER) cycle
            if (self%adam%maps%local_map_bc_crown(c, 9, 1) > 6_I4P) cycle

            i = self%adam%maps%local_map_bc_crown(c, 2, 1)
            j = self%adam%maps%local_map_bc_crown(c, 3, 1)
            k = self%adam%maps%local_map_bc_crown(c, 4, 1)
            face = fec_1_6_array(self%adam%maps%local_map_bc_crown(c, 9, 1))
            if (face /= face_stage) cycle
            if (.not. is_face_line_seed_cpu(face=face, i=i, j=j, k=k, ni=self%ni, nj=self%nj, nk=self%nk)) cycle

            call solve_silver_muller_normal_line_cpu(q=q, ngc=self%ngc, ni=self%ni, nj=self%nj, nk=self%nk,      &
                                                     b=b, face=face, i_seed=i, j_seed=j, k_seed=k, hs=hs,         &
                                                     dxyz=self%adam%field%dxyz(:,b), has_rho=has_rho, var_rho=var_rho)
         enddo
      enddo
   enddo
   endsubroutine enforce_silver_muller_normal_bc_cpu

   logical pure function is_face_line_seed_cpu(face, i, j, k, ni, nj, nk) result(is_seed)
   integer(I4P), intent(in) :: face, i, j, k, ni, nj, nk

   is_seed = .false.
   select case(face)
   case(1)
      is_seed = i == 0_I4P      .and. j >= 1_I4P .and. j <= nj .and. k >= 1_I4P .and. k <= nk
   case(2)
      is_seed = i == ni + 1_I4P .and. j >= 1_I4P .and. j <= nj .and. k >= 1_I4P .and. k <= nk
   case(3)
      is_seed = j == 0_I4P      .and. i >= 1_I4P .and. i <= ni .and. k >= 1_I4P .and. k <= nk
   case(4)
      is_seed = j == nj + 1_I4P .and. i >= 1_I4P .and. i <= ni .and. k >= 1_I4P .and. k <= nk
   case(5)
      is_seed = k == 0_I4P      .and. i >= 1_I4P .and. i <= ni .and. j >= 1_I4P .and. j <= nj
   case(6)
      is_seed = k == nk + 1_I4P .and. i >= 1_I4P .and. i <= ni .and. j >= 1_I4P .and. j <= nj
   endselect
   endfunction is_face_line_seed_cpu

   subroutine solve_silver_muller_normal_line_cpu(q, ngc, ni, nj, nk, b, face, i_seed, j_seed, k_seed, hs, dxyz, has_rho, var_rho)
   integer(I4P), intent(in)    :: ngc, ni, nj, nk, b, face, i_seed, j_seed, k_seed, hs, var_rho
   real(R8P),    intent(inout) :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(in)    :: dxyz(3)
   logical,      intent(in)    :: has_rho
   integer(I4P)                :: dir_t1, dir_t2
   integer(I4P)                :: var_d_t1, var_d_t2, var_d_n
   integer(I4P)                :: var_b_t1, var_b_t2, var_b_n

   call silver_muller_face_metadata_cpu(face=face, dir_t1=dir_t1, dir_t2=dir_t2,                              &
                                        var_d_t1=var_d_t1, var_d_t2=var_d_t2, var_d_n=var_d_n,               &
                                        var_b_t1=var_b_t1, var_b_t2=var_b_t2, var_b_n=var_b_n)

   call solve_silver_muller_normal_field_cpu(q=q, ngc=ngc, ni=ni, nj=nj, nk=nk, b=b, face=face,            &
                                             i_seed=i_seed, j_seed=j_seed, k_seed=k_seed, hs=hs,            &
                                             dxyz=dxyz, dir_t1=dir_t1, dir_t2=dir_t2,                        &
                                             var_t1=var_d_t1, var_t2=var_d_t2, var_n=var_d_n,               &
                                             use_target_var=has_rho, target_var=var_rho)

   call solve_silver_muller_normal_field_cpu(q=q, ngc=ngc, ni=ni, nj=nj, nk=nk, b=b, face=face,            &
                                             i_seed=i_seed, j_seed=j_seed, k_seed=k_seed, hs=hs,            &
                                             dxyz=dxyz, dir_t1=dir_t1, dir_t2=dir_t2,                        &
                                             var_t1=var_b_t1, var_t2=var_b_t2, var_n=var_b_n,               &
                                             use_target_var=.false., target_var=0_I4P)
   endsubroutine solve_silver_muller_normal_line_cpu

   subroutine silver_muller_face_metadata_cpu(face, dir_t1, dir_t2, var_d_t1, var_d_t2, var_d_n, var_b_t1, var_b_t2, var_b_n)
   integer(I4P), intent(in)  :: face
   integer(I4P), intent(out) :: dir_t1, dir_t2
   integer(I4P), intent(out) :: var_d_t1, var_d_t2, var_d_n
   integer(I4P), intent(out) :: var_b_t1, var_b_t2, var_b_n

   select case(face)
   case(1, 2)
      dir_t1 = 2_I4P ; dir_t2 = 3_I4P
      var_d_t1 = VAR_DY ; var_d_t2 = VAR_DZ ; var_d_n = VAR_DX
      var_b_t1 = VAR_BY ; var_b_t2 = VAR_BZ ; var_b_n = VAR_BX
   case(3, 4)
      dir_t1 = 3_I4P ; dir_t2 = 1_I4P
      var_d_t1 = VAR_DZ ; var_d_t2 = VAR_DX ; var_d_n = VAR_DY
      var_b_t1 = VAR_BZ ; var_b_t2 = VAR_BX ; var_b_n = VAR_BY
   case(5, 6)
      dir_t1 = 1_I4P ; dir_t2 = 2_I4P
      var_d_t1 = VAR_DX ; var_d_t2 = VAR_DY ; var_d_n = VAR_DZ
      var_b_t1 = VAR_BX ; var_b_t2 = VAR_BY ; var_b_n = VAR_BZ
   case default
      call mpih%error_stop(msg='silver_muller_face_metadata_cpu: invalid face index '//trim(str(face)))
   endselect
   endsubroutine silver_muller_face_metadata_cpu

   subroutine solve_silver_muller_normal_field_cpu(q, ngc, ni, nj, nk, b, face, i_seed, j_seed, k_seed, hs, dxyz, &
                                                   dir_t1, dir_t2, var_t1, var_t2, var_n, use_target_var, target_var)
   integer(I4P), intent(in)    :: ngc, ni, nj, nk, b, face, i_seed, j_seed, k_seed, hs
   integer(I4P), intent(in)    :: dir_t1, dir_t2, var_t1, var_t2, var_n, target_var
   real(R8P),    intent(inout) :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(in)    :: dxyz(3)
   logical,      intent(in)    :: use_target_var
   real(R8P)                   :: unknown(FDV_S_MAX)
   real(R8P)                   :: rhs, coeff, target
   integer(I4P)                :: eq, u, u_new, m, dir_n
   integer(I4P)                :: ic, jc, kc

   if (hs > FDV_S_MAX) call mpih%error_stop(msg='Silver_Muller hs exceeds FDV_S_MAX')

   unknown = 0._R8P
   dir_n = face_normal_direction_cpu(face)

   do eq=hs, 1, -1
      call interior_cell_from_face_cpu(face=face, eq=eq, ni=ni, nj=nj, nk=nk,                                 &
                                       i_seed=i_seed, j_seed=j_seed, k_seed=k_seed, ic=ic, jc=jc, kc=kc)

      if (use_target_var) then
         target = q(target_var, ic, jc, kc, b)
      else
         target = 0._R8P
      endif

      rhs = target
      rhs = rhs - tangential_divergence_at_cell_cpu(q=q, ngc=ngc, b=b, i=ic, j=jc, k=kc, hs=hs, dxyz=dxyz, &
                                                    dir_t1=dir_t1, dir_t2=dir_t2, var_t1=var_t1, var_t2=var_t2)
      rhs = rhs - normal_known_contribution_cpu(q=q, ngc=ngc, ni=ni, nj=nj, nk=nk, b=b, face=face,          &
                                                i=ic, j=jc, k=kc, hs=hs, dxyz=dxyz, var_n=var_n)

      u_new = hs - eq + 1_I4P
      do u=1, u_new-1
         m = eq + u - 1_I4P
         coeff = silver_muller_unknown_coefficient_cpu(face=face, m=m, hs=hs, ds=dxyz(dir_n))
         rhs = rhs - coeff * unknown(u)
      enddo

      coeff = silver_muller_unknown_coefficient_cpu(face=face, m=hs, hs=hs, ds=dxyz(dir_n))
      unknown(u_new) = rhs / coeff
   enddo

   do u=1, hs
      call ghost_cell_from_face_cpu(face=face, ghost_layer=u, ni=ni, nj=nj, nk=nk,                             &
                                    i_seed=i_seed, j_seed=j_seed, k_seed=k_seed, ic=ic, jc=jc, kc=kc)
      q(var_n, ic, jc, kc, b) = unknown(u)
   enddo
   endsubroutine solve_silver_muller_normal_field_cpu

   integer(I4P) pure function face_normal_direction_cpu(face) result(dir_n)
   integer(I4P), intent(in) :: face

   select case(face)
   case(1, 2)
      dir_n = 1_I4P
   case(3, 4)
      dir_n = 2_I4P
   case(5, 6)
      dir_n = 3_I4P
   case default
      dir_n = 0_I4P
   endselect
   endfunction face_normal_direction_cpu

   real(R8P) pure function silver_muller_unknown_coefficient_cpu(face, m, hs, ds) result(coeff)
   integer(I4P), intent(in) :: face, m, hs
   real(R8P),    intent(in) :: ds

   coeff = FD1_CC(m, hs) / ds
   if (face == 1_I4P .or. face == 3_I4P .or. face == 5_I4P) coeff = -coeff
   endfunction silver_muller_unknown_coefficient_cpu

   subroutine interior_cell_from_face_cpu(face, eq, ni, nj, nk, i_seed, j_seed, k_seed, ic, jc, kc)
   integer(I4P), intent(in)  :: face, eq, ni, nj, nk, i_seed, j_seed, k_seed
   integer(I4P), intent(out) :: ic, jc, kc

   select case(face)
   case(1)
      ic = eq
      jc = j_seed
      kc = k_seed
   case(2)
      ic = ni - eq + 1_I4P
      jc = j_seed
      kc = k_seed
   case(3)
      ic = i_seed
      jc = eq
      kc = k_seed
   case(4)
      ic = i_seed
      jc = nj - eq + 1_I4P
      kc = k_seed
   case(5)
      ic = i_seed
      jc = j_seed
      kc = eq
   case(6)
      ic = i_seed
      jc = j_seed
      kc = nk - eq + 1_I4P
   case default
      ic = 0_I4P
      jc = 0_I4P
      kc = 0_I4P
   endselect
   endsubroutine interior_cell_from_face_cpu

   subroutine ghost_cell_from_face_cpu(face, ghost_layer, ni, nj, nk, i_seed, j_seed, k_seed, ic, jc, kc)
   integer(I4P), intent(in)  :: face, ghost_layer, ni, nj, nk, i_seed, j_seed, k_seed
   integer(I4P), intent(out) :: ic, jc, kc

   select case(face)
   case(1)
      ic = 1_I4P - ghost_layer
      jc = j_seed
      kc = k_seed
   case(2)
      ic = ni + ghost_layer
      jc = j_seed
      kc = k_seed
   case(3)
      ic = i_seed
      jc = 1_I4P - ghost_layer
      kc = k_seed
   case(4)
      ic = i_seed
      jc = nj + ghost_layer
      kc = k_seed
   case(5)
      ic = i_seed
      jc = j_seed
      kc = 1_I4P - ghost_layer
   case(6)
      ic = i_seed
      jc = j_seed
      kc = nk + ghost_layer
   case default
      ic = 0_I4P
      jc = 0_I4P
      kc = 0_I4P
   endselect
   endsubroutine ghost_cell_from_face_cpu

   real(R8P) pure function tangential_divergence_at_cell_cpu(q, ngc, b, i, j, k, hs, dxyz, dir_t1, dir_t2, var_t1, var_t2) result(div_t)
   integer(I4P), intent(in) :: ngc, b, i, j, k, hs, dir_t1, dir_t2, var_t1, var_t2
   real(R8P),    intent(in) :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(in) :: dxyz(3)
   integer(I4P)             :: m

   div_t = 0._R8P
   do m=1, hs
      div_t = div_t + FD1_CC(m, hs) * (component_along_direction_cpu(q, ngc, b, var_t1, dir_t1, i, j, k,  m) - &
                                       component_along_direction_cpu(q, ngc, b, var_t1, dir_t1, i, j, k, -m)) / dxyz(dir_t1)
      div_t = div_t + FD1_CC(m, hs) * (component_along_direction_cpu(q, ngc, b, var_t2, dir_t2, i, j, k,  m) - &
                                       component_along_direction_cpu(q, ngc, b, var_t2, dir_t2, i, j, k, -m)) / dxyz(dir_t2)
   enddo
   endfunction tangential_divergence_at_cell_cpu

   real(R8P) pure function normal_known_contribution_cpu(q, ngc, ni, nj, nk, b, face, i, j, k, hs, dxyz, var_n) result(div_n_known)
   integer(I4P), intent(in) :: ngc, ni, nj, nk, b, face, i, j, k, hs, var_n
   real(R8P),    intent(in) :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(in) :: dxyz(3)
   integer(I4P)             :: dir_n, m, depth

   dir_n = face_normal_direction_cpu(face)
   depth = normal_depth_from_cell_cpu(face=face, i=i, j=j, k=k, ni=ni, nj=nj, nk=nk)
   div_n_known = 0._R8P

   select case(face)
   case(1, 3, 5)
      do m=1, hs
         div_n_known = div_n_known + FD1_CC(m, hs) * component_along_direction_cpu(q, ngc, b, var_n, dir_n, i, j, k,  m) / dxyz(dir_n)
      enddo
      do m=1, depth-1
         div_n_known = div_n_known - FD1_CC(m, hs) * component_along_direction_cpu(q, ngc, b, var_n, dir_n, i, j, k, -m) / dxyz(dir_n)
      enddo
   case(2, 4, 6)
      do m=1, depth-1
         div_n_known = div_n_known + FD1_CC(m, hs) * component_along_direction_cpu(q, ngc, b, var_n, dir_n, i, j, k,  m) / dxyz(dir_n)
      enddo
      do m=1, hs
         div_n_known = div_n_known - FD1_CC(m, hs) * component_along_direction_cpu(q, ngc, b, var_n, dir_n, i, j, k, -m) / dxyz(dir_n)
      enddo
   endselect
   endfunction normal_known_contribution_cpu

   integer(I4P) pure function normal_depth_from_cell_cpu(face, i, j, k, ni, nj, nk) result(depth)
   integer(I4P), intent(in) :: face, i, j, k, ni, nj, nk

   select case(face)
   case(1)
      depth = i
   case(2)
      depth = ni - i + 1_I4P
   case(3)
      depth = j
   case(4)
      depth = nj - j + 1_I4P
   case(5)
      depth = k
   case(6)
      depth = nk - k + 1_I4P
   case default
      depth = 0_I4P
   endselect
   endfunction normal_depth_from_cell_cpu

   real(R8P) pure function component_along_direction_cpu(q, ngc, b, var, dir, i, j, k, offset) result(value_)
   integer(I4P), intent(in) :: ngc, b, var, dir, i, j, k, offset
   real(R8P),    intent(in) :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)

   select case(dir)
   case(1)
      value_ = q(var, i + offset, j, k, b)
   case(2)
      value_ = q(var, i, j + offset, k, b)
   case(3)
      value_ = q(var, i, j, k + offset, b)
   case default
      value_ = 0._R8P
   endselect
   endfunction component_along_direction_cpu

   !subroutine compute_residuals_BC(self,s)
   !!< Compute residuals BCs.
   !!< La sua scrittura si lega all'ordine di interpolazione dell'operatore spaziale. Al momento
   !!< e' scritto per operatore di secondo ordine (1 punto ghost).
   !class(prism_cpu_object), intent(inout) :: self                    !< The equation.
   !integer(I4P),            intent(in)    :: s                       !< Stage counter.
   !real(R8P)                              :: ds                      !< Distanza tra le celle in x, y o z.
   !integer(I4P)                           :: i,j,k,b,c,v             !< Counter
   !integer(I4P)                           :: idelta,jdelta,kdelta    !< IJK delta step for extrapolation.
   !integer(I4P)                           :: bc_type                 !< Boundary condition type.
   !integer(I4P)                           :: crown                   !< Crown counter.
   !integer(I4P)                           :: fec                     !< Boundary fec (1 to 26).
   !integer(I4P)                           :: fec_1_6                 !< Boundary fec (1 to 6).
   !associate(local_map_bc_crown=>self%adam%maps%local_map_bc_crown,                                                          &
   !          nv=>self%nv, ngc=>self%ngc, q_bc_vars=>self%bc%q, dx=>self%adam%field%dxyz(1,:), dy=>self%adam%field%dxyz(2,:), &
   !          dz=>self%adam%field%dxyz(3,:), ni=>self%ni, nj=>self%nj, nk=>self%nk, dt=>self%time%dt, chi=>self%physics%chi,  &
   !          nv_c=>self%physics%nv_c, nv_cl=>self%physics%nv_cl, div_corr_var=>self%numerics%div_corr_var,                   &
   !          constrained_transport_B=>self%numerics%constrained_transport_B,                                                 &
   !          constrained_transport_D=>self%numerics%constrained_transport_D, q_rk=>self%rk%q_rk,                             &
   !          q_bc_rk=>self%rk_bc%q_bc_rk,dq_bc_rk=>self%rk_bc%dq_bc_rk)
   !if (allocated(self%adam%maps%local_map_bc_crown)) then
   !   do crown=1, ngc
   !      do c=1, size(local_map_bc_crown, dim=1)
   !         b = local_map_bc_crown(c, 1 ,crown)
   !         if (b>0) then
   !            bc_type = local_map_bc_crown(c, 8 ,crown)
   !            i       = local_map_bc_crown(c, 2 ,crown)
   !            j       = local_map_bc_crown(c, 3 ,crown)
   !            k       = local_map_bc_crown(c, 4 ,crown)
   !            idelta  = local_map_bc_crown(c, 5 ,crown)
   !            jdelta  = local_map_bc_crown(c, 6 ,crown)
   !            kdelta  = local_map_bc_crown(c, 7 ,crown)
   !            fec     = local_map_bc_crown(c, 9 ,crown) !da qua la faccia e quindi la normale
   !            fec_1_6 = fec_1_6_array(fec)
   !            if (fec <= 6) then
   !               select case(fec)
   !               case(1) !xmin
   !                  ds = dx(b)
   !               case(2) !xmax
   !                  ds = dx(b)
   !               case(3) !ymin
   !                  ds = dy(b)
   !               case(4) !ymax
   !                  ds = dy(b)
   !               case(5) !zmin
   !                  ds = dz(b)
   !               case(6) !zmax
   !                  ds = dz(b)
   !               end select
   !               do v = 1,6
   !                  dq_bc_rk(v,i,j,k,b) = -C0*(q_bc_rk(v,i,j,k,b,s)-q_rk(v,i-idelta,j-jdelta,k-kdelta,b,s))/ds
   !               enddo
   !               if (nv_cl == 1_I4P) then
   !                  dq_bc_rk(nv_c,i,j,k,b) = -chi*C0*(q_bc_rk(nv_c,i,j,k,b,s)- &
   !                                             q_rk(nv_c,i-idelta,j-jdelta,k-kdelta,b,s))/ds
   !               elseif (nv_cl == 2_I4P) then
   !                  dq_bc_rk(nv_c-1,i,j,k,b) = -chi*C0*(q_bc_rk(nv_c-1,i,j,k,b,s)- &
   !                                             q_rk(nv_c-1,i-idelta,j-jdelta,k-kdelta,b,s))/ds
   !                  dq_bc_rk(nv_c,i,j,k,b) = -chi*C0*(q_bc_rk(nv_c,i,j,k,b,s)- &
   !                                             q_rk(nv_c,i-idelta,j-jdelta,k-kdelta,b,s))/ds
   !               endif
   !               !if (div_corr_var == DIV_CORR_VAR_HYPER) then
   !               !   if (constrained_transport_D .and. .not.constrained_transport_B) &
   !               !      dq_bc_rk(nv_c,i,j,k,b) = -chi*C0*(q_bc_rk(nv_c,i,j,k,b,s)- &
   !               !                                 q_rk(nv_c,i-idelta,j-jdelta,k-kdelta,b,s))/ds
   !               !   if (.not.constrained_transport_D .and. constrained_transport_B) &
   !               !      dq_bc_rk(nv_c,i,j,k,b) = -chi*C0*(q_bc_rk(nv_c,i,j,k,b,s)- &
   !               !                                 q_rk(nv_c,i-idelta,j-jdelta,k-kdelta,b,s))/ds
   !               !   if (constrained_transport_D .and. constrained_transport_B) &
   !               !      dq_bc_rk(nv_c-1,i,j,k,b) = -chi*C0*(q_bc_rk(nv_c-1,i,j,k,b,s)- &
   !               !                                 q_rk(nv_c-1,i-idelta,j-jdelta,k-kdelta,b,s))/ds
   !               !   if (constrained_transport_D .and. constrained_transport_B) &
   !               !      dq_bc_rk(nv_c,i,j,k,b) = -chi*C0*(q_bc_rk(nv_c,i,j,k,b,s)- &
   !               !                                 q_rk(nv_c,i-idelta,j-jdelta,k-kdelta,b,s))/ds
   !               !endif
   !            endif
   !         endif
   !      enddo
   !   enddo
   !endif
   !endassociate
   !endsubroutine compute_residuals_BC

   subroutine set_initial_conditions(self, is_restart) !DA CORREGGERE CON NV_PIC QUANDO SERVE PER BC CARICA SE MODELLO PIC ATTIVO
   !< Set initial conditions and coils on field.
   class(prism_cpu_object), intent(inout) :: self       !< The equation.
   logical,                 intent(in)    :: is_restart !< Branching sentinel for restart/non restart path.
   real(R8P)                              :: max_div_D
   real(R8P)                              :: max_div_B
   real(R8P)                              :: max_div_J

   if (.not.is_restart) call self%ic%set_initial_conditions(physics=self%physics, field=self%adam%field, &
                                                            grid=self%adam%grid, q=self%q)
   call self%initialize_pic_time_zero()

   call self%initialize_coils
   call self%compute_coils_current_time_zero()

   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) &
                                 call self%impose_pic_fields_time_zero(ivar=VAR_DX)
   if (maxval(abs(self%q(self%physics%var_Jx:self%physics%var_Jz,:,:,:,:))) > 0.0_R8P) &
                                 call self%impose_pic_fields_time_zero(ivar=VAR_BX)

   call self%apply_fWL_correction(q=self%q)
   if (self%external_fields%ef_type/=EF_TYPE_NONE) &
      call self%external_fields%add_external_fields(field=self%adam%field, grid=self%adam%grid, &
                                                    time=0._R8P, dt=0._R8P, q=self%q)
   call self%weight_pic_fields_time_zero()

   call self%compute_divergence(hs=self%fdv_half_stencils(1), ivar=1_I4P, q=self%q(VAR_DX:VAR_DZ,:,:,:,:), &
                                   divergence=self%divergence(1,:,:,:,:))
   call self%compute_divergence(hs=self%fdv_half_stencils(1), ivar=1_I4P, q=self%q(VAR_BX:VAR_BZ,:,:,:,:), &
                                   divergence=self%divergence(2,:,:,:,:))
   call compute_max_divergence_outside_absorbing_layers(self=self, hs=self%fdv_half_stencils(1), max_div_D=max_div_D, &
                                                        max_div_B=max_div_B, max_div_J=max_div_J)
   call mpih%print_message('Initial conditions setting completed')
   if (self%physics%physical_model == EM_PHYSICAL_MODEL .or. self%physics%physical_model == ADIM_EM_PHYSICAL_MODEL) then
      call mpih%print_message('   max div(D) outside absorbing layers at t0='//trim(str(max_div_D)))
   elseif (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
      call mpih%print_message('   max div(D)-rho outside absorbing layers at t0='//trim(str(max_div_D)))
   endif
   call mpih%print_message('   max div(B) outside absorbing layers at t0='//trim(str(max_div_B)))
   endsubroutine set_initial_conditions

   subroutine update_ghost(self, q, step, s)
   !< Update ghost cells (intra-realm only).
   !< If not specified all steps are performed, synchronous computation.
   !<
   !< Inter-realm seam ghost cells are NO LONGER filled here: the forest's
   !< Phase 2 seam fill (between begin_stage and end_stage) populates them
   !< via the `fill_seam_from_peer_forest` TBP operating on the realm's
   !< currently-active buffer (selected internally via `stage_active`).
   !< The `realm(:)` optional dummy that used to thread through this
   !< signature has been retired.
   class(prism_cpu_object), intent(inout)                   :: self            !< The equation.
   real(R8P),               intent(inout)                   :: q(1:,         &
                                                                 1-self%ngc:,&
                                                                 1-self%ngc:,&
                                                                 1-self%ngc:,&
                                                                 1:)           !< Conservative variables.
   integer(I4P),            intent(in),    optional         :: step            !< Step to be perfordmed in asyncronous comp.
   integer(I4P),            intent(in),    optional         :: s               !< Stage counter.
   logical                                                  :: do_local_update !< Flag for triggering local update.
   logical                                                  :: do_set_bc       !< Flag for triggering setting bc.
   ! perform local update if step is not specified or if first step is selected
   do_local_update = .false.
   do_set_bc       = .false.
   if (.not.present(step)) then
      do_local_update = .true.
      do_set_bc       = .true.
   else
      if (step==1) do_local_update = .true.
      if (step==3) do_set_bc       = .true.
   endif
   if (do_local_update) then
      call self%adam%field%update_ghost_local(grid=self%adam%grid, maps=self%adam%maps, q=q)
   endif
   call self%adam%field%update_ghost_mpi(grid=self%adam%grid, maps=self%adam%maps, q=q, step=step)
   if (do_set_bc) then
      call self%set_boundary_conditions(q=q, s=s)
   endif
   if (present(s)) then
      call self%compute_coils_current(q=q, gamma=self%rk%gamm(s))
   else
      call self%compute_coils_current(q=q)
   endif
   endsubroutine update_ghost

   ! forest orchestrator contract methods overridings
   subroutine initialize_forest(self, filename, realms_number, memory_avail, nv, verbose)
   !< Initialize this realm from scratch: PRISM init, IC injection (or restart load), initial ghost update,
   !< initial diagnostics dump, IO files open, plus PIC/leapfrog priming if those schemes are active.
   !<
   !< Invoked by forest%initialize. The forest writes `self%realm_index = is`
   !< BEFORE calling this routine, so the body can already read the 1-based
   !< forest position from `self%realm_index` if needed.
   class(prism_cpu_object), intent(inout)           :: self          !< The realm.
   character(*),            intent(in)              :: filename      !< Input parameters file name.
   integer(I4P),            intent(in),    optional :: realms_number !< Realm count; divides the per-process budget.
   real(R8P),               intent(in),    optional :: memory_avail  !< Per-process memory budget override.
   integer(I4P),            intent(in),    optional :: nv            !< Number of field variables override.
   logical,                 intent(in),    optional :: verbose       !< Trigger verbose output.
   real(R8P)                                        :: max_div_D     !< Maximum of divergence of D field.
   real(R8P)                                        :: max_div_B     !< Maximum of divergence of B
   real(R8P)                                        :: max_div_J     !< Maximum of divergence of J field.
   integer(I4P)                                     :: r             !< Auxiliary variable to identify fWL presence
   real(R8P)                                        :: F_l(3)        !< Lorentz force for leapfrog preliminary integration.
   integer(I4P)                                     :: i, n, b, ind  !< Counters.

   call self%initialize_prism(filename=filename, realms_number=realms_number)
   if (self%io%restart) then
      call mpih%print_message('restart simulation from "'//trim(self%io%restart_basename)//'" files')
      call self%load_restart_files(t=self%time%it, time=self%time%time)
      call mpih%print_message('restart [t, time]: '//trim(str(self%time%it))//', '//trim(str(self%time%time)))
      call self%set_initial_conditions(is_restart=self%io%restart)
   else
      call mpih%print_message('impose initial conditions start')
      do i=1, self%ic%amr_iterations
         call mpih%print_message('  AMR/set IC iteration:'//trim(str(i,.true.)))
         call self%set_initial_conditions(is_restart=self%io%restart)
         call self%amr_update
      enddo
      call self%set_initial_conditions(is_restart=self%io%restart)
      call self%adam%make_comm_local_maps_ghost_bc
      self%time%time = 0._R8P
      self%time%it = 0
      call mpih%print_message('impose initial conditions finish')
   endif
   call self%update_ghost(q=self%q)
   !Il controllo sulla divergenza dei campi a t = 0 lo metterai poi a valle di questo update ghost
   call mpih%print_message('Coils initialization values')
   do n=1, self%coil%total_coils_number
      call self%compute_divergence(hs=self%fdv_half_stencils(1), ivar=1_I4P, q=self%coil%J_vec(1:3,:,:,:,:,n), &
                                   divergence=self%divergence(3,:,:,:,:))
      call mpih%print_message('Coil n='//trim(str(n,.true.)))
      call mpih%print_message('   max div(J)='//trim(str(maxval(abs(self%divergence(3,:,:,:,:)))*self%coil%coil_amplitude(n))))
   enddo

   call mpih%print_message('assigned block number: '//trim(str(self%adam%field%blocks_number,.true.)))
   do b = 1, self%adam%field%blocks_number
      call mpih%print_message('  b='//trim(str(b,.true.))//' code='//trim(str(self%adam%field%code(b))))
   enddo

   if (.not.self%io%restart .and. self%pic%problem_type == SINGLE_PARTICLE_TYPE_PROBLEM) then
      call initialize_single_particle_output(filename='single_particle_output.dat')
   endif

   associate(hs => self%fdv_half_stencil)
   call self%compute_divergence(hs=hs, ivar=1, q=self%q, divergence=self%divergence(1,:,:,:,:))
   call self%compute_divergence(hs=hs, ivar=4, q=self%q, divergence=self%divergence(2,:,:,:,:))
   call self%compute_divergence(hs=hs, ivar=self%physics%var_Jx, q=self%q, divergence=self%divergence(3,:,:,:,:))
   endassociate

   call compute_max_divergence_outside_absorbing_layers(self=self, hs=self%fdv_half_stencils(1), max_div_D=max_div_D, &
                                                        max_div_B=max_div_B, max_div_J=max_div_J)
   if (self%physics%physical_model == EM_PHYSICAL_MODEL .or. self%physics%physical_model == ADIM_EM_PHYSICAL_MODEL) then
      call mpih%print_message('   max div(D) outside absorbing layers at t0 after update_ghost='//trim(str(max_div_D)))
   elseif (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
      call mpih%print_message('   max div(D)-rho outside absorbing layers at t0 after update ghost='//trim(str(max_div_D)))
   endif
   call mpih%print_message('   max div(B) outside absorbing layers at t0 after update_ghost='//trim(str(max_div_B)))

   call self%save_simulation_data
   call self%compute_energy
   if (self%grms%do_save_history) call self%compute_grms
   !call self%save_energy_error(is_to_open=.true.)
   call self%save_energy_history(is_to_open=.true.)
   call self%save_grms_history(is_to_open=.true.)
   call self%save_divergence_history(is_to_open=.true., div_D=max_div_D, div_B=max_div_B, div_J=max_div_J)
   call self%io%open_file_residuals(nv=self%nv)

   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
      if(self%pic%scheme_time==NUM_SCHEME_TIME_PIC_LEAPFROG) then
         ! first time integration done apart with explicit euler scheme to iniziale leapfrog
         call self%leapfrog_pic%assign_step(grid=self%adam%grid, s=1, q_pic=self%q_pic)
         call self%compute_dt
         !< Pic residual computation
         !Qua ci va il calcolo dei campi nelle posizioni delle particelle se PIC
         !ma devi metterlo nell'inizializzazione per coerenza
         !< Integration of equations
         self%q_pic(1,:) = self%q_pic(1,:) + self%time%dt * self%q_pic(4,:)
         self%q_pic(2,:) = self%q_pic(2,:) + self%time%dt * self%q_pic(5,:)
         self%q_pic(3,:) = self%q_pic(3,:) + self%time%dt * self%q_pic(6,:)
         do i = 1, self%pic%particle_number
            F_l = crossproduct(self%q_pic(4:6,i), self%pic_fields(4:6,i))
            self%q_pic(4,i) = self%q_pic(4,i)+self%time%dt*self%q_pic(7,i)/self%q_pic(8,i)*(self%pic_fields(1,i)+F_l(1))
            self%q_pic(5,i) = self%q_pic(5,i)+self%time%dt*self%q_pic(7,i)/self%q_pic(8,i)*(self%pic_fields(2,i)+F_l(2))
            self%q_pic(6,i) = self%q_pic(6,i)+self%time%dt*self%q_pic(7,i)/self%q_pic(8,i)*(self%pic_fields(3,i)+F_l(3))
            !self%q_pic(4:6,i) = self%q_pic(4:6,i) + self%time%dt * self%q_pic(8,i) / self%q_pic(7,i) * &
            !                  (pic_fields(1:3,i) + crossproduct(self%q_pic(4:6,i), pic_fields(4:6,i)))
         enddo
      endif
   endif

   if (self%numerics%scheme_time==NUM_SCHEME_TIME_LEAPFROG) then
      ! first time integration done apart with explicit euler scheme to iniziale leapfrog
      call self%leapfrog%assign_step(field=self%adam%field, s=1, q=self%q)
      call self%compute_dt
      !Qua c'era il calcolo delle correnti delle particelle se PIC
      !Ora è nelle condizioni iniziali per coerenza con if legato a se ho pic o meno
      call self%compute_residuals(q=self%q, dq=self%dq)
      self%q = self%q + self%time%dt * self%dq
   endif
   endsubroutine initialize_forest

   subroutine compute_local_dt_forest(self, dt_local)
   !< Compute this realm's local stability-limited dt (no MPI reduction).
   !<
   !< Invoked by forest%compute_global_dt during the min reduction across all realms in the forest. The reduction
   !< itself is the orchestrator's job; this method computes only the value local to `self`.
   class(prism_cpu_object), intent(in)  :: self     !< The realm.
   real(R8P),               intent(out) :: dt_local !< Local stability-limited dt.
   real(R8P)                            :: umax     !< Maximum speed of waves propagation (light speed).
   real(R8P)                            :: dxyz_min !< Minimal space step.

   dxyz_min = huge(1._R8P)
   associate(blocks_number=>self%blocks_number, dxyz=>self%adam%field%dxyz, chi=>self%physics%chi, evmax=>self%physics%evmax)
   call compute_dxyz_min(blocks_number=blocks_number, dxyz=dxyz, dxyz_min=dxyz_min)
   umax = evmax
   dt_local = self%time%CFL*dxyz_min / umax
   endassociate
   endsubroutine compute_local_dt_forest

   subroutine advance_one_step_forest(self, dt)
   !< Advance this realm by one full timestep of size `dt`.
   !<
   !< Invoked by forest%evolve_one_step once per realm per timestep. Owns the integration itself, i.e. everything that turns
   !< `q` at time `t` into `q` at time `t + dt`.
   class(prism_cpu_object), intent(inout) :: self    !< The realm.
   real(R8P),               intent(in)    :: dt      !< Timestep size from the forest's global reduction.
   real(R8P)                              :: dt_step !< Local copy, possibly capped for time_max.

   self%time%it = self%time%it + 1
   dt_step = dt
   if ((self%time%it_max <= 0).and.(self%time%time+dt_step > self%time%time_max)) dt_step = self%time%time_max - self%time%time
   self%time%dt = dt_step
   call self%integrate
   self%time%time = self%time%time + dt_step
   call self%time%print_progress(nodes_number=self%adam%tree%nodes_number)
   endsubroutine advance_one_step_forest

   function stages_per_step_forest(self) result(K)
   !< Number of integrator stages this realm exposes per step.
   !<
   !< For the multi-realm path the forest drives the stage loop, so it
   !< needs to know `K` up front. ONLY `runge-kutta-ssp-*` is split into
   !< per-stage TBPs (`begin_stage_forest` / `end_stage_forest`): the
   !< staged protocol reads `gamm(k)` per stage and applies `beta(:)` in
   !< `close_step_forest` — coefficients the low-storage schemes do not
   !< even allocate (issue #25: LS schemes used to pass this gate silently
   !< and crash/corrupt far downstream on both backends). The forest
   !< queries this TBP only on the STAGED branch, so refusing here leaves
   !< the fused N=1/no-seam fast path — where LS schemes legitimately run —
   !< untouched. Other integrators (Yoshida, Leapfrog, Blanes-Moan, CFM)
   !< must error-stop here as well when they become selectable.
   class(prism_cpu_object), intent(in) :: self !< The realm.
   integer(I4P)                        :: K    !< Number of integrator stages per step.

   select case(self%rk%scheme)
   case(RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54)
      K = self%rk%nrk
   case default
      K = 0
      call mpih%error_stop(msg=': RK scheme "'//trim(adjustl(self%rk%scheme))//'" is not stage-splittable: '// &
                               'the staged forest path (multi-realm, or intra-realm AMR seam faces) requires '// &
                               'an SSP scheme (runge-kutta-ssp-*); low-storage schemes run only on the fused '// &
                               'single-realm/no-seam fast path')
   endselect
   endfunction stages_per_step_forest

   subroutine open_step_forest(self, dt)
   !< Per-step prologue on the multi-realm path: set dt, init RK stages.
   !<
   !< Mirrors the head of `integrate_rk_ssp`: external-field prelude (if active), `rk%initialize_stages(q=self%q)`, and
   !< the per-step time bookkeeping `advance_one_step_forest` does inline (it increment, dt cap for time_max, self%time%dt update).
   !< Time advance and progress print run in `close_step_forest`, mirroring the legacy ordering.
   class(prism_cpu_object), intent(inout) :: self    !< The realm.
   real(R8P),               intent(in)    :: dt      !< Timestep size from the forest.
   real(R8P)                              :: dt_step !< Local copy, possibly capped for time_max.

   self%time%it = self%time%it + 1
   dt_step = dt
   if ((self%time%it_max <= 0).and.(self%time%time+dt_step > self%time%time_max)) dt_step = self%time%time_max - self%time%time
   self%time%dt = dt_step
   if (self%external_fields%ef_type/=EF_TYPE_NONE) &
      call self%external_fields%sub_external_fields(field=self%adam%field, grid=self%adam%grid, &
                                                    time=self%time%time, dt=self%time%dt, q=self%q)
   call self%rk%initialize_stages(field=self%adam%field, q=self%q)
   endsubroutine open_step_forest

   subroutine begin_stage_forest(self, k, K_total, dt, realm)
   !< Begin integrator stage `k` on the multi-realm path.
   !<
   !< Mirrors the FIRST line of `integrate_rk_ssp`'s substage body:
   !< `rk%compute_stage(s=k, dt=self%time%dt [, phi=ib%phi])`. This populates
   !< `rk%q_rk(:,:,:,:,:,k)` from previously computed substages and from
   !< `self%q`, and publishes the stage buffer to peers by setting
   !< `self%stage_active = k`. No ghost reads, no peer-realm access — peer
   !< realms may not yet have opened their stage-`k` buffer when this fires.
   class(prism_cpu_object), intent(inout)                   :: self     !< The realm.
   integer(I4P),            intent(in)                      :: k        !< Stage index (1..K_total).
   integer(I4P),            intent(in)                      :: K_total  !< Forest-wide stage count for this step.
   real(R8P),               intent(in)                      :: dt       !< Timestep size from the forest.
   class(realm_object),     intent(inout), optional, target :: realm(:) !< Sibling realms (contract parity).

   self%stage_active = k
   if (self%ib%solids_number>0) then
      call self%rk%compute_stage(field=self%adam%field, s=k, dt=self%time%dt, phi=self%ib%phi)
   else
      call self%rk%compute_stage(field=self%adam%field, s=k, dt=self%time%dt)
   endif
   endsubroutine begin_stage_forest

   subroutine end_stage_forest(self, k, K_total, dt, realm, flux_register)
   !< Close integrator stage `k`: residuals + stage assignment in one sweep.
   class(prism_cpu_object),     intent(inout)                   :: self          !< The realm.
   integer(I4P),                intent(in)                      :: k             !< Stage index (1..K_total).
   integer(I4P),                intent(in)                      :: K_total       !< Forest-wide stage count for this step.
   real(R8P),                   intent(in)                      :: dt            !< Timestep size from the forest.
   class(realm_object),         intent(inout), optional, target :: realm(:)      !< Sibling realms (FNL parity only).
   class(flux_register_object), intent(inout), optional         :: flux_register !< Forest's flux register for FV reflux.

   if (present(realm)) continue
   if (present(flux_register)) then
      call self%compute_residuals(q=self%rk%q_rk(:,:,:,:,:,k), dq=self%dq, s=k, flux_register=flux_register)
   else
      call self%compute_residuals(q=self%rk%q_rk(:,:,:,:,:,k), dq=self%dq, s=k)
   endif
   if (self%ib%solids_number>0) then
      call self%rk%assign_stage(field=self%adam%field, s=k, q=self%dq, phi=self%ib%phi)
   else
      call self%rk%assign_stage(field=self%adam%field, s=k, q=self%dq)
   endif
   endsubroutine end_stage_forest

   subroutine close_step_forest(self, dt)
   !< Per-step epilogue on the multi-realm path: q assembly, BC, div-clean,
   !< residual save, coil source refresh, time advance, progress print.
   !<
   !< Mirrors the post-substage tail of `integrate_rk_ssp`:
   !< `rk%update_q + update_q_BC + save_residuals + compute_coils_current +
   !< impose_div_free + external_fields add`, plus the post-loop time
   !< advance + progress print that `advance_one_step_forest` does inline.
   class(prism_cpu_object), intent(inout) :: self !< The realm.
   real(R8P),               intent(in)    :: dt   !< Timestep size from the forest.

   associate(dt_unused => dt) ! dt is the global reduction; self%time%dt is the local capped value
   end associate
   if (self%ib%solids_number>0) then
      call self%rk%update_q(field=self%adam%field, dt=self%time%dt, phi=self%ib%phi, q=self%q)
      !call self%update_q_BC(dt=self%time%dt, phi=self%ib%phi)
   else
      call self%rk%update_q(field=self%adam%field, dt=self%time%dt, q=self%q, dq=self%dq)
      !call self%update_q_BC(dt=self%time%dt)
      call self%save_residuals
   endif
   call self%apply_fWL_correction(q=self%q)
   call self%compute_coils_current(q=self%q)
   call self%impose_div_free
   if (self%external_fields%ef_type/=EF_TYPE_NONE) &
      call self%external_fields%add_external_fields(field=self%adam%field, grid=self%adam%grid, &
                                                      time=self%time%time, dt=self%time%dt, q=self%q)
   self%time%time = self%time%time + self%time%dt
   call self%time%print_progress(nodes_number=self%adam%tree%nodes_number)
   ! Clear active-stage marker so any `fill_seam_from_peer_forest`
   ! query between steps reads from self%q (committed state) rather than
   ! a stale q_rk slot from the previous step's last stage.
   self%stage_active = 0_I4P
   endsubroutine close_step_forest

   subroutine fill_seam_from_peer_forest(self, peer, p_idx)
   !< Fill THIS realm's seam ghost cells for peer slot `p_idx` from
   !< `peer`'s interior cells. Single-TBP roundtrip: walks the receiver's
   !< sorted seam-map slice for the peer, and for each row copies
   !<     peer's INTERIOR(b_send, i_send, j_send, k_send, v)
   !<  →  self's GHOST (b_recv, i_recv, j_recv, k_recv, v)
   !< on the active conservative buffer (`q` when `stage_active == 0`, or
   !< `rk%q_rk(:,...,stage_active)` when active). The PEER's `stage_active`
   !< controls its read source; SELF's `stage_active` controls its write
   !< target.
   !<
   !< `select type(peer)` is confined here — the only place that needs
   !< to know the peer's concrete backend type. Cross-backend mixing
   !< (CPU realm reading from FNL peer) is not supported and triggers
   !< `error_stop`.
   class(prism_cpu_object), intent(inout)         :: self
   class(realm_object),     intent(in),    target :: peer
   integer(I4P),            intent(in)            :: p_idx
   integer(I4P)                                   :: c, row_start, row_count, row
   integer(I4P)                                   :: b_send, i_send, j_send, k_send
   integer(I4P)                                   :: b_recv, i_recv, j_recv, k_recv

   row_start = self%adam%maps%seam_local_peer_row_start(p_idx)
   row_count = self%adam%maps%seam_local_peer_row_count(p_idx)
   select type (peer)
   class is (prism_cpu_object)
      associate(rows => self%adam%maps%seam_local_map_ghost_cell)
      do c = 1_I4P, row_count
         row = row_start + c - 1_I4P
         b_send = rows(row, 2) ; i_send = rows(row, 4) ; j_send = rows(row, 5) ; k_send = rows(row, 6)
         b_recv = rows(row, 3) ; i_recv = rows(row, 7) ; j_recv = rows(row, 8) ; k_recv = rows(row, 9)
         if (self%stage_active > 0_I4P) then
            if (peer%stage_active > 0_I4P) then
               self%rk%q_rk(:, i_recv, j_recv, k_recv, b_recv, self%stage_active) = &
                  peer%rk%q_rk(:, i_send, j_send, k_send, b_send, peer%stage_active)
            else
               self%rk%q_rk(:, i_recv, j_recv, k_recv, b_recv, self%stage_active) = &
                  peer%q(:, i_send, j_send, k_send, b_send)
            endif
         else
            if (peer%stage_active > 0_I4P) then
               self%q(:, i_recv, j_recv, k_recv, b_recv) = &
                  peer%rk%q_rk(:, i_send, j_send, k_send, b_send, peer%stage_active)
            else
               self%q(:, i_recv, j_recv, k_recv, b_recv) = peer%q(:, i_send, j_send, k_send, b_send)
            endif
         endif
      enddo
      endassociate
   class default
      call mpih%error_stop(msg='prism_cpu_object%fill_seam_from_peer_forest: peer realm is not prism_cpu_object')
   end select
   endsubroutine fill_seam_from_peer_forest

   subroutine apply_reflux_to_stage_forest(self, stage, dt, flux_register)
   !< PRISM-CPU override of the Berger-Colella reflux correction TBP.
   !<
   !< **α.r1 cadence: end-of-step barrier.** The body fires exactly once per
   !< realm per step, at `stage == self%rk%nrk` (the realm's own final RK
   !< substage). Earlier stages return immediately as a no-op. This pairs
   !< with M4 (`accumulate_seam_fluxes_fv` gated on `stage_idx == rk%nrk`):
   !< the realm accumulates its end-of-step face flux into
   !< `F_coarse(:,:,1)` / `F_fine_sum(:,:,1)` at its final substage, and the
   !< same final substage immediately applies the Berger-Colella correction
   !< from those accumulators. This is the AMReX `Reflux` cadence; same-K
   !< and asymmetric-K both reduce to the same expression because the
   !< accumulators hold the final, committed face flux.
   !<
   !< For each face in `flux_register` where `face%coarse_realm == self%realm_index`,
   !< writes the per-cell correction
   !<     q_rk(:, seam_cell, b_coarse, stage) +=
   !<         ark(stage) * sign * (dt / dx_coarse) * (F_coarse(:,:,1) - F_fine_sum(:,:,1))
   !< on the matching one-cell-thick seam slice (i ∈ {1, ni} for x-normal
   !< faces, etc.). The third-axis index is **hardcoded to 1** under α.r1
   !< (register collapsed by M2); the `stage` parameter is retained on the
   !< signature for API stability and used only for the cadence gate and
   !< as the write target into `self%rk%q_rk(:, ..., stage)`.
   !<
   !< Empty-register / out-of-range guards short-circuit to a no-op, so a
   !< realm with no seam faces returns cleanly.
   class(prism_cpu_object),     intent(inout) :: self          !< The realm.
   integer(I4P),                intent(in)    :: stage         !< Integrator stage 1..K_total.
   real(R8P),                   intent(in)    :: dt            !< Time step.
   class(flux_register_object), intent(in)    :: flux_register !< Forest's flux register.
   integer(I4P)                               :: f, c, c0
   integer(I4P)                               :: axis, sgn
   integer(I4P)                               :: i_coarse, j_coarse, k_coarse
   integer(I4P)                               :: ni_, nj_, nk_
   real(R8P)                                  :: dx_coarse, scale_

   if (.not. flux_register%is_initialized_)    return
   if (flux_register%nfaces == 0_I4P)          return
   if (.not. allocated(flux_register%face))    return
   if (stage < 1_I4P .or. stage > self%rk%nrk) return
   ! α.r1 end-of-step gate: reflux fires once per realm per step at the
   ! realm's own final RK substage. Earlier stages no-op.
   if (stage /= self%rk%nrk)         return

   ni_ = self%adam%grid%ni
   nj_ = self%adam%grid%nj
   nk_ = self%adam%grid%nk

   do f = 1_I4P, flux_register%nfaces
      associate(face_f => flux_register%face(f))
      if (face_f%coarse_realm /= self%realm_index) cycle
      ! Issue #28 D4: only the coarse block's OWNER rank applies (and prints)
      ! this face — `coarse_block` is an owner-rank-LOCAL field slot, and on
      ! any other rank it aliases an unrelated local block. The accumulators
      ! are complete on every rank after reduce_fine_sums (#28 D3), so the
      ! owner applies the full correction exactly once.
      if (face_f%coarse_rank /= mpih%myrank)      cycle
      if (.not. allocated(face_f%F_coarse))       cycle
      if (.not. allocated(face_f%F_fine_sum))     cycle

      call face_axis_sign(face_f%coarse_face, axis, sgn)
      if (axis == 0_I4P) cycle  ! malformed face_code; defensive.

      dx_coarse = self%adam%field%dxyz(axis, face_f%coarse_block)
      if (dx_coarse <= 0._R8P) cycle  ! defensive (uninitialised block geometry)
      ! Register-level diagnostic (issue #23 R3): the flux mismatch this face is about
      ! to correct with. Format matched with the FNL twin so the two backends'
      ! register contents are directly comparable from the logs.
      call mpih%print_message('reflux face '//trim(str(f, .true.))//' coarse_block '//                     &
                              trim(str(face_f%coarse_block, .true.))//' max|F_coarse-F_fine_sum| = '//     &
                              trim(str(maxval(abs(face_f%F_coarse(:,:,1) - face_f%F_fine_sum(:,:,1))))))
      ! α.r1 end-of-step Berger-Colella correction applied DIRECTLY to the
      ! committed solution `self%q` (this TBP runs AFTER close_step_forest's
      ! update_q). The full step weight is `dt/dx_coarse` — NOT a stage RK
      ! coefficient: the AMReX convention reflux corrects the committed q once
      ! per step with the end-of-step flux mismatch, decoupled from the stage
      ! integration. (The earlier `ark(stage)*dt/dx` form was wrong on two
      ! counts: `ark` is allocated only for the low-storage RK family, never for
      ! the SSP family the staged forest path uses — so it silently no-op'd for
      ! every SSP run — and writing q_rk(stage) pre-update_q would entangle the
      ! correction with the stage beta weight.)
      scale_ = real(sgn, R8P) * dt / dx_coarse

      do c = 1_I4P, face_f%nface_cells
         c0 = c - 1_I4P
         select case (axis)
         case (1_I4P)  ! x-normal face: i fixed; tangentials (j, k) walk (mod nj, div nj).
            i_coarse = merge(ni_, 1_I4P, sgn > 0_I4P)
            j_coarse = 1_I4P + mod(c0, nj_)
            k_coarse = 1_I4P + c0 / nj_
         case (2_I4P)  ! y-normal face: j fixed; tangentials (i, k) walk (mod ni, div ni).
            i_coarse = 1_I4P + mod(c0, ni_)
            j_coarse = merge(nj_, 1_I4P, sgn > 0_I4P)
            k_coarse = 1_I4P + c0 / ni_
         case (3_I4P)  ! z-normal face: k fixed; tangentials (i, j) walk (mod ni, div ni).
            i_coarse = 1_I4P + mod(c0, ni_)
            j_coarse = 1_I4P + c0 / ni_
            k_coarse = merge(nk_, 1_I4P, sgn > 0_I4P)
         case default
            cycle
         end select

         ! α.r1: third axis hardcoded to 1 (register collapsed by M2).
         self%q(:, i_coarse, j_coarse, k_coarse, face_f%coarse_block) = &
            self%q(:, i_coarse, j_coarse, k_coarse, face_f%coarse_block) &
            + scale_ * (face_f%F_coarse(:, c, 1) - face_f%F_fine_sum(:, c, 1))
      enddo
      end associate
   enddo
   endsubroutine apply_reflux_to_stage_forest

   subroutine post_step_forest(self, dt, t, it, do_save_state, do_save_residuals, do_save_restart, do_amr, realm)
   !< Run PRISM-CPU's per-timestep post-step work: state IO, energy
   !< diagnostics, divergence diagnostics.
   !<
   !< Invoked by forest%post_step. v1 implementation is the verbatim post-
   !< step block formerly inline in `simulate` — every action runs every
   !< step, since today's cadence is enforced inside the save_* routines
   !< themselves (e.g. save_simulation_data honours `io%it_save`). The
   !< `do_*` flags are signature-only for now: when the forest takes over
   !< cadence the flags will gate the individual calls. For now they are
   !< accepted but unused, preserving present-day behavior bit-for-bit.
   !<
   !< `dt`, `t`, `it` are not consumed by the current body; they are on
   !< the contract so the forest can supply them once it owns time-state
   !< (today they are still read from the `time` module singleton).
   !<
   !< `realm` is accepted on contract parity with FNL; unused on CPU.
   class(prism_cpu_object), intent(inout)                   :: self              !< The realm.
   real(R8P),               intent(in)                      :: dt                !< Timestep size just advanced.
   real(R8P),               intent(in)                      :: t                 !< Simulation time after the advance.
   integer(I4P),            intent(in)                      :: it                !< Iteration index after the advance.
   logical,                 intent(in),    optional         :: do_save_state     !< Save state output this step.
   logical,                 intent(in),    optional         :: do_save_residuals !< Save residuals output this step.
   logical,                 intent(in),    optional         :: do_save_restart   !< Save restart dump this step.
   logical,                 intent(in),    optional         :: do_amr            !< Run AMR update this step.
   class(realm_object),     intent(inout), optional, target :: realm(:)          !< Sibling realms.
   real(R8P)                                                :: max_div_D         !< Maximum of divergence of D field.
   real(R8P)                                                :: max_div_B         !< Maximum of divergence of B
   real(R8P)                                                :: max_div_J         !< Maximum of divergence of J field.

   if (self%io%save_memory_status) then
      call save_memory_status(file_name='memory_cpu-'//mpih%myrankstr//'.dat', tag=str(self%time%it,.true.))
   endif
   if (mod(self%time%it,self%amr%frequency)==0) then
      ! call mpih%barrier(tictoc=.true.)
      !call self%amr_update
      ! call mpih%barrier(tictoc=.true.)
   endif
   associate(hs => self%fdv_half_stencil)
   call self%save_simulation_data
   call self%update_ghost(q=self%q)
   ! Issue #31: self%update_ghost fills the intra-realm ghosts + physical BCs but NOT the
   ! inter-realm seam (the forest owns that). The divergence diagnostic below reads an
   ! s1-deep stencil that, at seam-adjacent cells, needs the seam ghosts — which
   ! update_ghost leaves stale/BC-filled, producing a spurious div(B) at the seam skin
   ! (the evolved field is div-free; verified). Re-establish the inter-realm seam from
   ! peers on the committed q (stage_active==0 → fill writes self%q) before the diagnostic.
   if (present(realm) .and. allocated(self%adam%maps%seam_local_map_ghost_cell) .and. &
       allocated(self%adam%maps%seam_local_peer_realm)) then
      block
         integer(I4P) :: p_s
         do p_s = 1_I4P, int(size(self%adam%maps%seam_local_peer_realm), I4P)
            call self%fill_seam_from_peer_forest(peer=realm(self%adam%maps%seam_local_peer_realm(p_s)), p_idx=p_s)
         enddo
      endblock
   endif
   call self%compute_energy
   if (self%grms%do_save_history) call self%compute_grms
   !call self%save_energy_error
   call self%save_energy_history
   call self%save_grms_history
   call self%compute_divergence(hs=hs, ivar=1, q=self%q, divergence=self%divergence(1,:,:,:,:))
   call self%compute_divergence(hs=hs, ivar=4, q=self%q, divergence=self%divergence(2,:,:,:,:))
   call self%compute_divergence(hs=hs, ivar=self%physics%var_Jx, q=self%q, divergence=self%divergence(3,:,:,:,:))
   endassociate
   call compute_max_divergence_outside_absorbing_layers(self=self, hs=self%fdv_half_stencils(1), max_div_D=max_div_D, &
                                                        max_div_B=max_div_B, max_div_J=max_div_J)
   call self%save_divergence_history(div_D=max_div_D, div_B=max_div_B, div_J=max_div_J)
   endsubroutine post_step_forest

   subroutine is_done_forest(self, done)
   !< Decide whether this realm has reached its local termination criterion.
   !<
   !< Invoked by forest%is_done. PRISM-CPU override: matches the legacy
   !< condition inline in simulate — terminate when either the simulated
   !< time has reached self%time%time_max (time-driven mode, it_max <= 0) or
   !< the iteration count has reached self%time%it_max (iteration-driven mode).
   !< Today the test reads time-state from the `time` module singleton,
   !< so `self` is unused; once the forest takes over time bookkeeping
   !< the body will consume self%time%* instead.
   class(prism_cpu_object), intent(in)  :: self !< The realm.
   logical,                 intent(out) :: done !< True if this realm is done evolving.

   associate(self_unused => self)
   end associate
   done = ((self%time%it_max <= 0).and.(self%time%time >= self%time%time_max)).or.&
          ((self%time%it >= self%time%it_max).and.(self%time%it_max > 0))
   endsubroutine is_done_forest

   subroutine finalize_forest(self)
   !< Shut this realm down: final state dump, close residuals/energy/
   !< divergence history files, post-loop divergence diagnostics, finalize
   !< MPI handler.
   !<
   !< Invoked by forest%finalize. v1 implementation is the verbatim post-
   !< loop block formerly inline in `simulate`. Behavior unchanged.
   class(prism_cpu_object), intent(inout) :: self !< The realm.
   logical                                :: is_open !< Guard close() against unopened units.
   !call self%compute_energy_error
   call self%save_simulation_data
   call self%io%close_file_residuals
   !call self%save_energy_error(is_to_close=.true.)
   !call mpih%print_message('Initial/final energy of D field: '//trim(str(sqrt(self%energy_D(1))))//' '//&
   !                                                                  trim(str(sqrt(self%energy_D(size(self%energy_D))))))
   !call mpih%print_message('Initial/final energy of B field: '//trim(str(sqrt(self%energy_B(1))))//' '//&
   !                                                                  trim(str(sqrt(self%energy_D(size(self%energy_B))))))
   !call mpih%print_message('RMS Error of D field: '//trim(str(self%rms_energy_error_D)))
   !call mpih%print_message('RMS Error of B field: '//trim(str(self%rms_energy_error_B)))
   ! The final post_step_forest already writes the terminal energy/divergence row when
   ! the stop criterion is hit. finalize_forest must only close the files, otherwise the
   ! last iteration is appended twice in the CPU backend.
   if (mpih%myrank == 0) then
      inquire(unit=self%io%energy_history_unit, opened=is_open)
      if (is_open) close(self%io%energy_history_unit)
      if (self%grms%history_unit > 0_I4P) then
         inquire(unit=self%grms%history_unit, opened=is_open)
         if (is_open) close(self%grms%history_unit)
      endif
      inquire(unit=self%io%divergence_history_unit, opened=is_open)
      if (is_open) close(self%io%divergence_history_unit)
   endif

   ! NB: MPI_FINALIZE is NOT called here — it is process-global and runs once via
   ! forest%finalize -> finalize_mpi_forest after ALL realms finish.
   endsubroutine finalize_forest

   ! numerical methods
   subroutine compute_dt(self)
   !< Compute the global stability-limited dt and store it on `self%time%dt`.
   !<
   !< Body delegates the local computation to compute_local_dt_forest
   !< (orchestrator contract method), then performs the legacy
   !< MPI_ALLREDUCE on MPI_COMM_WORLD for backward compatibility with
   !< simulate. The forest's `compute_global_dt`
   !< will perform its own reduction, possibly on a per-realm sub-comm;
   !< the redundancy disappears once the legacy compute_dt is retired.
   class(prism_cpu_object), intent(inout) :: self !< The equation.

   call self%compute_local_dt_forest(dt_local=self%time%dt)
   call MPI_ALLREDUCE(MPI_IN_PLACE, self%time%dt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, mpih%error)
   endsubroutine compute_dt

   subroutine compute_energy(self)
   !< Compute energy.
   class(prism_cpu_object), intent(inout) :: self          !< The equation.
   real(R8P)                              :: energy_D      !< Energy of D field.
   real(R8P)                              :: energy_B      !< Energy of B field.
   real(R8P)                              :: coil_power    !< Coil power.
   real(R8P)                              :: poynting_flux !< Total Poynting flux from boundary.

   call compute_e(ivar=VAR_DX, energy=energy_D)
   call compute_e(ivar=VAR_BX, energy=energy_B)
   if (self%coil%total_coils_number > 0_I4P) then
      call compute_coil_power(ivar=self%physics%var_Jx, coil_power=coil_power)
   else
      coil_power = 0.0_R8P
   endif
   call compute_Poynting_flux(Poynting_flux=Poynting_flux)
   call MPI_ALLREDUCE(MPI_IN_PLACE, energy_D, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, energy_B, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, coil_power, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, Poynting_flux, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   if (allocated(self%energy_D).and.allocated(self%energy_B) &
      .and.allocated(self%coil_power).and.allocated(self%Poynting_flux)) then
      self%energy_D      = [self%energy_D,      energy_D     ]
      self%energy_B      = [self%energy_B,      energy_B     ]
      self%coil_power    = [self%coil_power,    coil_power   ]
      self%poynting_flux = [self%poynting_flux, poynting_flux]
   else
      allocate(self%energy_D(     1:self%time%it))
      allocate(self%energy_B(     1:self%time%it))
      allocate(self%coil_power(   1:self%time%it))
      allocate(self%poynting_flux(1:self%time%it))
      self%energy_D      = energy_D
      self%energy_B      = energy_B
      self%coil_power    = coil_power
      self%poynting_flux = poynting_flux
   endif
   contains
      subroutine compute_e(ivar, energy)
      !< Compute electric/magnetic energy of vector field starting from ivar.
      integer(I4P), intent(in)  :: ivar    !< Starting position of vector field.
      real(R8P),    intent(out) :: energy  !< Energy of the vector field starting from ivar.
      integer(I4P)              :: i,j,k,b !< Counter.
      real(R8P)                 :: const   !< Costant for the energy computation.

      if (ivar==VAR_DX) then
         const = EPS0
      elseif (ivar==VAR_BX) then
         const = MU0
      endif
      energy = 0.0_R8P
      associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number, &
                dx=>self%adam%field%dxyz(1,:), dy=>self%adam%field%dxyz(2,:), dz=>self%adam%field%dxyz(3,:))
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         energy = energy + 0.5_R8P * (self%q(ivar  ,i,j,k,b)*self%q(ivar  ,i,j,k,b) + &
                                      self%q(ivar+1,i,j,k,b)*self%q(ivar+1,i,j,k,b) + &
                                      self%q(ivar+2,i,j,k,b)*self%q(ivar+2,i,j,k,b))/const*(dx(b)*dy(b)*dz(b))
      enddo
      enddo
      enddo
      enddo
      endassociate
      endsubroutine compute_e

      subroutine compute_coil_power(ivar, coil_power)
      !< Compute coil power of vector field starting from ivar.
      integer(I4P), intent(in)  :: ivar        !< Starting position of vector field.
      real(R8P),    intent(out) :: coil_power  !< Coil power of the vector field.
      integer(I4P)              :: i,j,k,b     !< Counter.

      coil_power = 0.0_R8P
      associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,         &
                dx=>self%adam%field%dxyz(1,:), dy=>self%adam%field%dxyz(2,:), dz=>self%adam%field%dxyz(3,:), &
                J_vec=>self%coil%J_vec)
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         if (abs(minval(J_vec(:,i,j,k,b,:))) >= 1E-12_R8P) then
            coil_power = coil_power - (self%q(VAR_DX  ,i,j,k,b)*self%q(ivar  ,i,j,k,b) + &
                                       self%q(VAR_DX+1,i,j,k,b)*self%q(ivar+1,i,j,k,b) + &
                                       self%q(VAR_DX+2,i,j,k,b)*self%q(ivar+2,i,j,k,b))/EPS0*(dx(b)*dy(b)*dz(b))
         endif
      enddo
      enddo
      enddo
      enddo
      endassociate
      endsubroutine compute_coil_power

      subroutine compute_Poynting_flux(poynting_flux)
      !< Compute Poynting flux.
      real(R8P),    intent(out) :: poynting_flux  !< Power irradiated outside computational domain.
      integer(I4P)              :: i,j,k,b,v      !< Counter.
      real(R8P)                 :: q_boundary(6)  !< Variables at boundary for the Poynting flux computation.
      real(R8P)                 :: n(3)           !< Boundary normal direction

      associate(ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc,blocks_number=>self%blocks_number,        &
                s=>self%fdv_half_stencils(1), dx=>self%adam%field%dxyz(1,:), dy=>self%adam%field%dxyz(2,:), &
                dz=>self%adam%field%dxyz(3,:))
      poynting_flux = 0.0_R8P
      !Faccia -x
      n = [-1.0_R8P, 0.0_R8P, 0.0_R8P]
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
         do v=1, 6
            call compute_reconstruction_r_fd_centered(s=s,q=self%q(v,1-s:s,j,k,b),qr=q_boundary(v))
         enddo
         Poynting_flux = Poynting_flux + dotproduct(crossproduct(q_boundary(1:3), q_boundary(4:6))/MU0,n)*(dy(b)*dz(b))
      enddo
      enddo
      enddo
      !Faccia +x
      n = [1.0_R8P, 0.0_R8P, 0.0_R8P]
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
         do v=1, 6
            call compute_reconstruction_r_fd_centered(s=s,q=self%q(v,ni+1-s:ni+s,j,k,b),qr=q_boundary(v))
         enddo
         Poynting_flux = Poynting_flux + dotproduct(crossproduct(q_boundary(1:3),q_boundary(4:6))/MU0,n)*(dy(b)*dz(b))
      enddo
      enddo
      enddo
      !Faccia -y
      n = [0.0_R8P, -1.0_R8P, 0.0_R8P]
      do b=1, blocks_number
      do k=1, nk
      do i=1, ni
         do v=1, 6
            call compute_reconstruction_r_fd_centered(s=s,q=self%q(v,i,1-s:s,k,b),qr=q_boundary(v))
         enddo
         Poynting_flux = Poynting_flux + dotproduct(crossproduct(q_boundary(1:3),q_boundary(4:6))/MU0,n)*(dx(b)*dz(b))
      enddo
      enddo
      enddo
      !Faccia +y
      n = [0.0_R8P, 1.0_R8P, 0.0_R8P]
      do b=1, blocks_number
      do k=1, nk
      do i=1, ni
         do v=1, 6
            call compute_reconstruction_r_fd_centered(s=s,q=self%q(v,i,nj+1-s:nj+s,k,b),qr=q_boundary(v))
         enddo
         Poynting_flux = Poynting_flux + dotproduct(crossproduct(q_boundary(1:3),q_boundary(4:6))/MU0,n)*(dx(b)*dz(b))
      enddo
      enddo
      enddo
      !Faccia -z
      n = [0.0_R8P, 0.0_R8P, -1.0_R8P]
      do b=1, blocks_number
      do j=1, nj
      do i=1, ni
         do v=1, 6
            call compute_reconstruction_r_fd_centered(s=s,q=self%q(v,i,j,1-s:s,b),qr=q_boundary(v))
         enddo
         Poynting_flux = Poynting_flux + dotproduct(crossproduct(q_boundary(1:3),q_boundary(4:6))/MU0,n)*(dx(b)*dy(b))
      enddo
      enddo
      enddo
      !Faccia +z
      n = [0.0_R8P, 0.0_R8P, 1.0_R8P]
      do b=1, blocks_number
      do j=1, nj
      do i=1, ni
         do v=1, 6
            call compute_reconstruction_r_fd_centered(s=s,q=self%q(v,i,j,nk+1-s:nk+s,b),qr=q_boundary(v))
         enddo
         Poynting_flux = Poynting_flux + dotproduct(crossproduct(q_boundary(1:3),q_boundary(4:6))/MU0,n)*(dx(b)*dy(b))
      enddo
      enddo
      enddo
      endassociate
   endsubroutine compute_Poynting_flux
   endsubroutine compute_energy

   subroutine compute_energy_error(self)
   !< Compute energy error.
   class(prism_cpu_object), intent(inout) :: self       !< The equation.
   real(R8P)                              :: energy_D0  !< Initial energy of D field.
   real(R8P)                              :: energy_B0  !< Initial energy of B field.
   real(R8P), allocatable                 :: error_D(:) !< Error (normalized) of energy of D field.
   real(R8P), allocatable                 :: error_B(:) !< Error (normalized) of energy of B field.

   if (allocated(self%energy_D).and.allocated(self%energy_B)) then
      energy_D0 = self%energy_D(1)
      energy_B0 = self%energy_B(1)

      error_D = abs(self%energy_D - energy_D0) / abs(energy_D0)
      error_B = abs(self%energy_B - energy_B0) / abs(energy_B0)

      self%rms_energy_error_D = sqrt(sum(error_D)/size(error_D))
      self%rms_energy_error_B = sqrt(sum(error_B)/size(error_B))
   endif
   endsubroutine compute_energy_error

   subroutine compute_grms(self)
   !< Compute the RMS gradient of the rotating magnetic-field amplitude over the selected domain
   !< and over its -3 dB subset.
   class(prism_cpu_object), intent(inout) :: self
   integer(I4P)                          :: b, i, j, k, s
   integer(I4P)                          :: lo_i, hi_i, lo_j, hi_j, lo_k, hi_k
   integer(I4P)                          :: b_ref, i_ref, j_ref, k_ref
   integer(I8P)                          :: cells_local_3db, cells_global_3db
   integer(I8P)                          :: cells_local_domain, cells_global_domain
   real(R8P)                             :: b_amp
   real(R8P)                             :: b_minus, b_plus
   real(R8P)                             :: best_r2_local, best_r2_global
   real(R8P)                             :: bref_candidate, bref_local, bref_global
   real(R8P)                             :: dBdx, dBdy, dBdz
   real(R8P)                             :: domain_center(3)
   real(R8P)                             :: grad2
   real(R8P)                             :: half_length, radius2
   real(R8P)                             :: measure_local_3db, measure_global_3db
   real(R8P)                             :: measure_local_domain, measure_global_domain
   real(R8P)                             :: threshold
   real(R8P)                             :: weighted_sum_local_3db, weighted_sum_global_3db
   real(R8P)                             :: weighted_sum_local_domain, weighted_sum_global_domain
   real(R8P)                             :: x, y, z
   logical                               :: use_cylinder

   use_cylinder = self%grms%use_cylindrical_region
   domain_center = [0._R8P, 0._R8P, 0._R8P]
   if (use_cylinder) domain_center = self%grms%center
   half_length = 0.5_R8P * self%grms%length
   radius2 = self%grms%radius * self%grms%radius
   best_r2_local = huge(1.0_R8P)
   bref_local = 0.0_R8P
   b_ref = 0_I4P ; i_ref = 1_I4P ; j_ref = 1_I4P ; k_ref = 1_I4P
   do b = 1, self%blocks_number
      call get_valid_window(self=self, hs=self%fdv_half_stencils(1), b=b, lo_i=lo_i, hi_i=hi_i, lo_j=lo_j, hi_j=hi_j, &
                            lo_k=lo_k, hi_k=hi_k)
      if (lo_i > hi_i .or. lo_j > hi_j .or. lo_k > hi_k) cycle
      do k = lo_k, hi_k
         z = self%adam%field%z_cell(k,b)
         do j = lo_j, hi_j
            y = self%adam%field%y_cell(j,b)
            do i = lo_i, hi_i
               x = self%adam%field%x_cell(i,b)
               if (.not. is_inside_selected_domain(x=x, y=y, z=z)) cycle
               if (reference_distance2(x=x, y=y, z=z) < best_r2_local) then
                  best_r2_local = reference_distance2(x=x, y=y, z=z)
                  b_ref = b ; i_ref = i ; j_ref = j ; k_ref = k
               endif
            enddo
         enddo
      enddo
   enddo

   if (best_r2_local < huge(1.0_R8P)) then
      bref_local = sqrt(self%q(VAR_BX,i_ref,j_ref,k_ref,b_ref)**2 + self%q(VAR_BY,i_ref,j_ref,k_ref,b_ref)**2)
   endif
   call MPI_ALLREDUCE(best_r2_local, best_r2_global, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, mpih%error)
   bref_candidate = 0.0_R8P
   if (abs(best_r2_local - best_r2_global) <= 1.0E-12_R8P * max(1.0_R8P, abs(best_r2_global))) bref_candidate = bref_local
   bref_global = bref_candidate
   call MPI_ALLREDUCE(MPI_IN_PLACE, bref_global, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, mpih%error)

   threshold = GRMS_3DB_RATIO * bref_global
   weighted_sum_local_domain = 0.0_R8P
   weighted_sum_local_3db = 0.0_R8P
   measure_local_domain = 0.0_R8P
   measure_local_3db = 0.0_R8P
   cells_local_domain = 0_I8P
   cells_local_3db = 0_I8P
   do b = 1, self%blocks_number
      call get_valid_window(self=self, hs=self%fdv_half_stencils(1), b=b, lo_i=lo_i, hi_i=hi_i, lo_j=lo_j, hi_j=hi_j, &
                            lo_k=lo_k, hi_k=hi_k)
      if (lo_i > hi_i .or. lo_j > hi_j .or. lo_k > hi_k) cycle
      do k = lo_k, hi_k
         z = self%adam%field%z_cell(k,b)
         do j = lo_j, hi_j
            y = self%adam%field%y_cell(j,b)
            do i = lo_i, hi_i
               x = self%adam%field%x_cell(i,b)
               if (.not. is_inside_selected_domain(x=x, y=y, z=z)) cycle
               b_amp = sqrt(self%q(VAR_BX,i,j,k,b)**2 + self%q(VAR_BY,i,j,k,b)**2)
               dBdx = 0.0_R8P
               dBdy = 0.0_R8P
               dBdz = 0.0_R8P
               do s = 1, self%fdv_half_stencils(1)
                  b_plus = sqrt(self%q(VAR_BX,i+s,j,k,b)**2 + self%q(VAR_BY,i+s,j,k,b)**2)
                  b_minus = sqrt(self%q(VAR_BX,i-s,j,k,b)**2 + self%q(VAR_BY,i-s,j,k,b)**2)
                  dBdx = dBdx + FD1_CC(s,self%fdv_half_stencils(1)) * (b_plus - b_minus) / self%adam%field%dxyz(1,b)
                  b_plus = sqrt(self%q(VAR_BX,i,j+s,k,b)**2 + self%q(VAR_BY,i,j+s,k,b)**2)
                  b_minus = sqrt(self%q(VAR_BX,i,j-s,k,b)**2 + self%q(VAR_BY,i,j-s,k,b)**2)
                  dBdy = dBdy + FD1_CC(s,self%fdv_half_stencils(1)) * (b_plus - b_minus) / self%adam%field%dxyz(2,b)
                  b_plus = sqrt(self%q(VAR_BX,i,j,k+s,b)**2 + self%q(VAR_BY,i,j,k+s,b)**2)
                  b_minus = sqrt(self%q(VAR_BX,i,j,k-s,b)**2 + self%q(VAR_BY,i,j,k-s,b)**2)
                  dBdz = dBdz + FD1_CC(s,self%fdv_half_stencils(1)) * (b_plus - b_minus) / self%adam%field%dxyz(3,b)
               enddo
               grad2 = dBdx*dBdx + dBdy*dBdy + dBdz*dBdz
               weighted_sum_local_domain = weighted_sum_local_domain + grad2 * product(self%adam%field%dxyz(:,b))
               measure_local_domain = measure_local_domain + product(self%adam%field%dxyz(:,b))
               cells_local_domain = cells_local_domain + 1_I8P
               if (b_amp >= threshold) then
                  weighted_sum_local_3db = weighted_sum_local_3db + grad2 * product(self%adam%field%dxyz(:,b))
                  measure_local_3db = measure_local_3db + product(self%adam%field%dxyz(:,b))
                  cells_local_3db = cells_local_3db + 1_I8P
               endif
            enddo
         enddo
      enddo
   enddo

   weighted_sum_global_domain = weighted_sum_local_domain
   weighted_sum_global_3db = weighted_sum_local_3db
   measure_global_domain = measure_local_domain
   measure_global_3db = measure_local_3db
   cells_global_domain = cells_local_domain
   cells_global_3db = cells_local_3db
   call MPI_ALLREDUCE(MPI_IN_PLACE, weighted_sum_global_domain, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, weighted_sum_global_3db, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, measure_global_domain, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, measure_global_3db, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, cells_global_domain, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, cells_global_3db, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, mpih%error)

   self%grms%reference_B = bref_global
   self%grms%threshold_B = threshold
   self%grms%domain_measure = measure_global_domain
   self%grms%measure_3db = measure_global_3db
   self%grms%domain_cells_number = cells_global_domain
   self%grms%cells_number_3db = cells_global_3db
   self%grms%grms_domain_B = 0.0_R8P
   self%grms%grms_3db_B = 0.0_R8P
   if (measure_global_domain > 0.0_R8P) self%grms%grms_domain_B = sqrt(weighted_sum_global_domain / measure_global_domain)
   if (measure_global_3db > 0.0_R8P) self%grms%grms_3db_B = sqrt(weighted_sum_global_3db / measure_global_3db)
   contains
      subroutine get_valid_window(self, hs, b, lo_i, hi_i, lo_j, hi_j, lo_k, hi_k)
      class(prism_cpu_object), intent(in)  :: self
      integer(I4P),            intent(in)  :: hs
      integer(I4P),            intent(in)  :: b
      integer(I4P),            intent(out) :: lo_i, hi_i, lo_j, hi_j, lo_k, hi_k

      lo_i = 1_I4P ; hi_i = self%ni
      lo_j = 1_I4P ; hi_j = self%nj
      lo_k = 1_I4P ; hi_k = self%nk
      if (allocated(self%fWLayer%ni_fWL)) then
         call exclude_stencil_contaminated_face(lo=lo_i, hi=hi_i, face_first=self%fWLayer%ni_fWL(1,b,PML_FACE_X_M), &
                                                face_last=self%fWLayer%ni_fWL(2,b,PML_FACE_X_M), is_minus=.true., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_i, hi=hi_i, face_first=self%fWLayer%ni_fWL(1,b,PML_FACE_X_P), &
                                                face_last=self%fWLayer%ni_fWL(2,b,PML_FACE_X_P), is_minus=.false., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_j, hi=hi_j, face_first=self%fWLayer%nj_fWL(1,b,PML_FACE_Y_M), &
                                                face_last=self%fWLayer%nj_fWL(2,b,PML_FACE_Y_M), is_minus=.true., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_j, hi=hi_j, face_first=self%fWLayer%nj_fWL(1,b,PML_FACE_Y_P), &
                                                face_last=self%fWLayer%nj_fWL(2,b,PML_FACE_Y_P), is_minus=.false., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_k, hi=hi_k, face_first=self%fWLayer%nk_fWL(1,b,PML_FACE_Z_M), &
                                                face_last=self%fWLayer%nk_fWL(2,b,PML_FACE_Z_M), is_minus=.true., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_k, hi=hi_k, face_first=self%fWLayer%nk_fWL(1,b,PML_FACE_Z_P), &
                                                face_last=self%fWLayer%nk_fWL(2,b,PML_FACE_Z_P), is_minus=.false., hs=hs)
      endif
      if (allocated(self%pml%ni_pml)) then
         call exclude_stencil_contaminated_face(lo=lo_i, hi=hi_i, face_first=self%pml%ni_pml(1,b,PML_FACE_X_M), &
                                                face_last=self%pml%ni_pml(2,b,PML_FACE_X_M), is_minus=.true., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_i, hi=hi_i, face_first=self%pml%ni_pml(1,b,PML_FACE_X_P), &
                                                face_last=self%pml%ni_pml(2,b,PML_FACE_X_P), is_minus=.false., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_j, hi=hi_j, face_first=self%pml%nj_pml(1,b,PML_FACE_Y_M), &
                                                face_last=self%pml%nj_pml(2,b,PML_FACE_Y_M), is_minus=.true., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_j, hi=hi_j, face_first=self%pml%nj_pml(1,b,PML_FACE_Y_P), &
                                                face_last=self%pml%nj_pml(2,b,PML_FACE_Y_P), is_minus=.false., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_k, hi=hi_k, face_first=self%pml%nk_pml(1,b,PML_FACE_Z_M), &
                                                face_last=self%pml%nk_pml(2,b,PML_FACE_Z_M), is_minus=.true., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_k, hi=hi_k, face_first=self%pml%nk_pml(1,b,PML_FACE_Z_P), &
                                                face_last=self%pml%nk_pml(2,b,PML_FACE_Z_P), is_minus=.false., hs=hs)
      endif
      endsubroutine get_valid_window

      logical function is_inside_selected_domain(x, y, z)
      real(R8P), intent(in) :: x, y, z
      real(R8P)             :: axial_distance, radial2_local

      if (.not. use_cylinder) then
         is_inside_selected_domain = .true.
         return
      endif
      axial_distance = (x - self%grms%center(1)) * self%grms%axis(1) + &
                       (y - self%grms%center(2)) * self%grms%axis(2) + &
                       (z - self%grms%center(3)) * self%grms%axis(3)
      if (abs(axial_distance) > half_length) then
         is_inside_selected_domain = .false.
         return
      endif
      radial2_local = (x - self%grms%center(1))**2 + (y - self%grms%center(2))**2 + (z - self%grms%center(3))**2 - &
                      axial_distance * axial_distance
      is_inside_selected_domain = radial2_local <= radius2
      endfunction is_inside_selected_domain

      real(R8P) function reference_distance2(x, y, z)
      real(R8P), intent(in) :: x, y, z

      reference_distance2 = (x - domain_center(1))**2 + (y - domain_center(2))**2 + (z - domain_center(3))**2
      endfunction reference_distance2
   endsubroutine compute_grms

   subroutine compute_max_divergence_outside_absorbing_layers(self, hs, max_div_D, max_div_B, max_div_J)
   !< Compute divergence maxima excluding absorbing layers and any cell whose
   !< centered-divergence stencil intersects them.
   class(prism_cpu_object), intent(in)  :: self
   integer(I4P),            intent(in)  :: hs
   real(R8P),               intent(out) :: max_div_D
   real(R8P),               intent(out) :: max_div_B
   real(R8P),               intent(out) :: max_div_J
   integer(I4P)                         :: b
   integer(I4P)                         :: lo_i, hi_i, lo_j, hi_j, lo_k, hi_k
   integer(I4P)                         :: rho_ivar

   max_div_D = 0._R8P
   max_div_B = 0._R8P
   max_div_J = 0._R8P
   rho_ivar = self%nv
   do b = 1, self%blocks_number
      lo_i = 1_I4P ; hi_i = self%ni
      lo_j = 1_I4P ; hi_j = self%nj
      lo_k = 1_I4P ; hi_k = self%nk
      if (allocated(self%fWLayer%ni_fWL)) then
         call exclude_stencil_contaminated_face(lo=lo_i, hi=hi_i, face_first=self%fWLayer%ni_fWL(1,b,PML_FACE_X_M), &
                                                face_last=self%fWLayer%ni_fWL(2,b,PML_FACE_X_M), is_minus=.true., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_i, hi=hi_i, face_first=self%fWLayer%ni_fWL(1,b,PML_FACE_X_P), &
                                                face_last=self%fWLayer%ni_fWL(2,b,PML_FACE_X_P), is_minus=.false., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_j, hi=hi_j, face_first=self%fWLayer%nj_fWL(1,b,PML_FACE_Y_M), &
                                                face_last=self%fWLayer%nj_fWL(2,b,PML_FACE_Y_M), is_minus=.true., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_j, hi=hi_j, face_first=self%fWLayer%nj_fWL(1,b,PML_FACE_Y_P), &
                                                face_last=self%fWLayer%nj_fWL(2,b,PML_FACE_Y_P), is_minus=.false., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_k, hi=hi_k, face_first=self%fWLayer%nk_fWL(1,b,PML_FACE_Z_M), &
                                                face_last=self%fWLayer%nk_fWL(2,b,PML_FACE_Z_M), is_minus=.true., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_k, hi=hi_k, face_first=self%fWLayer%nk_fWL(1,b,PML_FACE_Z_P), &
                                                face_last=self%fWLayer%nk_fWL(2,b,PML_FACE_Z_P), is_minus=.false., hs=hs)
      endif
      if (allocated(self%pml%ni_pml)) then
         call exclude_stencil_contaminated_face(lo=lo_i, hi=hi_i, face_first=self%pml%ni_pml(1,b,PML_FACE_X_M), &
                                                face_last=self%pml%ni_pml(2,b,PML_FACE_X_M), is_minus=.true., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_i, hi=hi_i, face_first=self%pml%ni_pml(1,b,PML_FACE_X_P), &
                                                face_last=self%pml%ni_pml(2,b,PML_FACE_X_P), is_minus=.false., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_j, hi=hi_j, face_first=self%pml%nj_pml(1,b,PML_FACE_Y_M), &
                                                face_last=self%pml%nj_pml(2,b,PML_FACE_Y_M), is_minus=.true., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_j, hi=hi_j, face_first=self%pml%nj_pml(1,b,PML_FACE_Y_P), &
                                                face_last=self%pml%nj_pml(2,b,PML_FACE_Y_P), is_minus=.false., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_k, hi=hi_k, face_first=self%pml%nk_pml(1,b,PML_FACE_Z_M), &
                                                face_last=self%pml%nk_pml(2,b,PML_FACE_Z_M), is_minus=.true., hs=hs)
         call exclude_stencil_contaminated_face(lo=lo_k, hi=hi_k, face_first=self%pml%nk_pml(1,b,PML_FACE_Z_P), &
                                                face_last=self%pml%nk_pml(2,b,PML_FACE_Z_P), is_minus=.false., hs=hs)
      endif
      if (lo_i > hi_i .or. lo_j > hi_j .or. lo_k > hi_k) cycle
      if (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
         ! In PIC runs rho is appended as the last state variable, so the most useful
         ! electric-field constraint monitor is max|div(D) - rho|.
         max_div_D = max(max_div_D, maxval(abs(self%divergence(1,lo_i:hi_i,lo_j:hi_j,lo_k:hi_k,b:b) - &
                                               self%q(rho_ivar,lo_i:hi_i,lo_j:hi_j,lo_k:hi_k,b:b))))
      else
         max_div_D = max(max_div_D, maxval(abs(self%divergence(1,lo_i:hi_i,lo_j:hi_j,lo_k:hi_k,b:b))))
      endif
      max_div_B = max(max_div_B, maxval(abs(self%divergence(2,lo_i:hi_i,lo_j:hi_j,lo_k:hi_k,b:b))))
      max_div_J = max(max_div_J, maxval(abs(self%divergence(3,lo_i:hi_i,lo_j:hi_j,lo_k:hi_k,b:b))))
   enddo
   endsubroutine compute_max_divergence_outside_absorbing_layers

   pure subroutine exclude_stencil_contaminated_face(lo, hi, face_first, face_last, is_minus, hs)
   !< Shrink a 1D diagnostic window so that no retained cell uses a centered
   !< divergence stencil intersecting the selected absorbing layer face.
   integer(I4P), intent(inout) :: lo
   integer(I4P), intent(inout) :: hi
   integer(I4P), intent(in)    :: face_first
   integer(I4P), intent(in)    :: face_last
   logical,      intent(in)    :: is_minus
   integer(I4P), intent(in)    :: hs

   if (face_first <= 0_I4P .or. face_last <= 0_I4P) return

   if (is_minus) then
      lo = max(lo, face_last + hs + 1_I4P)
   else
      hi = min(hi, face_first - hs - 1_I4P)
   endif
   endsubroutine exclude_stencil_contaminated_face

   subroutine impose_div_free(self)
   !< Impose divergence-free property.
   class(prism_cpu_object), intent(inout) :: self !< The equation.

   associate(constrained_transport_D=>self%numerics%constrained_transport_D,&
             constrained_transport_B=>self%numerics%constrained_transport_B,div_corr_var=>self%numerics%div_corr_var)
   if (constrained_transport_D.and.div_corr_var==DIV_CORR_VAR_POISS) call self%impose_ct_correction(ivar=1_I4P)
   if (constrained_transport_B.and.div_corr_var==DIV_CORR_VAR_POISS) call self%impose_ct_correction(ivar=4_I4P)
   ! here should go also other corrections...
   endassociate
   endsubroutine impose_div_free

   subroutine print_elliptic_progress(label, progress_done, progress_total, last_percent, converged)
   character(*),  intent(in)    :: label        !< Progress label.
   integer(I4P),  intent(in)    :: progress_done  !< Completed smoothing sweeps.
   integer(I4P),  intent(in)    :: progress_total !< Planned smoothing sweeps.
   integer(I4P),  intent(inout) :: last_percent !< Last printed percent.
   logical,       intent(in), optional :: converged !< Force progress to 100% when converged.
   integer(I4P)                  :: progress_total_ !< Safe planned smoothing sweeps.
   integer(I4P)                  :: progress_done_  !< Safe completed smoothing sweeps.
   integer(I4P)                  :: shown_done      !< Printed completed smoothing sweeps.
   integer(I4P)                  :: progress      !< Current progress percentage.
   logical                       :: converged_    !< Convergence flag.

   progress_total_ = max(1_I4P, progress_total)
   progress_done_  = min(max(progress_done, 0_I4P), progress_total_)
   converged_ = .false.
   if (present(converged)) converged_ = converged

   shown_done = progress_done_
   if (converged_) shown_done = progress_total_
   progress = int(100._R8P * real(shown_done, R8P) / real(progress_total_, R8P), kind=I4P)
   if (converged_) progress = 100_I4P
   progress = min(100_I4P, max(0_I4P, progress))

   if (progress > last_percent) then
      call mpih%print_message(trim(label)//' progress: '//trim(str(progress,.true.))//'% ('//trim(str(shown_done,.true.))// &
                              '/'//trim(str(progress_total_,.true.))//')')
      last_percent = progress
   endif
   endsubroutine print_elliptic_progress

   subroutine simulate(self, filename)
   !< Perform the simulation: legacy single-realm entry point.
   class(prism_cpu_object), intent(inout) :: self      !< The equation.
   character(*),            intent(in)    :: filename  !< Input file name.
   logical                                :: loop_done !< Termination predicate.

   call self%initialize_forest(filename=filename)
   integration: do
      call self%compute_dt
      call self%advance_one_step_forest(dt=self%time%dt)
      call self%post_step_forest(dt=self%time%dt, t=self%time%time, it=self%time%it)
      call self%is_done_forest(done=loop_done)
      if (loop_done) exit integration
   enddo integration
   call self%finalize_forest
   call self%finalize_mpi_forest ! finalize_forest no longer finalizes MPI; single-realm legacy path
   endsubroutine simulate

   ! pointer TBP concrete implementations
   subroutine compute_residuals_fd_centered(self, q, dq, s, flux_register)
   !< Compute residuals of equation, space operator, centered finite difference schemes.
   class(prism_cpu_object),     intent(inout)                   :: self                     !< The equation.
   real(R8P),                   intent(inout)                   :: q(1:,          &
                                                                     1-self%ngc:, &
                                                                     1-self%ngc:, &
                                                                     1-self%ngc:, &
                                                                     1:)                    !< Conservative variables.
   real(R8P),                   intent(inout)                   :: dq(1:,          &
                                                                      1-self%ngc:, &
                                                                      1-self%ngc:, &
                                                                      1-self%ngc:, &
                                                                      1:)                   !< Residuals.
   integer(I4P),                intent(in),    optional         :: s                        !< Stage counter.
   class(flux_register_object), intent(inout), optional         :: flux_register            !< Flux register.
   integer(I4P)                                                 :: i,j,k,b                  !< Counter
   real(R8P)                                                    :: curlD(3), curlB(3)       !< Curl of D and B.
   real(R8P)                                                    :: gradphi(3), gradpsi(3)   !< Graident of phi and psi.
   real(R8P)                                                    :: divergenceD, divergenceB !< Divergence of D and B.
   real(R8P)                                                    :: KO_Dx_x,KO_Dx_y,KO_Dx_z  !< Buffer for KO correction, D x.
   real(R8P)                                                    :: KO_Dy_x,KO_Dy_y,KO_Dy_z  !< Buffer for KO correction, D y.
   real(R8P)                                                    :: KO_Dz_x,KO_Dz_y,KO_Dz_z  !< Buffer for KO correction, D z.
   real(R8P)                                                    :: KO_Bx_x,KO_Bx_y,KO_Bx_z  !< Buffer for KO correction, B x.
   real(R8P)                                                    :: KO_By_x,KO_By_y,KO_By_z  !< Buffer for KO correction, B y.
   real(R8P)                                                    :: KO_Bz_x,KO_Bz_y,KO_Bz_z  !< Buffer for KO correction, B z.
   real(R8P), parameter                                         :: sigma = 1000.01_R8P
   real(R8P)                                                    :: min_curlD,max_curlD
   real(R8P)                                                    :: damping_coeff            !< Optional GLM parabolic damping coefficient.

   min_curlD =  huge(1._R8P)
   max_curlD = -huge(1._R8P)

   !call self%apply_fWL_correction(q=q)
   if (present(s)) then
      call self%update_ghost(q=q, s=s)
   else
      call self%update_ghost(q=q)
   endif
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv_c=>self%nv_c,blocks_number=>self%blocks_number, &
             dxyz=>self%adam%field%dxyz,                                                                              &
             s1=>self%fdv_half_stencils(1),                                                                           &
             s4=>self%fdv_half_stencils(4),                                                                           &
             chi =>self%physics%chi, c_r=>self%physics%c_r,                                                          &
             constrained_transport_D=>self%numerics%constrained_transport_D,                                         &
             constrained_transport_B=>self%numerics%constrained_transport_B,                                          &
             var_Jx=>self%physics%var_Jx, var_Jy=>self%physics%var_Jy, var_Jz=>self%physics%var_Jz)
   if (blocks_number > 0) then
      if (self%physics%physical_model == ADIM_EM_PHYSICAL_MODEL) then !Adimensional equations
         if (self%numerics%div_corr_var == DIV_CORR_VAR_HYPER .and. constrained_transport_D .and. &
            .not.constrained_transport_B) then
            ! RHS:
            ! dD/dt = curl(B) - grad(phi) - J
            ! dB/dt = -curl(D)
            ! dphi/dt = -ch^2*div(D)

            do b=1,blocks_number
            do k=1,nk
            do j=1,nj
            do i=1,ni
               call compute_curl_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                   &
                                           q=q(VAR_DX:VAR_DZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),        &
                                           curl=curlD)
               call compute_curl_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                   &
                                           q=q(VAR_BX:VAR_BZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),        &
                                           curl=curlB)
               call compute_gradient_fd_centered(s=s1,dxyz=dxyz(1:3,b),                               &
                                                 q=q(nv_c,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),           &
                                                 gradient=gradphi)
               call compute_divergence_fd_centered(s=s1,dxyz=dxyz(1:3,b),                             &
                                                   q=q(VAR_DX:VAR_DZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),&
                                                   divergence = divergenceD)
               dq(VAR_DX,i,j,k,b) =  curlB(1) - gradphi(1) - q(var_Jx,i,j,k,b)
               dq(VAR_DY,i,j,k,b) =  curlB(2) - gradphi(2) - q(var_Jy,i,j,k,b)
               dq(VAR_DZ,i,j,k,b) =  curlB(3) - gradphi(3) - q(var_Jz,i,j,k,b)
               dq(VAR_BX,i,j,k,b) = -curlD(1)                                 
               dq(VAR_BY,i,j,k,b) = -curlD(2)                                 
               dq(VAR_BZ,i,j,k,b) = -curlD(3)                                 
               dq(nv_c,i,j,k,b)   = -(chi)**2*divergenceD
               if (c_r > 0._R8P) then
                  damping_coeff = chi / (c_r * minval(dxyz(1:3,b)))
                  dq(nv_c,i,j,k,b) = dq(nv_c,i,j,k,b) - damping_coeff * q(nv_c,i,j,k,b)
               endif
            enddo
            enddo
            enddo
            enddo
         elseif (self%numerics%div_corr_var == DIV_CORR_VAR_HYPER .and. .not.constrained_transport_D .and. &
                 constrained_transport_B) then
            ! RHS:
            ! dD/dt = curl(B) - J
            ! dB/dt = -curl(D) -grad(psi)
            ! dpsi/dt = -ch^2*div(B)

            do b=1,blocks_number
            do k=1,nk
            do j=1,nj
            do i=1,ni
               call compute_curl_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                    &
                                           q=q(VAR_DX:VAR_DZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),         &
                                           curl=curlD)
               call compute_curl_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                    &
                                           q=q(VAR_BX:VAR_BZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),         &
                                           curl=curlB)
               call compute_gradient_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                &
                                                 q=q(nv_c,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),            &
                                                 gradient=gradpsi)
               call compute_divergence_fd_centered(s=s1,dxyz=dxyz(1:3,b),                              &
                                                   q=q(VAR_BX:VAR_BZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b), &
                                                   divergence = divergenceB)
               dq(VAR_DX,i,j,k,b) =  curlB(1) - q(var_Jx,i,j,k,b)
               dq(VAR_DY,i,j,k,b) =  curlB(2) - q(var_Jy,i,j,k,b)
               dq(VAR_DZ,i,j,k,b) =  curlB(3) - q(var_Jz,i,j,k,b)
               dq(VAR_BX,i,j,k,b) = -curlD(1) - gradpsi(1)       
               dq(VAR_BY,i,j,k,b) = -curlD(2) - gradpsi(2)       
               dq(VAR_BZ,i,j,k,b) = -curlD(3) - gradpsi(3)       
               dq(nv_c,i,j,k,b)   = -(chi)**2*divergenceB
               if (c_r > 0._R8P) then
                  damping_coeff = chi / (c_r * minval(dxyz(1:3,b)))
                  dq(nv_c,i,j,k,b) = dq(nv_c,i,j,k,b) - damping_coeff * q(nv_c,i,j,k,b)
               endif
            enddo
            enddo
            enddo
            enddo
         elseif (self%numerics%div_corr_var == DIV_CORR_VAR_HYPER .and. constrained_transport_D .and. &
                  constrained_transport_B) then
            ! RHS:
            ! dD/dt = curl(B) - grad(phi) - J
            ! dB/dt = -curl(D) -grad(psi)
            ! dphi/dt = -ch^2*div(D)
            ! dpsi/dt = -ch^2*div(B)
            do b=1,blocks_number
            do k=1,nk
            do j=1,nj
            do i=1,ni
               call compute_curl_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                        &
                                           q=q(VAR_DX:VAR_DZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),             &
                                           curl=curlD)
               call compute_curl_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                        &
                                           q=q(VAR_BX:VAR_BZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),             &
                                           curl=curlB)
               call compute_gradient_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                    &
                                                 q=q(nv_c-1_I4P,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),          &
                                                 gradient=gradphi)
               call compute_gradient_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                    &
                                                 q=q(nv_c,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),                &
                                                 gradient=gradpsi)
               call compute_divergence_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                  &
                                                   q=q(VAR_DX:VAR_DZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),     &
                                                   divergence = divergenceD)
               call compute_divergence_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                  &
                                                   q=q(VAR_BX:VAR_BZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),     &
                                                   divergence = divergenceB)
               dq(VAR_DX,i,j,k,b)      =  curlB(1)  - gradphi(1) - q(var_Jx,i,j,k,b)
               dq(VAR_DY,i,j,k,b)      =  curlB(2)  - gradphi(2) - q(var_Jy,i,j,k,b)
               dq(VAR_DZ,i,j,k,b)      =  curlB(3)  - gradphi(3) - q(var_Jz,i,j,k,b)
               dq(VAR_BX,i,j,k,b)      = -curlD(1) - gradpsi(1)
               dq(VAR_BY,i,j,k,b)      = -curlD(2) - gradpsi(2)
               dq(VAR_BZ,i,j,k,b)      = -curlD(3) - gradpsi(3)
               dq(nv_c-1_I4P,i,j,k,b)  = -(chi)**2*divergenceD
               dq(nv_c,i,j,k,b)        = -(chi)**2*divergenceB
               if (c_r > 0._R8P) then
                  damping_coeff = chi / (c_r * minval(dxyz(1:3,b)))
                  dq(nv_c-1_I4P,i,j,k,b) = dq(nv_c-1_I4P,i,j,k,b) - damping_coeff * q(nv_c-1_I4P,i,j,k,b)
                  dq(nv_c,i,j,k,b)       = dq(nv_c,i,j,k,b)       - damping_coeff * q(nv_c,i,j,k,b)
               endif
            enddo
            enddo
            enddo
            enddo
         else
            ! RHS:
            ! dD/dt = curl(B) - J
            ! dB/dt = -curl(D)
            do b=1,blocks_number
            do k=1,nk
            do j=1,nj
            do i=1,ni
               call compute_curl_fd_centered(s=s1,dxyz=dxyz(1:3,b),                             &
                                          q=q(VAR_DX:VAR_DZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),   &
                                          curl=curlD)
               call compute_curl_fd_centered(s=s1,dxyz=dxyz(1:3,b),                             &
                                          q=q(VAR_BX:VAR_BZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),   &
                                          curl=curlB)
               dq(VAR_DX,i,j,k,b) =  curlB(1) - q(var_Jx,i,j,k,b)
               dq(VAR_DY,i,j,k,b) =  curlB(2) - q(var_Jy,i,j,k,b)
               dq(VAR_DZ,i,j,k,b) =  curlB(3) - q(var_Jz,i,j,k,b)
               dq(VAR_BX,i,j,k,b) = -curlD(1)                    
               dq(VAR_BY,i,j,k,b) = -curlD(2)                    
               dq(VAR_BZ,i,j,k,b) = -curlD(3)                    
            enddo
            enddo
            enddo
            enddo
         endif
      else !Dimensional equations
         if (self%numerics%div_corr_var == DIV_CORR_VAR_HYPER .and. constrained_transport_D .and. &
            .not.constrained_transport_B) then
            ! RHS:
            ! dD/dt = curl(B/MU0) - grad(phi) - J
            ! dB/dt = -curl(D/EPS0)
            ! dphi/dt = -ch^2*div(D)
            do b=1,blocks_number
            do k=1,nk
            do j=1,nj
            do i=1,ni
               call compute_curl_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                   &
                                           q=q(VAR_DX:VAR_DZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),        &
                                           curl=curlD)
               call compute_curl_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                   &
                                           q=q(VAR_BX:VAR_BZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),        &
                                           curl=curlB)
               call compute_gradient_fd_centered(s=s1,dxyz=dxyz(1:3,b),                               &
                                                 q=q(nv_c,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),           &
                                                 gradient=gradphi)
               call compute_divergence_fd_centered(s=s1,dxyz=dxyz(1:3,b),                             &
                                                   q=q(VAR_DX:VAR_DZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),&
                                                   divergence = divergenceD)
               dq(VAR_DX,i,j,k,b) =  curlB(1)/MU0 - gradphi(1) - q(var_Jx,i,j,k,b) !- sigma*C0*(KO_Dx_x+KO_Dx_y+KO_Dx_z)/16._R8P
               dq(VAR_DY,i,j,k,b) =  curlB(2)/MU0 - gradphi(2) - q(var_Jy,i,j,k,b) !- sigma*C0*(KO_Dy_x+KO_Dy_y+KO_Dy_z)/16._R8P
               dq(VAR_DZ,i,j,k,b) =  curlB(3)/MU0 - gradphi(3) - q(var_Jz,i,j,k,b) !- sigma*C0*(KO_Dz_x+KO_Dz_y+KO_Dz_z)/16._R8P
               dq(VAR_BX,i,j,k,b) = -curlD(1)/EPS0                                 !- sigma*C0*(KO_Bx_x+KO_Bx_y+KO_Bx_z)/16._R8P
               dq(VAR_BY,i,j,k,b) = -curlD(2)/EPS0                                 !- sigma*C0*(KO_By_x+KO_By_y+KO_By_z)/16._R8P
               dq(VAR_BZ,i,j,k,b) = -curlD(3)/EPS0                                 !- sigma*C0*(KO_Bz_x+KO_Bz_y+KO_Bz_z)/16._R8P
               dq(nv_c,  i,j,k,b) = -(chi*C0)**2*divergenceD
               if (c_r > 0._R8P) then
                  damping_coeff = chi * C0 / (c_r * minval(dxyz(1:3,b)))
                  dq(nv_c,i,j,k,b) = dq(nv_c,i,j,k,b) - damping_coeff * q(nv_c,i,j,k,b)
               endif
            enddo
            enddo
            enddo
            enddo
         elseif (self%numerics%div_corr_var == DIV_CORR_VAR_HYPER .and. .not.constrained_transport_D .and. &
                 constrained_transport_B) then
            ! RHS:
            ! dD/dt = curl(B/MU0)  - J
            ! dB/dt = -curl(D/EPS0) - grad(psi)
            ! dpsi/dt = -ch^2*div(B)
            do b=1,blocks_number
            do k=1,nk
            do j=1,nj
            do i=1,ni
               call compute_curl_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                    &
                                           q=q(VAR_DX:VAR_DZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),         &
                                           curl=curlD)
               call compute_curl_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                    &
                                           q=q(VAR_BX:VAR_BZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),         &
                                           curl=curlB)
               call compute_gradient_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                &
                                                 q=q(nv_c,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),            &
                                                 gradient=gradpsi)
               call compute_divergence_fd_centered(s=s1,dxyz=dxyz(1:3,b),                              &
                                                   q=q(VAR_BX:VAR_BZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b), &
                                                   divergence = divergenceB)
               dq(VAR_DX,i,j,k,b) =  curlB(1)/MU0 - q(var_Jx,i,j,k,b) !- sigma*C0*(KO_Dx_x+KO_Dx_y+KO_Dx_z)/16._R8P
               dq(VAR_DY,i,j,k,b) =  curlB(2)/MU0 - q(var_Jy,i,j,k,b) !- sigma*C0*(KO_Dy_x+KO_Dy_y+KO_Dy_z)/16._R8P
               dq(VAR_DZ,i,j,k,b) =  curlB(3)/MU0 - q(var_Jz,i,j,k,b) !- sigma*C0*(KO_Dz_x+KO_Dz_y+KO_Dz_z)/16._R8P
               dq(VAR_BX,i,j,k,b) = -curlD(1)/EPS0 - gradpsi(1)       !- sigma*C0*(KO_Bx_x+KO_Bx_y+KO_Bx_z)/16._R8P
               dq(VAR_BY,i,j,k,b) = -curlD(2)/EPS0 - gradpsi(2)       !- sigma*C0*(KO_By_x+KO_By_y+KO_By_z)/16._R8P
               dq(VAR_BZ,i,j,k,b) = -curlD(3)/EPS0 - gradpsi(3)       !- sigma*C0*(KO_Bz_x+KO_Bz_y+KO_Bz_z)/16._R8P
               dq(nv_c,  i,j,k,b) = -(chi*C0)**2*divergenceB
               if (c_r > 0._R8P) then
                  damping_coeff = chi * C0 / (c_r * minval(dxyz(1:3,b)))
                  dq(nv_c,i,j,k,b) = dq(nv_c,i,j,k,b) - damping_coeff * q(nv_c,i,j,k,b)
               endif
            enddo
            enddo
            enddo
            enddo
         elseif (self%numerics%div_corr_var == DIV_CORR_VAR_HYPER .and. constrained_transport_D .and. &
                 constrained_transport_B) then
            ! RHS:
            ! dD/dt = curl(B/MU0) - grad(phi) - J
            ! dB/dt = -curl(D/EPS0) - grad(psi)
            ! dphi/dt = -ch^2*div(D)
            ! dpsi/dt = -ch^2*div(B)
            do b=1,blocks_number
            do k=1,nk
            do j=1,nj
            do i=1,ni
               call compute_curl_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                    &
                                             q=q(VAR_DX:VAR_DZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),       &
                                             curl=curlD)
               call compute_curl_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                    &
                                             q=q(VAR_BX:VAR_BZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),       &
                                             curl=curlB)
               call compute_gradient_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                &
                                                 q=q(nv_c-1_I4P,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),      &
                                                 gradient=gradphi)
               call compute_gradient_fd_centered(s=s1,dxyz=dxyz(1:3,b),                                &
                                                 q=q(nv_c,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),            &
                                                 gradient=gradpsi)
               call compute_divergence_fd_centered(s=s1,dxyz=dxyz(1:3,b),                              &
                                                   q=q(VAR_DX:VAR_DZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b), &
                                                   divergence = divergenceD)
               call compute_divergence_fd_centered(s=s1,dxyz=dxyz(1:3,b),                              &
                                                   q=q(VAR_BX:VAR_BZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b), &
                                                   divergence = divergenceB)
               dq(VAR_DX,    i,j,k,b) =  curlB(1)/MU0  - gradphi(1) - q(var_Jx,i,j,k,b)
               dq(VAR_DY,    i,j,k,b) =  curlB(2)/MU0  - gradphi(2) - q(var_Jy,i,j,k,b)
               dq(VAR_DZ,    i,j,k,b) =  curlB(3)/MU0  - gradphi(3) - q(var_Jz,i,j,k,b)
               dq(VAR_BX,    i,j,k,b) = -curlD(1)/EPS0 - gradpsi(1)                    
               dq(VAR_BY,    i,j,k,b) = -curlD(2)/EPS0 - gradpsi(2)                    
               dq(VAR_BZ,    i,j,k,b) = -curlD(3)/EPS0 - gradpsi(3)                    
               dq(nv_c-1_I4P,i,j,k,b) = -(chi*C0)**2*divergenceD
               dq(nv_c,      i,j,k,b) = -(chi*C0)**2*divergenceB
               if (c_r > 0._R8P) then
                  damping_coeff = chi * C0 / (c_r * minval(dxyz(1:3,b)))
                  dq(nv_c-1_I4P,i,j,k,b) = dq(nv_c-1_I4P,i,j,k,b) - damping_coeff * q(nv_c-1_I4P,i,j,k,b)
                  dq(nv_c,i,j,k,b)       = dq(nv_c,i,j,k,b)       - damping_coeff * q(nv_c,i,j,k,b)
               endif
            enddo
            enddo
            enddo
            enddo
         else
            ! RHS:
            ! dD/dt = curl(B/MU0) - J
            ! dB/dt = -curl(D/EPS0)
            do b=1,blocks_number
            do k=1,nk
            do j=1,nj
            do i=1,ni
               call compute_curl_fd_centered(s=s1,dxyz=dxyz(1:3,b),                             &
                                          q=q(VAR_DX:VAR_DZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),   &
                                          curl=curlD)
               call compute_curl_fd_centered(s=s1,dxyz=dxyz(1:3,b),                             &
                                          q=q(VAR_BX:VAR_BZ,i-s1:i+s1,j-s1:j+s1,k-s1:k+s1,b),   &
                                          curl=curlB)
               dq(VAR_DX,i,j,k,b) =  curlB(1)/MU0 - q(var_Jx,i,j,k,b)
               dq(VAR_DY,i,j,k,b) =  curlB(2)/MU0 - q(var_Jy,i,j,k,b)
               dq(VAR_DZ,i,j,k,b) =  curlB(3)/MU0 - q(var_Jz,i,j,k,b)
               dq(VAR_BX,i,j,k,b) = -curlD(1)/EPS0                   
               dq(VAR_BY,i,j,k,b) = -curlD(2)/EPS0                   
               dq(VAR_BZ,i,j,k,b) = -curlD(3)/EPS0                   
            enddo
            enddo
            enddo
            enddo
         endif
      endif
   endif
   ! Issue #29: the matched-difference-operator 2:1-seam correction (E3) was removed
   ! here — measured eigenvalue-unstable (both-operator horizon diverged to 1e18),
   ! the non-SBP construction Mattsson-Carpenter 2010 / Ranocha 2019 predict. The
   ! prototype is archived (untracked/c3-fine-side-prototype-*.bak) and its reusable
   ! seam derivative primitives live in adam_fdv_operators_library. The standing seam
   ! div(B) fix is fix-class A (SBP-norm-compatible seam + weak SAT), still to build.
   endassociate
   if (self%pml%enabled) call apply_pml_fd_centered_ade(self=self, q=q, dq=dq, s=s)
   endsubroutine compute_residuals_fd_centered

   subroutine apply_pml_fd_centered_ade(self, q, dq, s)
   !< Add ADE-PML source terms to the FD-centered Maxwell residuals and build the 12 auxiliary RHSs.
   class(prism_cpu_object), intent(inout)           :: self
   real(R8P),               intent(in)              :: q(1:,          &
                                                         1-self%ngc:, &
                                                         1-self%ngc:, &
                                                         1-self%ngc:, &
                                                         1:)
   real(R8P),               intent(inout)           :: dq(1:,          &
                                                          1-self%ngc:, &
                                                          1-self%ngc:, &
                                                          1-self%ngc:, &
                                                          1:)
   integer(I4P),            intent(in), optional    :: s
   real(R8P)                                      :: inv_eps_scale
   real(R8P)                                      :: inv_mu_scale

   if (.not. self%pml%enabled) return

   call select_pml_field_scales(self=self, inv_eps_scale=inv_eps_scale, inv_mu_scale=inv_mu_scale)
   call self%rk_pml%reset_rhs()

   if (present(s)) then
      if (allocated(self%rk_pml%q_pml_x_m_rk)) then
         call apply_x_face_field_correction(self=self, q=q, q_face=self%rk_pml%q_pml_x_m_rk(:,:,:,:,:,s), dq=dq, &
                                            block_ids=self%pml%blocks_x_m, face=PML_FACE_X_M, inv_eps_scale=inv_eps_scale, &
                                            inv_mu_scale=inv_mu_scale)
         call compute_x_face_pml_rhs(self=self, q=q, q_face=self%rk_pml%q_pml_x_m_rk(:,:,:,:,:,s), &
                                     dq_face=self%rk_pml%dq_pml_x_m, block_ids=self%pml%blocks_x_m, face=PML_FACE_X_M, &
                                     inv_eps_scale=inv_eps_scale, inv_mu_scale=inv_mu_scale)
      endif
      if (allocated(self%rk_pml%q_pml_x_p_rk)) then
         call apply_x_face_field_correction(self=self, q=q, q_face=self%rk_pml%q_pml_x_p_rk(:,:,:,:,:,s), dq=dq, &
                                            block_ids=self%pml%blocks_x_p, face=PML_FACE_X_P, inv_eps_scale=inv_eps_scale, &
                                            inv_mu_scale=inv_mu_scale)
         call compute_x_face_pml_rhs(self=self, q=q, q_face=self%rk_pml%q_pml_x_p_rk(:,:,:,:,:,s), &
                                     dq_face=self%rk_pml%dq_pml_x_p, block_ids=self%pml%blocks_x_p, face=PML_FACE_X_P, &
                                     inv_eps_scale=inv_eps_scale, inv_mu_scale=inv_mu_scale)
      endif
      if (allocated(self%rk_pml%q_pml_y_m_rk)) then
         call apply_y_face_field_correction(self=self, q=q, q_face=self%rk_pml%q_pml_y_m_rk(:,:,:,:,:,s), dq=dq, &
                                            block_ids=self%pml%blocks_y_m, face=PML_FACE_Y_M, inv_eps_scale=inv_eps_scale, &
                                            inv_mu_scale=inv_mu_scale)
         call compute_y_face_pml_rhs(self=self, q=q, q_face=self%rk_pml%q_pml_y_m_rk(:,:,:,:,:,s), &
                                     dq_face=self%rk_pml%dq_pml_y_m, block_ids=self%pml%blocks_y_m, face=PML_FACE_Y_M, &
                                     inv_eps_scale=inv_eps_scale, inv_mu_scale=inv_mu_scale)
      endif
      if (allocated(self%rk_pml%q_pml_y_p_rk)) then
         call apply_y_face_field_correction(self=self, q=q, q_face=self%rk_pml%q_pml_y_p_rk(:,:,:,:,:,s), dq=dq, &
                                            block_ids=self%pml%blocks_y_p, face=PML_FACE_Y_P, inv_eps_scale=inv_eps_scale, &
                                            inv_mu_scale=inv_mu_scale)
         call compute_y_face_pml_rhs(self=self, q=q, q_face=self%rk_pml%q_pml_y_p_rk(:,:,:,:,:,s), &
                                     dq_face=self%rk_pml%dq_pml_y_p, block_ids=self%pml%blocks_y_p, face=PML_FACE_Y_P, &
                                     inv_eps_scale=inv_eps_scale, inv_mu_scale=inv_mu_scale)
      endif
      if (allocated(self%rk_pml%q_pml_z_m_rk)) then
         call apply_z_face_field_correction(self=self, q=q, q_face=self%rk_pml%q_pml_z_m_rk(:,:,:,:,:,s), dq=dq, &
                                            block_ids=self%pml%blocks_z_m, face=PML_FACE_Z_M, inv_eps_scale=inv_eps_scale, &
                                            inv_mu_scale=inv_mu_scale)
         call compute_z_face_pml_rhs(self=self, q=q, q_face=self%rk_pml%q_pml_z_m_rk(:,:,:,:,:,s), &
                                     dq_face=self%rk_pml%dq_pml_z_m, block_ids=self%pml%blocks_z_m, face=PML_FACE_Z_M, &
                                     inv_eps_scale=inv_eps_scale, inv_mu_scale=inv_mu_scale)
      endif
      if (allocated(self%rk_pml%q_pml_z_p_rk)) then
         call apply_z_face_field_correction(self=self, q=q, q_face=self%rk_pml%q_pml_z_p_rk(:,:,:,:,:,s), dq=dq, &
                                            block_ids=self%pml%blocks_z_p, face=PML_FACE_Z_P, inv_eps_scale=inv_eps_scale, &
                                            inv_mu_scale=inv_mu_scale)
         call compute_z_face_pml_rhs(self=self, q=q, q_face=self%rk_pml%q_pml_z_p_rk(:,:,:,:,:,s), &
                                     dq_face=self%rk_pml%dq_pml_z_p, block_ids=self%pml%blocks_z_p, face=PML_FACE_Z_P, &
                                     inv_eps_scale=inv_eps_scale, inv_mu_scale=inv_mu_scale)
      endif
   else
      if (allocated(self%pml%q_pml_x_m)) then
         call apply_x_face_field_correction(self=self, q=q, q_face=self%pml%q_pml_x_m, dq=dq, block_ids=self%pml%blocks_x_m, &
                                            face=PML_FACE_X_M, inv_eps_scale=inv_eps_scale, inv_mu_scale=inv_mu_scale)
         call compute_x_face_pml_rhs(self=self, q=q, q_face=self%pml%q_pml_x_m, dq_face=self%rk_pml%dq_pml_x_m, &
                                     block_ids=self%pml%blocks_x_m, face=PML_FACE_X_M, inv_eps_scale=inv_eps_scale, &
                                     inv_mu_scale=inv_mu_scale)
      endif
      if (allocated(self%pml%q_pml_x_p)) then
         call apply_x_face_field_correction(self=self, q=q, q_face=self%pml%q_pml_x_p, dq=dq, block_ids=self%pml%blocks_x_p, &
                                            face=PML_FACE_X_P, inv_eps_scale=inv_eps_scale, inv_mu_scale=inv_mu_scale)
         call compute_x_face_pml_rhs(self=self, q=q, q_face=self%pml%q_pml_x_p, dq_face=self%rk_pml%dq_pml_x_p, &
                                     block_ids=self%pml%blocks_x_p, face=PML_FACE_X_P, inv_eps_scale=inv_eps_scale, &
                                     inv_mu_scale=inv_mu_scale)
      endif
      if (allocated(self%pml%q_pml_y_m)) then
         call apply_y_face_field_correction(self=self, q=q, q_face=self%pml%q_pml_y_m, dq=dq, block_ids=self%pml%blocks_y_m, &
                                            face=PML_FACE_Y_M, inv_eps_scale=inv_eps_scale, inv_mu_scale=inv_mu_scale)
         call compute_y_face_pml_rhs(self=self, q=q, q_face=self%pml%q_pml_y_m, dq_face=self%rk_pml%dq_pml_y_m, &
                                     block_ids=self%pml%blocks_y_m, face=PML_FACE_Y_M, inv_eps_scale=inv_eps_scale, &
                                     inv_mu_scale=inv_mu_scale)
      endif
      if (allocated(self%pml%q_pml_y_p)) then
         call apply_y_face_field_correction(self=self, q=q, q_face=self%pml%q_pml_y_p, dq=dq, block_ids=self%pml%blocks_y_p, &
                                            face=PML_FACE_Y_P, inv_eps_scale=inv_eps_scale, inv_mu_scale=inv_mu_scale)
         call compute_y_face_pml_rhs(self=self, q=q, q_face=self%pml%q_pml_y_p, dq_face=self%rk_pml%dq_pml_y_p, &
                                     block_ids=self%pml%blocks_y_p, face=PML_FACE_Y_P, inv_eps_scale=inv_eps_scale, &
                                     inv_mu_scale=inv_mu_scale)
      endif
      if (allocated(self%pml%q_pml_z_m)) then
         call apply_z_face_field_correction(self=self, q=q, q_face=self%pml%q_pml_z_m, dq=dq, block_ids=self%pml%blocks_z_m, &
                                            face=PML_FACE_Z_M, inv_eps_scale=inv_eps_scale, inv_mu_scale=inv_mu_scale)
         call compute_z_face_pml_rhs(self=self, q=q, q_face=self%pml%q_pml_z_m, dq_face=self%rk_pml%dq_pml_z_m, &
                                     block_ids=self%pml%blocks_z_m, face=PML_FACE_Z_M, inv_eps_scale=inv_eps_scale, &
                                     inv_mu_scale=inv_mu_scale)
      endif
      if (allocated(self%pml%q_pml_z_p)) then
         call apply_z_face_field_correction(self=self, q=q, q_face=self%pml%q_pml_z_p, dq=dq, block_ids=self%pml%blocks_z_p, &
                                            face=PML_FACE_Z_P, inv_eps_scale=inv_eps_scale, inv_mu_scale=inv_mu_scale)
         call compute_z_face_pml_rhs(self=self, q=q, q_face=self%pml%q_pml_z_p, dq_face=self%rk_pml%dq_pml_z_p, &
                                     block_ids=self%pml%blocks_z_p, face=PML_FACE_Z_P, inv_eps_scale=inv_eps_scale, &
                                     inv_mu_scale=inv_mu_scale)
      endif
   endif
   endsubroutine apply_pml_fd_centered_ade

   subroutine select_pml_field_scales(self, inv_eps_scale, inv_mu_scale)
   !< Convert D/B derivatives into E/H derivatives for the current physical scaling.
   class(prism_cpu_object), intent(in)  :: self
   real(R8P),               intent(out) :: inv_eps_scale
   real(R8P),               intent(out) :: inv_mu_scale

   select case (self%physics%physical_model)
   case (ADIM_EM_PHYSICAL_MODEL)
      inv_eps_scale = 1._R8P
      inv_mu_scale  = 1._R8P
   case default
      inv_eps_scale = 1._R8P / EPS0
      inv_mu_scale  = 1._R8P / MU0
   endselect
   endsubroutine select_pml_field_scales

   subroutine apply_x_face_field_correction(self, q, q_face, dq, block_ids, face, inv_eps_scale, inv_mu_scale)
   class(prism_cpu_object), intent(in)    :: self
   real(R8P),               intent(in)    :: q(1:,          &
                                               1-self%ngc:, &
                                               1-self%ngc:, &
                                               1-self%ngc:, &
                                               1:)
   real(R8P),               intent(in)    :: q_face(1:,1:,1:,1:,1:)
   real(R8P),               intent(inout) :: dq(1:,          &
                                                 1-self%ngc:, &
                                                 1-self%ngc:, &
                                                 1-self%ngc:, &
                                                 1:)
   integer(I4P),            intent(in)    :: block_ids(1:)
   integer(I4P),            intent(in)    :: face
   real(R8P),               intent(in)    :: inv_eps_scale
   real(R8P),               intent(in)    :: inv_mu_scale
   integer(I4P)                           :: b, i, j, k, lid, li
   integer(I4P)                           :: cells
   integer(I4P)                           :: i0
   integer(I4P)                           :: s1
   real(R8P)                              :: center_distance
   real(R8P)                              :: d_field
   real(R8P)                              :: kappa
   real(R8P)                              :: span

   s1 = self%fdv_half_stencils(1)
   span = self%pml%profile_span(face)
   do lid=1, size(block_ids) !itero sul numero di blocchi che hanno PML sul lato x
      b     = block_ids(lid) !ottengo l'id globale del blocco lid-esimo di pml
      i0    = self%pml%ni_pml(1,b,face) !ottengo l'indice di cella in cui il pml inizia nel blocco lid-esimo
      cells = self%pml%ni_pml(2,b,face) - i0 + 1_I4P !calcolo il numero di celle del pml nel blocco lid-esimo
      do k=1, self%nk
         do j=1, self%nj
            do li=1, cells
               i = i0 + li - 1_I4P !ottengo l'indice globale della cella li-esima del pml nel blocco lid-esimo.
                                   ! Ho così ottenuto l'indice globale della cella in cui applicare la correzione del pml
               center_distance = merge(self%adam%field%emin(1,b) - self%adam%grid%domain_emin(1) + real(i - 1_I4P, R8P) * &
                                       self%adam%field%dxyz(1,b),                                                          &
                                       self%adam%grid%domain_emax(1) - self%adam%field%emax(1,b) + real(self%ni - i, R8P) * &
                                       self%adam%field%dxyz(1,b), face == PML_FACE_X_M)
               kappa = compute_pml_kappa(self=self, center_distance=center_distance, span=span)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(1,b), q=q(VAR_BZ,i-s1:i+s1,j,k,b), dq_ds=d_field)
               dq(VAR_DY,i,j,k,b) = dq(VAR_DY,i,j,k,b) + (1._R8P / kappa - 1._R8P) * (-d_field * inv_mu_scale)
               dq(VAR_DY,i,j,k,b) = dq(VAR_DY,i,j,k,b) + q_face(4,li,j,k,lid)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(1,b), q=q(VAR_BY,i-s1:i+s1,j,k,b), dq_ds=d_field)
               dq(VAR_DZ,i,j,k,b) = dq(VAR_DZ,i,j,k,b) + (1._R8P / kappa - 1._R8P) * (d_field * inv_mu_scale)
               dq(VAR_DZ,i,j,k,b) = dq(VAR_DZ,i,j,k,b) - q_face(3,li,j,k,lid)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(1,b), q=q(VAR_DZ,i-s1:i+s1,j,k,b), dq_ds=d_field)
               dq(VAR_BY,i,j,k,b) = dq(VAR_BY,i,j,k,b) + (1._R8P / kappa - 1._R8P) * (d_field * inv_eps_scale)
               dq(VAR_BY,i,j,k,b) = dq(VAR_BY,i,j,k,b) - q_face(2,li,j,k,lid)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(1,b), q=q(VAR_DY,i-s1:i+s1,j,k,b), dq_ds=d_field)
               dq(VAR_BZ,i,j,k,b) = dq(VAR_BZ,i,j,k,b) + (1._R8P / kappa - 1._R8P) * (-d_field * inv_eps_scale)
               dq(VAR_BZ,i,j,k,b) = dq(VAR_BZ,i,j,k,b) + q_face(1,li,j,k,lid)
            enddo
         enddo
      enddo
   enddo
   endsubroutine apply_x_face_field_correction

   subroutine apply_y_face_field_correction(self, q, q_face, dq, block_ids, face, inv_eps_scale, inv_mu_scale)
   class(prism_cpu_object), intent(in)    :: self
   real(R8P),               intent(in)    :: q(1:,          &
                                               1-self%ngc:, &
                                               1-self%ngc:, &
                                               1-self%ngc:, &
                                               1:)
   real(R8P),               intent(in)    :: q_face(1:,1:,1:,1:,1:)
   real(R8P),               intent(inout) :: dq(1:,          &
                                                 1-self%ngc:, &
                                                 1-self%ngc:, &
                                                 1-self%ngc:, &
                                                 1:)
   integer(I4P),            intent(in)    :: block_ids(1:)
   integer(I4P),            intent(in)    :: face
   real(R8P),               intent(in)    :: inv_eps_scale
   real(R8P),               intent(in)    :: inv_mu_scale
   integer(I4P)                           :: b, i, j, k, lid, lj
   integer(I4P)                           :: cells
   integer(I4P)                           :: j0
   integer(I4P)                           :: s1
   real(R8P)                              :: center_distance
   real(R8P)                              :: d_field
   real(R8P)                              :: kappa
   real(R8P)                              :: span

   s1 = self%fdv_half_stencils(1)
   span = self%pml%profile_span(face)
   do lid=1, size(block_ids)
      b     = block_ids(lid)
      j0    = self%pml%nj_pml(1,b,face)
      cells = self%pml%nj_pml(2,b,face) - j0 + 1_I4P
      do k=1, self%nk
         do lj=1, cells
            j = j0 + lj - 1_I4P
            center_distance = merge(self%adam%field%emin(2,b) - self%adam%grid%domain_emin(2) + real(j - 1_I4P, R8P) * &
                                    self%adam%field%dxyz(2,b),                                                          &
                                    self%adam%grid%domain_emax(2) - self%adam%field%emax(2,b) + real(self%nj - j, R8P) * &
                                    self%adam%field%dxyz(2,b), face == PML_FACE_Y_M)
            kappa = compute_pml_kappa(self=self, center_distance=center_distance, span=span)
            do i=1, self%ni
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(2,b), q=q(VAR_BZ,i,j-s1:j+s1,k,b), dq_ds=d_field)
               dq(VAR_DX,i,j,k,b) = dq(VAR_DX,i,j,k,b) + (1._R8P / kappa - 1._R8P) * (d_field * inv_mu_scale)
               dq(VAR_DX,i,j,k,b) = dq(VAR_DX,i,j,k,b) - q_face(4,i,lj,k,lid)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(2,b), q=q(VAR_BX,i,j-s1:j+s1,k,b), dq_ds=d_field)
               dq(VAR_DZ,i,j,k,b) = dq(VAR_DZ,i,j,k,b) + (1._R8P / kappa - 1._R8P) * (-d_field * inv_mu_scale)
               dq(VAR_DZ,i,j,k,b) = dq(VAR_DZ,i,j,k,b) + q_face(3,i,lj,k,lid)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(2,b), q=q(VAR_DZ,i,j-s1:j+s1,k,b), dq_ds=d_field)
               dq(VAR_BX,i,j,k,b) = dq(VAR_BX,i,j,k,b) + (1._R8P / kappa - 1._R8P) * (-d_field * inv_eps_scale)
               dq(VAR_BX,i,j,k,b) = dq(VAR_BX,i,j,k,b) + q_face(2,i,lj,k,lid)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(2,b), q=q(VAR_DX,i,j-s1:j+s1,k,b), dq_ds=d_field)
               dq(VAR_BZ,i,j,k,b) = dq(VAR_BZ,i,j,k,b) + (1._R8P / kappa - 1._R8P) * (d_field * inv_eps_scale)
               dq(VAR_BZ,i,j,k,b) = dq(VAR_BZ,i,j,k,b) - q_face(1,i,lj,k,lid)
            enddo
         enddo
      enddo
   enddo
   endsubroutine apply_y_face_field_correction

   subroutine apply_z_face_field_correction(self, q, q_face, dq, block_ids, face, inv_eps_scale, inv_mu_scale)
   class(prism_cpu_object), intent(in)    :: self
   real(R8P),               intent(in)    :: q(1:,          &
                                               1-self%ngc:, &
                                               1-self%ngc:, &
                                               1-self%ngc:, &
                                               1:)
   real(R8P),               intent(in)    :: q_face(1:,1:,1:,1:,1:)
   real(R8P),               intent(inout) :: dq(1:,          &
                                                 1-self%ngc:, &
                                                 1-self%ngc:, &
                                                 1-self%ngc:, &
                                                 1:)
   integer(I4P),            intent(in)    :: block_ids(1:)
   integer(I4P),            intent(in)    :: face
   real(R8P),               intent(in)    :: inv_eps_scale
   real(R8P),               intent(in)    :: inv_mu_scale
   integer(I4P)                           :: b, i, j, k, lid, lk
   integer(I4P)                           :: cells
   integer(I4P)                           :: k0
   integer(I4P)                           :: s1
   real(R8P)                              :: center_distance
   real(R8P)                              :: d_field
   real(R8P)                              :: kappa
   real(R8P)                              :: span

   s1 = self%fdv_half_stencils(1)
   span = self%pml%profile_span(face)
   do lid=1, size(block_ids)
      b     = block_ids(lid)
      k0    = self%pml%nk_pml(1,b,face)
      cells = self%pml%nk_pml(2,b,face) - k0 + 1_I4P
      do lk=1, cells
         k = k0 + lk - 1_I4P
         center_distance = merge(self%adam%field%emin(3,b) - self%adam%grid%domain_emin(3) + real(k - 1_I4P, R8P) * &
                                 self%adam%field%dxyz(3,b),                                                          &
                                 self%adam%grid%domain_emax(3) - self%adam%field%emax(3,b) + real(self%nk - k, R8P) * &
                                 self%adam%field%dxyz(3,b), face == PML_FACE_Z_M)
         kappa = compute_pml_kappa(self=self, center_distance=center_distance, span=span)
         do j=1, self%nj
            do i=1, self%ni
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(3,b), q=q(VAR_BY,i,j,k-s1:k+s1,b), dq_ds=d_field)
               dq(VAR_DX,i,j,k,b) = dq(VAR_DX,i,j,k,b) + (1._R8P / kappa - 1._R8P) * (-d_field * inv_mu_scale)
               dq(VAR_DX,i,j,k,b) = dq(VAR_DX,i,j,k,b) + q_face(4,i,j,lk,lid)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(3,b), q=q(VAR_BX,i,j,k-s1:k+s1,b), dq_ds=d_field)
               dq(VAR_DY,i,j,k,b) = dq(VAR_DY,i,j,k,b) + (1._R8P / kappa - 1._R8P) * (d_field * inv_mu_scale)
               dq(VAR_DY,i,j,k,b) = dq(VAR_DY,i,j,k,b) - q_face(3,i,j,lk,lid)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(3,b), q=q(VAR_DY,i,j,k-s1:k+s1,b), dq_ds=d_field)
               dq(VAR_BX,i,j,k,b) = dq(VAR_BX,i,j,k,b) + (1._R8P / kappa - 1._R8P) * (d_field * inv_eps_scale)
               dq(VAR_BX,i,j,k,b) = dq(VAR_BX,i,j,k,b) - q_face(2,i,j,lk,lid)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(3,b), q=q(VAR_DX,i,j,k-s1:k+s1,b), dq_ds=d_field)
               dq(VAR_BY,i,j,k,b) = dq(VAR_BY,i,j,k,b) + (1._R8P / kappa - 1._R8P) * (-d_field * inv_eps_scale)
               dq(VAR_BY,i,j,k,b) = dq(VAR_BY,i,j,k,b) + q_face(1,i,j,lk,lid)
            enddo
         enddo
      enddo
   enddo
   endsubroutine apply_z_face_field_correction

   subroutine compute_x_face_pml_rhs(self, q, q_face, dq_face, block_ids, face, inv_eps_scale, inv_mu_scale)
   class(prism_cpu_object), intent(in)    :: self
   real(R8P),               intent(in)    :: q(1:,          &
                                               1-self%ngc:, &
                                               1-self%ngc:, &
                                               1-self%ngc:, &
                                               1:)
   real(R8P),               intent(in)    :: q_face(1:,1:,1:,1:,1:)
   real(R8P),               intent(inout) :: dq_face(1:,1:,1:,1:,1:)
   integer(I4P),            intent(in)    :: block_ids(1:)
   integer(I4P),            intent(in)    :: face
   real(R8P),               intent(in)    :: inv_eps_scale
   real(R8P),               intent(in)    :: inv_mu_scale
   integer(I4P)                           :: b, i, j, k, lid, li
   integer(I4P)                           :: cells
   integer(I4P)                           :: i0
   integer(I4P)                           :: s1
   real(R8P)                              :: center_distance
   real(R8P)                              :: gamma
   real(R8P)                              :: alpha
   real(R8P)                              :: kappa
   real(R8P)                              :: d_field
   real(R8P)                              :: span

   s1 = self%fdv_half_stencils(1)
   span = self%pml%profile_span(face)
   do lid=1, size(block_ids)
      b     = block_ids(lid)
      i0    = self%pml%ni_pml(1,b,face)
      cells = self%pml%ni_pml(2,b,face) - i0 + 1_I4P
      do k=1, self%nk
         do j=1, self%nj
            do li=1, cells
               i     = i0 + li - 1_I4P
               center_distance = merge(self%adam%field%emin(1,b) - self%adam%grid%domain_emin(1) + real(i - 1_I4P, R8P) * &
                                       self%adam%field%dxyz(1,b),                                                          &
                                       self%adam%grid%domain_emax(1) - self%adam%field%emax(1,b) + real(self%ni - i, R8P) * &
                                       self%adam%field%dxyz(1,b), face == PML_FACE_X_M)
               call compute_pml_coefficients(self=self, center_distance=center_distance, span=span, gamma=gamma, alpha=alpha, &
                                             kappa=kappa)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(1,b), q=q(VAR_DY,i-s1:i+s1,j,k,b), dq_ds=d_field)
               dq_face(1,li,j,k,lid) = gamma / kappa**2 * d_field * inv_eps_scale - (alpha + gamma / kappa) * q_face(1,li,j,k,lid)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(1,b), q=q(VAR_DZ,i-s1:i+s1,j,k,b), dq_ds=d_field)
               dq_face(2,li,j,k,lid) = gamma / kappa**2 * d_field * inv_eps_scale - (alpha + gamma / kappa) * q_face(2,li,j,k,lid)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(1,b), q=q(VAR_BY,i-s1:i+s1,j,k,b), dq_ds=d_field)
               dq_face(3,li,j,k,lid) = gamma / kappa**2 * d_field * inv_mu_scale - (alpha + gamma / kappa) * q_face(3,li,j,k,lid)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(1,b), q=q(VAR_BZ,i-s1:i+s1,j,k,b), dq_ds=d_field)
               dq_face(4,li,j,k,lid) = gamma / kappa**2 * d_field * inv_mu_scale - (alpha + gamma / kappa) * q_face(4,li,j,k,lid)
            enddo
         enddo
      enddo
   enddo
   endsubroutine compute_x_face_pml_rhs

   subroutine compute_y_face_pml_rhs(self, q, q_face, dq_face, block_ids, face, inv_eps_scale, inv_mu_scale)
   class(prism_cpu_object), intent(in)    :: self
   real(R8P),               intent(in)    :: q(1:,          &
                                               1-self%ngc:, &
                                               1-self%ngc:, &
                                               1-self%ngc:, &
                                               1:)
   real(R8P),               intent(in)    :: q_face(1:,1:,1:,1:,1:)
   real(R8P),               intent(inout) :: dq_face(1:,1:,1:,1:,1:)
   integer(I4P),            intent(in)    :: block_ids(1:)
   integer(I4P),            intent(in)    :: face
   real(R8P),               intent(in)    :: inv_eps_scale
   real(R8P),               intent(in)    :: inv_mu_scale
   integer(I4P)                           :: b, i, j, k, lid, lj
   integer(I4P)                           :: cells
   integer(I4P)                           :: j0
   integer(I4P)                           :: s1
   real(R8P)                              :: center_distance
   real(R8P)                              :: gamma
   real(R8P)                              :: alpha
   real(R8P)                              :: kappa
   real(R8P)                              :: d_field
   real(R8P)                              :: span

   s1 = self%fdv_half_stencils(1)
   span = self%pml%profile_span(face)
   do lid=1, size(block_ids)
      b     = block_ids(lid)
      j0    = self%pml%nj_pml(1,b,face)
      cells = self%pml%nj_pml(2,b,face) - j0 + 1_I4P
      do k=1, self%nk
         do lj=1, cells
            j     = j0 + lj - 1_I4P
            center_distance = merge(self%adam%field%emin(2,b) - self%adam%grid%domain_emin(2) + real(j - 1_I4P, R8P) * &
                                    self%adam%field%dxyz(2,b),                                                          &
                                    self%adam%grid%domain_emax(2) - self%adam%field%emax(2,b) + real(self%nj - j, R8P) * &
                                    self%adam%field%dxyz(2,b), face == PML_FACE_Y_M)
            call compute_pml_coefficients(self=self, center_distance=center_distance, span=span, gamma=gamma, alpha=alpha, &
                                          kappa=kappa)
            do i=1, self%ni
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(2,b), q=q(VAR_DX,i,j-s1:j+s1,k,b), dq_ds=d_field)
               dq_face(1,i,lj,k,lid) = gamma / kappa**2 * d_field * inv_eps_scale - (alpha + gamma / kappa) * q_face(1,i,lj,k,lid)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(2,b), q=q(VAR_DZ,i,j-s1:j+s1,k,b), dq_ds=d_field)
               dq_face(2,i,lj,k,lid) = gamma / kappa**2 * d_field * inv_eps_scale - (alpha + gamma / kappa) * q_face(2,i,lj,k,lid)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(2,b), q=q(VAR_BX,i,j-s1:j+s1,k,b), dq_ds=d_field)
               dq_face(3,i,lj,k,lid) = gamma / kappa**2 * d_field * inv_mu_scale - (alpha + gamma / kappa) * q_face(3,i,lj,k,lid)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(2,b), q=q(VAR_BZ,i,j-s1:j+s1,k,b), dq_ds=d_field)
               dq_face(4,i,lj,k,lid) = gamma / kappa**2 * d_field * inv_mu_scale - (alpha + gamma / kappa) * q_face(4,i,lj,k,lid)
            enddo
         enddo
      enddo
   enddo
   endsubroutine compute_y_face_pml_rhs

   subroutine compute_z_face_pml_rhs(self, q, q_face, dq_face, block_ids, face, inv_eps_scale, inv_mu_scale)
   class(prism_cpu_object), intent(in)    :: self
   real(R8P),               intent(in)    :: q(1:,          &
                                               1-self%ngc:, &
                                               1-self%ngc:, &
                                               1-self%ngc:, &
                                               1:)
   real(R8P),               intent(in)    :: q_face(1:,1:,1:,1:,1:)
   real(R8P),               intent(inout) :: dq_face(1:,1:,1:,1:,1:)
   integer(I4P),            intent(in)    :: block_ids(1:)
   integer(I4P),            intent(in)    :: face
   real(R8P),               intent(in)    :: inv_eps_scale
   real(R8P),               intent(in)    :: inv_mu_scale
   integer(I4P)                           :: b, i, j, k, lid, lk
   integer(I4P)                           :: cells
   integer(I4P)                           :: k0
   integer(I4P)                           :: s1
   real(R8P)                              :: center_distance
   real(R8P)                              :: gamma
   real(R8P)                              :: alpha
   real(R8P)                              :: kappa
   real(R8P)                              :: d_field
   real(R8P)                              :: span

   s1 = self%fdv_half_stencils(1)
   span = self%pml%profile_span(face)
   do lid=1, size(block_ids)
      b     = block_ids(lid)
      k0    = self%pml%nk_pml(1,b,face)
      cells = self%pml%nk_pml(2,b,face) - k0 + 1_I4P
      do lk=1, cells
         k     = k0 + lk - 1_I4P
         center_distance = merge(self%adam%field%emin(3,b) - self%adam%grid%domain_emin(3) + real(k - 1_I4P, R8P) * &
                                 self%adam%field%dxyz(3,b),                                                          &
                                 self%adam%grid%domain_emax(3) - self%adam%field%emax(3,b) + real(self%nk - k, R8P) * &
                                 self%adam%field%dxyz(3,b), face == PML_FACE_Z_M)
         call compute_pml_coefficients(self=self, center_distance=center_distance, span=span, gamma=gamma, alpha=alpha, &
                                       kappa=kappa)
         do j=1, self%nj
            do i=1, self%ni
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(3,b), q=q(VAR_DX,i,j,k-s1:k+s1,b), dq_ds=d_field)
               dq_face(1,i,j,lk,lid) = gamma / kappa**2 * d_field * inv_eps_scale - (alpha + gamma / kappa) * q_face(1,i,j,lk,lid)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(3,b), q=q(VAR_DY,i,j,k-s1:k+s1,b), dq_ds=d_field)
               dq_face(2,i,j,lk,lid) = gamma / kappa**2 * d_field * inv_eps_scale - (alpha + gamma / kappa) * q_face(2,i,j,lk,lid)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(3,b), q=q(VAR_BX,i,j,k-s1:k+s1,b), dq_ds=d_field)
               dq_face(3,i,j,lk,lid) = gamma / kappa**2 * d_field * inv_mu_scale - (alpha + gamma / kappa) * q_face(3,i,j,lk,lid)
               call compute_derivative1_fd_centered(s=s1, ds=self%adam%field%dxyz(3,b), q=q(VAR_BY,i,j,k-s1:k+s1,b), dq_ds=d_field)
               dq_face(4,i,j,lk,lid) = gamma / kappa**2 * d_field * inv_mu_scale - (alpha + gamma / kappa) * q_face(4,i,j,lk,lid)
            enddo
         enddo
      enddo
   enddo
   endsubroutine compute_z_face_pml_rhs

   pure subroutine compute_pml_coefficients(self, center_distance, span, gamma, alpha, kappa)
   !< Return the face-local ADE-PML coefficients for the supported variants.
   class(prism_cpu_object), intent(in) :: self
   real(R8P),               intent(in) :: center_distance
   real(R8P),               intent(in) :: span
   real(R8P),               intent(out):: gamma
   real(R8P),               intent(out):: alpha
   real(R8P),               intent(out):: kappa
   real(R8P)                           :: depth
   real(R8P)                           :: distance_to_outer
   real(R8P), parameter                :: BERMUDEZ_EPS = 1.e-6_R8P
   integer(I4P), parameter             :: CFS_PROFILE_EXPONENT = 2_I4P

   if (span < 0._R8P) then
      gamma = 0._R8P
      alpha = 0._R8P
      kappa = 1._R8P
      return
   endif

   gamma = 0._R8P
   alpha = 0._R8P
   kappa = 1._R8P

   select case (trim(self%pml%pml_type))
   case ('CLASSIC')
      if (span > 0._R8P) then
         depth = 1._R8P - center_distance / span
      else
         depth = 1._R8P
      endif
      depth = max(0._R8P, min(1._R8P, depth))
      gamma = self%pml%gamma_max * depth**self%pml%gamma_exponent
   case ('BERMUDEZ')
      if (span > 0._R8P) then
         depth = 1._R8P - center_distance / span
         distance_to_outer = self%pml%width * center_distance / span + BERMUDEZ_EPS
      else
         depth = 1._R8P
         distance_to_outer = BERMUDEZ_EPS
      endif
      depth = max(0._R8P, min(1._R8P, depth))
      gamma = self%pml%beta / distance_to_outer**self%pml%gamma_exponent
   case ('CFS')
      if (span > 0._R8P) then
         depth = 1._R8P - center_distance / span
      else
         depth = 1._R8P
      endif
      depth = max(0._R8P, min(1._R8P, depth))
      gamma = self%pml%gamma_max * depth**CFS_PROFILE_EXPONENT
      alpha = self%pml%alpha_max * (1._R8P - depth)**CFS_PROFILE_EXPONENT
      kappa = 1._R8P + (self%pml%k_max - 1._R8P) * depth**CFS_PROFILE_EXPONENT
   case default
      gamma = 0._R8P
   endselect
   endsubroutine compute_pml_coefficients

   pure real(R8P) function compute_pml_kappa(self, center_distance, span) result(kappa)
   class(prism_cpu_object), intent(in) :: self
   real(R8P),               intent(in) :: center_distance
   real(R8P),               intent(in) :: span
   real(R8P)                           :: gamma
   real(R8P)                           :: alpha

   call compute_pml_coefficients(self=self, center_distance=center_distance, span=span, gamma=gamma, alpha=alpha, kappa=kappa)
   endfunction compute_pml_kappa

   subroutine compute_residuals_fv_centered(self, q, dq, s, flux_register)
   !< Compute residuals of equation, space operator, centered finite volume schemes.
   class(prism_cpu_object),     intent(inout)                    :: self                      !< The equation.
   real(R8P),                   intent(inout)                    :: q(1:,       &
                                                                    1-self%ngc:,&
                                                                    1-self%ngc:,&
                                                                    1-self%ngc:,&
                                                                    1:)                       !< Conservative variables.
   real(R8P),                   intent(inout)                    :: dq(1:,       &
                                                                     1-self%ngc:,&
                                                                     1-self%ngc:,&
                                                                     1-self%ngc:,&
                                                                     1:)                      !< Residuals.
   integer(I4P),                intent(in),    optional         :: s                          !< Stage counter.
   class(flux_register_object), intent(inout), optional         :: flux_register              !< Forest's flux register; FV reflux.
   integer(I4P)                                                 :: i,j,k,b,d,v                !< Counter
   integer(I4P)                                                 :: stage_idx                  !< Integrator stage index.
   real(R8P), parameter                                         :: sir(3,3) = reshape([1._R8P,0._R8P,0._R8P,  &
                                                                                       0._R8P,1._R8P,0._R8P,  &
                                                                                       0._R8P,0._R8P,1._R8P], &
                                                                                       [3,3]) !< Direction versor, real.
   real(R8P)                                                    :: damping_coeff

   ! Capture stage BEFORE the associate block rebinds `s` to the
   ! reconstruction stencil half-width. The inter-realm FV reflux hook
   ! below needs the stage index to address `forest_flux_register`'s
   ! per-stage accumulators.
   if (present(s)) then
      stage_idx = s
   else
      stage_idx = 0_I4P
   endif

   if (present(s)) then
      call self%update_ghost(q=q, s=s)
   else
      call self%update_ghost(q=q)
   endif
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv_c=>self%nv_c,blocks_number=>self%blocks_number, &
             dxyz=>self%adam%field%dxyz, flxyz_c=>self%flxyz_c, flx_f=>self%flx_f, fly_f=>self%fly_f, flz_f=>self%flz_f,        &
             s=>self%fdv_half_stencils(1),                                                                            &
             var_Jx=>self%physics%var_Jx, var_Jy=>self%physics%var_Jy, var_Jz=>self%physics%var_Jz, chi=>self%physics%chi,      &
             c_r=>self%physics%c_r)
   if (blocks_number > 0) then
      ! compute fluxes at cell centers
      do b=1, blocks_number
      do k=1-ngc, nk+ngc
      do j=1-ngc, nj+ngc
      do i=1-ngc, ni+ngc
         do d=1,3 ! x, y, z
            call compute_fluxes_maxwell(sir=sir(:,d),q=q(:,i,j,k,b),f=flxyz_c(:,1,d,i,j,k,b),chi=chi)
         enddo
      enddo
      enddo
      enddo
      enddo
      ! reconstruct fluxes at cell faces
      d = 1
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=0, ni
         do v=1, nv_c
            call compute_reconstruction_r_fv_centered(s=s,q=flxyz_c(v,1,d,1+i-s:i+s,j,k,b),qr=flx_f(v,i,j,k,b))
         enddo
      enddo
      enddo
      enddo
      enddo
      d = 2
      do b=1, blocks_number
      do k=1, nk
      do j=0, nj
      do i=1, ni
         do v=1, nv_c
            call compute_reconstruction_r_fv_centered(s=s,q=flxyz_c(v,1,d,i,1+j-s:j+s,k,b),qr=fly_f(v,i,j,k,b))
         enddo
      enddo
      enddo
      enddo
      enddo
      d = 3
      do b=1, blocks_number
      do k=0, nk
      do j=1, nj
      do i=1, ni
         do v=1, nv_c
            call compute_reconstruction_r_fv_centered(s=s,q=flxyz_c(v,1,d,i,j,1+k-s:k+s,b),qr=flz_f(v,i,j,k,b))
         enddo
      enddo
      enddo
      enddo
      enddo
      ! Inter-realm seam flux accumulation.
      !
      ! After face fluxes are reconstructed, before the conservative
      ! update consumes them, walk every (block, fec) pair and check the
      ! per-realm signed register lookup. Non-zero entries name a seam
      ! face: pack the corresponding flux skin from `flx_f`/`fly_f`/
      ! `flz_f` into the register's `(nv, nface_cells)` accumulator
      ! shape, then dispatch by sign (+ = coarse side → F_coarse, − =
      ! fine side → F_fine_sum). Symmetric mirror seams produce equal
      ! coarse/fine fluxes → end-of-step reflux is round-off zero by
      ! expectation and the FD-centered case is unreachable here (the
      ! hook only fires for fv_centered).
      !
      ! α.r1 end-of-step gate: accumulate ONLY at the realm's final RK
      ! substage. The register holds a single end-of-step bucket and the
      ! reflux correction consumes it at the same final substage.
      ! Earlier substages skip the accumulation entirely — their face
      ! fluxes are intermediate-stage values not used by the
      ! Berger-Colella end-of-step correction.
      if (present(flux_register) .and. stage_idx == self%rk%nrk &
                                  .and. allocated(self%adam%maps%inter_realm_face_register_index)) then
         if (flux_register%nfaces > 0_I4P) then
            call accumulate_seam_fluxes_fv(self, ni, nj, nk, nv_c, blocks_number, flux_register)
         endif
      endif
      ! compute fluxes difference for RHS, dF/ds = (F(i+1/2)-F(i-1/2))/Ds
      do b=1,blocks_number
      do k=1,nk
      do j=1,nj
      do i=1,ni
         do v=1, nv_c
            dq(v,i,j,k,b) = - (flx_f(v,i,j,k,b) - flx_f(v,i-1,j,k,b)) / dxyz(1,b) &
                            - (fly_f(v,i,j,k,b) - fly_f(v,i,j-1,k,b)) / dxyz(2,b) &
                            - (flz_f(v,i,j,k,b) - flz_f(v,i,j,k-1,b)) / dxyz(3,b)
         enddo
         ! J sources
         dq(VAR_DX,i,j,k,b) = dq(VAR_DX,i,j,k,b) - q(var_Jx,i,j,k,b)
         dq(VAR_DY,i,j,k,b) = dq(VAR_DY,i,j,k,b) - q(var_Jy,i,j,k,b)
         dq(VAR_DZ,i,j,k,b) = dq(VAR_DZ,i,j,k,b) - q(var_Jz,i,j,k,b)
         if (self%fv_add_phi_damping) then
            if (c_r > 0._R8P) then
               if (self%physics%physical_model == ADIM_EM_PHYSICAL_MODEL) then
                  damping_coeff = chi / (c_r * minval(dxyz(1:3,b)))
               else
                  damping_coeff = chi * C0 / (c_r * minval(dxyz(1:3,b)))
               endif
               dq(self%fv_ivar_phi,i,j,k,b) = dq(self%fv_ivar_phi,i,j,k,b) - &
                                              damping_coeff * q(self%fv_ivar_phi,i,j,k,b)
            endif
         endif
         if (self%fv_add_psi_damping) then
            if (c_r > 0._R8P) then
               if (self%physics%physical_model == ADIM_EM_PHYSICAL_MODEL) then
                  damping_coeff = chi / (c_r * minval(dxyz(1:3,b)))
               else
                  damping_coeff = chi * C0 / (c_r * minval(dxyz(1:3,b)))
               endif
               dq(self%fv_ivar_psi,i,j,k,b) = dq(self%fv_ivar_psi,i,j,k,b) - &
                                              damping_coeff * q(self%fv_ivar_psi,i,j,k,b)
            endif
         endif
      enddo
      enddo
      enddo
      enddo
   endif
   endassociate
   endsubroutine compute_residuals_fv_centered

   subroutine accumulate_seam_fluxes_fv(self, ni, nj, nk, nv_c, blocks_number, flux_register)
   !< Accumulate end-of-step face fluxes on inter-realm seam faces of
   !< the FV-centered scheme. Gated to the realm's final RK substage
   !< under α.r1.
   !<
   !< Walks every (block, fec) pair and, for non-zero entries in
   !< `self%adam%maps%inter_realm_face_register_index(b, fec)`, packs
   !< the appropriate `flx_f`/`fly_f`/`flz_f` face skin into a
   !< `(nv, nface_cells)` slab and routes it to the right register
   !< accumulator: positive index → coarse side → `accumulate_coarse_flux`,
   !< negative index → fine side → `accumulate_fine_flux`. The call site
   !< above gates on `stage_idx == self%rk%nrk`, so this routine fires
   !< exactly once per realm per step at the realm's end-of-step.
   !<
   !< The register accumulators are sized `(nv, nface_cells, 1)` under
   !< α.r1 (M2 reshape) with the full state-vector width `nv` (set at
   !< `register_face` time from `realm%adam%field%nv`). The FV scheme
   !< only fills the conservative slice `1:nv_c` of the flux arrays;
   !< the remaining rows are left at zero in `flux_slab` so the
   !< register's `F_coarse - F_fine_sum` correction is identically zero
   !< outside the conservative variables. The third-axis index passed
   !< to `accumulate_*_flux` is hardcoded to `1` (collapsed register).
   !<
   !< For a same-resolution mirror seam (rmf-2realm), each coarse and fine side
   !< carries an identical end-of-step flux; `F_coarse - F_fine_sum` is round-off
   !< zero in expectation and the downstream q-correction is a no-op. For a true
   !< 2:1 intra-realm AMR jump (#13 §7.5 M3) the coarse and fine sides differ:
   !<
   !<   * Coarse side (`sgn_idx > 0`): packs its full `nface_cells`-cell face skin
   !<     directly — the linear index `c` runs (inner axis fastest, outer axis
   !<     slowest) over the coarse cells, exactly as before.
   !<   * Fine side (`sgn_idx < 0`): each of the `ratio/2` (4 in 3D) fine blocks
   !<     covers ONE quadrant of the coarse face. Its `nj*nk` fine-face cells are
   !<     2:1-restricted (2x2 conservative average) onto the `(nj/2)*(nk/2)` coarse
   !<     cells of that quadrant, written into the coarse-skin-shaped slab at the
   !<     quadrant offset (zero elsewhere), then accumulated. The quadrant offset
   !<     is read from `maps%amr_seam_quadrant`, precomputed at registration from
   !<     the two blocks' Morton codes (issue #28 D2 — the coarse partner may
   !<     live on another rank, so its geometry is not readable here). The fine
   !<     blocks' `accumulate_fine_flux` calls sum into disjoint quadrants, so
   !<     `F_fine_sum` ends up holding the area-averaged fine flux over the WHOLE
   !<     coarse face, matching `F_coarse`'s shape for the Berger-Colella delta.
   !<
   !< Conservative 2:1 averaging (arithmetic mean of the 2x2 fine sub-faces) is the
   !< correct restriction for a face FLUX (Berger-Colella 1989 §4; Olivares 2019
   !< Eq. 26-27): the coarse-face flux equals the average of the fine-face fluxes
   !< covering it, so the telescoping `∮F·n` across the seam cancels to round-off.
   !<
   !< The FD-centered case does not reach this routine (different residual target).
   class(prism_cpu_object),     intent(inout) :: self            !< The realm.
   integer(I4P),                intent(in)    :: ni, nj, nk      !< Interior cell counts.
   integer(I4P),                intent(in)    :: nv_c            !< Number of conservative variables.
   integer(I4P),                intent(in)    :: blocks_number   !< Number of local blocks.
   class(flux_register_object), intent(inout) :: flux_register   !< Forest's flux register.
   integer(I4P)                               :: sgn_idx         !< Signed register index for (b, fec).
   integer(I4P)                               :: face_idx        !< |sgn_idx| → register face.
   integer(I4P)                               :: nv_reg          !< Register state-vector width.
   integer(I4P)                               :: nface_cells     !< Coarse-face skin cell count.
   integer(I4P)                               :: fec, b          !< Face, block counters.
   integer(I4P)                               :: inner_n, outer_n   !< Coarse cell counts along tangential axes.
   integer(I4P)                               :: ioff, joff      !< Fine-block quadrant offset (inner,outer) ∈ {0,1}.
   real(R8P), allocatable                     :: flux_slab(:,:)  !< Coarse-skin-shaped contribution.

   do b=1, blocks_number
      do fec=1, 6
         sgn_idx = self%adam%maps%inter_realm_face_register_index(b, fec)
         if (sgn_idx == 0_I4P) cycle
         face_idx = abs(sgn_idx)
         if (face_idx > flux_register%nfaces) cycle
         if (.not. allocated(flux_register%face(face_idx)%F_coarse)) cycle
         nv_reg      = size(flux_register%face(face_idx)%F_coarse, dim=1)
         nface_cells = flux_register%face(face_idx)%nface_cells
         if (allocated(flux_slab)) deallocate(flux_slab)
         allocate(flux_slab(1:nv_reg, 1:nface_cells))
         flux_slab = 0._R8P
         ! Tangential axes per face (inner fastest in the linear index c):
         ! x-faces (1,2): inner=y(2), outer=z(3); y-faces (3,4): inner=x(1),
         ! outer=z(3); z-faces (5,6): inner=x(1), outer=y(2).
         select case (fec)
         case (1_I4P, 2_I4P) ; inner_n = nj ; outer_n = nk
         case (3_I4P, 4_I4P) ; inner_n = ni ; outer_n = nk
         case (5_I4P, 6_I4P) ; inner_n = ni ; outer_n = nj
         case default ; cycle
         end select

         if (sgn_idx > 0_I4P) then
            ! Coarse side: pack the full face skin directly.
            call pack_coarse_face(self, fec, ni, nj, nk, nv_c, b, inner_n, outer_n, flux_slab)
            call flux_register%accumulate_coarse_flux(face_index=face_idx, stage=1_I4P, flux_face=flux_slab)
         else
            ! Fine side: 2:1-restrict this fine block's face into its quadrant of
            ! the coarse skin. Quadrant offsets are PRECOMPUTED at registration
            ! from the two blocks' Morton codes (issue #28 D2): the register's
            ! `coarse_block` is an owner-rank-LOCAL index, so deriving the
            ! quadrant from that block's emin/emax here would read an unrelated
            ! local block's geometry whenever the coarse partner lives on
            ! another rank. The quadrant table is allocated only by the
            ! intra-realm AMR registration pass; an inter-realm mirror seam
            ! (table unallocated) has no 2:1 quadrant — offsets are zero.
            if (allocated(self%adam%maps%amr_seam_quadrant)) then
               ioff = self%adam%maps%amr_seam_quadrant(1, b, fec)
               joff = self%adam%maps%amr_seam_quadrant(2, b, fec)
            else
               ioff = 0_I4P ; joff = 0_I4P
            endif
            call restrict_fine_face(self, fec, ni, nj, nk, nv_c, b, inner_n, outer_n, ioff, joff, flux_slab)
            call flux_register%accumulate_fine_flux(face_index=face_idx, stage=1_I4P, flux_face=flux_slab)
         endif
      enddo
   enddo
   if (allocated(flux_slab)) deallocate(flux_slab)

   contains
      subroutine pack_coarse_face(self, fec, ni, nj, nk, nv_c, b, inner_n, outer_n, slab)
      !< Pack a coarse block's face skin into `slab(1:nv, 1:inner_n*outer_n)`,
      !< inner axis fastest. Conservative slice `1:nv_c` only; rest left zero.
      class(prism_cpu_object), intent(in)    :: self
      integer(I4P),            intent(in)    :: fec, ni, nj, nk, nv_c, b
      integer(I4P),            intent(in)    :: inner_n, outer_n
      real(R8P),               intent(inout) :: slab(:,:)
      integer(I4P)                           :: i, j, k, v, c

      c = 0_I4P
      select case (fec)
      case (1_I4P) ; do k=1,nk ; do j=1,nj ; c=c+1 ; do v=1,nv_c ; slab(v,c)=self%flx_f(v,0_I4P,j,k,b) ; enddo ; enddo ; enddo
      case (2_I4P) ; do k=1,nk ; do j=1,nj ; c=c+1 ; do v=1,nv_c ; slab(v,c)=self%flx_f(v,ni,j,k,b)   ; enddo ; enddo ; enddo
      case (3_I4P) ; do k=1,nk ; do i=1,ni ; c=c+1 ; do v=1,nv_c ; slab(v,c)=self%fly_f(v,i,0_I4P,k,b) ; enddo ; enddo ; enddo
      case (4_I4P) ; do k=1,nk ; do i=1,ni ; c=c+1 ; do v=1,nv_c ; slab(v,c)=self%fly_f(v,i,nj,k,b)   ; enddo ; enddo ; enddo
      case (5_I4P) ; do j=1,nj ; do i=1,ni ; c=c+1 ; do v=1,nv_c ; slab(v,c)=self%flz_f(v,i,j,0_I4P,b) ; enddo ; enddo ; enddo
      case (6_I4P) ; do j=1,nj ; do i=1,ni ; c=c+1 ; do v=1,nv_c ; slab(v,c)=self%flz_f(v,i,j,nk,b)   ; enddo ; enddo ; enddo
      end select
      if (c /= inner_n*outer_n) &
         call mpih%error_stop(msg='prism_cpu_object: FV coarse seam-flux pack count != nface_cells')
      endsubroutine pack_coarse_face

      subroutine restrict_fine_face(self, fec, ni, nj, nk, nv_c, b, inner_n, outer_n, ioff, joff, slab)
      !< Pack fine block `b`'s tangential face flux into a raw (nv_c, inner_n,
      !< outer_n) array, then delegate the 2:1 restriction to the pure module
      !< helper `restrict_fine_face_to_quadrant`. `slab` is the coarse-skin slab;
      !< only this block's (ioff,joff) quadrant is written.
      class(prism_cpu_object), intent(in)    :: self
      integer(I4P),            intent(in)    :: fec, ni, nj, nk, nv_c, b
      integer(I4P),            intent(in)    :: inner_n, outer_n, ioff, joff
      real(R8P),               intent(inout) :: slab(:,:)
      real(R8P)                              :: fine_face(1:nv_c, 1:inner_n, 1:outer_n)
      integer(I4P)                           :: fi, fo

      do fo = 1_I4P, outer_n
         do fi = 1_I4P, inner_n
            fine_face(1:nv_c, fi, fo) = fine_face_cell(self, fec, ni, nj, nk, nv_c, b, fi, fo)
         enddo
      enddo
      call restrict_fine_face_to_quadrant(fine_face=fine_face, inner_n=inner_n, outer_n=outer_n, &
                                          ioff=ioff, joff=joff, slab=slab)
      endsubroutine restrict_fine_face

      function fine_face_cell(self, fec, ni, nj, nk, nv_c, b, fi, fo) result(val)
      !< Conservative-slice face flux at fine-block cell (inner=fi, outer=fo) on
      !< face `fec` of block `b`. Maps (fi,fo) onto the face's (i,j,k) per fec.
      class(prism_cpu_object), intent(in) :: self
      integer(I4P),            intent(in) :: fec, ni, nj, nk, nv_c, b, fi, fo
      real(R8P)                           :: val(1:nv_c)

      select case (fec)
      case (1_I4P) ; val(1:nv_c) = self%flx_f(1:nv_c, 0_I4P, fi, fo, b)   ! -x: inner=j, outer=k
      case (2_I4P) ; val(1:nv_c) = self%flx_f(1:nv_c, ni,    fi, fo, b)   ! +x
      case (3_I4P) ; val(1:nv_c) = self%fly_f(1:nv_c, fi, 0_I4P, fo, b)   ! -y: inner=i, outer=k
      case (4_I4P) ; val(1:nv_c) = self%fly_f(1:nv_c, fi, nj,    fo, b)   ! +y
      case (5_I4P) ; val(1:nv_c) = self%flz_f(1:nv_c, fi, fo, 0_I4P, b)   ! -z: inner=i, outer=j
      case (6_I4P) ; val(1:nv_c) = self%flz_f(1:nv_c, fi, fo, nk,    b)   ! +z
      case default ; val(1:nv_c) = 0._R8P
      end select
      endfunction fine_face_cell
   endsubroutine accumulate_seam_fluxes_fv

   subroutine compute_residuals_weno(self, q, dq, s, flux_register)
   !< Compute residuals of equation, space operator, WENO schemes.
   class(prism_cpu_object),     intent(inout)                   :: self          !< The equation.
   real(R8P),                   intent(inout)                   :: q(1:,         &
                                                                     1-self%ngc:,&
                                                                     1-self%ngc:,&
                                                                     1-self%ngc:,&
                                                                     1:)         !< Conservative variables.
   real(R8P),                   intent(inout)                   :: dq(1:,         &
                                                                      1-self%ngc:,&
                                                                      1-self%ngc:,&
                                                                      1-self%ngc:,&
                                                                      1:)        !< Residuals.
   integer(I4P),                intent(in),    optional         :: s             !< Stage counter.
   class(flux_register_object), intent(inout), optional         :: flux_register !< Flux register.

   if (present(s)) then
      call self%update_ghost(q=q, s=s)
   else
      call self%update_ghost(q=q)
   endif
   !call self%integrate_eikonal_coils(q=q)
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv, nv_c=>self%nv_c,blocks_number=>self%blocks_number,&
             dx=>self%adam%field%dxyz(1,:), dy=>self%adam%field%dxyz(2,:), dz=>self%adam%field%dxyz(3,:),                         &
             flx=>self%flxyz_c(:,1,1,:,:,:,:), fly=>self%flxyz_c(:,1,2,:,:,:,:), flz=>self%flxyz_c(:,1,3,:,:,:,:),                &
             weno_s=>self%weno%S, weno_zeps=>self%weno%zeps,                                                                      &
             weno_a=>self%weno%a, weno_p=>self%weno%p, weno_d=>self%weno%d, weno_c=>self%weno%c,                                  &
             var_Jx=>self%physics%var_Jx, var_Jy=>self%physics%var_Jy, var_Jz=>self%physics%var_Jz, chi=>self%physics%chi,        &
             evmax=>self%physics%evmax, erw=>self%physics%erw, elw=>self%physics%elw)

   if (blocks_number > 0) then
      call compute_fluxes_convective_weno(dir=1,blocks_number=blocks_number,ni=ni,nj=nj,nk=nk,ngc=ngc,nv_c=nv_c,               &
                                     weno_s=weno_S,weno_a=weno_a,weno_p=weno_p,weno_d=weno_d,weno_c=weno_c,weno_zeps=weno_zeps,&
                                     evmax=evmax,erw=erw,elw=elw,chi=chi,                                                      &
                                     q=q,fluxes=flx)
      call compute_fluxes_convective_weno(dir=2,blocks_number=blocks_number,ni=ni,nj=nj,nk=nk,ngc=ngc,nv_c=nv_c,               &
                                     weno_s=weno_S,weno_a=weno_a,weno_p=weno_p,weno_d=weno_d,weno_c=weno_c,weno_zeps=weno_zeps,&
                                     evmax=evmax,erw=erw,elw=elw,chi=chi,                                                      &
                                     q=q,fluxes=fly)
      call compute_fluxes_convective_weno(dir=3,blocks_number=blocks_number,ni=ni,nj=nj,nk=nk,ngc=ngc,nv_c=nv_c,               &
                                     weno_s=weno_S,weno_a=weno_a,weno_p=weno_p,weno_d=weno_d,weno_c=weno_c,weno_zeps=weno_zeps,&
                                     evmax=evmax,erw=erw,elw=elw,chi=chi,                                                      &
                                     q=q,fluxes=flz)
      call compute_fluxes_difference(blocks_number=blocks_number, ni=ni, nj=nj, nk=nk, ngc=ngc, nv_c=nv_c, &
                                     var_Jx=var_Jx, var_Jy=var_Jy, var_Jz=var_Jz,                          &
                                     dx=dx, dy=dy, dz=dz, flx=flx, fly=fly, flz=flz, dq=dq, q=q)
   endif
   endassociate
   endsubroutine compute_residuals_weno

   subroutine integrate_blanesmoan(self)
   !< Integrate equation, time operator, Yoshida RK scheme.
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   associate(nc=>self%blanesmoan%nc,a=>self%blanesmoan%a,b=>self%blanesmoan%b)
   call self%compute_coils_current(q=self%q)
   do s=1, nc
      call self%compute_residuals(q=self%q, dq=self%dq)
      if (s==1) call self%save_residuals
      self%q(VAR_BX,:,:,:,:) = self%q(VAR_BX,:,:,:,:) + b(s) * self%time%dt * self%dq(VAR_BX,:,:,:,:)
      self%q(VAR_BY,:,:,:,:) = self%q(VAR_BY,:,:,:,:) + b(s) * self%time%dt * self%dq(VAR_BY,:,:,:,:)
      self%q(VAR_BZ,:,:,:,:) = self%q(VAR_BZ,:,:,:,:) + b(s) * self%time%dt * self%dq(VAR_BZ,:,:,:,:)
      call self%compute_residuals(q=self%q, dq=self%dq)
      self%q(VAR_DX,:,:,:,:) = self%q(VAR_DX,:,:,:,:) + a(s) * self%time%dt * self%dq(VAR_DX,:,:,:,:)
      self%q(VAR_DY,:,:,:,:) = self%q(VAR_DY,:,:,:,:) + a(s) * self%time%dt * self%dq(VAR_DY,:,:,:,:)
      self%q(VAR_DZ,:,:,:,:) = self%q(VAR_DZ,:,:,:,:) + a(s) * self%time%dt * self%dq(VAR_DZ,:,:,:,:)
   enddo
   call self%apply_fWL_correction(q=self%q)
   call self%compute_coils_current(q=self%q)
   call self%impose_div_free
   endassociate
   endsubroutine integrate_blanesmoan

   subroutine integrate_cfm(self)
   !< Integrate equation, time operator, Commutator-Free Magnus integrator.
   class(prism_cpu_object), intent(inout) :: self             !< The equation.
   real(R8P), parameter                   :: toll=1.0e-14_R8P !< CFM coefficients tollerance.
   integer(I4P)                           :: s,ss             !< Counter.

   call self%compute_coils_current(q=self%q)
   associate(dt=>self%time%dt,s_coeffs=>self%cfm%s_coeffs,e_coeffs=>self%cfm%e_coeffs)
   self%cfm%q = self%q
   call self%compute_residuals(q=self%cfm%q, dq=self%cfm%dq(:,:,:,:,:,1))
   do s=2, self%cfm%n_stages
      do ss=1, s-1
         if (abs(s_coeffs(s,ss))>toll) &
            call self%cfm%compute_exponential_update(alpha=dt*s_coeffs(s,ss),dq=self%cfm%dq(:,:,:,:,:,ss),q=self%cfm%q)
      enddo
      call self%compute_residuals(q=self%cfm%q, dq=self%cfm%dq(:,:,:,:,:,s))
   enddo
   self%cfm%q = self%q
   do s=1, self%cfm%n_stages
      if (abs(self%cfm%e_coeffs(s))>toll) &
         call self%cfm%compute_exponential_update(alpha=dt*e_coeffs(s),dq=self%cfm%dq(:,:,:,:,:,s),q=self%cfm%q)
   enddo
   self%q = self%cfm%q
   endassociate
   call self%apply_fWL_correction(q=self%q)
   call self%compute_coils_current(q=self%q)
   call self%impose_div_free
   endsubroutine integrate_cfm

   subroutine integrate_leapfrog(self)
   !< Integrate equation, time operator, leapfrog scheme.
   class(prism_cpu_object), intent(inout) :: self !< The equation.

   call self%compute_coils_current(q=self%q)
   call self%compute_residuals(q=self%q, dq=self%dq)
   call self%save_residuals
   call self%leapfrog%integrate(field=self%adam%field, dt=self%time%dt, q=self%q, dq=self%dq)
   call self%apply_fWL_correction(q=self%q)
   call self%compute_coils_current(q=self%q)
   call self%impose_div_free
   endsubroutine integrate_leapfrog

   subroutine integrate_leapfrog_pic(self)
   !< Integrate equation, time operator, leapfrog scheme for particle in cell
   class(prism_cpu_object), intent(inout) :: self !< The equation.

   !< Maxwell source terms computation: particles and coils
   call self%pic%particle_cartesian_grid_index(field=self%adam%field, grid=self%adam%grid, q_pic=self%q_pic)
   call self%pic%current_weighting(field=self%adam%field, grid=self%adam%grid, q=self%q, q_pic=self%q_pic, nv=self%nv)
   call self%verify_no_pic_deposition_on_coils(q=self%q, check_current=.true., &
                                               context='integrate_leapfrog_pic(current)')
   call self%compute_coils_current(q=self%q)
   !< Maxwell residuals computation
   call self%compute_residuals(q=self%q, dq=self%dq)
   call self%save_residuals
   !< Pic residual computation
   call self%pic%field_weighting(field=self%adam%field, grid=self%adam%grid, q=self%q, q_pic=self%q_pic, &
                                 pic_fields=self%pic_fields, nv=self%nv)
   !< Integration of equations
   call self%leapfrog%integrate(field=self%adam%field, dt=self%time%dt, q=self%q, dq=self%dq)
   call self%leapfrog_pic%integrate(dt=self%time%dt, q_pic=self%q_pic, pic_fields=self%pic_fields)
   call self%apply_fWL_correction(q=self%q)
   call self%impose_div_free
   endsubroutine integrate_leapfrog_pic

   subroutine integrate_rk_ls(self)
   !< Integrate equation, time operator, RK classical low storage schemes.
   !< Low storage RK working on q_rk(:,:,:,:,:,1)/q as stages, update q in place.
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   call self%compute_coils_current(q=self%q) !da modificare per avere i tempi corretti
   call self%rk%initialize_stages(field=self%adam%field, q=self%q)
   do s=1, self%rk%nrk
      call self%compute_residuals(q=self%q, dq=self%dq)
      if (s==1) call self%save_residuals
      if (self%ib%solids_number>0) then
         call self%rk%compute_stage_ls(field=self%adam%field, s=s,dt=self%time%dt,phi=self%ib%phi,dq=self%dq,q=self%q)
      else
         call self%rk%compute_stage_ls(field=self%adam%field, s=s,dt=self%time%dt,dq=self%dq,q=self%q)
      endif
   enddo
   call self%apply_fWL_correction(q=self%q)
   call self%compute_coils_current(q=self%q)
   call self%impose_div_free
   endsubroutine integrate_rk_ls

   subroutine integrate_rk_ssp(self)
   !< Integrate equation, time operator, SSP RK schemes.
   !< SSP RK working on q_rk as stages.
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   if (self%external_fields%ef_type/=EF_TYPE_NONE) &
      call self%external_fields%sub_external_fields(field=self%adam%field, grid=self%adam%grid, &
                                                      time=self%time%time, dt=self%time%dt, q=self%q)
   call self%rk%initialize_stages(field=self%adam%field, q=self%q)
   if (self%pml%enabled) call self%rk_pml%initialize_stages(pml=self%pml)
   do s=1, self%rk%nrk
      if (self%ib%solids_number>0) then
         call self%rk%compute_stage(field=self%adam%field, s=s, dt=self%time%dt, phi=self%ib%phi)
      else
         call self%rk%compute_stage(field=self%adam%field, s=s, dt=self%time%dt)
      endif
      if (self%pml%enabled) call self%rk_pml%compute_stage(s=s, dt=self%time%dt)
      !call self%compute_coils_current(q=rk%q_rk(:,:,:,:,:,s), gamma=rk%gamm(s)) !Spostato in update_ghost
      call self%compute_residuals(q=self%rk%q_rk(:,:,:,:,:,s), dq=self%dq, s=s)
      !if (s==1) call self%save_residuals
      if (self%ib%solids_number>0) then
         call self%rk%assign_stage(field=self%adam%field, s=s, q=self%dq, phi=self%ib%phi)
      else
         call self%rk%assign_stage(field=self%adam%field, s=s, q=self%dq)
      endif
      if (self%pml%enabled) call self%rk_pml%assign_stage(s=s)
   enddo
   if (self%ib%solids_number>0) then
      call self%rk%update_q(field=self%adam%field, dt=self%time%dt, phi=self%ib%phi, q=self%q)
      !call self%update_q_BC(dt=self%time%dt, phi=self%ib%phi)
   else
      call self%rk%update_q(field=self%adam%field, dt=self%time%dt, q=self%q, dq=self%dq)
      !call self%update_q_BC(dt=self%time%dt)
      call self%save_residuals
   endif
   if (self%pml%enabled) call self%rk_pml%update_q_pml(dt=self%time%dt, pml=self%pml)
   call self%apply_fWL_correction(q=self%q)
   call self%compute_coils_current(q=self%q)
   call self%impose_div_free
   if (self%external_fields%ef_type/=EF_TYPE_NONE) &
      call self%external_fields%add_external_fields(field=self%adam%field, grid=self%adam%grid, &
                                                      time=self%time%time, dt=self%time%dt, q=self%q)
   endsubroutine integrate_rk_ssp

   subroutine integrate_rk_ssp_pic(self)
   !< Integrate equation, time operator, SSP RK schemes.
   !< SSP RK working on q_rk as stages.
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.
   real(R8P), allocatable                 :: q_stage(:,:,:,:,:) !< Contiguous stage scratch to avoid huge slice temporaries.

   !call sub_external_fields(self = self%external_fields, field = field, &
   !                        time = self%time%time, dt = self%time%dt, q = self%q)

   !Inizializzo stadi RK per campi e PIC
   call self%rk%initialize_stages(field=self%adam%field, q=self%q)
   call self%rk_pic%initialize_stages(q_pic=self%q_pic)
   if (self%pml%enabled) call self%rk_pml%initialize_stages(pml=self%pml)
   call allocate_variable(var=q_stage,                              &
                          ulb=reshape([1,self%nv,                   &
                                       1-self%ngc,self%ni+self%ngc, &
                                       1-self%ngc,self%nj+self%ngc, &
                                       1-self%ngc,self%nk+self%ngc, &
                                       1,self%nb],[2,5]),           &
                          msg=mpih%myrankstr//'integrate_rk_ssp_pic allocate q_stage')

   do s=1, self%rk%nrk
      !Calcolo stadio RK per campi e PIC
      if (self%ib%solids_number>0) then
         call self%rk%compute_stage(field=self%adam%field, s=s, dt=self%time%dt, phi=self%ib%phi)
      else
         call self%rk%compute_stage(field=self%adam%field, s=s, dt=self%time%dt)
      endif
      call self%rk_pic%compute_stage(s=s, dt=self%time%dt)
      if (self%pml%enabled) call self%rk_pml%compute_stage(s=s, dt=self%time%dt)
      q_stage = self%rk%q_rk(:,:,:,:,:,s)
      !Calcolo termini sorgente Maxwell da particelle e bobine
      call self%pic%particle_cartesian_grid_index(field=self%adam%field, grid=self%adam%grid, q_pic=self%rk_pic%q_pic_rk(:,:,s))
      call self%pic%current_weighting(field=self%adam%field, grid=self%adam%grid, q=q_stage, &
                                       q_pic=self%rk_pic%q_pic_rk(:,:,s), nv=self%nv)
      call self%verify_no_pic_deposition_on_coils(q=q_stage, check_current=.true., &
                                                  context='integrate_rk_ssp_pic(stage current)')
      call self%compute_coils_current(q=q_stage, gamma=self%rk%gamm(s))
      !Calcolo residui Maxwell
      call self%compute_residuals(q=q_stage, dq=self%dq, s=s)
      if (s==1) call self%save_residuals
      !Calcolo residui PIC: calcolati direttamente nell'assegnazione dello stadio RK
      !Interpolo quindi i campi (probabilmente è qui che ti conviene sommare e sottrarre i campi esterni)
      call self%pic%field_weighting(field=self%adam%field, grid=self%adam%grid, q=q_stage, &
                                    q_pic=self%rk_pic%q_pic_rk(:,:,s), pic_fields=self%pic_fields, nv=self%nv)
      !Assegno lo stadio RK per campi e PIC
      if (self%ib%solids_number>0) then
         call self%rk%assign_stage(field=self%adam%field, s=s, q=self%dq, phi=self%ib%phi)
      else
         call self%rk%assign_stage(field=self%adam%field, s=s, q=self%dq)
      endif
      call self%rk_pic%assign_stage(s=s, pic_fields=self%pic_fields)
      if (self%pml%enabled) call self%rk_pml%assign_stage(s=s)
   enddo
   ! Completo l'integrazione temporale
   if (self%ib%solids_number>0) then
      call self%rk%update_q(field=self%adam%field, dt=self%time%dt, phi=self%ib%phi, q=self%q)
      !call self%update_q_BC(dt=self%time%dt, phi=self%ib%phi)
   else
      call self%rk%update_q(field=self%adam%field, dt=self%time%dt, q=self%q)
      !call self%update_q_BC(dt=self%time%dt)
   endif
   call self%rk_pic%update_q_pic(dt=self%time%dt, q_pic=self%q_pic)
   if (self%pml%enabled) call self%rk_pml%update_q_pml(dt=self%time%dt, pml=self%pml)
   !Aggiorno i termini sorgente di Maxwell al tempo in cui andrò a plottare i risultati
   call self%apply_fWL_correction(q=self%q)
   call self%impose_div_free
   call self%pic%particle_cartesian_grid_index(field=self%adam%field, grid=self%adam%grid, q_pic=self%q_pic)
   call self%pic%current_weighting(field=self%adam%field, grid=self%adam%grid, q=self%q, q_pic=self%q_pic, nv=self%nv)
   call self%pic%particle_weighting(field=self%adam%field, grid=self%adam%grid, q=self%q, q_pic=self%q_pic, nv=self%nv)
   call self%verify_no_pic_deposition_on_coils(q=self%q, check_current=.true., check_charge=.true., &
                                               context='integrate_rk_ssp_pic(final deposition)')
   call self%compute_coils_current(q=self%q)
   !call add_external_fields(self = self%external_fields, field = field, &
   !                        time = self%time%time, dt = self%time%dt, q = self%q)
   endsubroutine integrate_rk_ssp_pic

   !subroutine update_q_BC(self, dt, phi)
   !!< Update RK q ghost cells.
   !class(prism_cpu_object), intent(inout) :: self             !< RK object.
   !real(R8P),        intent(in)           :: dt               !< Current time step.
   !real(R8P),        intent(in), optional :: phi(1:,          &
   !                                              1-self%ngc:, &
   !                                              1-self%ngc:, &
   !                                              1-self%ngc:, &
   !                                              1:)          !< IB distance.
   !integer(I4P)                           :: all_solids       !< Last phi index, all solids summary.
   !integer(I4P)                           :: i, j, k, b, v, s, c !< Counter.
!
   !integer(I4P)                           :: idelta,jdelta,kdelta    !< IJK delta step for extrapolation.
   !integer(I4P)                           :: bc_type                 !< Boundary condition type.
   !integer(I4P)                           :: crown                   !< Crown counter.
   !associate(local_map_bc_crown=>self%adam%maps%local_map_bc_crown,                                                             &
   !             nv=>self%nv, ngc=>self%ngc, q_bc_vars=>self%bc%q, dx=>self%adam%field%dxyz(1,:), dy=>self%adam%field%dxyz(2,:), &
   !             dz=>self%adam%field%dxyz(3,:), ni=>self%ni, nj=>self%nj, nk=>self%nk, dt=>self%time%dt, chi=>self%physics%chi,  &
   !             nv_c=>self%physics%nv_c, nv_cl=>self%physics%nv_cl,                                                             &
   !             constrained_transport_B=>self%numerics%constrained_transport_B,                                                 &
   !             constrained_transport_D=>self%numerics%constrained_transport_D, nrk=>self%rk_bc%nrk,                            &
   !             q_bc_rk=>self%rk_bc%q_bc_rk, blocks_number=>self%blocks_number, beta=>self%rk_bc%beta)
!
   !if (present(phi)) then
   !   all_solids = ubound(phi, dim=1)
   !   !$omp parallel do collapse(6) default(firstprivate) shared(phi,self)
   !   do s=1, nrk
   !      do b=1, blocks_number
   !         do k=1-ngc, nk+ngc
   !            do j=1-ngc, nj+ngc
   !               do i=1-ngc, ni+ngc
   !                  !(O cambiare i cicli da 1:nv_c a nv_c-nv_cl+1:nv_c)
   !                  do v=1, nv_c
   !                     if (phi(all_solids,i,j,k,b) < 0._R8P) then
   !                        q_bc_rk(v,i,j,k,b,nrk+1) = q_bc_rk(v,i,j,k,b,nrk+1) + dt * beta(s) * q_bc_rk(1,i,j,k,b,s)
   !                     endif
   !                  enddo
   !               enddo
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !   !$omp end parallel do
   !else
   !   !$omp parallel do collapse(6) default(firstprivate) shared(self)
   !   do s=1, nrk
   !      do b=1, blocks_number
   !         do k=1-ngc, nk+ngc
   !            do j=1-ngc, nj+ngc
   !               do i=1-ngc, ni+ngc
   !                  do v=1, nv_c
   !                     q_bc_rk(v,i,j,k,b,nrk+1) = q_bc_rk(v,i,j,k,b,nrk+1) + dt * beta(s) * q_bc_rk(1,i,j,k,b,s)
   !                  enddo
   !               enddo
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   !   !$omp end parallel do
   !endif
   !if (allocated(self%adam%maps%local_map_bc_crown)) then
   !   do crown=1, ngc
   !      do c=1, size(local_map_bc_crown, dim=1)
   !         b = local_map_bc_crown(c, 1 ,crown)
   !         if (b>0) then
   !            bc_type = local_map_bc_crown(c, 8 ,crown)
   !            if (bc_type == BC_radiative) then
   !               i       = local_map_bc_crown(c, 2 ,crown)
   !               j       = local_map_bc_crown(c, 3 ,crown)
   !               k       = local_map_bc_crown(c, 4 ,crown)
   !               idelta  = local_map_bc_crown(c, 5 ,crown)
   !               jdelta  = local_map_bc_crown(c, 6 ,crown)
   !               kdelta  = local_map_bc_crown(c, 7 ,crown)
   !               do v=1, nv_c
   !                  ! this seems bugged, loop over v and using always nv_c or 1...
   !                  ! self%q(nv_c,i,j,k,b) = 2*q_bc_rk(1,i,j,k,b,nrk+1)-self%q(nv_c,i-idelta,j-jdelta,k-kdelta,b)
   !                  ! probable fix below
   !                  self%q(v,i,j,k,b) = 2*q_bc_rk(v,i,j,k,b,nrk+1)-self%q(v,i-idelta,j-jdelta,k-kdelta,b)
   !               enddo
   !               !print *, 'Updating BC SM', b, ' cell (', i, ',', j, ',', k, ')'
   !            endif
   !         endif
   !         !Qua ci aggiungi gli altri if a seconda elle variabili su cui vuoi implementare questa BC
   !      enddo
   !   enddo
   !endif
   !endassociate
   !endsubroutine update_q_BC

   subroutine integrate_rk_yoshida(self)
   !< Integrate equation, time operator, Yoshida RK scheme.
   class(prism_cpu_object), intent(inout) :: self !< The equation.
   integer(I4P)                           :: s    !< Counter.

   call self%compute_coils_current(q=self%q)
   do s=1, self%rk%nrk - 1
      call self%compute_residuals(q=self%q, dq=self%dq)
      if (s==1) call self%save_residuals
      self%q(VAR_BX,:,:,:,:) = self%q(VAR_BX,:,:,:,:) + self%rk%ssa(s) * self%time%dt * self%dq(VAR_BX,:,:,:,:)
      self%q(VAR_BY,:,:,:,:) = self%q(VAR_BY,:,:,:,:) + self%rk%ssa(s) * self%time%dt * self%dq(VAR_BY,:,:,:,:)
      self%q(VAR_BZ,:,:,:,:) = self%q(VAR_BZ,:,:,:,:) + self%rk%ssa(s) * self%time%dt * self%dq(VAR_BZ,:,:,:,:)
      call self%compute_residuals(q=self%q, dq=self%dq)
      self%q(VAR_DX,:,:,:,:) = self%q(VAR_DX,:,:,:,:) + self%rk%ssb(s) * self%time%dt * self%dq(VAR_DX,:,:,:,:)
      self%q(VAR_DY,:,:,:,:) = self%q(VAR_DY,:,:,:,:) + self%rk%ssb(s) * self%time%dt * self%dq(VAR_DY,:,:,:,:)
      self%q(VAR_DZ,:,:,:,:) = self%q(VAR_DZ,:,:,:,:) + self%rk%ssb(s) * self%time%dt * self%dq(VAR_DZ,:,:,:,:)
   enddo
   call self%compute_residuals(q=self%q, dq=self%dq)
   self%q(VAR_BX,:,:,:,:) = self%q(VAR_BX,:,:,:,:) + self%rk%ssa(self%rk%nrk) * self%time%dt * self%dq(VAR_BX,:,:,:,:)
   self%q(VAR_BY,:,:,:,:) = self%q(VAR_BY,:,:,:,:) + self%rk%ssa(self%rk%nrk) * self%time%dt * self%dq(VAR_BY,:,:,:,:)
   self%q(VAR_BZ,:,:,:,:) = self%q(VAR_BZ,:,:,:,:) + self%rk%ssa(self%rk%nrk) * self%time%dt * self%dq(VAR_BZ,:,:,:,:)
   call self%apply_fWL_correction(q=self%q)
   call self%compute_coils_current(q=self%q)
   call self%impose_div_free
   endsubroutine integrate_rk_yoshida

   ! non TBP
   subroutine compute_dxyz_min(blocks_number, dxyz, dxyz_min)
   !< Compute minimum dxyz space step.
   integer(I4P), intent(in)  :: blocks_number !< Number of blocks.
   real(R8P),    intent(in)  :: dxyz(:,:)     !< XYZ space steps.
   real(R8P),    intent(out) :: dxyz_min      !< Minimum space step.
   integer(I4P)              :: b             !< Counter.

   dxyz_min = huge(0._R8P)
   !$omp parallel do shared(dxyz) reduction(min:dxyz_min)
   do b=1, blocks_number
      dxyz_min = min(dxyz_min, dxyz(1,b), dxyz(2,b), dxyz(3,b))
   enddo
   dxyz_min = dxyz_min * 0.5_R8P
   endsubroutine compute_dxyz_min

   subroutine compute_fluxes_convective_weno(dir,blocks_number,ni,nj,nk,ngc,nv_c,weno_s,weno_a,weno_p,weno_d,weno_c,weno_zeps,&
                                             evmax,erw,elw,chi,q,fluxes)
   !< Compute convective fluxes along direction `dir`, WENO scheme for space operator.
   integer(I4P), intent(in)    :: dir                                !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)    :: blocks_number                      !< Number of blocks.
   integer(I4P), intent(in)    :: ni                                 !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                 !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                 !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                !< Ghost cells number.
   integer(I4P), intent(in)    :: nv_c                               !< Number of conservative varibales.
   integer(I4P), intent(in)    :: weno_s                             !< Weno stencils number/dimension.
   real(R8P),    intent(in)    :: weno_a(1:,0:,1:)                   !< Optimal weights.
   real(R8P),    intent(in)    :: weno_p(1:,0:,0:,1:)                !< Polinomials coefficients.
   real(R8P),    intent(in)    :: weno_d(0:,0:,0:,1:)                !< Smoothness indicators coefficients.
   real(R8P),    intent(in)    :: weno_c(1-weno_s:,1:)               !< Centered polinomials coefficients.
   real(R8P),    intent(in)    :: weno_zeps                          !< Parameter for avoiding division by zero in computing IS.
   real(R8P),    intent(in)    :: evmax                              !< Maximum waves speed estimation.
   real(R8P),    intent(in)    :: erw(1:,1:,1:)                      !< Right eigenvectors for WENO reconstruction.
   real(R8P),    intent(in)    :: elw(1:,1:,1:)                      !< Left  eigenvectors for WENO reconstruction.
   real(R8P),    intent(in)    :: chi                                !< Speed parameter for divergence cleaning.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Field variables.
   real(R8P),    intent(inout) :: fluxes(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Fluxes.
   integer(I4P)                :: si(3), si_i, si_j, si_k            !< Directional (1=x,2=y,3=z) increment.
   real(R8P)                   :: sir(3)                             !< Directional (1=x,2=y,3=z) increment, real.
   integer(I4P)                :: b, i, j, k                         !< Counter.

   select case(dir)
   case(1)
      si = [1,0,0]
   case(2)
      si = [0,1,0]
   case(3)
      si = [0,0,1]
   endselect
   sir = real(si,R8P)
   si_i = 1-si(1)
   si_j = 1-si(2)
   si_k = 1-si(3)

   do b=1, blocks_number
   do k=si_k, nk
   do j=si_j, nj
   do i=si_i, ni
      call compute_fluxes_convective_ri_weno(dir=dir,b=b,i=i,j=j,k=k,ngc=ngc,nv_c=nv_c,                 &
                                             weno_s=weno_s, weno_zeps=weno_zeps,                        &
                                             weno_a=weno_a, weno_p=weno_p, weno_d=weno_d, weno_c=weno_c,&
                                             evmax=evmax,erw=erw,elw=elw,chi=chi,                       &
                                             si=si,sir=sir,q=q,fluxes=fluxes)
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_fluxes_convective_weno

   subroutine compute_fluxes_convective_ri_weno(dir,b,i,j,k,ngc,nv_c,                         &
                                                weno_s,weno_zeps,weno_a,weno_p,weno_d,weno_c, &
                                                evmax,erw,elw,chi,si,sir,q,fluxes)
   !< Compute convective fluxes at right interface of b,i,j,k.
   integer(I4P), intent(in)    :: dir                                 !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)    :: b, i, j, k                          !< Counter.
   integer(I4P), intent(in)    :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in)    :: nv_c                                !< Number of conservative varibales in q vector.
   integer(I4P), intent(in)    :: weno_s                              !< Weno stencils number/dimension.
   real(R8P),    intent(in)    :: weno_zeps                           !< Parameter to avoid division by zero.
   real(R8P),    intent(in)    :: weno_a(1:,0:,1:)                    !< Optimal weights.
   real(R8P),    intent(in)    :: weno_p(1:,0:,0:,1:)                 !< Polinomials coefficients.
   real(R8P),    intent(in)    :: weno_d(0:,0:,0:,1:)                 !< Smoothness indicators coefficients.
   real(R8P),    intent(in)    :: weno_c(1-weno_s:,1:)                !< Centered polinomials coefficients.
   real(R8P),    intent(in)    :: evmax                               !< Maximum waves speed estimation.
   real(R8P),    intent(in)    :: erw(1:,1:,1:)                       !< Right eigenvectors for WENO reconstruction.
   real(R8P),    intent(in)    :: elw(1:,1:,1:)                       !< Left  eigenvectors for WENO reconstruction.
   real(R8P),    intent(in)    :: chi                                 !< Speed parameter for divergence cleaning.
   integer(I4P), intent(in)    :: si(3)                               !< Stencil increment.
   real(R8P),    intent(in)    :: sir(3)                              !< Stencil increment, real cast.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)       !< Fields variables.
   real(R8P),    intent(inout) :: fluxes(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Fluxes.
   real(R8P)                   :: fmpc(1:2,1-S_MAX:-1+S_MAX,1:NV_MAX) !< Fluxes -+ decomposition in c. space.
   real(R8P)                   :: fpmr(1:2,1:NV_MAX)                  !< Fluxes +- reconstructed.
   integer(I4P)                :: v, vv                               !< Counter.

   call decompose_fluxes_convective(dir=dir, si=si, sir=sir,                     &
                                    b=b, i=i, j=j, k=k, ngc=ngc, nv_c=nv_c,      &
                                    weno_s=weno_s, evmax=evmax, elw=elw, chi=chi,&
                                    q=q, fmpc=fmpc(1:2,1-weno_s:-1+weno_s,1:NV_MAX))
   do v=1, nv_c
      call weno_reconstruct_upwind(S=weno_s, weno_a=weno_a, weno_p=weno_p, weno_d=weno_d,&
                                   weno_zeps=weno_zeps, V=fmpc(1:2,1-weno_s:-1+weno_s,v), VR=fpmr(1:2,v))
   enddo
   ! back projection in conservative variables space
   do v=1, nv_c
      fluxes(v,i,j,k,b) = 0._R8P
      do vv=1,nv_c
         fluxes(v,i,j,k,b) = fluxes(v,i,j,k,b) + erw(vv,v,dir) * (fpmr(1,vv) + fpmr(2,vv))
      enddo
   enddo
   endsubroutine compute_fluxes_convective_ri_weno

   subroutine compute_fluxes_difference(blocks_number, ni, nj, nk, ngc, nv_c, var_Jx, var_Jy, var_Jz, &
                                       dx, dy, dz, flx, fly, flz, q, dq)
   !< Compute fluxes difference.
   integer(I4P), intent(in)    :: blocks_number                   !< Number of blocks.
   integer(I4P), intent(in)    :: ni                              !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                              !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                              !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                             !< Ghost cells number.
   integer(I4P), intent(in)    :: nv_c                            !< Number of conservative varibales in q.
   integer(I4P), intent(in)    :: var_Jx, var_Jy, var_Jz          !< Current variable indices.
   real(R8P),    intent(in)    :: dx(1:), dy(1:), dz(1:)          !< Space steps.
   real(R8P),    intent(in)    :: flx(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< X direction fluxes.
   real(R8P),    intent(in)    :: fly(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Y direction fluxes.
   real(R8P),    intent(in)    :: flz(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Z direction fluxes.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Fields.
   real(R8P),    intent(inout) :: dq(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Fluxes differences.
   integer(I4P)                :: b, i, j, k, v                   !< Counter.

   do b=1,blocks_number
   do k=1,nk
   do j=1,nj
   do i=1,ni
      do v=1, nv_c
         dq(v,i,j,k,b) = - (flx(v,i,j,k,b) - flx(v,i-1,j,k,b)) / dx(b) &
                         - (fly(v,i,j,k,b) - fly(v,i,j-1,k,b)) / dy(b) &
                         - (flz(v,i,j,k,b) - flz(v,i,j,k-1,b)) / dz(b)
      enddo
      ! J sources
      dq(VAR_DX,i,j,k,b) = dq(VAR_DX,i,j,k,b) - q(var_Jx,i,j,k,b)
      dq(VAR_DY,i,j,k,b) = dq(VAR_DY,i,j,k,b) - q(var_Jy,i,j,k,b)
      dq(VAR_DZ,i,j,k,b) = dq(VAR_DZ,i,j,k,b) - q(var_Jz,i,j,k,b)
      ! corrections
      ! if (d_divergence_cleaner .and. .not.b_divergence_cleaner .and. eta>0._R8P) then
      !    dq(10,i,j,k,b) = dq(10,i,j,k,b) - chi/eta*chi/eta*q(10,i,j,k,b)
      ! elseif (D_divergence_cleaner .and. B_divergence_cleaner .and. eta>0._R8P) then
      !    dq(10,i,j,k,b) = dq(10,i,j,k,b) - chi/eta*chi/eta*q(10,i,j,k,b)
      !    dq(11,i,j,k,b) = dq(11,i,j,k,b) - chi/eta*chi/eta*q(11,i,j,k,b)
      ! endif
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_fluxes_difference

   subroutine decompose_fluxes_convective(dir,si,sir,b,i,j,k,ngc,nv_c,weno_s,evmax,elw,q,chi,fmpc)
   !< Decompose convective fluxes.
   integer(I4P), intent(in)    :: dir                           !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)    :: si(3)                         !< Stencil increment.
   real(R8P),    intent(in)    :: sir(3)                        !< Stencil increment, real cast.
   integer(I4P), intent(in)    :: b, i, j, k                    !< Counter.
   integer(I4P), intent(in)    :: ngc                           !< Ghost cells number.
   integer(I4P), intent(in)    :: nv_c                          !< Number of conservative varibales in q vector.
   integer(I4P), intent(in)    :: weno_s                        !< Weno stencils number/dimension.
   real(R8P),    intent(in)    :: evmax                         !< Maximum eigenvalue.
   real(R8P),    intent(in)    :: elw(1:,1:,1:)                 !< Left eigenvectors for WENO reconstruction.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(in)    :: chi                           !< Speed coefficient for D & B div-cleaning.
   real(R8P),    intent(inout) :: fmpc(1:,1-weno_s:,1:)         !< Fluxes -+ decomposition in characteristics space.
   real(R8P)                   :: fmp(2)                        !< Fluxes -+ decomposition in each cell stencils.
   real(R8P)                   :: gc, wc                        !< Increments for fluxes decomposition.
   integer(I4P)                :: v, vv, s, is, js, ks          !< Counter.
   real(R8P)                   :: f(NV_MAX)                     !< Conservative fluxes.
   real(R8P)                   :: p(NV_MAX)                     !< Eigenvalues.

   do s=1-weno_s, weno_s
      is = i + (s) * si(1) ; js = j + (s) * si(2) ; ks = k + (s) * si(3)
      call compute_fluxes_maxwell(sir=sir,q=q(:,is,js,ks,b),f=f,chi=chi)
      !call compute_eigenvalues_vector(sir=sir, p=p)
      do v=1, nv_c
         wc = 0._R8P
         gc = 0._R8P
         do vv=1, nv_c
            wc = wc + elw(vv,v,dir) * q(vv,is,js,ks,b)
            gc = gc + elw(vv,v,dir) * f(vv)
         enddo
         fmp(2) = 0.5_R8P * (gc + evmax * wc)
         !fmp(2) = 0.5_R8P * (gc + evmax * p(v) * wc)
         fmp(1) = gc - fmp(2)
         if (s<weno_s)   fmpc(2,s  ,v) = fmp(2)
         if (s>1-weno_s) fmpc(1,s-1,v) = fmp(1)
      enddo
   enddo
   endsubroutine decompose_fluxes_convective

   subroutine write_current_behavior_tab(filename, current_density, time)
   character(len=1), parameter  :: TAB = achar(9)
   character(len=*), intent(in) :: filename
   real(R8P),        intent(in) :: current_density
   real(R8P),        intent(in) :: time
   integer(I4P)                 :: iu, ios

   open(newunit=iu, file=trim(filename), status='unknown', action='write', &
        form='formatted', position='append', iostat=ios)
   if (ios /= 0) then
     write(*,'(a,i0)') 'write_current_tab: errore open(), iostat=', ios
     error stop
   end if
   write(iu,'(ES24.16,a,ES24.16)') time, TAB, current_density
   close(iu)
   endsubroutine write_current_behavior_tab

   subroutine compute_field_mean_value(self, q, n_x, n_y, n_z, n_b, mean_value)
   !< Compute mean value of the field in a certain region of the domain.
   class(prism_cpu_object), intent(in)  :: self                                      !< The object
   real(R8P), intent(in)                :: q(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Field variables.
   integer(I4P), intent(in)             :: n_x(2), n_y(2), n_z(2), n_b(2)
                                                                                     !< Number of cells in each direction and number
                                                                                     !< of blocks
   real(R8P), intent(out)               :: mean_value                                !< Mean value of the field out of the fWLayer
   integer(I4P)                         :: i, j, k, b                                !< Counter.

   mean_value = 0._R8P
   do b=n_b(1), n_b(2)
      do k=n_z(1), n_z(2)
         do j=n_y(1), n_y(2)
            do i=n_x(1), n_x(2)
               mean_value = mean_value + q(i,j,k,b)
            enddo
         enddo
      enddo
   enddo
   mean_value = mean_value / real((n_x(2)-n_x(1)+1)*(n_y(2)-n_y(1)+1)*(n_z(2)-n_z(1)+1)*(n_b(2)-n_b(1)+1), R8P)
   endsubroutine compute_field_mean_value

   function dotproduct(a, b) result(dot)
   !< Compute the scalar (dot) product.
   real(R8P), intent(in) :: a(3) !< Left hand side.
   real(R8P), intent(in) :: b(3) !< Left hand side.
   real(R8P)             :: dot  !< Dot product.

   dot = (a(1) * b(1)) + (a(2) * b(2)) + (a(3) * b(3))
   endfunction dotproduct

   function crossproduct(a, b) result(cross)
   real(R8P), intent(in) :: a(3)     !< Left hand side.
   real(R8P), intent(in) :: b(3)     !< Left hand side.
   real(R8P)             :: cross(3) !< Cross product.

   cross(1) = (a(2) * b(3)) - (a(3) * b(2))
   cross(2) = (a(3) * b(1)) - (a(1) * b(3))
   cross(3) = (a(1) * b(2)) - (a(2) * b(1))
   endfunction crossproduct
endmodule adam_prism_cpu_object
