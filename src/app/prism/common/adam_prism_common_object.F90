!< ADAM, Maxwell equations system class definition, common data to all backends.
module adam_prism_common_object

! ADAM modules
use adam_adam_object
use adam_amr_object
use adam_blanes_moan_object
use adam_cfm_object
use adam_field_object
use adam_flail_object
use adam_grid_object
use adam_ib_object
use adam_leapfrog_object
use adam_mpih_object
use adam_rk_object
use adam_slices_object
use adam_weno_object
! PRISM modules
use adam_prism_bc_object
use adam_prism_coil_object
use adam_prism_external_fields_object
use adam_prism_fWLayer_object
use adam_prism_ic_object
use adam_prism_io_object
use adam_prism_numerics_object
use adam_prism_physics_object
use adam_prism_pic_object
use adam_prism_rk_bc_object
use adam_prism_time_object
! third party modules
use penf
! use ISO_C_BINDING

implicit none
private
public :: prism_common_object

type :: prism_common_object
   !< Maxwell equations system class definition, common data to all backends.
   ! ADAM library objects
   type(mpih_object)           :: mpih          !< MPI handler.
   type(adam_object)           :: adam          !< ADAM.
   type(field_object), pointer :: field=>null() !< The field.
   type(grid_object),  pointer :: grid=>null()  !< The grid.
   type(amr_object)            :: amr           !< AMR marker handler.
   type(ib_object)             :: ib            !< Immersed Boundary (IB) handler.
   type(slices_object)         :: slices        !< Slices handler.
   type(blanesmoan_object)     :: blanesmoan    !< Blanes-Moan integrator.
   type(cfm_object)            :: cfm           !< Commutator-Free Magnus integrator.
   type(leapfrog_object)       :: leapfrog      !< Leapfrog integrator.
   type(rk_object)             :: rk            !< RK integrator.
   type(weno_object)           :: weno          !< WENO reconstructor.
   type(flail_object)          :: flail         !< Linear algebra methods handler.
   ! PRISM library objects
   type(prism_io_object)              :: io              !< IO handler.
   type(prism_numerics_object)        :: numerics        !< Numerics handler.
   type(prism_physics_object)         :: physics         !< Fluids physiscs handler.
   type(prism_ic_object)              :: ic              !< Initial Conditions (IC) handler.
   type(prism_bc_object)              :: bc              !< Boundary Conditions (BC) handler.
   type(prism_rk_bc_object)           :: rk_bc           !< RK integrator for BC.
   type(prism_time_object)            :: time            !< Time handler.
   type(prism_fWLayer_object)         :: fWLayer         !< fWLayer handler.
   type(prism_coil_object)            :: coil            !< Coils handler.
   type(prism_external_fields_object) :: external_fields !< External fields handler.
   type(prism_pic_object)             :: pic             !< Particle-in-Cell (PIC) handler.
   ! grid/field data replica for easy handling
   integer(I4P), pointer :: ngc=>null()           !< Number of ghost cells.
   integer(I4P), pointer :: ni=>null()            !< Number of cells in i direction.
   integer(I4P), pointer :: nj=>null()            !< Number of cells in j direction.
   integer(I4P), pointer :: nk=>null()            !< Number of cells in k direction.
   integer(I4P), pointer :: nb=>null()            !< Total blocks number for MPI.
   integer(I4P), pointer :: blocks_number=>null() !< Actual blocks number.
   integer(I4P), pointer :: nv=>null()            !< Number of variables in q vector.
   integer(I4P), pointer :: nv_c=>null()          !< Number of conservative variables in q vector.
   integer(I4P), pointer :: nv_s=>null()          !< Number of source variables in q vector.
   integer(I4P), pointer :: nv_cl=>null()         !< Number of divergence cleaning variables in q vector.
   ! fields data [1:nv,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb].
   real(R8P),    allocatable :: q(         :,:,:,:,:) !< Conservative cell centered variables.
   real(R8P),    allocatable :: dq(        :,:,:,:,:) !< Residuals right hand side.
   real(R8P),    allocatable :: q_pic(           :,:) !< PIC variables.
   real(R8P),    allocatable :: dq_pic(          :,:) !< PIC variables derivatives.
   real(R8P),    allocatable :: curl(      :,:,:,:,:) !< Curl fields.
   real(R8P),    allocatable :: divergence(:,:,:,:,:) !< Divergence fields.
   character(3), allocatable :: q_name(:)             !< Fields names [1:nv].
   ! auxiliary data
   real(R8P), allocatable :: energy_D(:)                !< Energy of field D, time history.
   real(R8P), allocatable :: energy_B(:)                !< Energy of field B, time history.
   real(R8P)              :: rms_energy_error_D=0.0_R8P !< RMS energy error of D field.
   real(R8P)              :: rms_energy_error_B=0.0_R8P !< RMS energy error of B field.
   contains
      procedure, pass(self) :: allocate_common   !< Allocate common data.
      procedure, pass(self) :: initialize_common !< Initialize the equation common data.
endtype prism_common_object
contains
   subroutine allocate_common(self)
   !< Allocate common data.
   class(prism_common_object), intent(inout) :: self !< The equation.

   associate(nv=>self%nv, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, nb=>self%nb)
   call allocate_variable(var=self%dq,               &
                          ulb=reshape([1,nv,         &
                                       1-ngc,ni+ngc, &
                                       1-ngc,nj+ngc, &
                                       1-ngc,nk+ngc, &
                                       1,nb],[2,5]), &
                          msg=self%mpih%myrankstr//'prism_common_object%allocate_common(dq) ', verbose=.true.)
   self%dq = 0._R8P
   call allocate_variable(var=self%divergence,        &
                          ulb=reshape([1,self%nv,    &
                                       1-ngc,ni+ngc, &
                                       1-ngc,nj+ngc, &
                                       1-ngc,nk+ngc, &
                                       1,nb],[2,5]), &
                          msg=self%mpih%myrankstr//'prism_common_object%allocate_common(divergence) ', verbose=.true.)
   self%divergence = 0._R8P
   self%divergence(4,:,:,:,:) = self%fWLayer%f(1,:,:,:,:)
   self%divergence(5,:,:,:,:) = self%fWLayer%f(2,:,:,:,:)
   self%divergence(6,:,:,:,:) = self%fWLayer%f(3,:,:,:,:)
   call allocate_variable(var=self%curl,             &
                          ulb=reshape([1,self%nv,    &
                                       1-ngc,ni+ngc, &
                                       1-ngc,nj+ngc, &
                                       1-ngc,nk+ngc, &
                                       1,nb],[2,5]), &
                          msg=self%mpih%myrankstr//'prism_common_object%allocate_common(curl) ', verbose=.true.)
   self%curl = 0._R8P

   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) then
      call allocate_variable(var=self%q_pic,           &
                             ulb=reshape([1,self%pic%particle_number,  &
                                          1,8],[2,2]), &
                             msg=self%mpih%myrankstr//'prism_common_object%allocate_common(q_pic) ', verbose=.true.)
      self%q_pic = 0._R8P
      call allocate_variable(var=self%dq_pic,          &
                             ulb=reshape([1,self%pic%particle_number,  &
                                          1,8],[2,2]), &
                             msg=self%mpih%myrankstr//'prism_common_object%allocate_common(dq_pic) ', verbose=.true.)
      self%dq_pic = 0._R8P
   endif
   endassociate
   endsubroutine allocate_common

   subroutine initialize_common(self, field, filename, memory_avail, do_mpi_init, verbose)
   !< Initialize the equation common data.
   class(prism_common_object), intent(inout), target :: self         !< The equation.
   type(field_object),         intent(inout)         :: field        !< The field.
   character(*),               intent(in)            :: filename     !< Input file name.
   real(R8P),                  intent(in),value      :: memory_avail !< Memory available for single MPI process.
   logical,                    intent(in), optional  :: do_mpi_init  !< Flag to activate MPI init call.
   logical,                    intent(in), optional  :: verbose      !< Trigger verbose output.
   logical                                           :: verbose_     !< Trigger verbose output, local variable.
   integer(I8P)                                      :: nodes_number !< Allocated nodes on tree.
   integer(I4P)                                      :: nb           !< Number of allocated blocks.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   call self%mpih%initialize(do_mpi_init=do_mpi_init, verbose=verbose_)
   if (verbose_) call self%mpih%print_message('prism_common_object%initialize start')
   call self%io%initialize(filename=trim(filename))
   associate(file_parameters=>self%io%file_parameters)
   call self%bc%initialize(file_parameters=file_parameters)
   call self%numerics%initialize(file_parameters=file_parameters)
   call self%physics%initialize(file_parameters=file_parameters, reconstruction_vars=self%numerics%reconstruction_vars, &
                                 div_corr_var=self%numerics%div_corr_var,                                               &
                                 constrained_transport_D=self%numerics%constrained_transport_D,                         &
                                 constrained_transport_B=self%numerics%constrained_transport_B)
   if (self%physics%physical_model == PIC_PHYSICAL_MODEL) call self%pic%initialize(file_parameters=file_parameters)
   call self%adam%grid%initialize(file_parameters=file_parameters,bc_type=self%bc%bc_type, verbose=.true.)
   call self%adam%compute_blocks_number(memory_avail=memory_avail, fields_number=80, nb=nb, nodes_number=nodes_number)
   call self%adam%initialize(file_parameters=file_parameters, &
                             do_tree_init=.true.,             &
                             do_maps_init=.true.,             &
                             do_field_init=.true.,            &
                             nv=self%physics%nv, nb=1, nodes_number=11_I8P, q=self%q) !nb = nb !nodes_number = nodes_number
   call associate_adam_data(grid=self%adam%grid, field=self%adam%field, physics=self%physics)
   call self%adam%refine_uniform(refinement_levels=self%adam%tree%iu_ref_levels, do_blocks_reorder=.false.,q=self%q)
   call self%adam%prune(ijkl_prune=self%adam%tree%ijkl_prune, do_blocks_reorder=.false.,q=self%q)
   call self%amr%initialize(file_parameters=file_parameters)
   call self%time%initialize(file_parameters=file_parameters)
   call self%ic%initialize(file_parameters=file_parameters)
   call self%fWLayer%initialize(file_parameters=file_parameters, physics=self%physics, field=self%field)
   call self%coil%initialize(file_parameters=file_parameters, field=self%field)
   call self%ib%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field)
   call self%slices%initialize(file_parameters=file_parameters)
   call self%external_fields%initialize(file_parameters=file_parameters)
   if (self%numerics%scheme_time==NUM_SCHEME_TIME_RUNGE_KUTTA) &
   call self%rk%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field)
   if (self%numerics%scheme_time==NUM_SCHEME_TIME_RUNGE_KUTTA) &
   call self%rk_bc%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field, rk=self%rk, &
                              physics=self%physics)
   if (self%numerics%scheme_time==NUM_SCHEME_TIME_LEAPFROG) &
   call self%leapfrog%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field)
   if (self%numerics%scheme_time==NUM_SCHEME_TIME_BLANES_MOAN) &
   call self%blanesmoan%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field)
   if (self%numerics%scheme_time==NUM_SCHEME_TIME_CFM) &
   call self%cfm%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field)
   if (self%numerics%scheme_space==NUM_SCHEME_SPACE_WENO) &
   call self%weno%initialize(file_parameters=file_parameters, nb=self%nb, ngc=self%ngc, ni=self%ni, nj=self%nj, nk=self%nk)
   call self%flail%initialize(file_parameters=file_parameters)
   call check_ngc_number
   call self%allocate_common
   call io_initialize
   endassociate
   if (verbose_) call self%mpih%print_message('prism_common_object%initialize finish')
   contains
      subroutine associate_adam_data(grid, field, physics)
      !< Associate objects data to equation for easy handling.
      type(grid_object),          intent(in), target :: grid    !< The grid.
      type(field_object),         intent(in), target :: field   !< The field.
      type(prism_physics_object), intent(in), target :: physics !< The physics.

      self%grid          => grid
      self%field         => field
      self%blocks_number => field%blocks_number
      self%ni            => field%grid%ni
      self%nj            => field%grid%nj
      self%nk            => field%grid%nk
      self%ngc           => field%grid%ngc
      self%nb            => field%nb
      self%nv            => physics%nv
      self%nv_c          => physics%nv_c
      self%nv_s          => physics%nv_s
      self%nv_cl         => physics%nv_cl
      !self%nv_pic        => physics%nv_pic
      endsubroutine associate_adam_data

      subroutine check_ngc_number
      !< Check if the number of ghost cells is consistent with the numerical schemes used, if not an error is echoed and
      !< the simulation is stop.

      if (self%numerics%scheme_space==NUM_SCHEME_SPACE_WENO) then
         if (self%weno%S > self%grid%ngc) &
            call self%mpih%error_stop(msg=': ghost cells number (ngc) must be >= of weno stencil number (weno%S):'//&
                                      ' ngc='//trim(str(self%grid%ngc))//' weno%S='//trim(str(self%weno%S)))
      endif
      if (self%numerics%fdv_half_stencil > self%grid%ngc) &
         call self%mpih%error_stop(msg=': ghost cells number (ngc) must be >= of FDV half stencil number (numerics%fdv_hs):'//&
                                  ' ngc='//trim(str(self%grid%ngc))//' numerics%fdv_hs='//trim(str(self%numerics%fdv_half_stencil)))
      endsubroutine check_ngc_number

      subroutine io_initialize
      !< Initialize IO data.
      character(6),  allocatable :: q1_R8P_name(:) !< Variables names buffer.
      character(5),  allocatable :: q2_R8P_name(:) !< Variables names buffer.
      character(7),  allocatable :: q3_R8P_name(:) !< Variables names buffer.
      character(12), allocatable :: q4_R8P_name(:) !< Variables names buffer.
      character(7),  allocatable :: q5_R8P_name(:) !< Variables names buffer.
      integer(I4P)               :: c              !< Counter.

      call self%adam%io%initialize(grid=self%adam%grid, field=self%adam%field)
      if (self%physics%physical_model == EM_PHYSICAL_MODEL) then
         select case(self%numerics%div_corr_var)
         case(DIV_CORR_VAR_POISS)
            self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ']
            q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_Jx','res_Jy','res_Jz']
            q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09']
         case(DIV_CORR_VAR_HYPER)
            if (self%numerics%constrained_transport_D .and. .not.self%numerics%constrained_transport_B) then
               self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','phi','Jx ','Jy ','Jz ']
               q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_ph','res_Jx','res_Jy','res_Jz']
               q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09','div10']
            elseif (.not.self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
               self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','psi','Jx ','Jy ','Jz ']
               q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_ps','res_Jx','res_Jy','res_Jz']
               q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09','div10']
            elseif (self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
               self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','phi','psi','Jx ','Jy ','Jz ']
               q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_ph','res_ps','res_Jx','res_Jy','res_Jz']
               q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09','div10','div11']
            endif
         case default
            self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ']
            q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_Jx','res_Jy','res_Jz']
            q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09']
         endselect
         q5_R8P_name = ['curlD_x','curlD_y','curlD_z','curlB_x','curlB_y','curlB_z','curlJ_x','curlJ_y','curlJ_z']
         if (self%io%save_residual_fields) call self%adam%io%register_aux_field(q1_R8P=self%dq,q1_R8P_name=q1_R8P_name)
         if (self%io%save_divergence_fields) call self%adam%io%register_aux_field(q2_R8P=self%divergence,q2_R8P_name=q2_R8P_name)
         if (self%coil%total_coils_number>0) then
            q3_R8P_name = ['j_vec_1','j_vec_2','j_vec_3','f_Gauss']
            q4_R8P_name = [('coil_phi_'//trim(strz(c,3)),c=1,self%coil%total_coils_number)]
            call self%adam%io%register_aux_field(q3_R8P=self%coil%j_vec,    q3_R8P_name=q3_R8P_name)
            call self%adam%io%register_aux_field(q4_R8P=self%coil%phi,      q4_R8P_name=q4_R8P_name)
            call self%adam%io%register_aux_field(s1_I4P=self%coil%coil_flag,s1_I4P_name='coil_flag')
         endif
         if (self%io%save_curl_fields) call self%adam%io%register_aux_field(q5_R8P=self%curl,q5_R8P_name=q5_R8P_name)
      elseif(self%physics%physical_model == PIC_PHYSICAL_MODEL) then
         select case(self%numerics%div_corr_var)
         case(DIV_CORR_VAR_POISS)
            self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ','rho']
            q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_Jx','res_Jy','res_Jz','res_rh']
            q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09','div10']
         case(DIV_CORR_VAR_HYPER)
            if (self%numerics%constrained_transport_D .and. .not.self%numerics%constrained_transport_B) then
               self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','phi','Jx ','Jy ','Jz ','rho']
               q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_ph','res_Jx','res_Jy','res_Jz','res_rh']
               q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09','div10','div11']
            elseif (.not.self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
               self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','psi','Jx ','Jy ','Jz ','rho']
               q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_ps','res_Jx','res_Jy','res_Jz','res_rh']
               q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09','div10','div11']
            elseif (self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
               self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','phi','psi','Jx ','Jy ','Jz ','rho']
               q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_ph','res_ps','res_Jx','res_Jy', &
                              'res_Jz','res_rh']
               q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09','div10','div11','div12']
            endif
         case default
            self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ','rho']
            q1_R8P_name = ['res_Dx','res_Dy','res_Dz','res_Bx','res_By','res_Bz','res_Jx','res_Jy','res_Jz','res_rh']
            q2_R8P_name = ['div_D','div_B','div_J','fWL_x','fWL_y','fWL_z','div07','div08','div09','div10']
         endselect
         q5_R8P_name = ['curlD_x','curlD_y','curlD_z','curlB_x','curlB_y','curlB_z','curlJ_x','curlJ_y','curlJ_z']
         if (self%io%save_residual_fields) call self%adam%io%register_aux_field(q1_R8P=self%dq,q1_R8P_name=q1_R8P_name)
         if (self%io%save_divergence_fields) call self%adam%io%register_aux_field(q2_R8P=self%divergence,q2_R8P_name=q2_R8P_name)
         if (self%coil%total_coils_number>0) then
            q3_R8P_name = ['j_vec_1','j_vec_2','j_vec_3','f_Gauss']
            q4_R8P_name = [('coil_phi_'//trim(strz(c,3)),c=1,self%coil%total_coils_number)]
            call self%adam%io%register_aux_field(q3_R8P=self%coil%j_vec,    q3_R8P_name=q3_R8P_name)
            call self%adam%io%register_aux_field(q4_R8P=self%coil%phi,      q4_R8P_name=q4_R8P_name)
            call self%adam%io%register_aux_field(s1_I4P=self%coil%coil_flag,s1_I4P_name='coil_flag')
         endif
         if (self%io%save_curl_fields) call self%adam%io%register_aux_field(q5_R8P=self%curl,q5_R8P_name=q5_R8P_name)
      endif
      endsubroutine io_initialize
   endsubroutine initialize_common
endmodule adam_prism_common_object
