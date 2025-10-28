!< ADAM, Maxwell equations system class definition, common data to all backends.
module adam_prism_common_object

! ADAM modules
use adam_adam_object
use adam_amr_object
use adam_field_object
use adam_flail_object
use adam_grid_object
use adam_ib_object
use adam_mpih_object
use adam_rk_object
use adam_slices_object
use adam_weno_object
! PRISM modules
use adam_prism_bc_object
use adam_prism_coil_object
use adam_prism_ic_object
use adam_prism_io_object
use adam_prism_numerics_object
use adam_prism_physics_object
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
   type(rk_object)             :: rk            !< RK integrator.
   type(weno_object)           :: weno          !< WENO reconstructor.
   type(flail_object)          :: flail         !< Linear algebra methods handler.
   ! PRISM library objects
   type(prism_io_object)       :: io       !< IO handler.
   type(prism_numerics_object) :: numerics !< Numerics handler.
   type(prism_physics_object)  :: physics  !< Fluids physiscs handler.
   type(prism_ic_object)       :: ic       !< Initial Conditions (IC) handler.
   type(prism_bc_object)       :: bc       !< Boundary Conditions (BC) handler.
   type(prism_time_object)     :: time     !< Time handler.
   type(prism_coil_object)     :: coil     !< Oggetto con informazioni su spire.
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
   ! fields data
   real(R8P), allocatable    :: field_div(:,:,:,:,:) !< Field divergence.
   real(R8P), allocatable    :: q(:,:,:,:,:)         !< Conservative cell centered variables.
   real(R8P), allocatable    :: dq(:,:,:,:,:)        !< Residuals right hand side.
   character(3), allocatable :: q_name(:)            !< Fields names.
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
   call allocate_variable(var=self%field_div,        &
                          ulb=reshape([1,self%nv,    &
                                       1-ngc,ni+ngc, &
                                       1-ngc,nj+ngc, &
                                       1-ngc,nk+ngc, &
                                       1,nb],[2,5]), &
                          msg=self%mpih%myrankstr//'prism_common_object%allocate_common(field_div) ', verbose=.true.)
   self%field_div = 0._R8P
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
   call self%coil%initialize(file_parameters=file_parameters, field=self%field)
   call self%ib%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field)
   call self%slices%initialize(file_parameters=file_parameters)
   if (self%numerics%scheme_time==NUM_SCHEME_TIME_RUNGE_KUTTA) &
   call self%rk%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field)
   call self%weno%initialize(file_parameters=file_parameters, nb=self%nb, ngc=self%ngc, ni=self%ni, nj=self%nj, nk=self%nk)
   call self%flail%initialize(file_parameters=file_parameters)
   call self%allocate_common

   select case(self%numerics%div_corr_var)
   case(DIV_CORR_VAR_POISS)
      self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ']
      call self%adam%io%initialize(grid=self%adam%grid, field=self%adam%field,                                       &
                             q1_R8P=self%field_div,      q1_R8P_name=['DivD_d','DivB_d','DivJ_d','DivG0_',           &
                                                                      'DivG1_','DivG2_','DivG3_','DivG4_','DivG5_'], &
                             q2_R8P=self%coil%j_vec,     q2_R8P_name=['j_vec_1','j_vec_2','j_vec_3','f_Gauss'],      &
                             !q4_R8P=self%dq,             q4_R8P_name=['res_Dx','res_Dy','res_Dz',                    &
                             !                                         'res_Bx','res_By','res_Bz',                    &
                             !                                         'res_Jx','res_Jy','res_Jz'],                   &
                             s1_I4P=self%coil%coil_flag, s1_I4P_name='coil_flag',                                    &
                             s1_R8P=self%coil%phi(1,:,:,:,:),       s1_R8P_name='coil_phi')
   case(DIV_CORR_VAR_HYPER)
      if (self%numerics%constrained_transport_D .and. .not.self%numerics%constrained_transport_B) then
         self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ', 'phi']
         call self%adam%io%initialize(grid=self%adam%grid, field=self%adam%field,                                       &
                                q1_R8P=self%field_div,      q1_R8P_name=['DivD_d','DivB_d','DivJ_d','DivG0_',           &
                                                                         'DivG1_','DivG2_','DivG3_','DivG4_','DivG5_'], &
                                q2_R8P=self%coil%j_vec,     q2_R8P_name=['j_vec_1','j_vec_2','j_vec_3','f_Gauss'],      &
                                q4_R8P=self%dq,             q4_R8P_name=['res_Dx','res_Dy','res_Dz',                    &
                                                                         'res_Bx','res_By','res_Bz','res_ph',           &
                                                                         'res_Jx','res_Jy','res_Jz'],                   &
                                s1_I4P=self%coil%coil_flag, s1_I4P_name='coil_flag',                                    &
                                s1_R8P=self%coil%phi(1,:,:,:,:),       s1_R8P_name='coil_phi')
      elseif (.not.self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
         self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ', 'psi']
         call self%adam%io%initialize(grid=self%adam%grid, field=self%adam%field,                                       &
                                q1_R8P=self%field_div,      q1_R8P_name=['DivD_d','DivB_d','DivJ_d','DivG0_',           &
                                                                         'DivG1_','DivG2_','DivG3_','DivG4_','DivG5_'], &
                                q2_R8P=self%coil%j_vec,     q2_R8P_name=['j_vec_1','j_vec_2','j_vec_3','f_Gauss'],      &
                                q4_R8P=self%dq,             q4_R8P_name=['res_Dx','res_Dy','res_Dz',                    &
                                                                         'res_Bx','res_By','res_Bz','res_ps',           &
                                                                         'res_Jx','res_Jy','res_Jz'],                   &
                                s1_I4P=self%coil%coil_flag, s1_I4P_name='coil_flag',                                    &
                                s1_R8P=self%coil%phi(1,:,:,:,:),       s1_R8P_name='coil_phi')
      elseif (self%numerics%constrained_transport_D .and. self%numerics%constrained_transport_B) then
         self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ', 'phi', 'psi']
               call self%adam%io%initialize(grid=self%adam%grid, field=self%adam%field,                                 &
                                q1_R8P=self%field_div,      q1_R8P_name=['DivD_d','DivB_d','DivJ_d','DivG0_',           &
                                                                         'DivG1_','DivG2_','DivG3_','DivG4_','DivG5_'], &
                                q2_R8P=self%coil%j_vec,     q2_R8P_name=['j_vec_1','j_vec_2','j_vec_3','f_Gauss'],      &
                                q4_R8P=self%dq,             q4_R8P_name=['res_Dx','res_Dy','res_Dz',                    &
                                                                         'res_Bx','res_By','res_Bz','res_ph',           &
                                                                         'res_ps','res_Jx','res_Jy','res_Jz'],          &
                                s1_I4P=self%coil%coil_flag, s1_I4P_name='coil_flag',                                    &
                                s1_R8P=self%coil%phi(1,:,:,:,:),       s1_R8P_name='coil_phi')
      endif
   case default
      self%q_name = ['Dx ','Dy ','Dz ','Bx ','By ','Bz ','Jx ','Jy ','Jz ']
      call self%adam%io%initialize(grid=self%adam%grid, field=self%adam%field,                                    &
                          q1_R8P=self%field_div,      q1_R8P_name=['DivD_d','DivB_d','DivJ_d','DivG0_',           &
                                                                   'DivG1_','DivG2_','DivG3_','DivG4_','DivG5_'], &
                          q2_R8P=self%coil%j_vec,     q2_R8P_name=['j_vec_1','j_vec_2','j_vec_3','f_Gauss'],      &
                          q4_R8P=self%dq,             q4_R8P_name=['res_Dx','res_Dy','res_Dz',                    &
                                                                   'res_Bx','res_By','res_Bz',                    &
                                                                   'res_Jx','res_Jy','res_Jz'],                   &
                          s1_I4P=self%coil%coil_flag, s1_I4P_name='coil_flag',                                    &
                          s1_R8P=self%coil%phi(1,:,:,:,:),       s1_R8P_name='coil_phi')
   endselect
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
      endsubroutine associate_adam_data
   endsubroutine initialize_common
endmodule adam_prism_common_object
