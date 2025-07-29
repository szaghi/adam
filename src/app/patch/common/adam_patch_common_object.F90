!< ADAM, PATCH Poisson solver with Adpative mesh Refinement for HPC computing class definition, CPU backend.
module adam_patch_common_object
!< ADAM, PATCH Poisson solver with Adpative mesh Refinement for HPC computing class definition, CPU backend.

! ADAM modules
use adam_adam_object
use adam_amr_object
use adam_field_object
use adam_flail_object
use adam_grid_object
use adam_ib_object
use adam_mpih_object
! PATCH modules
use adam_patch_ic_object
use adam_patch_io_object
use adam_patch_bc_object
use adam_patch_time_object
! third party modules
use penf
! use ISO_C_BINDING

implicit none
private
public :: patch_common_object

type :: patch_common_object
   !< Maxwell equations system class definition, common data to all backends.
   ! ADAM library objects
   type(mpih_object)           :: mpih          !< MPI handler.
   type(adam_object)           :: adam          !< ADAM.
   type(field_object), pointer :: field=>null() !< The field.
   type(grid_object),  pointer :: grid=>null()  !< The grid.
   type(amr_object)            :: amr           !< AMR marker handler.
   type(ib_object)             :: ib            !< Immersed Boundary (IB) handler.
   type(flail_object)          :: flail         !< Linear algebra methods handler.
   ! PATCH library objects
   type(patch_io_object)      :: io   !< IO handler.
   type(patch_ic_object)      :: ic   !< Initial Conditions (IC) handler.
   type(patch_bc_object)      :: bc   !< Boundary Conditions (BC) handler.
   type(patch_time_object)    :: time !< Time handler.
   ! grid/field data replica for easy handling
   integer(I4P), pointer :: ngc=>null()           !< Number of ghost cells.
   integer(I4P), pointer :: ni=>null()            !< Number of cells in i direction.
   integer(I4P), pointer :: nj=>null()            !< Number of cells in j direction.
   integer(I4P), pointer :: nk=>null()            !< Number of cells in k direction.
   integer(I4P), pointer :: nb=>null()            !< Total blocks number for MPI.
   integer(I4P), pointer :: blocks_number=>null() !< Actual blocks number.
   integer(I4P), pointer :: nv=>null()            !< Number of conservative/primitive variables.
   ! fields data
   real(R8P), allocatable :: q(:,:,:,:,:)    !< Potential field cell centered variable.
   real(R8P), allocatable :: r(:,:,:,:,:)    !< Rho function.
   character(3)           :: q_name(1)='phi' !< Potential field name.
   contains
      procedure, pass(self) :: allocate_common   !< Allocate common data.
      procedure, pass(self) :: initialize_common !< Initialize the equation common data.
endtype patch_common_object
contains
   subroutine allocate_common(self)
   !< Allocate common data.
   class(patch_common_object), intent(inout) :: self !< The equation.

   associate(ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, nb=>self%nb)
   call allocate_variable(var=self%r,                &
                          ulb=reshape([1,self%nv,    &
                                       1-ngc,ni+ngc, &
                                       1-ngc,nj+ngc, &
                                       1-ngc,nk+ngc, &
                                       1,nb],[2,5]), &
                          msg=self%mpih%myrankstr//'prsim_common_object%allocate_common(r) ', verbose=.true.)
   self%r = 0._R8P
   endassociate
   endsubroutine allocate_common

   subroutine initialize_common(self, field, filename, memory_avail, do_mpi_init, verbose)
   !< Initialize the equation common data.
   class(patch_common_object), intent(inout), target :: self         !< The equation.
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
   if (verbose_) call self%mpih%print_message('patch_common_object%initialize start')
   call self%io%initialize(filename=trim(filename))
   associate(file_parameters=>self%io%file_parameters)
   call self%bc%initialize(file_parameters=file_parameters)
   call self%adam%grid%initialize(file_parameters=file_parameters,bc_type=self%bc%bc_type, verbose=.true.)
   call self%adam%compute_blocks_number(memory_avail=memory_avail, fields_number=80, nb=nb, nodes_number=nodes_number)
   call self%adam%initialize(file_parameters=file_parameters, &
                             do_tree_init=.true.,             &
                             do_maps_init=.true.,             &
                             do_field_init=.true.,            &
                             nv=1, nb=1, nodes_number=11_I8P, q=self%q) !nb = nb !nodes_number = nodes_number
   call associate_adam_data(grid=self%adam%grid, field=self%adam%field)
   call self%adam%refine_uniform(refinement_levels=self%adam%tree%iu_ref_levels, do_blocks_reorder=.false.,q=self%q)
   call self%adam%prune(ijkl_prune=self%adam%tree%ijkl_prune, do_blocks_reorder=.false.,q=self%q)
   call self%amr%initialize(file_parameters=file_parameters)
   call self%time%initialize(file_parameters=file_parameters)
   call self%ic%initialize(file_parameters=file_parameters)
   call self%ib%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field)
   call self%flail%initialize(file_parameters=file_parameters)
   call self%allocate_common
   call self%adam%io%initialize(grid=self%adam%grid, field=self%adam%field, q1_R8P=self%r, q1_R8P_name=['rho'])
   endassociate
   if (verbose_) call self%mpih%print_message('patch_common_object%initialize finish')
   contains
      subroutine associate_adam_data(grid, field)
      !< Associate objects data to equation for easy handling.
      type(grid_object),          intent(in), target :: grid    !< The grid.
      type(field_object),         intent(in), target :: field   !< The field.

      self%grid          => grid
      self%field         => field
      self%blocks_number => field%blocks_number
      self%ni            => field%grid%ni
      self%nj            => field%grid%nj
      self%nk            => field%grid%nk
      self%ngc           => field%grid%ngc
      self%nb            => field%nb
      self%nv            => field%nv
      endsubroutine associate_adam_data
   endsubroutine initialize_common
endmodule adam_patch_common_object
