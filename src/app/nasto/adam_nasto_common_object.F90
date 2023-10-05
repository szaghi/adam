!< ADAM, Navier-Stokes equations system class definition, common data to all backends.
module adam_nasto_common_object
!< ADAM, Navier-Stokes equations system class definition, common data to all backends.

use adam_adam_object
use adam_amr_object
use adam_field_object
use adam_grid_object
use adam_ib_object
use adam_mpih_object
use adam_slices_object
use adam_nasto_ic_object
use adam_nasto_io_object
use adam_nasto_bc_object
use adam_nasto_physics_object
use adam_nasto_parameters
use adam_nasto_schemesh_object
use adam_nasto_timeh_object
use finer
use penf
use ISO_C_BINDING

implicit none
private
public :: nasto_common_object

type :: nasto_common_object
   !< Navier-Stokes equations system class definition, common data to all backends.
   ! ADAM library objects
   type(mpih_object)           :: mpih          !< MPI handler.
   type(adam_object)           :: adam          !< ADAM.
   type(field_object), pointer :: field=>null() !< The field.
   type(grid_object),  pointer :: grid=>null()  !< The grid.
   type(amr_object)            :: amr           !< AMR marker handler.
   type(ib_object)             :: ib            !< Immersed Boundary (IB) handler.
   type(slices_object)         :: slices        !< Slices handler.
   ! NASTO library objects
   type(nasto_io_object)       :: io       !< IO handler.
   type(nasto_physics_object)  :: physics  !< Fluids physiscs handler.
   type(nasto_ic_object)       :: ic       !< Initial Conditions (IC) handler.
   type(nasto_bc_object)       :: bc       !< Boundary Conditions (BC) handler.
   type(nasto_timeh_object)    :: timeh    !< Time handler.
   type(nasto_schemesh_object) :: schemesh !< Schemes handler.
   ! grid/field data replica for easy handling
   integer(I4P), pointer :: ngc=>null()           !< Number of ghost cells.
   integer(I4P), pointer :: ni=>null()            !< Number of cells in i direction.
   integer(I4P), pointer :: nj=>null()            !< Number of cells in j direction.
   integer(I4P), pointer :: nk=>null()            !< Number of cells in k direction.
   integer(I4P), pointer :: nb=>null()            !< Total blocks number for MPI.
   integer(I4P), pointer :: blocks_number=>null() !< Actual blocks number.
   integer(I4P), pointer :: ns=>null()            !< Number of fluids specie.
   integer(I4P), pointer :: nv=>null()            !< Number of conservative variables.
   integer(I4P), pointer :: nv_aux=>null()        !< Number of auxiliary variables.
   ! auxiliary fields data: see nasto parameters definition for the arrangement of conservative and auxiliary variables
   real(R8P), allocatable :: q_aux(:,:,:,:,:) !< Auxiliary cell centered variables.
   real(R8P), allocatable ::   phi(:,:,:,:,:) !< Distance function cell centered.

   type(c_ptr), allocatable :: ptree(:) !< CGAL trees for solids.

   contains
      procedure, pass(self) :: allocate_common   !< Allocate common data.
      procedure, pass(self) :: initialize_common !< Initialize the equation common data.
endtype nasto_common_object

contains
   subroutine allocate_common(self)
   !< Allocate common data.
   class(nasto_common_object), intent(inout) :: self !< The equation.

   associate(nv_aux=>self%nv_aux, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, nb=>self%nb, &
             solids_number=>self%ib%solids_number)
   allocate(self%q_aux(1:nv_aux, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
   self%q_aux = 0._R8P
   if (solids_number > 0) then
      allocate(self%phi(1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:solids_number))
      self%phi = -1._R8P
   endif
   endassociate
   endsubroutine allocate_common

   subroutine initialize_common(self, filename, nb, nodes_number, do_mpi_init)
   !< Initialize the equation common data.
   class(nasto_common_object), intent(inout)        :: self         !< The equation.
   character(*),               intent(in)           :: filename     !< Input file name.
   integer(I8P),               intent(in)           :: nodes_number !< Allocated nodes on tree.
   integer(I4P),               intent(in)           :: nb           !< Number of allocated blocks.
   logical,                    intent(in), optional :: do_mpi_init  !< Flag to activate MPI init call.

   call self%mpih%initialize(do_mpi_init=do_mpi_init)
   print '(A)', self%mpih%myrankstr//'nasto_common_object%initialize start'
   call self%io%initialize(filename=trim(filename))
   associate(file_parameters=>self%io%file_parameters)
   call self%physics%initialize(file_parameters=file_parameters)
   call self%adam%grid%initialize(file_parameters=file_parameters, verbose=.true.)
   call self%adam%initialize(file_parameters=file_parameters, &
                             do_tree_init=.true.,                  &
                             do_field_init=.true.,                 &
                             nv=self%physics%nv, nb=nb, nodes_number=nodes_number)
   call associate_adam_data(grid=self%adam%grid, field=self%adam%field, physics=self%physics)
   call self%adam%refine_uniform(refinement_levels=self%adam%tree%iu_ref_levels, do_blocks_reorder=.false.)
   call self%adam%prune(ijkl_prune=self%adam%tree%ijkl_prune, do_blocks_reorder=.false.)
   call self%amr%initialize(file_parameters=file_parameters)
   call self%timeh%initialize(file_parameters=file_parameters)
   call self%ic%initialize(file_parameters=file_parameters)
   call self%bc%initialize(file_parameters=file_parameters, grid=self%grid)
   call self%ib%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field)
   call self%slices%initialize(file_parameters=file_parameters)
   call self%schemesh%initialize(file_parameters=file_parameters, nb=self%nb, ngc=self%ngc, ni=self%ni, nj=self%nj, nk=self%nk)
   endassociate
   call self%allocate_common
   print '(A)', self%mpih%myrankstr//'nasto_common_object%initialize finish'
   contains
      subroutine associate_adam_data(grid, field, physics)
      !< Associate objects data to equation for easy handling.
      type(grid_object),          intent(in), target :: grid    !< The grid.
      type(field_object),         intent(in), target :: field   !< The field.
      type(nasto_physics_object), intent(in), target :: physics !< The physics.

      self%grid          => grid
      self%field         => field
      self%blocks_number => field%blocks_number
      self%ni            => field%grid%ni
      self%nj            => field%grid%nj
      self%nk            => field%grid%nk
      self%ngc           => field%grid%ngc
      self%nb            => field%nb
      self%ns            => physics%ns
      self%nv            => physics%nv
      self%nv_aux        => physics%nv_aux
      endsubroutine associate_adam_data
   endsubroutine initialize_common
endmodule adam_nasto_common_object
