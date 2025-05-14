!< ADAM, Maxwell equations system class definition, common data to all backends.
module adam_prism_common_object

use adam_adam_object
use adam_amr_object
use adam_field_object
use adam_grid_object
use adam_ib_object
use adam_mpih_object
use adam_rk_object
use adam_slices_object
use adam_weno_object
use adam_prism_ic_object
use adam_prism_coil_object !aggiunto coil
use adam_prism_io_object
use adam_prism_bc_object
use adam_prism_physics_object
use adam_prism_time_object
use penf
use ISO_C_BINDING

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
   ! PRISM library objects
   type(prism_io_object)      :: io      !< IO handler.
   type(prism_physics_object) :: physics !< Fluids physiscs handler.
   type(prism_ic_object)      :: ic      !< Initial Conditions (IC) handler.
   type(prism_bc_object)      :: bc      !< Boundary Conditions (BC) handler.
   type(prism_time_object)    :: time    !< Time handler.
   type(prism_coil_object)    :: coil    !< Oggetto con informazioni su spire.
   ! grid/field data replica for easy handling
   integer(I4P), pointer :: ngc=>null()           !< Number of ghost cells.
   integer(I4P), pointer :: ni=>null()            !< Number of cells in i direction.
   integer(I4P), pointer :: nj=>null()            !< Number of cells in j direction.
   integer(I4P), pointer :: nk=>null()            !< Number of cells in k direction.
   integer(I4P), pointer :: nb=>null()            !< Total blocks number for MPI.
   integer(I4P), pointer :: blocks_number=>null() !< Actual blocks number.
!   integer(I4P), pointer :: ns=>null()            !< Number of fluids specie.
   integer(I4P), pointer :: nv=>null()            !< Number of conservative variables.
   ! auxiliary fields data: see nasto parameters definition for the arrangement of conservative and auxiliary variables
   real(R8P), pointer    :: q(:,:,:,:,:)=>null()  !< Conservative cell centered variables.

   type(c_ptr), allocatable :: ptree(:) !< CGAL trees for solids.

   contains
      procedure, pass(self) :: allocate_common   !< Allocate common data.
      procedure, pass(self) :: initialize_common !< Initialize the equation common data.
endtype prism_common_object
contains
   subroutine allocate_common(self) !se così commentata questa subroutine diviene di fatto inutile e commentabile in initialize_common
   !< Allocate common data.
   class(nasto_common_object), intent(inout) :: self !< The equation.

   associate(nv=>self%nv, ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, nb=>self%nb, &
             solids_number=>self%ib%solids_number)
   !allocate(self%q(1:nv, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
   !self%q = 0._R8P !da chiarire vedendo field ma, avendo messo il puntatore in initialize_common, non dovrebbe servirmi a nulla
             !allocare o dichiarare
   endassociate
   endsubroutine allocate_common
   

   subroutine initialize_common(self, filename, memory_avail, do_mpi_init, verbose)
   !< Initialize the equation common data.
   class(prism_common_object), intent(inout)        :: self         !< The equation.
   character(*),               intent(in)           :: filename     !< Input file name.
   real(R8P),                  intent(in)           :: memory_avail !< Memory available for single MPI process.
   logical,                    intent(in), optional :: do_mpi_init  !< Flag to activate MPI init call.
   logical,                    intent(in), optional :: verbose      !< Trigger verbose output.
   logical                                          :: verbose_     !< Trigger verbose output, local variable.
   integer(I8P)                                     :: nodes_number !< Allocated nodes on tree.
   integer(I4P)                                     :: nb           !< Number of allocated blocks.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose

   call self%mpih%initialize(do_mpi_init=do_mpi_init, verbose=verbose_)
   if (verbose_) call self%mpih%print_message('prism_common_object%initialize start')
   call self%io%initialize(filename=trim(filename))
   associate(file_parameters=>self%io%file_parameters)
   call self%bc%initialize(file_parameters=file_parameters)
   call self%physics%initialize(file_parameters=file_parameters)
   call self%adam%grid%initialize(file_parameters=file_parameters,bc_type=self%bc%bc_type, verbose=.true.)
   call self%adam%compute_blocks_number(memory_avail=memory_avail, fields_number=80, nb=nb, nodes_number=nodes_number)
   call self%adam%initialize(file_parameters=file_parameters, &
                             do_tree_init=.true.,             &
                             do_maps_init=.true.,             &
                             do_field_init=.true.,            &
                             nv=self%physics%nv, nb=nb, nodes_number=nodes_number)
   call associate_adam_data(grid=self%adam%grid, field=self%adam%field, physics=self%physics)
   call self%adam%refine_uniform(refinement_levels=self%adam%tree%iu_ref_levels, do_blocks_reorder=.false.)
   call self%adam%prune(ijkl_prune=self%adam%tree%ijkl_prune, do_blocks_reorder=.false.)
   call self%amr%initialize(file_parameters=file_parameters)
   call self%time%initialize(file_parameters=file_parameters)
   call self%ic%initialize(file_parameters=file_parameters)
   call self%coil%initialize(file_parameters=file_parameters)
   call self%ib%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field)
   call self%slices%initialize(file_parameters=file_parameters)
   call self%rk%initialize(file_parameters=file_parameters, grid=self%grid, field=self%field)
   call self%weno%initialize(file_parameters=file_parameters, nb=self%nb, ngc=self%ngc, ni=self%ni, nj=self%nj, nk=self%nk)
   endassociate
   call self%allocate_common
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
!      self%ns            => physics%ns
      self%nv            => physics%nv
!      self%nv_aux        => physics%nv_aux
      self%q             => field%q !Ho cambiato modo di scrivere q: è necessario cambiarlo in field (identificandolo come un target)?
      endsubroutine associate_adam_data
   endsubroutine initialize_common
endmodule adam_prism_common_object