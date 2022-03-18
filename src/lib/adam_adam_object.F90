!< ADAM, ADAM class definition.
module adam_adam_object
!< ADAM, ADAM class definition.

use adam_base_mpi_object, only : base_mpi_object
use adam_field_object, only : field_object
use adam_grid_object, only : grid_object
use adam_parameters
use adam_tree_node_object, only : tree_node_object
use adam_tree_bucket_object, only : tree_bucket_object
use adam_tree_object, only : tree_object
use FINER, only : file_ini
use PENF
use STRINGIFOR
use VTK_FORTRAN
use HDF5
use MPI
use memorysaver

implicit none
private
public :: adam_object

type :: adam_object
   !< ADAM class definition.
   type(base_mpi_object) :: base_mpi !< The MPI backend.
   type(grid_object)     :: grid     !< The grid.
   type(tree_object)     :: tree     !< The tree.
   type(field_object)    :: field    !< The field.
   contains
      ! public methods
      procedure, pass(self) :: adapt                         !< Adapt tree/field accordingly to refine/derefine necessity.
      procedure, pass(self) :: amr_update                    !< Update AMR status.
      procedure, pass(self) :: blocks_reorder                !< Reorder blocks (for asyncrhonous MPI)
      procedure, pass(self) :: check_blocks_number           !< Check if blocks number is groving too much.
      procedure, pass(self) :: initialize                    !< Initialize ADAM.
      procedure, pass(self) :: interpolate_at_point          !< Interpolate a scalar variable at a given point.
      procedure, pass(self) :: load_restart_files            !< Load restart files.
      procedure, pass(self) :: make_comm_local_maps_ghost_bc !< Make communication/local maps of ghost cells.
      procedure, pass(self) :: mpi_gather_refinement_needed  !< Gather refinement needed.
      procedure, pass(self) :: mpi_redistribute              !< Redistribute nodes/blocks to processes, load balancing.
      procedure, pass(self) :: print_status                  !< Print status of main data.
      procedure, pass(self) :: prune                         !< Prune nodes/blocks.
      procedure, pass(self) :: refine_uniform                !< Refine all blocks uniformly.
      procedure, pass(self) :: save_hdf5                     !< Save ADAM in HDF5 format.
      procedure, pass(self) :: save_restart_files            !< Save restart files.
      procedure, pass(self) :: save_slice                    !< Save slice.
      procedure, pass(self) :: save_vtk                      !< Save ADAM in VTK  format.
endtype adam_object

contains
   ! public methods
   subroutine adapt(self)
   !< Adapt tree/field accordingly to refine/derefine necessity.
   class(adam_object), intent(inout) :: self !< ADAM.

   call self%tree%adapt
   call self%check_blocks_number
   call self%field%adapt(ratio=self%tree%ratio,                                                            &
                         block_to_refine=self%tree%block_to_refine, block_refined=self%tree%block_refined, &
                         block_to_derefine=self%tree%block_to_derefine, block_derefined=self%tree%block_derefined)
   endsubroutine adapt

   subroutine amr_update(self, is_marked_by_field, is_marked_by_tree, do_mpi_redistribute, do_blocks_reorder, &
                         print_mpi_stats, is_grid_changed)
   !< Update AMR status.
   !<
   !< Note: AMR update can be safely called only *after* update_ghost has been called for *q* variables, otherwise
   !< refine is not well done.
   !< Note: only if the AMR is UNIFORM and GLOBALLY made by tree, i.e. using mark_all_nodes, the mpi_redistribute can be avoided,
   !< otherwise mpi_gather_refinement_nedeed is not safe (having wrong nodes number counters).
   class(adam_object), intent(inout)         :: self                 !< ADAM.
   logical,            intent(in),  optional :: is_marked_by_field   !< Flag to check if marker is field.
   logical,            intent(in),  optional :: is_marked_by_tree    !< Flag to check if marker is tree.
   logical,            intent(in),  optional :: do_mpi_redistribute  !< Flag to activate MPI redistribute.
   logical,            intent(in),  optional :: do_blocks_reorder    !< Flag to activate blocks reorder.
   logical,            intent(in),  optional :: print_mpi_stats      !< Flag to activate MPI statistics print.
   logical,            intent(out), optional :: is_grid_changed      !< Flag to check if grid is changed.
   logical                                   :: do_mpi_redistribute_ !< Flag to activate MPI redistribute, local var.
   logical                                   :: do_blocks_reorder_   !< Flag to activate blocks reorder, local var.

   do_mpi_redistribute_ = .true. ; if (present(do_mpi_redistribute )) do_mpi_redistribute_ = do_mpi_redistribute
   do_blocks_reorder_ = .true. ; if (present(do_blocks_reorder)) do_blocks_reorder_ = do_blocks_reorder

   call self%mpi_gather_refinement_needed(is_marked_by_field=is_marked_by_field, is_marked_by_tree=is_marked_by_tree)

   call self%adapt

   if (present(is_grid_changed)) is_grid_changed = (size(self%tree%node_to_refine,   dim=1)>0_I4P).or.&
                                                   (size(self%tree%node_to_derefine, dim=1)>0_I4P)

   if (do_mpi_redistribute_) call self%mpi_redistribute(print_mpi_stats=print_mpi_stats)

   if (do_blocks_reorder_) call self%blocks_reorder

   call self%make_comm_local_maps_ghost_bc
   endsubroutine amr_update

   subroutine blocks_reorder(self)
   !< Reorder blocks (for asyncrhonous MPI)
   class(adam_object), intent(inout) :: self !< ADAM.

   call self%tree%blocks_reorder
   call self%field%blocks_reorder(inner_outer_block_map=self%tree%inner_outer_block_map, &
                                  inner_blocks_number=self%tree%inner_blocks_number)
   endsubroutine blocks_reorder

   subroutine check_blocks_number(self)
   !< Check if blocks number is groving too much.
   class(adam_object), intent(inout) :: self     !< ADAM.
   type(tree_node_object), pointer   :: node_ptr !< Pointer to current node.
   integer(I4P)                      :: max_nb   !< Maximum number of blocks desidered.

   max_nb = 0
   do while(self%tree%loop(node_ptr=node_ptr))
      max_nb = max(max_nb, node_ptr%block_index)
   enddo
   if (max_nb > self%field%nb) then
      print '(A)', self%base_mpi%myrankstr//'ERROR: the number of new blocks after AMR is greater than Nb'
      print '(A)', self%base_mpi%myrankstr//'max blocks numer available [Nb]: '//trim(str(self%field%nb))
      print '(A)', self%base_mpi%myrankstr//'blocks numer required after AMR: '//trim(str(node_ptr%block_index))
      call MPI_ABORT(MPI_COMM_WORLD, -101, self%base_mpi%error)
   endif
   print '(A)', self%base_mpi%myrankstr//'maximum number of blocks created after AMR update: '//str(max_nb)//'/'//str(self%field%nb)
   endsubroutine check_blocks_number

   subroutine load_restart_files(self, basename, t, time)
   !< Load restart files.
   class(adam_object), intent(inout) :: self      !< ADAM.
   character(*),       intent(in)    :: basename  !< Base name of output files.
   integer(I4P),       intent(out)   :: t         !< Time iteration.
   real(R8P),          intent(out)   :: time      !< Time.
   integer(I4P)                      :: file_unit !< Output file unit.

   open(newunit=file_unit, file=trim(adjustl(basename))//'.time', form='UNFORMATTED', access='STREAM')
   read(unit=file_unit) t, time
   close(file_unit)
   call self%tree%load_nodes(file_name=trim(adjustl(basename))//'.tnd')
   call self%field%load_blocks(basename=basename)
   endsubroutine load_restart_files

   subroutine initialize(self, file_parameters,                                              &
                         ni, nj, nk, ngc, emin, emax, bc_type, do_grid_init,                 &
                         max_load, nodes_number, buckets_number, ratio, max_level, add_adam, &
                         iu_ref_levels, i_prune, j_prune, k_prune, l_prune, do_tree_init,    &
                         nv, nb, do_field_init)
   !< Initialize ADAM.
   class(adam_object), intent(inout)           :: self               !< ADAM.
   type(file_ini),     intent(inout), optional :: file_parameters    !< INI file handler.
   ! grid options
   integer(I4P),       intent(in),    optional :: ni                 !< Number of cells in X direction.
   integer(I4P),       intent(in),    optional :: nj                 !< Number of cells in Y direction.
   integer(I4P),       intent(in),    optional :: nk                 !< Number of cells in Z direction.
   integer(I4P),       intent(in),    optional :: ngc                !< Number of ghost cells.
   real(R8P),          intent(in),    optional :: emin(3)            !< Coordinates of minium abscissa.
   real(R8P),          intent(in),    optional :: emax(3)            !< Coordinates of maxium abscissa.
   integer(I4P),       intent(in),    optional :: bc_type(6)         !< Type of boundary conditions in the 6 faces of grid.
   logical,            intent(in),    optional :: do_grid_init       !< Flag to activate grid initialize.
   ! tree options
   real(R8P),          intent(in),    optional :: max_load           !< Maximum load of tree buckets.
   integer(I8P),       intent(in),    optional :: nodes_number       !< Nodes number to be stored in the tree.
   integer(I8P),       intent(in),    optional :: buckets_number     !< Number of buckets for initialize the tree.
   integer(I4P),       intent(in),    optional :: ratio              !< Refinement ratio.
   integer(I4P),       intent(in),    optional :: max_level          !< Maximum refinement level.
   logical,            intent(in),    optional :: add_adam           !< Add ADAM node, the ancestor of all nodes.
   integer(I4P),       intent(in),    optional :: iu_ref_levels      !< Uniform initial refinement.
   integer(I4P),       intent(in),    optional :: i_prune            !< Pruning along x.
   integer(I4P),       intent(in),    optional :: j_prune            !< Pruning along y.
   integer(I4P),       intent(in),    optional :: k_prune            !< Pruning along z.
   integer(I4P),       intent(in),    optional :: l_prune            !< Pruning level.
   logical,            intent(in),    optional :: do_tree_init       !< Flag to activate tree initialize.
   ! field options
   integer(I4P),       intent(in),    optional :: nv                 !< Number of field variables.
   integer(I4P),       intent(in),    optional :: nb                 !< Number of all blocks that can be stored in field.
   logical,            intent(in),    optional :: do_field_init      !< Flag to activate field initialize.
   ! local var
   logical                                     :: do_grid_init_      !< Flag to activate grid initialize, local var.
   logical                                     :: do_tree_init_      !< Flag to activate tree initialize, local var.
   logical                                     :: do_field_init_     !< Flag to activate field initialize, local var.

   do_grid_init_  = .false. ; if (present(do_grid_init))  do_grid_init_  = do_grid_init
   do_tree_init_  = .false. ; if (present(do_tree_init))  do_tree_init_  = do_tree_init
   do_field_init_ = .false. ; if (present(do_field_init)) do_field_init_ = do_field_init
   call self%base_mpi%initialize
   print '(A)', self%base_mpi%myrankstr//'adam%initialize start'
   if (do_grid_init_) &
      call self%grid%initialize(file_parameters=file_parameters, &
                                ni=ni,                           &
                                nj=nj,                           &
                                nk=nk,                           &
                                ngc=ngc,                         &
                                emin=emin,                       &
                                emax=emax,                       &
                                bc_type=bc_type)
   if (do_tree_init_) &
      call self%tree%initialize(grid=self%grid, &
                                file_parameters=file_parameters, &
                                max_load=max_load,               &
                                nodes_number=nodes_number,       &
                                buckets_number=buckets_number,   &
                                ratio=ratio,                     &
                                max_level=max_level,             &
                                add_adam=add_adam,               &
                                iu_ref_levels=iu_ref_levels,     &
                                i_prune=i_prune,                 &
                                j_prune=j_prune,                 &
                                k_prune=k_prune,                 &
                                l_prune=l_prune)
   if (do_field_init_) &
      call self%field%initialize(grid=self%grid, file_parameters=file_parameters, nv=nv, nb=nb)
   call self%amr_update
   print '(A)', self%base_mpi%myrankstr//'adam%initialize finish'
   endsubroutine initialize

   subroutine interpolate_at_point(self, itype, point, q, qp, is_mine, p, qc, ijk, xyz, code, v)
   !< Interpolate a scalar variable at a given point.
   class(adam_object), intent(in)            :: self      !< ADAM.
   character(*),       intent(in)            :: itype     !< Type of interpolation.
   real(R8P),          intent(in)            :: point(3)  !< Interpolation point xyz coordinates.
   real(R8P),          intent(in)            :: q(1:,              &
                                                  1-self%grid%ngc:,&
                                                  1-self%grid%ngc:,&
                                                  1-self%grid%ngc:,&
                                                  1:)     !< Q variables to be interpolated.
   real(R8P),          intent(out)           :: qp(1:)    !< Q variables interpolated at given point.
   logical,            intent(out)           :: is_mine   !< Flag to check if point interpolation belongs to myrank.
   integer(I4P),       intent(in),  optional :: p         !< Power parameter.
   real(R8P),          intent(out), optional :: qc(1:,1:) !< Closest cells q-variable values.
   integer(I4P),       intent(out), optional :: ijk(3,8)  !< Closest cells indexes.
   real(R8P),          intent(out), optional :: xyz(3,8)  !< Closest cells center-coordinates.
   integer(I8P),       intent(out), optional :: code      !< Closest block Morton code.
   integer(I4P),       intent(out), optional :: v         !< Closest vertex index.
   integer(I8P)                              :: code_     !< Closest block Morton code, local var.
   type(tree_node_object), pointer           :: node      !< Pointer to node.
   integer(I4P)                              :: ijk_(3,8) !< Closest cells indexes, local var.
   real(R8P)                                 :: xyz_(3,8) !< Closest cells center-coordinates, local var.
   integer(I4P)                              :: i, j, k   !< Counter.

   code_ = self%tree%get_closest_block(point=point)
   node => self%tree%node(code=code_)
   if (node%myrank == self%base_mpi%myrank) then
      is_mine = .true.
      call self%tree%get_closest_cells(point=point, code=code_, ijk=ijk_, xyz=xyz_)
      select case(trim(itype))
      case('inverse_distance')
         call inverse_distance_interpolation
      case('trilinear')
         call trilinear_interpolation
      endselect
      if (present(qc)) then
         do j=1,8
           do i=1,size(q, dim=1)
              qc(i,j) = q(i,ijk_(1,j),ijk_(2,j),ijk_(3,j),node%block_index)
           enddo
         enddo
      endif
      if (present(ijk))  ijk  = ijk_
      if (present(xyz))  xyz  = xyz_
      if (present(code)) code = code_
   else
      is_mine = .false.
   endif
   contains
      subroutine inverse_distance_interpolation
      !< Compute inverse distance interpolation.
      real(R8P)    :: w       !< Interpolation weights of closest cells.
      real(R8P)    :: ws      !< Sum of interpolation weights of closest cells.
      integer(I4P) :: p_      !< Power parameter, local var.

      p_ = 3 ; if (present(p)) p_ = p
      qp = 0._R8P
      ws = 0._R8P
      do j=1, 8
         w = 0._R8P
         do i=1, 3
            w = w + (xyz_(i,j) - point(i)) ** 2 ! square distance from point
         enddo
         w = 1._R8P / w ** p_
         ws = ws + w
         do k=1,size(q, dim=1)
            qp(k) = qp(k) + w * q(k,ijk_(1,j),ijk_(2,j),ijk_(3,j),node%block_index)
         enddo
      enddo
      qp = qp / ws
      endsubroutine inverse_distance_interpolation

      subroutine trilinear_interpolation
      !< Compute trilinear interpolation.
      real(R8P)    :: q1, q2, q3 !< Linear distances.
      real(R8P)    :: p1, p2, p3 !< Linear distances complements.
      real(R8P)    :: qx( 4)     !< X linear interpolations.
      real(R8P)    :: qxy(2)     !< XY linear interpolations.

      q1 = (point(1)-xyz_(1,1))/(xyz_(1,8)-xyz_(1,1))
      q2 = (point(2)-xyz_(2,1))/(xyz_(2,8)-xyz_(2,1))
      q3 = (point(3)-xyz_(3,1))/(xyz_(3,8)-xyz_(3,1))
      p1 = 1._R8P - q1
      p2 = 1._R8P - q2
      p3 = 1._R8P - q3

      do i=1, size(q, dim=1)
         qx( 1) = p1*q(i,ijk_(1,1),ijk_(2,1),ijk_(3,1),node%block_index) + q1*q(i,ijk_(1,2),ijk_(2,2),ijk_(3,2),node%block_index)
         qx( 2) = p1*q(i,ijk_(1,3),ijk_(2,3),ijk_(3,3),node%block_index) + q1*q(i,ijk_(1,4),ijk_(2,4),ijk_(3,4),node%block_index)
         qx( 3) = p1*q(i,ijk_(1,5),ijk_(2,5),ijk_(3,5),node%block_index) + q1*q(i,ijk_(1,6),ijk_(2,6),ijk_(3,6),node%block_index)
         qx( 4) = p1*q(i,ijk_(1,7),ijk_(2,7),ijk_(3,7),node%block_index) + q1*q(i,ijk_(1,8),ijk_(2,8),ijk_(3,8),node%block_index)
         qxy(1) = p2*qx(1) + q2*qx(2)
         qxy(2) = p2*qx(3) + q2*qx(4)

         qp(i) = p3*qxy(1) + q3*qxy(2)
      enddo
      endsubroutine trilinear_interpolation
   endsubroutine interpolate_at_point

   subroutine make_comm_local_maps_ghost_bc(self)
   !< Make communication/local maps of ghost cells and boundary conditions.
   class(adam_object), intent(inout) :: self !< ADAM.

   call self%tree%make_comm_local_maps_ghost
   call self%tree%make_local_maps_bc
   call self%field%prepare_comm_local_ghost(local_map_ghost         = self%tree%local_map_ghost,         &
                                            comm_map_n_send_ghost   = self%tree%comm_map_n_send_ghost,   &
                                            comm_map_n_recv_ghost   = self%tree%comm_map_n_recv_ghost,   &
                                            comm_map_send_ptr_ghost = self%tree%comm_map_send_ptr_ghost, &
                                            comm_map_recv_ptr_ghost = self%tree%comm_map_recv_ptr_ghost, &
                                            comm_map_send_ghost     = self%tree%comm_map_send_ghost,     &
                                            comm_map_recv_ghost     = self%tree%comm_map_recv_ghost)
   call self%field%prepare_local_bc(local_map_bc_face   = self%tree%local_map_bc_face, &
                                    local_map_bc_edge   = self%tree%local_map_bc_edge, &
                                    local_map_bc_corner = self%tree%local_map_bc_corner)
   endsubroutine make_comm_local_maps_ghost_bc

   subroutine mpi_gather_refinement_needed(self, is_marked_by_field, is_marked_by_tree)
   !< Gather refinement needed.
   class(adam_object), intent(inout)         :: self                !< ADAM.
   logical,            intent(in),  optional :: is_marked_by_field  !< Flag to check if marker is field.
   logical,            intent(in),  optional :: is_marked_by_tree   !< Flag to check if marker is tree.
   logical                                   :: is_marked_by_field_ !< Flag to check if marker is field, local var.
   logical                                   :: is_marked_by_tree_  !< Flag to check if marker is tree, local var.

   is_marked_by_field_ = .false. ; if (present(is_marked_by_field)) is_marked_by_field_ = is_marked_by_field
   is_marked_by_tree_  = .false. ; if (present(is_marked_by_tree )) is_marked_by_tree_  = is_marked_by_tree

   if (is_marked_by_field_) then
      call self%field%mpi_gather_refinements_needed
      call self%tree%import_refinements_needed(refinements_needed_all=self%field%refinements_needed_all, &
                                               disp_count=self%field%disp_count)
   endif

   if (is_marked_by_tree_) then
      call self%tree%mpi_gather_nodes_data(node_member='refinement_needed')
   endif
   endsubroutine mpi_gather_refinement_needed

   subroutine mpi_redistribute(self, print_mpi_stats)
   !< Redistribute nodes/blocks to processes, load balancing.
   class(adam_object), intent(inout)         :: self             !< ADAM.
   logical,            intent(in),  optional :: print_mpi_stats  !< Flag to activate MPI statistics print.
   logical                                   :: print_mpi_stats_ !< Flag to activate MPI statistics print, local var.

   print_mpi_stats_ = .false. ; if (present(print_mpi_stats)) print_mpi_stats_ = print_mpi_stats
   call self%tree%mpi_redistribute
   if (print_mpi_stats_) call self%tree%mpi_print_stats
   call self%field%mpi_redistribute(comm_map_send=self%tree%comm_map_send,         &
                                    comm_map_recv=self%tree%comm_map_recv,         &
                                    comm_map_send_ptr=self%tree%comm_map_send_ptr, &
                                    comm_map_recv_ptr=self%tree%comm_map_recv_ptr, &
                                    local_map=self%tree%local_map,                 &
                                    coordinates=self%tree%block_coordinates,       &
                                    code=self%tree%block_code)
   endsubroutine mpi_redistribute

   subroutine print_status(self)
   !< Print status of main data.
   class(adam_object), intent(in) :: self !< Adam.

   call self%grid%print_status
   call self%tree%print_status
   call self%field%print_status
   endsubroutine print_status

   subroutine prune(self, ijkl_prune, print_mpi_stats, do_blocks_reorder)
   !< Prune nodes/blocks.
   class(adam_object), intent(inout)        :: self               !< Adam.
   integer(I4P),       intent(inout)        :: ijkl_prune(4)      !< Maximum coordinates after which the prune operates.
   logical,            intent(in), optional :: print_mpi_stats    !< Flag to activate MPI statistics print.
   logical,            intent(in), optional :: do_blocks_reorder  !< Flag to activate blocks reorder.
   logical                                  :: do_blocks_reorder_ !< Flag to activate blocks reorder, local var.

   do_blocks_reorder_ = .true.  ; if (present(do_blocks_reorder)) do_blocks_reorder_ = do_blocks_reorder

   call self%tree%prune(ijkl_prune=ijkl_prune)

   call self%mpi_redistribute(print_mpi_stats=print_mpi_stats)

   if (do_blocks_reorder_) call self%blocks_reorder

   call self%make_comm_local_maps_ghost_bc
   endsubroutine prune

   subroutine refine_uniform(self, refinement_levels, do_mpi_redistribute, do_blocks_reorder)
   !< Refine all blocks uniformly.
   class(adam_object), intent(inout)        :: self                 !< Adam.
   integer(I4P),       intent(in)           :: refinement_levels    !< Number of refinement to be performed.
   logical,            intent(in), optional :: do_mpi_redistribute  !< Flag to activate MPI redistribute.
   logical,            intent(in), optional :: do_blocks_reorder    !< Flag to activate blocks reorder.
   integer(I4P)                             :: l                    !< Counter.

   print '(A)', self%base_mpi%myrankstr//'uniformly refine mesh with '//trim(str(refinement_levels))//' levels'
   do l=1, refinement_levels
      call self%tree%mark_all_nodes(mark=TO_BE_REFINED)
      call self%amr_update(do_mpi_redistribute=do_mpi_redistribute, do_blocks_reorder=do_blocks_reorder)
   enddo
   endsubroutine refine_uniform

   subroutine save_restart_files(self, basename, t, time)
   !< Save restart files.
   class(adam_object), intent(in) :: self      !< ADAM.
   character(*),       intent(in) :: basename  !< Base name of output files.
   integer(I4P),       intent(in) :: t         !< Time iteration.
   real(R8P),          intent(in) :: time      !< Time.
   integer(I4P)                   :: file_unit !< Output file unit.

   if (self%base_mpi%myrank==0) then
      open(newunit=file_unit, file=trim(adjustl(basename))//'.time', form='UNFORMATTED', access='STREAM')
      write(unit=file_unit) t, time
      close(file_unit)
      call self%tree%save_nodes(file_name=trim(adjustl(basename))//'.tnd')
   endif
   call self%field%save_blocks(basename=basename)
   endsubroutine save_restart_files

   subroutine save_hdf5(self, basename, q, q_aux, q_name, q_aux_name, directory, with_ghost, with_cell_morton, t, time)
   !< Save ADAM in HDF5 format.
   class(adam_object), intent(inout)        :: self                    !< ADAM.
   character(*),       intent(in)           :: basename                !< Base name of output files.
   real(R8P),          intent(in)           :: q(1:,              &
                                                 1-self%grid%ngc:,&
                                                 1-self%grid%ngc:,&
                                                 1-self%grid%ngc:,&
                                                 1:)                   !< Q variables to be saved.
   real(R8P),          intent(in), optional :: q_aux(1:,              &
                                                     1-self%grid%ngc:,&
                                                     1-self%grid%ngc:,&
                                                     1-self%grid%ngc:,&
                                                     1:)               !< Q auxiliary variables to be saved.
   character(*),       intent(in), optional :: q_name(:)               !< Q variables names.
   character(*),       intent(in), optional :: q_aux_name(:)           !< Q auxiliary variables names.
   character(*),       intent(in), optional :: directory               !< Directory name of output files.
   logical,            intent(in), optional :: with_ghost              !< Flag to save ghost cells.
   logical,            intent(in), optional :: with_cell_morton        !< Flag to save Morton code also in cells.
   integer(I4P),       intent(in), optional :: t                       !< Time iteration.
   real(R8P),          intent(in), optional :: time                    !< Time.
   character(:), allocatable                :: q_name_(:)              !< Q variables names, local var.
   character(:), allocatable                :: q_aux_name_(:)          !< Q auxiliary variables names, local var.
   character(:), allocatable                :: directory_              !< Directory name of output files, local var.
   logical                                  :: with_ghost_             !< Flag to save ghost cells, local var.
   logical                                  :: with_cell_morton_       !< Flag to save Morton code also in cells, local var.
   type(tree_node_object), pointer          :: node                    !< Pointer to node.
   real(R8P)                                :: emin(3)                 !< Minimum abscissa of current block.
   real(R8P)                                :: emax(3)                 !< Maximum abscissa of current block.
   integer(I8P)                             :: b                       !< Counter.
   integer(I4P)                             :: v                       !< Counter.
   integer(I4P)                             :: xdmf_unit               !< XDMF file handler.
   character(len=:), allocatable            :: h5_file_name            !< H5 Dataset name.
   integer(HID_T)                           :: h5_file_id              !< H5 File identifier.
   integer(HID_T)                           :: h5_dspace_id            !< H5 Dataspace identifier.
   real(R8P)                                :: dxyz(3)                 !< Space steps.
   integer(I4P)                             :: ngc                     !< Ghost cells saved.
   integer(I8P), allocatable                :: codes(:)                !< Codes list, sorted by level.
   integer(I8P)                             :: c                       !< Codes counter.
   integer(I4P)                             :: node_level              !< Node level counter.
   integer(I4P)                             :: i, j, k, l              !< Counter.
   integer(I4P)                             :: ijk(3,2)                !< Blocks extents.

   if (present(q_name)) then
      allocate(character(len(q_name(1))):: q_name_(size(q, dim=1)))
      q_name_ = q_name
   else
      allocate(character(4):: q_name_(size(q, dim=1)))
      do v=1, size(q, dim=1)
         q_name_(v) = 'q-'//trim(strz(v,2))
      enddo
   endif
   if (present(q_aux_name).and.present(q_aux)) then
      allocate(character(len(q_aux_name(1))):: q_aux_name_(size(q_aux, dim=1)))
      q_aux_name_ = q_aux_name
   elseif (present(q_aux)) then
      allocate(character(7):: q_aux_name_(size(q_aux, dim=1)))
      do v=1, size(q_aux, dim=1)
         q_aux_name_(v) = 'qaux-'//trim(strz(v,2))
      enddo
   endif
   directory_ = '' ; if (present(directory)) directory_ = trim(directory)
   with_ghost_ = .false. ; if (present(with_ghost)) with_ghost_ = with_ghost
   with_cell_morton_ = .false. ; if (present(with_cell_morton)) with_cell_morton_ = with_cell_morton
   if (with_ghost_) then
      ngc = self%grid%ngc
   else
      ngc = 0_I4P
   endif

   associate(ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk)
   ijk(1,:) = [1-ngc,ni+ngc]
   ijk(2,:) = [1-ngc,nj+ngc]
   ijk(3,:) = [1-ngc,nk+ngc]
   endassociate

   ! save H5 file (one for each process)
   call open_hdf5(h5_file_name=directory_//trim(basename)//'-proc'//trim(strz(self%base_mpi%myrank,6))//'.h5', &
                  ni=int(ijk(1,2)-ijk(1,1)+1,I8P),                                                             &
                  nj=int(ijk(2,2)-ijk(2,1)+1,I8P),                                                             &
                  nk=int(ijk(3,2)-ijk(3,1)+1,I8P),                                                             &
                  h5_file_id=h5_file_id,                                                                       &
                  h5_dspace_id=h5_dspace_id)
   ! save all blocks in process
   do b=1, self%field%blocks_number
      call save_hdf5_block(h5_file_id=h5_file_id,                                          &
                           h5_dspace_id=h5_dspace_id,                                      &
                           myrank=self%base_mpi%myrank,                                    &
                           code=self%field%code(b),                                        &
                           block_index=b,                                                  &
                           ii=ijk(1,:),                                                    &
                           jj=ijk(2,:),                                                    &
                           kk=ijk(3,:),                                                    &
                           q=q(:,ijk(1,1):ijk(1,2),ijk(2,1):ijk(2,2),ijk(3,1):ijk(3,2),b), &
                           q_name=q_name_,                                                 &
                           with_cell_morton=with_cell_morton_,                             &
                           q_aux_name=q_aux_name_,                                         &
                           q_aux=q_aux(:,ijk(1,1):ijk(1,2),ijk(2,1):ijk(2,2),ijk(3,1):ijk(3,2),b))
   enddo
   call close_hdf5(h5_file_id=h5_file_id, h5_dspace_id=h5_dspace_id)

   ! save XDMF file (only master process does)
   if (self%base_mpi%myrank == 0_I4P) then
      call open_xdmf(file_name=directory_//trim(basename)//'.xdmf', file_unit=xdmf_unit)
      codes = self%tree%codes(sort_by_level=.true.)
      node_level = self%tree%level(code=codes(1))
      write(xdmf_unit, '(A)') '      <Grid Name="level-'//trim(str(node_level, .true.))//'" GridType="Collection">'
      do c=1, size(codes, dim=1)
         node => self%tree%node(code=codes(c))
         if (self%tree%level(code=node%code) > node_level) then
            write(xdmf_unit, '(A)') '      </Grid>'
            node_level = self%tree%level(code=node%code)
            write(xdmf_unit, '(A)') '      <Grid Name="level-'//trim(str(node_level, .true.))//'" GridType="Collection">'
         endif
         h5_file_name = directory_//trim(basename)//'-proc'//trim(strz(node%myrank,6))//'.h5'
         call self%tree%morton_to_coordinates(code=node%code, i=i, j=j, k=k, l=l)
         call self%grid%compute_metrics(coordinates=[i,j,k,l], emin=emin, emax=emax, dx=dxyz(1), dy=dxyz(2), dz=dxyz(3))
         emin = [emin(1)-ngc*dxyz(1),emin(2)-ngc*dxyz(2),emin(3)-ngc*dxyz(3)]
         call save_xdmf_block(file_unit=xdmf_unit,                                                &
                              h5_file_name=h5_file_name,                                          &
                              myrank=node%myrank,                                                 &
                              code=node%code,                                                     &
                              block_index=node%block_index,                                       &
                              emin=emin,                                                          &
                              dxyz=dxyz,                                                          &
                              nijk=[ijk(1,2)-ijk(1,1)+2,ijk(2,2)-ijk(2,1)+2,ijk(3,2)-ijk(3,1)+2], &
                              q_name=q_name_,                                                     &
                              with_cell_morton=with_cell_morton_,                                 &
                              q_aux_name=q_aux_name_,                                             &
                              t=t,                                                                &
                              time=time)
      enddo
      write(xdmf_unit, '(A)') '      </Grid>'
      call close_xdmf(file_unit=xdmf_unit)
   endif
   endsubroutine save_hdf5

   subroutine save_slice(self, itype, points, basename, q, q_name, phi, t, time)
   !< Save slice.
   class(adam_object), intent(inout)        :: self                  !< ADAM.
   character(*),       intent(in)           :: itype                 !< Type of interpolation.
   real(R8P),          intent(in)           :: points(1:,1:,1:,1:)   !< Interpolation points coordinates [1:3,1:ni,1:nj,1:nk].
   character(*),       intent(in)           :: basename              !< Base name of output files.
   real(R8P),          intent(in)           :: q(1:,              &
                                                 1-self%grid%ngc:,&
                                                 1-self%grid%ngc:,&
                                                 1-self%grid%ngc:,&
                                                 1:)                 !< Q variables to be saved.
   character(*),       intent(in), optional :: q_name(:)             !< Variables names.
   real(R8P),          intent(in), optional :: phi(1:,              &
                                                   1-self%grid%ngc:,&
                                                   1-self%grid%ngc:,&
                                                   1-self%grid%ngc:) !< Distance function.
   integer(I4P),       intent(in), optional :: t                     !< Time iteration.
   real(R8P),          intent(in), optional :: time                  !< Time.
   character(:), allocatable                :: q_name_(:)            !< Variables names, local var.
   character(:), allocatable                :: q_aux_name_(:)        !< Q auxiliary variables names, l. var.
   real(R8P)                                :: qp(1:size(q,dim=1))   !< Q variables interpolated at given point.
   logical                                  :: is_mine               !< Flag to check if point interpolation belongs to myrank.
   integer(I4P)                             :: nijkv(4)              !< Points grid dimensions.
   integer(I4P)                             :: i, j, k, v, p         !< Counter.
   integer(I4P)                             :: MPI_IO_FILE_unit      !< Unit file.
   integer(MPI_OFFSET_KIND)                 :: offset, offset_head   !< MPI record offset.
   integer(I4P)                             :: ijkc(3,8)             !< IJK indexes of cells containing interpolation point.
   integer(I8P)                             :: code                  !< Morton code of block containing interpolation point.
   integer(I4P)                             :: ijk                   !< IJK counter (linearized).
   integer(I4P)                             :: error                 !< MPI error traping flag.
   type(tree_node_object), pointer          :: node                  !< Pointer to node.

   if (present(q_name)) then
      allocate(character(len(q_name(1))):: q_name_(size(q, dim=1)))
      q_name_ = q_name
   else
      allocate(character(4):: q_name_(size(q, dim=1)))
      do v=1, size(q, dim=1)
         q_name_(v) = 'q-'//trim(strz(v,2))
      enddo
   endif
   p = 0 ; if (present(phi)) p = 1
   nijkv(1) = size(points, dim=2)
   nijkv(2) = size(points, dim=3)
   nijkv(3) = size(points, dim=4)
   nijkv(4) = 3 + size(q, dim=1) + p

   call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(adjustl(basename))//'.mat', &
                      MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, MPI_IO_FILE_unit, error)
   offset_head = 0
   if (self%base_mpi%myrank==0) then
      ! call MPI_FILE_WRITE_AT_ALL(MPI_IO_FILE_unit, offset_head, nijkv, 4, MPI_INT, MPI_STATUS_IGNORE, error)
      call MPI_FILE_WRITE_AT(MPI_IO_FILE_unit, offset_head, nijkv, 4, MPI_INT, MPI_STATUS_IGNORE, error)
   endif
   offset_head = offset_head + 4 * 4
   ijk = 0
   do k=lbound(points, dim=4), ubound(points, dim=4)
      do j=lbound(points, dim=3), ubound(points, dim=3)
         do i=lbound(points, dim=2), ubound(points, dim=2)
            call self%interpolate_at_point(itype=itype, point=points(:,i,j,k), q=q, qp=qp, is_mine=is_mine, ijk=ijkc, code=code)
            if (is_mine) then
               offset = offset_head + ijk * 8 * (3 + size(q, dim=1) + p)
               ! call MPI_FILE_WRITE_AT_ALL(MPI_IO_FILE_unit, offset, points(:,i,j,k), 3, MPI_REAL8, MPI_STATUS_IGNORE, error)
               call MPI_FILE_WRITE_AT(MPI_IO_FILE_unit, offset, points(:,i,j,k), 3, MPI_REAL8, MPI_STATUS_IGNORE, error)
               offset = offset + 8 * 3
               ! call MPI_FILE_WRITE_AT_ALL(MPI_IO_FILE_unit, offset, qp, size(q, dim=1), MPI_REAL8, MPI_STATUS_IGNORE, error)
               call MPI_FILE_WRITE_AT(MPI_IO_FILE_unit, offset, qp, size(q, dim=1), MPI_REAL8, MPI_STATUS_IGNORE, error)
               if (present(phi)) then
                  node => self%tree%node(code=code)
                  offset = offset + 8 * size(q, dim=1)
                  ! call MPI_FILE_WRITE_AT_ALL(MPI_IO_FILE_unit, offset, phi(node%block_index,ijkc(1,1),ijkc(2,1),ijkc(3,1)), 1, &
                  !                            MPI_REAL8, MPI_STATUS_IGNORE, error)
                  call MPI_FILE_WRITE_AT(MPI_IO_FILE_unit, offset, phi(node%block_index,ijkc(1,1),ijkc(2,1),ijkc(3,1)), 1, &
                                         MPI_REAL8, MPI_STATUS_IGNORE, error)
               endif
            endif
            ijk = ijk + 1
         enddo
      enddo
   enddo
   call MPI_FILE_CLOSE(MPI_IO_FILE_unit, error)
   endsubroutine save_slice

   subroutine save_vtk(self, basename, q, q_aux, q_name, q_aux_name, directory, with_ghost, with_cell_morton, t, time)
   !< Save ADAM in VTK files.
   class(adam_object), intent(inout)        :: self                                          !< ADAM.
   character(*),       intent(in)           :: basename                                      !< Base name of output files.
   real(R8P),          intent(in)           :: q(1:,              &
                                                 1-self%grid%ngc:,&
                                                 1-self%grid%ngc:,&
                                                 1-self%grid%ngc:,&
                                                 1:)                                         !< Q variables to be saved.
   real(R8P),          intent(in), optional :: q_aux(1:,              &
                                                     1-self%grid%ngc:,&
                                                     1-self%grid%ngc:,&
                                                     1-self%grid%ngc:,&
                                                     1:)                                     !< Q auxiliary variables to be saved.
   character(*),       intent(in), optional :: directory                                     !< Output directory name.
   character(*),       intent(in), optional :: q_name(:)                                     !< Variables names.
   character(*),       intent(in), optional :: q_aux_name(:)                                 !< Q auxiliary variables names.
   logical,            intent(in), optional :: with_ghost                                    !< Flag to save ghost cells.
   logical,            intent(in), optional :: with_cell_morton                              !< Flag to save Morton code in cells.
   integer(I4P),       intent(in), optional :: t                                             !< Time iteration.
   real(R8P),          intent(in), optional :: time                                          !< Time.
   character(:), allocatable                :: directory_                                    !< Output directory name, local var.
   character(:), allocatable                :: q_name_(:)                                    !< Variables names, local var.
   character(:), allocatable                :: q_aux_name_(:)                                !< Q auxiliary variables names, l. var.
   logical                                  :: with_ghost_                                   !< Flag to save ghost cells, local var.
   logical                                  :: with_cell_morton_                             !< Flag to save Morton code in cells.
   type(vtk_file)                           :: vtk                                           !< VTK file handler.
   type(vtm_file)                           :: vtm                                           !< VTM file handler.
   type(tree_node_object), pointer          :: node                                          !< Pointer to node.
   integer(I4P)                             :: b, l, v                                       !< Counter.
   integer(I4P)                             :: i, j, k                                       !< Counter.
   integer(I4P)                             :: max_level                                     !< Maximum level.
   integer(I4P)                             :: ngc                                           !< Ghost cells saved.
   real(R8P)                                :: x(0-self%grid%ngc:self%grid%ni+self%grid%ngc) !< X coordinates.
   real(R8P)                                :: y(0-self%grid%ngc:self%grid%nj+self%grid%ngc) !< Y coordinates.
   real(R8P)                                :: z(0-self%grid%ngc:self%grid%nk+self%grid%ngc) !< Z coordinates.

   if (present(q_name)) then
      allocate(character(len(q_name(1))):: q_name_(size(q, dim=1)))
      q_name_ = q_name
   else
      allocate(character(4):: q_name_(size(q, dim=1)))
      do v=1, size(q, dim=1)
         q_name_(v) = 'q-'//trim(strz(v,2))
      enddo
   endif
   if (present(q_aux_name).and.present(q_aux)) then
      allocate(character(len(q_aux_name(1))):: q_aux_name_(size(q_aux, dim=1)))
      q_aux_name_ = q_aux_name
   endif
   directory_ = '' ; if (present(directory)) directory_ = trim(directory)
   with_ghost_ = .false. ; if (present(with_ghost)) with_ghost_ = with_ghost
   with_cell_morton_ = .false. ; if (present(with_cell_morton)) with_cell_morton_ = with_cell_morton
   if (with_ghost_) then
      ngc = self%grid%ngc
   else
      ngc = 0_I4P
   endif
   associate(ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, error=>self%base_mpi%error)
      max_level = 0_I4P
      vtr_loop : do b=1, self%field%blocks_number
         call self%grid%compute_metrics(coordinates=self%field%coordinates(:,b), x_node=x, y_node=y, z_node=z)
         max_level = max(max_level, self%field%coordinates(4,b))
         error = vtk%initialize(format='raw', filename=directory_//trim(basename)//                              &
                                                       '-morton-'//trim(str(self%field%code(b),.true.))//        &
                                                       '-block-'//trim(str(b,.true.))//                          &
                                                       '-proc-'//trim(str(self%base_mpi%myrank,.true.))//'.vtr', &
                                mesh_topology='RectilinearGrid',                                                 &
                                nx1=0-ngc, nx2=ni+ngc, ny1=0-ngc, ny2=nj+ngc, nz1=0-ngc, nz2=nk+ngc)
                            error = vtk%xml_writer%write_fielddata(action='open')
                            error = vtk%xml_writer%write_fielddata(data_name='Morton', x=self%field%code(b))
                            error = vtk%xml_writer%write_fielddata(data_name='myrank', x=self%base_mpi%myrank)
         if (present(t))    error = vtk%xml_writer%write_fielddata(data_name='t', x=t)
         if (present(time)) error = vtk%xml_writer%write_fielddata(data_name='time', x=time)
                            error = vtk%xml_writer%write_fielddata(action='close')
         error = vtk%xml_writer%write_piece(nx1=0-ngc, nx2=ni+ngc, ny1=0-ngc, ny2=nj+ngc, nz1=0-ngc, nz2=nk+ngc)
         error = vtk%xml_writer%write_geo(x=x(0-ngc:ni+ngc), y=y(0-ngc:nj+ngc), z=z(0-ngc:nk+ngc))
         error = vtk%xml_writer%write_dataarray(location='cell', action='open')
         do v=1, size(q, dim=1)
            error = vtk%xml_writer%write_dataarray(data_name=trim(q_name_(v)), &
                                                   x=[q(v,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,b)])
         enddo
         if (present(q_aux)) then
            do v=1, size(q_aux, dim=1)
               error = vtk%xml_writer%write_dataarray(data_name=trim(q_aux_name_(v)), &
                                                      x=[q_aux(v,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,b)])
            enddo
         endif
         if (with_cell_morton_) then
            error = vtk%xml_writer%write_dataarray(data_name='morton', &
                                                   x=reshape([(self%field%code(b),i=1,(ni+2*ngc)*(nj+2*ngc)*(nk+2*ngc))], &
                                                             [ni+2*ngc,nj+2*ngc,nk+2*ngc]))
         endif
         error = vtk%xml_writer%write_dataarray(location='cell', action='close')
         error = vtk%xml_writer%write_piece()
         error = vtk%finalize()
      enddo vtr_loop

      ! save VTM file (only master process does)
      if (self%base_mpi%myrank == 0_I4P) then
         error = vtm%initialize(filename=directory_//trim(basename)//'.vtm', scratch_units_number=max_level)
         vtm_group_loop : do l=1, max_level
            error = vtm%write_block(scratch=l, action='open', name='level-'//trim(str(l,.true.)))
         enddo vtm_group_loop
         vtm_filenames_loop : do while(self%tree%loop(node_ptr=node))
            b = node%block_index
            l = self%tree%level(code=node%code)
            error = vtm%write_block(scratch=l, action='write', filename=trim(basename)//                                   &
                                                                        '-morton-'//trim(str(self%field%code(b),.true.))// &
                                                                        '-block-'//trim(str(b,.true.))//                   &
                                                                        '-proc-'//trim(str(self%base_mpi%myrank,.true.))//'.vtr')
         enddo vtm_filenames_loop
         error = vtm%finalize()
      endif
   endassociate
   endsubroutine save_vtk

   ! non TBP
   ! HDF5
   subroutine close_hdf5(h5_file_id, h5_dspace_id)
   !< Close HDF5 file.
   integer(HID_T), intent(in) :: h5_file_id   !< H5 File identifier.
   integer(HID_T), intent(in) :: h5_dspace_id !< H5 Dataspace identifier.
   integer(I4P)               :: error        !< Error traping flag.

   ! terminate access to the data space
   call h5sclose_f(h5_dspace_id, error)
   ! close the file
   call h5fclose_f(h5_file_id, error)
   ! close FORTRAN interface
   call h5close_f(error)
   endsubroutine close_hdf5

   subroutine open_hdf5(h5_file_name, ni, nj, nk, h5_file_id, h5_dspace_id)
   !< Open HDF5 file.
   character(*),   intent(in)  :: h5_file_name     !< H5 file name.
   integer(I8P),   intent(in)  :: ni, nj, nk       !< Blocks dimensions.
   integer(HID_T), intent(out) :: h5_file_id       !< H5 File identifier.
   integer(HID_T), intent(out) :: h5_dspace_id     !< H5 Dataspace identifier.
   integer(I4P)                :: error            !< Error traping flag.

   ! open fortran interface
   call h5open_f(error)
   ! create a new file using default properties
   call h5fcreate_f(h5_file_name, H5F_ACC_TRUNC_F, h5_file_id, error)
   ! create the dataspace for 3D fields
   call h5screate_simple_f(3_I4P, [ni,nj,nk], h5_dspace_id, error)
   endsubroutine open_hdf5

   subroutine save_hdf5_block(h5_file_id, h5_dspace_id, &
                              myrank, code, block_index, ii, jj, kk, q, q_name, with_cell_morton, q_aux_name, q_aux)
   !< Save block into HDF5 file.
   integer(HID_T),            intent(in)           :: h5_file_id                     !< H5 File identifier.
   integer(HID_T),            intent(in)           :: h5_dspace_id                   !< H5 Dataspace identifier.
   integer(I4P),              intent(in)           :: myrank                         !< MPI rank.
   integer(I8P),              intent(in)           :: code                           !< Block Morton code.
   integer(I8P),              intent(in)           :: block_index                    !< Block index.
   integer(I4P),              intent(in)           :: ii(2)                          !< First and last i indexes.
   integer(I4P),              intent(in)           :: jj(2)                          !< First and last j indexes.
   integer(I4P),              intent(in)           :: kk(2)                          !< First and last k indexes.
   real(R8P),                 intent(in)           :: q(1:,ii(1):,jj(1):,kk(1):)     !< Q variables to be saved.
   character(*),              intent(in)           :: q_name(:)                      !< Q variables names.
   logical,                   intent(in)           :: with_cell_morton               !< Flag to save Morton code also in cells.
   character(*), allocatable, intent(in)           :: q_aux_name(:)                  !< Q auxiliary variables names.
   real(R8P),                 intent(in), optional :: q_aux(1:,ii(1):,jj(1):,kk(1):) !< Q auxiliary variables to be saved.
   character(len=:), allocatable                   :: h5_dset_name                   !< H5 Dataset name.
   integer(HID_T)                                  :: h5_dset_id                     !< H5 Dataset identifier.
   integer(I4P)                                    :: v, i                           !< Counter.
   integer(I4P)                                    :: error                          !< Error traping flag.

   do v=1, size(q, dim=1)
      h5_dset_name = trim(q_name(v))//'-'//trim(str(myrank,.true.))//'-'//trim(str(block_index,.true.))
      call h5dcreate_f(h5_file_id, h5_dset_name, H5T_NATIVE_DOUBLE, h5_dspace_id, h5_dset_id, error)
      call h5dwrite_f(h5_dset_id, H5T_NATIVE_DOUBLE, q(v,ii(1):ii(2),jj(1):jj(2),kk(1):kk(2)), &
                      [int(ii(2)-ii(1)+1,I8P),int(jj(2)-jj(1)+1,I8P),int(kk(2)-kk(1)+1,I8P)], error)
      call h5dclose_f(h5_dset_id, error)
   enddo
   if (present(q_aux)) then
      do v=1, size(q_aux, dim=1)
         h5_dset_name = trim(q_aux_name(v))//'-'//trim(str(myrank,.true.))//'-'//trim(str(block_index,.true.))
         call h5dcreate_f(h5_file_id, h5_dset_name, H5T_NATIVE_DOUBLE, h5_dspace_id, h5_dset_id, error)
         call h5dwrite_f(h5_dset_id, H5T_NATIVE_DOUBLE, q_aux(v,ii(1):ii(2),jj(1):jj(2),kk(1):kk(2)), &
                         [int(ii(2)-ii(1)+1,I8P),int(jj(2)-jj(1)+1,I8P),int(kk(2)-kk(1)+1,I8P)], error)
         call h5dclose_f(h5_dset_id, error)
      enddo
   endif
   if (with_cell_morton) then
      h5_dset_name = 'morton-'//trim(str(myrank,.true.))//'-'//trim(str(block_index,.true.))
      call h5dcreate_f(h5_file_id, h5_dset_name, H5T_NATIVE_DOUBLE, h5_dspace_id, h5_dset_id, error)
      call h5dwrite_f(h5_dset_id, H5T_NATIVE_DOUBLE,                                                  &
                      reshape([(real(code,R8P),i=1,(ii(2)-ii(1)+1)*(ii(2)-ii(1)+1)*(ii(2)-ii(1)+1))], &
                              [ii(2)-ii(1)+1,jj(2)-jj(1)+1,kk(2)-kk(1)+1]),                           &
                      [int(ii(2)-ii(1)+1,I8P),int(jj(2)-jj(1)+1,I8P),int(kk(2)-kk(1)+1,I8P)], error)
      call h5dclose_f(h5_dset_id, error)
   endif
   endsubroutine save_hdf5_block

   ! XDMF
   subroutine close_xdmf(file_unit)
   !< Close XDMF file.
   integer(I4P), intent(in) :: file_unit !< XDMF file unit.

   write(file_unit, '(A)') '    </Grid>'
   write(file_unit, '(A)') '  </Domain>'
   write(file_unit, '(A)') '</Xdmf>'
   close(file_unit)
   endsubroutine close_xdmf

   subroutine open_xdmf(file_name, file_unit)
   !< Open XDMF file.
   character(*), intent(in)  :: file_name !< XDMF file name.
   integer(I4P), intent(out) :: file_unit !< XDMF file unit.

   open(newunit=file_unit, file=trim(file_name))
   write(file_unit, '(A)') '<?xml version="1.0" encoding="utf-8"?>'
   write(file_unit, '(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="3.0">'
   write(file_unit, '(A)') '  <Domain>'
   write(file_unit, '(A)') '    <Grid Name="ADAM" GridType="Collection">'
   endsubroutine open_xdmf

   subroutine save_xdmf_block(file_unit, h5_file_name, myrank, code, block_index, emin, dxyz, nijk, &
                              q_name, with_cell_morton, q_aux_name, t, time)
   !< Save XDMF block.
   integer(I4P),              intent(in)           :: file_unit        !< XDMF file unit.
   character(*),              intent(in)           :: h5_file_name     !< H5 file name.
   integer(I4P),              intent(in)           :: myrank           !< MPI rank.
   integer(I8P),              intent(in)           :: code             !< Block Morton code.
   integer(I8P),              intent(in)           :: block_index      !< Block index.
   real(R8P),                 intent(in)           :: emin(3)          !< Block minimum extents.
   real(R8P),                 intent(in)           :: dxyz(3)          !< Block space steps.
   integer(I4P),              intent(in)           :: nijk(3)          !< Block dimensions.
   character(*),              intent(in)           :: q_name(:)        !< Q variables names.
   logical,                   intent(in)           :: with_cell_morton !< Flag to save Morton code also in cells.
   character(*), allocatable, intent(in)           :: q_aux_name(:)    !< Q auxiliary variables names.
   integer(I4P),              intent(in), optional :: t                !< Time iteration.
   real(R8P),                 intent(in), optional :: time             !< Time.
   character(:), allocatable                       :: h5_dset_name     !< Dataset name.
   integer(I4P)                                    :: v                !< Counter.

   write(file_unit, '(A)') '        <Grid Name="block-'//trim(str(code, .true.))//'">'
   write(file_unit, '(A)') '          <Geometry Origin="" Type="ORIGIN_DXDYDZ">'
   write(file_unit, '(A)') '            <DataItem DataType="Float" Dimensions="3" Format="XML" Precision="8">'//&
                      trim(str([emin(3),emin(2),emin(1)],separator=' '))//'</DataItem>'
   write(file_unit, '(A)') '            <DataItem DataType="Float" Dimensions="3" Format="XML" Precision="8">'//&
                      trim(str([dxyz(3),dxyz(2),dxyz(1)],separator=' '))//'</DataItem>'
   write(file_unit, '(A)') '          </Geometry>'
   write(file_unit, '(A)') '          <Topology Dimensions="'//trim(str([nijk(3),nijk(2),nijk(1)],separator=' '))//&
                           '" Type="3DCoRectMesh"/>'
   do v=1, size(q_name, dim=1)
      h5_dset_name = trim(q_name(v))//'-'//trim(str(myrank,.true.))//'-'//trim(str(block_index,.true.))
      write(file_unit, '(A)') '          <Attribute Name="'//trim(q_name(v))//&
                         '" Center="Cell" ElementDegree="0" Type="Scalar">'
      write(file_unit, '(A)') '            <DataItem DataType="Float" Dimensions="'//&
                              trim(str([nijk(3),nijk(2),nijk(1)],separator=' '))//   &
                              '" Format="HDF" Precision="8">'//h5_file_name//':'//h5_dset_name//'</DataItem>'
      write(file_unit, '(A)') '          </Attribute>'
   enddo
   if (allocated(q_aux_name)) then
      do v=1, size(q_aux_name, dim=1)
         h5_dset_name = trim(q_aux_name(v))//'-'//trim(str(myrank,.true.))//'-'//trim(str(block_index,.true.))
         write(file_unit, '(A)') '          <Attribute Name="'//trim(q_aux_name(v))//&
                            '" Center="Cell" ElementDegree="0" Type="Scalar">'
         write(file_unit, '(A)') '            <DataItem DataType="Float" Dimensions="'//&
                                 trim(str([nijk(3),nijk(2),nijk(1)],separator=' '))//   &
                                 '" Format="HDF" Precision="8">'//h5_file_name//':'//h5_dset_name//'</DataItem>'
         write(file_unit, '(A)') '          </Attribute>'
      enddo
   endif
   if (with_cell_morton) then
      h5_dset_name = 'morton-'//trim(str(myrank,.true.))//'-'//trim(str(block_index,.true.))
      write(file_unit, '(A)') '          <Attribute Name="morton" Center="Cell" ElementDegree="0" Type="Scalar">'
      write(file_unit, '(A)') '            <DataItem DataType="Float" Dimensions="'//&
                              trim(str([nijk(3),nijk(2),nijk(1)],separator=' '))//   &
                              '" Format="HDF" Precision="8">'//h5_file_name//':'//h5_dset_name//'</DataItem>'
      write(file_unit, '(A)') '          </Attribute>'
   endif
   write(file_unit, '(A)') '          <Attribute Name="Morton" Center="Grid">'
   write(file_unit, '(A)') '            <DataItem Dimensions="1" Format="XML" DataType="Int">'//trim(str(code))//'</DataItem>'
   write(file_unit, '(A)') '          </Attribute>'
   write(file_unit, '(A)') '          <Attribute Name="block-index" Center="Grid">'
   write(file_unit, '(A)') '            <DataItem Dimensions="1" Format="XML" DataType="Int">'//trim(str(block_index))//'</DataItem>'
   write(file_unit, '(A)') '          </Attribute>'
   write(file_unit, '(A)') '          <Attribute Name="myrank" Center="Grid">'
   write(file_unit, '(A)') '            <DataItem Dimensions="1" Format="XML" DataType="Int">'//trim(str(myrank))//'</DataItem>'
   write(file_unit, '(A)') '          </Attribute>'
   if (present(t)) then
      write(file_unit, '(A)') '          <Attribute Name="t" Center="Grid">'
      write(file_unit, '(A)') '            <DataItem Dimensions="1" Format="XML" DataType="Int">'//trim(str(t))//'</DataItem>'
      write(file_unit, '(A)') '          </Attribute>'
   endif
   if (present(time)) then
      write(file_unit, '(A)') '          <Attribute Name="time" Center="Grid">'
      write(file_unit, '(A)') '            <DataItem Dimensions="1" Format="XML" DataType="Float">'//trim(str(time))//'</DataItem>'
      write(file_unit, '(A)') '          </Attribute>'
   endif
   write(file_unit, '(A)') '        </Grid>'
   endsubroutine save_xdmf_block
endmodule adam_adam_object
