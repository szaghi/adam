!< ADAM, ADAM class definition.
module adam_adam_object
!< ADAM, ADAM class definition.

! ADAM modules
use :: adam_field_object
use :: adam_grid_object
use :: adam_maps_object
use :: adam_mpih_object
use :: adam_parameters
use :: adam_tree_node_object
use :: adam_tree_bucket_object
use :: adam_tree_object
! third party modules
use :: finer, only : file_ini
use :: motion
use :: penf
use :: stringifor
use :: vtk_fortran
use :: hdf5
use :: mpi

implicit none
private
public :: adam_object

type :: adam_object
   !< ADAM class definition.
   type(mpih_object)  :: mpih  !< The MPI handler.
   type(grid_object)  :: grid  !< The grid.
   type(tree_object)  :: tree  !< The tree.
   type(maps_object)  :: maps  !< The maps.
   type(field_object) :: field !< The field.
   contains
      ! public methods
      procedure, pass(self) :: adapt                         !< Adapt tree/field accordingly to refine/derefine necessity.
      procedure, pass(self) :: amr_update                    !< Update AMR status.
      procedure, pass(self) :: blocks_reorder                !< Reorder blocks (for asyncrhonous MPI)
      procedure, pass(self) :: check_blocks_number           !< Check if blocks number is groving too much.
      procedure, pass(self) :: compute_blocks_number         !< Compute maximum blocks number allocatable on memory available.
      procedure, pass(self) :: description                   !< Return pretty-printed object description.
      procedure, pass(self) :: initialize                    !< Initialize ADAM.
      procedure, pass(self) :: interpolate_at_point          !< Interpolate a scalar variable at a given point.
      procedure, pass(self) :: load_restart_files            !< Load restart files.
      procedure, pass(self) :: make_comm_local_maps_ghost_bc !< Make communication/local maps of ghost cells.
      procedure, pass(self) :: mpi_gather_refinement_needed  !< Gather refinement needed.
      procedure, pass(self) :: mpi_redistribute              !< Redistribute nodes/blocks to processes, load balancing.
      procedure, pass(self) :: prune                         !< Prune nodes/blocks.
      procedure, pass(self) :: refine_uniform                !< Refine all blocks uniformly.
      procedure, pass(self) :: save_restart_files            !< Save restart files.
      procedure, pass(self) :: save_slice                    !< Save slice.
      procedure, pass(self) :: save_vtk                      !< Save ADAM in VTK  format.
endtype adam_object

contains
   ! public methods
   subroutine adapt(self, q)
   !< Adapt tree/field accordingly to refine/derefine necessity.
   class(adam_object), intent(inout) :: self  !< ADAM.
   real(R8P),          intent(inout) :: q(1:,              &
                                          1-self%grid%ngc:,&
                                          1-self%grid%ngc:,&
                                          1-self%grid%ngc:,&
                                          1:) !< Field cell centered variables.

   call self%tree%adapt

   call self%check_blocks_number

   call self%field%adapt(ratio=self%tree%ratio,                                                            &
                         block_to_refine=self%tree%block_to_refine, block_refined=self%tree%block_refined, &
                         block_to_derefine=self%tree%block_to_derefine, block_derefined=self%tree%block_derefined, q=q)
   endsubroutine adapt

   subroutine amr_update(self, q, is_marked_by_field, is_marked_by_tree, do_mpi_redistribute, do_blocks_reorder, is_grid_changed)
   !< Update AMR status.
   !<
   !< Note: AMR update can be safely called only *after* update_ghost has been called for *q* variables, otherwise
   !< refine is not well done.
   !< Note: only if the AMR is UNIFORM and GLOBALLY made by tree, i.e. using mark_all_nodes, the mpi_redistribute can be avoided,
   !< otherwise mpi_gather_refinement_nedeed is not safe (having wrong nodes number counters).
   class(adam_object), intent(inout)         :: self                 !< ADAM.
   real(R8P),          intent(inout)         :: q(1:,              &
                                                  1-self%grid%ngc:,&
                                                  1-self%grid%ngc:,&
                                                  1-self%grid%ngc:,&
                                                  1:)                !< Field cell centered variables.
   logical,            intent(in),  optional :: is_marked_by_field   !< Flag to check if marker is field.
   logical,            intent(in),  optional :: is_marked_by_tree    !< Flag to check if marker is tree.
   logical,            intent(in),  optional :: do_mpi_redistribute  !< Flag to activate MPI redistribute.
   logical,            intent(in),  optional :: do_blocks_reorder    !< Flag to activate blocks reorder.
   logical,            intent(out), optional :: is_grid_changed      !< Flag to check if grid is changed.
   logical                                   :: do_mpi_redistribute_ !< Flag to activate MPI redistribute, local var.
   logical                                   :: do_blocks_reorder_   !< Flag to activate blocks reorder, local var.
   call self%mpih%print_message('adam_object%amr_update start')
   do_mpi_redistribute_ = .true. ; if (present(do_mpi_redistribute )) do_mpi_redistribute_ = do_mpi_redistribute
   do_blocks_reorder_ = .false. ; if (present(do_blocks_reorder)) do_blocks_reorder_ = do_blocks_reorder

   call self%mpi_gather_refinement_needed(is_marked_by_field=is_marked_by_field, is_marked_by_tree=is_marked_by_tree)

   call self%adapt(q=q)

   if (present(is_grid_changed)) is_grid_changed = (size(self%tree%node_to_refine,   dim=1)>0_I4P).or.&
                                                   (size(self%tree%node_to_derefine, dim=1)>0_I4P)

   if (do_mpi_redistribute_) call self%mpi_redistribute(q=q)

   if (do_blocks_reorder_) call self%blocks_reorder(q=q)

   call self%make_comm_local_maps_ghost_bc
   call self%mpih%print_message('adam_object%amr_update finish')
   endsubroutine amr_update

   subroutine blocks_reorder(self, q)
   !< Reorder blocks (for asyncrhonous MPI)
   class(adam_object), intent(inout) :: self !< ADAM.
   real(R8P),          intent(inout) :: q(1:,              &
                                          1-self%grid%ngc:,&
                                          1-self%grid%ngc:,&
                                          1-self%grid%ngc:,&
                                          1:)!< Field cell centered variables.

   call self%maps%blocks_reorder
   call self%field%blocks_reorder(q=q)
   endsubroutine blocks_reorder

   subroutine check_blocks_number(self)
   !< Check if blocks number is groving too much.
   class(adam_object), intent(inout) :: self             !< ADAM.
   type(tree_node_object), pointer   :: node_ptr         !< Pointer to current node.
   integer(I8P)                      :: max_nb           !< Maximum number of blocks desidered.
   character(len=1), parameter       :: NL=new_line('a') !< New line character.

   max_nb = 0
   do while(self%tree%loop(node_ptr=node_ptr))
      max_nb = max(max_nb, node_ptr%block_index)
   enddo
   if (max_nb > self%field%nb) then
      call self%mpih%abort(error_code=-101, msg='ERROR: the number of new blocks after AMR is greater than Nb'//NL//&
                                                'max blocks numer available [Nb]: '//trim(str(self%field%nb))//NL//&
                                                'blocks numer required after AMR: '//trim(str(max_nb)))
   endif
   call self%mpih%print_message('maximum number of blocks created after AMR update: '//str(max_nb)//'/'//str(self%field%nb))
   endsubroutine check_blocks_number

   subroutine compute_blocks_number(self, memory_avail, fields_number, nb, nodes_number)
   !< Compute maximum blocks number allocatable on memory available.
   class(adam_object), intent(in)           :: self          !< ADAM.
   real(R8P),          intent(in)           :: memory_avail  !< Memory available for single MPI process (GBytes).
   integer(I4P),       intent(in)           :: fields_number !< Fields number.
   integer(I4P),       intent(out)          :: nb            !< Maximum blocks number for single MPI process.
   integer(I8P),       intent(out)          :: nodes_number  !< Maximum blocks number for all MPI processes (nodes).
   integer(I4P)                             :: size_of_real  !< Size (bytes) of (one) real.
   real(R8P)                                :: save_factor   !< Factor to avoid memory completely full.

   size_of_real = storage_size(1._R8P)/8._R8P
   save_factor = 0.95_R8P
   nb = nint(save_factor * memory_avail*1e9 / (fields_number * self%grid%block_weight * size_of_real))
   nodes_number  = nb * self%mpih%procs_number
   endsubroutine compute_blocks_number

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(adam_object), intent(in) :: self             !< Adam.
   character(len=:), allocatable  :: desc             !< Description.
   character(len=1), parameter    :: NL=new_line('a') !< New line character.

   desc = self%grid%description()//NL//self%tree%description()//NL//self%field%description()
   endfunction description

   subroutine load_restart_files(self, basename, t, time, q)
   !< Load restart files.
   class(adam_object), intent(inout) :: self      !< ADAM.
   character(*),       intent(in)    :: basename  !< Base name of output files.
   integer(I4P),       intent(out)   :: t         !< Time iteration.
   real(R8P),          intent(out)   :: time      !< Time.
   real(R8P),          intent(inout) :: q(1:,              &
                                          1-self%grid%ngc:,&
                                          1-self%grid%ngc:,&
                                          1-self%grid%ngc:,&
                                          1:)              !< Field cell centered variables.
   integer(I4P)                      :: file_unit !< Output file unit.

   open(newunit=file_unit, file=trim(adjustl(basename))//'.time', form='UNFORMATTED', access='STREAM')
   read(unit=file_unit) t, time
   close(file_unit)
   call self%tree%load_nodes(file_name=trim(adjustl(basename))//'.tnd')
   call self%field%load_blocks(basename=basename, q=q)
   endsubroutine load_restart_files

   subroutine initialize(self, file_parameters, add_adam, verbose)
   !< Initialize ADAM.
   class(adam_object),     intent(inout)        :: self            !< ADAM.
   type(file_ini),         intent(inout)        :: file_parameters !< INI file handler.
   logical,                intent(in), optional :: add_adam        !< Add ADAM node, the ancestor of all nodes.
   logical,                intent(in), optional :: verbose         !< Trigger verbose output.
   logical                                      :: verbose_        !< Trigger verbose output, local variable.
   integer(I8P)                                 :: nodes_number    !< Nodes number to be stored in the tree.
   integer(I4P)                                 :: nb              !< Number of all blocks that can be stored in field.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   call self%mpih%initialize(verbose=verbose_)
   if (verbose_) call self%mpih%print_message('adam_object%initialize start')
  !call self%compute_blocks_number(memory_avail=memory_avail,&
  !                                fields_number=80,         & ! remember to change
  !                                nb=nb,                    &
  !                                nodes_number=nodes_number)
   nodes_number = 11
   self%field%blocks_number = 1
   call self%grid%initialize(file_parameters=file_parameters,verbose=verbose_) ! remember to call self%adam%grid%set_bc_type
   call self%tree%initialize(grid=self%grid,                 &
                             file_parameters=file_parameters,&
                             nodes_number=nodes_number,      &
                             add_adam=add_adam,              &
                             verbose=verbose_)
   call self%maps%initialize(grid=self%grid,tree=self%tree,verbose=verbose_)
   call self%field%initialize(grid=self%grid,                 &
                              maps=self%maps,                 &
                              file_parameters=file_parameters,&
                              nb=1,                           &! remember to change
                              verbose=verbose_)
   if (verbose_) call self%mpih%print_message('adam_object%initialize finish')
   if (verbose_) call self%mpih%print_message('blocks number (maximum) for single MPI [nb]: '//trim(str(self%field%nb)))
   if (verbose_) call self%mpih%print_message('blocks number for all MPI [nodes_number]: '//trim(str(self%tree%nodes_number)))
   endsubroutine initialize

   subroutine interpolate_at_point(self, itype, point, q, qp, is_mine, p, qc, ijk, xyz, code, v)
   !< Interpolate a scalar variable at a given point.
   class(adam_object), intent(inout)         :: self      !< ADAM.
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
   if (node%myrank == self%mpih%myrank) then
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

   call self%mpih%print_message('adam_object%make_comm_local_maps_ghost_bc start')
   call self%maps%make_comm_local_maps_ghost(nv=self%field%nv)
   call self%maps%make_local_maps_bc
   call self%mpih%print_message('adam_object%make_comm_local_maps_ghost_bc finish')
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
      call self%maps%mpi_gather_nodes_data(node_member='refinement_needed')
   endif
   endsubroutine mpi_gather_refinement_needed

   subroutine mpi_redistribute(self, q)
   !< Redistribute nodes/blocks to processes, load balancing.
   class(adam_object), intent(inout) :: self !< ADAM.
   real(R8P),          intent(inout) :: q(1:,              &
                                          1-self%grid%ngc:,&
                                          1-self%grid%ngc:,&
                                          1-self%grid%ngc:,&
                                          1:)!< Field cell centered variables.

   call self%tree%mpi_redistribute
   call self%maps%make_comm_local_maps
   call self%field%mpi_redistribute(q=q)

   endsubroutine mpi_redistribute

   subroutine prune(self, q, ijkl_prune, do_blocks_reorder)
   !< Prune nodes/blocks.
   class(adam_object), intent(inout)        :: self               !< Adam.
   real(R8P),          intent(inout)        :: q(1:,              &
                                                 1-self%grid%ngc:,&
                                                 1-self%grid%ngc:,&
                                                 1-self%grid%ngc:,&
                                                 1:)              !< Field cell centered variables.
   integer(I4P),       intent(inout)        :: ijkl_prune(4)      !< Maximum coordinates after which the prune operates.
   logical,            intent(in), optional :: do_blocks_reorder  !< Flag to activate blocks reorder.
   logical                                  :: do_blocks_reorder_ !< Flag to activate blocks reorder, local var.

   do_blocks_reorder_ = .true.  ; if (present(do_blocks_reorder)) do_blocks_reorder_ = do_blocks_reorder
   call self%tree%prune(ijkl_prune=ijkl_prune)
   call self%mpi_redistribute(q=q)
   if (do_blocks_reorder_) call self%blocks_reorder(q=q)
   call self%make_comm_local_maps_ghost_bc
   endsubroutine prune

   subroutine refine_uniform(self, refinement_levels, q, do_mpi_redistribute, do_blocks_reorder)
   !< Refine all blocks uniformly.
   class(adam_object), intent(inout)        :: self                 !< Adam.
   integer(I4P),       intent(in)           :: refinement_levels    !< Number of refinement to be performed.
   real(R8P),          intent(inout)        :: q(1:,              &
                                                 1-self%grid%ngc:,&
                                                 1-self%grid%ngc:,&
                                                 1-self%grid%ngc:,&
                                                 1:)                !< Field cell centered variables.
   logical,            intent(in), optional :: do_mpi_redistribute  !< Flag to activate MPI redistribute.
   logical,            intent(in), optional :: do_blocks_reorder    !< Flag to activate blocks reorder.
   integer(I4P)                             :: l                    !< Counter.

   call self%mpih%print_message('uniformly refine mesh with '//trim(str(refinement_levels))//' levels')
   do l=1, refinement_levels
      call self%tree%mark_all_nodes(mark=TO_BE_REFINED)
      call self%amr_update(q=q, do_mpi_redistribute=do_mpi_redistribute, do_blocks_reorder=do_blocks_reorder)
   enddo
   endsubroutine refine_uniform

   subroutine save_restart_files(self, basename, t, time, q)
   !< Save restart files.
   class(adam_object), intent(in) :: self      !< ADAM.
   character(*),       intent(in) :: basename  !< Base name of output files.
   integer(I4P),       intent(in) :: t         !< Time iteration.
   real(R8P),          intent(in) :: time      !< Time.
   real(R8P),          intent(in) :: q(1:,              &
                                       1-self%grid%ngc:,&
                                       1-self%grid%ngc:,&
                                       1-self%grid%ngc:,&
                                       1:)     !< Field cell centered variables.
   integer(I4P)                   :: file_unit !< Output file unit.

   if (self%mpih%myrank==0) then
      open(newunit=file_unit, file=trim(adjustl(basename))//'.time', form='UNFORMATTED', access='STREAM')
      write(unit=file_unit) t, time
      close(file_unit)
      call self%tree%save_nodes(file_name=trim(adjustl(basename))//'.tnd')
   endif
   call self%field%save_blocks(basename=basename, q=q)
   endsubroutine save_restart_files

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
                      MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, MPI_IO_FILE_unit, self%mpih%error)
   offset_head = 0
   if (self%mpih%myrank==0) then
      ! call MPI_FILE_WRITE_AT_ALL(MPI_IO_FILE_unit, offset_head, nijkv, 4, MPI_INT, MPI_STATUS_IGNORE, error)
      call MPI_FILE_WRITE_AT(MPI_IO_FILE_unit, offset_head, nijkv, 4, MPI_INT, MPI_STATUS_IGNORE, self%mpih%error)
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
               call MPI_FILE_WRITE_AT(MPI_IO_FILE_unit, offset, points(:,i,j,k), 3, MPI_REAL8, MPI_STATUS_IGNORE, self%mpih%error)
               offset = offset + 8 * 3
               ! call MPI_FILE_WRITE_AT_ALL(MPI_IO_FILE_unit, offset, qp, size(q, dim=1), MPI_REAL8, MPI_STATUS_IGNORE, error)
               call MPI_FILE_WRITE_AT(MPI_IO_FILE_unit, offset, qp, size(q, dim=1), MPI_REAL8, MPI_STATUS_IGNORE, self%mpih%error)
               if (present(phi)) then
                  node => self%tree%node(code=code)
                  offset = offset + 8 * size(q, dim=1)
                  ! call MPI_FILE_WRITE_AT_ALL(MPI_IO_FILE_unit, offset, phi(node%block_index,ijkc(1,1),ijkc(2,1),ijkc(3,1)), 1, &
                  !                            MPI_REAL8, MPI_STATUS_IGNORE, error)
                  call MPI_FILE_WRITE_AT(MPI_IO_FILE_unit, offset, phi(node%block_index,ijkc(1,1),ijkc(2,1),ijkc(3,1)), 1, &
                                         MPI_REAL8, MPI_STATUS_IGNORE, self%mpih%error)
               endif
            endif
            ijk = ijk + 1
         enddo
      enddo
   enddo
   call MPI_FILE_CLOSE(MPI_IO_FILE_unit, self%mpih%error)
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
   integer(I4P)                             :: i                                             !< Counter.
   integer(I4P)                             :: max_level                                     !< Maximum level.
   integer(I4P)                             :: ngc                                           !< Ghost cells saved.
   integer(I4P)                             :: error                                         !< Error traping flag.
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
   associate(ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk)
      max_level = 0_I4P
      vtr_loop : do b=1, self%field%blocks_number
         call self%grid%compute_metrics(coordinates=self%field%coordinates(:,b), x_node=x, y_node=y, z_node=z)
         max_level = max(max_level, self%field%coordinates(4,b))
         error = vtk%initialize(format='raw', filename=directory_//trim(basename)//                          &
                                                       '-morton-'//trim(str(self%field%code(b),.true.))//    &
                                                       '-block-'//trim(str(b,.true.))//                      &
                                                       '-proc-'//trim(str(self%mpih%myrank,.true.))//'.vtr', &
                                mesh_topology='RectilinearGrid',                                             &
                                nx1=0-ngc, nx2=ni+ngc, ny1=0-ngc, ny2=nj+ngc, nz1=0-ngc, nz2=nk+ngc)
                            error = vtk%xml_writer%write_fielddata(action='open')
                            error = vtk%xml_writer%write_fielddata(data_name='Morton', x=self%field%code(b))
                            error = vtk%xml_writer%write_fielddata(data_name='myrank', x=self%mpih%myrank)
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
      if (self%mpih%myrank == 0_I4P) then
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
                                                                        '-proc-'//trim(str(self%mpih%myrank,.true.))//'.vtr')
         enddo vtm_filenames_loop
         error = vtm%finalize()
      endif
   endassociate
   endsubroutine save_vtk
endmodule adam_adam_object
