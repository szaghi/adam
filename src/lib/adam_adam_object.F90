!< ADAM, ADAM class definition.
module adam_adam_object
!< ADAM, ADAM class definition.

use adam_field_object
use adam_grid_object
use adam_parameters
use adam_tree_node_object
use adam_tree_bucket_object
use adam_tree_object
use PENF
use stringifor
use vtk_fortran
use HDF5
#ifdef _MPI_
use MPI
#endif

implicit none
private
public :: adam_object

type :: adam_object
   !< ADAM class definition.
   type(grid_object)  :: grid        !< The grid.
   type(tree_object)  :: tree        !< The tree.
   type(field_object) :: field       !< The field.
   integer(I4P)       :: error=0_I4P !< Error traping flag.
   ! MPI data
   integer(I4P) :: procs_number=1_I4P !< MPI Number of processes.
   integer(I4P) :: myrank=0_I4P       !< MPI rank process.
   contains
      ! public methods
      procedure, pass(self) :: amr_update       !< Update AMR status.
      procedure, pass(self) :: destroy          !< Destroy ADAM.
      procedure, pass(self) :: finalize         !< Finalize ADAM.
      procedure, pass(self) :: initialize       !< Initialize ADAM.
      procedure, pass(self) :: mpi_redistribute !< Redistribute nodes/blocks to processes, load balancing.
      procedure, pass(self) :: save_hdf5        !< Save ADAM in HDF5 format.
      procedure, pass(self) :: save_vtk         !< Save ADAM in VTK  format.
      ! operators
      generic :: assignment(=) => adam_assign_adam      !< Overload `=`.
      procedure, pass(lhs), private :: adam_assign_adam !< Operator `=`.
endtype adam_object

contains
   ! public methods
   subroutine amr_update(self, is_marked_by_field, is_marked_by_tree, is_grid_changed)
   !< Update AMR status.
   class(adam_object), intent(inout)         :: self                      !< ADAM.
   logical,            intent(in),  optional :: is_marked_by_field        !< Flag to check if marker is field.
   logical,            intent(in),  optional :: is_marked_by_tree         !< Flag to check if marker is tree.
   logical,            intent(out), optional :: is_grid_changed           !< Flag to check if grid is changed.
   logical                                   :: is_marked_by_field_       !< Flag to check if marker is field, local var.
   logical                                   :: is_marked_by_tree_        !< Flag to check if marker is tree, local var.

   is_marked_by_field_ = .false. ; if (present(is_marked_by_field)) is_marked_by_field_ = is_marked_by_field
   is_marked_by_tree_  = .false. ; if (present(is_marked_by_tree )) is_marked_by_tree_  = is_marked_by_tree

   if (is_marked_by_field_) then
      call self%field%mpi_gather_refinements_needed
      call self%tree%import_refinements_needed(refinements_needed_all=self%field%refinements_needed_all, &
                                               disp_count=self%field%disp_count)
   endif
   if (is_marked_by_tree_) then
      call self%tree%mpi_gather_refinements_needed
   endif
   call self%tree%adapt
   call self%field%adapt(ratio=self%tree%ratio,                                                            &
                         block_to_refine=self%tree%block_to_refine, block_refined=self%tree%block_refined, &
                         block_to_derefine=self%tree%block_to_derefine, block_derefined=self%tree%block_derefined)
   if (present(is_grid_changed)) &
      is_grid_changed = (size(self%tree%node_to_refine, dim=1)>0_I4P).or.(size(self%tree%node_to_derefine, dim=1)>0_I4P)
   endsubroutine amr_update

   subroutine destroy(self)
   !< Destroy ADAM.
   class(adam_object), intent(inout) :: self  !< ADAM.
   type(adam_object)                 :: fresh !< Fresh ADAM.

   self = fresh
   endsubroutine destroy

   subroutine finalize(self)
   !< Finalize ADAM.
   class(adam_object), intent(inout) :: self !< ADAM.

#ifdef _MPI_
   call MPI_FINALIZE(self%error)
#endif
   stop
   endsubroutine finalize

   subroutine initialize(self,                                                               &
                         ni, nj, nk, gc, emin, emax,                                         &
                         max_load, nodes_number, buckets_number, ratio, max_level, add_adam, &
                         nv, nb)
   !< Initialize ADAM.
   class(adam_object), intent(inout)        :: self           !< ADAM.
   ! grid options
   integer(I4P),       intent(in), optional :: ni             !< Number of cells in X direction.
   integer(I4P),       intent(in), optional :: nj             !< Number of cells in Y direction.
   integer(I4P),       intent(in), optional :: nk             !< Number of cells in Z direction.
   integer(I4P),       intent(in), optional :: gc(6)          !< Number of ghost cells in each direction.
   real(R8P),          intent(in), optional :: emin(3)        !< Coordinates of minium abscissa.
   real(R8P),          intent(in), optional :: emax(3)        !< Coordinates of maxium abscissa.
   ! tree options
   real(R8P),          intent(in), optional :: max_load       !< Maximum load of tree buckets.
   integer(I8P),       intent(in), optional :: nodes_number   !< Nodes number to be stored in the tree.
   integer(I8P),       intent(in), optional :: buckets_number !< Number of buckets for initialize the tree.
   integer(I4P),       intent(in), optional :: ratio          !< Refinement ratio.
   integer(I4P),       intent(in), optional :: max_level      !< Maximum refinement level.
   logical,            intent(in), optional :: add_adam       !< Add ADAM node, the ancestor of all nodes.
   ! field options
   integer(I4P),       intent(in), optional :: nv             !< Number of field variables.
   integer(I4P),       intent(in), optional :: nb             !< Number of all blocks that can be stored in field.

#ifdef _MPI_
   call MPI_INIT(self%error)
#endif
   call self%destroy
#ifdef _MPI_
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
#endif
   call self%grid%initialize(ni=ni, nj=nj, nk=nk, gc=gc, emin=emin, emax=emax)
   call self%tree%initialize(grid=self%grid, max_load=max_load, nodes_number=nodes_number, buckets_number=buckets_number, &
                             ratio=ratio, max_level=max_level, add_adam=add_adam)
   call self%field%initialize(grid=self%grid, nv=nv, nb=nb)
   endsubroutine initialize

   subroutine mark_sphere_nodes(self, center, radius, threshold)
   !< Mark all nodes inside a sphere to be refined.
   class(adam_object), intent(inout)        :: self            !< ADAM.
   real(R8P),          intent(in)           :: center(3)       !< Sphere center coordinates [x,y,z].
   real(R8P),          intent(in)           :: radius          !< Sphere radius.
   real(R8P),          intent(in), optional :: threshold       !< Threshold for sphere proximity.
   real(R8P)                                :: threshold_      !< Threshold for sphere proximity, local var.
   type(tree_node_object), pointer          :: node            !< Pointer to current node.
   real(R8P)                                :: block_center(3) !< block center coordinates.
   real(R8P)                                :: block_diagonal  !< block diagonal.
   real(R8P)                                :: distance(0:8)   !< Distances between block and sphere.
   real(R8P)                                :: max_cell_delta  !< Max cell delta.

   threshold_ = 2.2_R8P ; if (present(threshold)) threshold_ = threshold
   do while(self%tree%loop(node=node))
      if (self%myrank /= node%myrank) cycle ! mark only my nodes
      block_center = (self%field%emax(:,node%block_index) + self%field%emin(:,node%block_index)) / 2._R8P
      block_diagonal = sqrt((self%field%emax(1,node%block_index) - self%field%emin(1,node%block_index))**2 + &
                            (self%field%emax(2,node%block_index) - self%field%emin(2,node%block_index))**2 + &
                            (self%field%emax(3,node%block_index) - self%field%emin(3,node%block_index))**2)

      associate (emin=>self%field%emin(:,node%block_index), emax=>self%field%emax(:,node%block_index), &
                 ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk)
      distance(0) = sphere_distance(point=block_center)
      distance(1) = sphere_distance(point=[emin(1), emin(2), emin(3)])
      distance(2) = sphere_distance(point=[emax(1), emin(2), emin(3)])
      distance(3) = sphere_distance(point=[emin(1), emax(2), emin(3)])
      distance(4) = sphere_distance(point=[emax(1), emax(2), emin(3)])
      distance(5) = sphere_distance(point=[emin(1), emin(2), emax(3)])
      distance(6) = sphere_distance(point=[emax(1), emin(2), emax(3)])
      distance(7) = sphere_distance(point=[emin(1), emax(2), emax(3)])
      distance(8) = sphere_distance(point=[emax(1), emax(2), emax(3)])
      if (maxval(distance(0:8),dim=1)*minval(distance(0:8),dim=1) < 0._R8P) then
         distance(0) = 0._R8P
      endif

      max_cell_delta = self%field%max_cell_delta(distance=distance(0))

      if (block_diagonal/min(ni,nj,nk) > max_cell_delta) then
         node%refinement_needed = TO_BE_REFINED
      elseif (block_diagonal/min(ni,nj,nk) * threshold_ < max_cell_delta) then
         node%refinement_needed = TO_BE_DEREFINED
      endif
      endassociate
   enddo
   contains
      pure function sphere_distance(point)
      !< Return the distance from a point to the sphere surface, with sign.
      real(R8P), intent(in) :: point(3)        !< Point coordinates.
      real(R8P)             :: sphere_distance !< Distance from sphere surface.

      sphere_distance = sqrt((center(1) - point(1))**2 + &
                             (center(2) - point(2))**2 + &
                             (center(3) - point(3))**2) - radius
      endfunction sphere_distance
   endsubroutine mark_sphere_nodes

   subroutine mpi_redistribute(self)
   !< Redistribute nodes/blocks to processes, load balancing.
   class(adam_object), intent(inout) :: self !< ADAM.

   call self%tree%mpi_redistribute
   call self%field%mpi_redistribute(comm_map_send=self%tree%comm_map_send,         &
                                    comm_map_recv=self%tree%comm_map_recv,         &
                                    comm_map_send_ptr=self%tree%comm_map_send_ptr, &
                                    comm_map_recv_ptr=self%tree%comm_map_recv_ptr, &
                                    local_map=self%tree%local_map,                 &
                                    coordinates=self%tree%block_coordinates)
   endsubroutine mpi_redistribute

   subroutine save_hdf5(self, basename, directory)
   !< Save ADAM in HDF5 format.
   class(adam_object), intent(in)           :: self             !< ADAM.
   character(*),       intent(in)           :: basename         !< Base name of output files.
   character(*),       intent(in), optional :: directory        !< Directory name of output files.
   character(:), allocatable                :: directory_       !< Directory name of output files, local var.
   type(tree_node_object), pointer          :: node             !< Pointer to node.
   real(R8P)                                :: emin(3), emax(3) !< Minimum/maximum abscissa of current block.
   integer(I4P)                             :: error            !< Error trapping flag.
   integer(I4P)                             :: b                !< Counter.
   integer(I4P)                             :: xdmf             !< XDMF file handler.
   character(len=:), allocatable            :: h5_file_name     !< H5 Dataset name.
   character(len=:), allocatable            :: h5_dset_name     !< H5 Dataset name.
   integer(HID_T)                           :: h5_file_id       !< H5 File identifier.
   integer(HID_T)                           :: h5_dspace_id     !< H5 Dataspace identifier.
   integer(HID_T)                           :: h5_dset_id       !< H5 Dataset identifier.
   real(R8P)                                :: dxl, dyl, dzl    !< Local delta space.
   integer(I4P)                             :: i, j, k, l       !< Counter.

   directory_ = '' ; if (present(directory)) directory_ = trim(directory)
   associate(grid=>self%grid, ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk)
   ! save H5 file (one for each process)
   ! open fortran interface
   call h5open_f(error)
   ! create a new file using default properties
   h5_file_name = directory_//trim(basename)//trim(str(self%myrank,.true.))//'.h5'
   call h5fcreate_f(h5_file_name, H5F_ACC_TRUNC_F, h5_file_id, error)
   ! create the dataspace
   call h5screate_simple_f(3_I4P, [int(ni,I8P),int(nj,I8P),int(nk,I8P)], h5_dspace_id, error)
   ! save all blocks in process
   do b=1, self%field%blocks_number
      h5_dset_name = 'u-'//trim(str(self%myrank,.true.))//'-'//trim(str(b,.true.))
      call h5dcreate_f(h5_file_id, h5_dset_name, H5T_NATIVE_DOUBLE, h5_dspace_id, h5_dset_id, error)
      call h5dwrite_f(h5_dset_id, H5T_NATIVE_DOUBLE, self%field%u(1:ni,1:nj,1:nk,b), [int(ni,I8P),int(nj,I8P),int(nk,I8P)], error)
      call h5dclose_f(h5_dset_id, error)
   enddo
   ! terminate access to the data space
   call h5sclose_f(h5_dspace_id, error)
   ! close the file
   call h5fclose_f(h5_file_id, error)
   ! close FORTRAN interface
   call h5close_f(error)

   ! save XDMF file (only master process does)
   if (self%myrank == 0_I4P) then
      open(newunit=xdmf, file=directory_//trim(basename)//'.xdmf')
      write(xdmf, '(A)') '<?xml version="1.0" encoding="utf-8"?>'
      write(xdmf, '(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="3.0">'
      write(xdmf, '(A)') '  <Domain>'
      write(xdmf, '(A)') '    <Grid Name="ADAM" GridType="Collection">'
      do while(self%tree%loop(node=node))
         b = node%block_index
         h5_file_name = trim(basename)//trim(str(node%myrank,.true.))//'.h5'
         h5_dset_name = 'u-'//trim(str(node%myrank,.true.))//'-'//trim(str(b,.true.))
         call self%tree%morton_to_coordinates(code=node%code, i=i, j=j, k=k, l=l)
         call grid%compute_emin_emax(coordinates=[i,j,k,l], emin=emin, emax=emax)
         dxl = grid%dxyz(emin=emin, emax=emax, axis='x')
         dyl = grid%dxyz(emin=emin, emax=emax, axis='y')
         dzl = grid%dxyz(emin=emin, emax=emax, axis='z')

         write(xdmf, '(A)') '      <Grid Name="'//trim(str(node%code))//'">'
         write(xdmf, '(A)') '        <Geometry Origin="" Type="ORIGIN_DXDYDZ">'
         write(xdmf, '(A)') '          <DataItem DataType="Float" Dimensions="3" Format="XML" Precision="8">'// &
                                        trim(str([emin(3),emin(2),emin(1)],separator=' '))//'</DataItem>'
         write(xdmf, '(A)') '          <DataItem DataType="Float" Dimensions="3" Format="XML" Precision="8">'// &
                                        trim(str([dxl,dyl,dzl],separator=' '))//'</DataItem>'
         write(xdmf, '(A)') '        </Geometry>'
         write(xdmf, '(A)') '        <Topology Dimensions="'//trim(str([ni+1,nj+1,nk+1],separator=' '))// &
                                      '" Type="3DCoRectMesh"/>'
         write(xdmf, '(A)') '        <Attribute Center="Cell" ElementCell="" ElementDegree="0" ElementFamily=""'// &
                                      ' ItemType="" Name="u" Type="Scalar"> '
         write(xdmf, '(A)') '          <DataItem DataType="Float" Dimensions="'//trim(str([ni+1,nj+1,nk+1],separator=' '))// &
                                        '" Format="HDF" Precision="8">'//h5_file_name//':'//h5_dset_name//'</DataItem>'
         write(xdmf, '(A)') '        </Attribute>'
         write(xdmf, '(A)') '      </Grid>'
      enddo
      write(xdmf, '(A)') '    </Grid>'
      write(xdmf, '(A)') '  </Domain>'
      write(xdmf, '(A)') '</Xdmf>'
      close(xdmf)
   endif
   endassociate
   endsubroutine save_hdf5

   subroutine save_vtk(self, basename, directory)
   !< Save ADAM in VTK files.
   class(adam_object), intent(in)           :: self       !< ADAM.
   character(*),       intent(in)           :: basename   !< Base name of output files.
   character(*),       intent(in), optional :: directory  !< Directory name of output files.
   character(:), allocatable                :: directory_ !< Directory name of output files, local var.
   integer(I4P)                             :: error      !< Error trapping flag.
   type(vtk_file)                           :: vtk        !< VTK file handler.
   type(vtm_file)                           :: vtm        !< VTM file handler.
   type(tree_node_object), pointer          :: node       !< Pointer to node.
   integer(I4P)                             :: b, l       !< Counter.
   integer(I4P)                             :: i, j, k    !< Counter.
   integer(I4P)                             :: max_level  !< Maximum level.

   directory_ = '' ; if (present(directory)) directory_ = trim(directory)
   associate(emin=>self%field%emin, emax=>self%field%emax, ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk)
      max_level = 0_I4P
      vtr_loop : do while(self%tree%loop(node=node))
         if (self%myrank /= node%myrank) cycle ! only the process having node can save VTR file
         b = node%block_index
         max_level = max(max_level, self%tree%level(code=node%code))
         error = vtk%initialize(format='raw', filename=directory_//trim(basename)//'-block-'//trim(str(b,.true.))//&
                                                       '-proc-'//trim(str(node%myrank,.true.))//'.vtr',            &
                                mesh_topology='RectilinearGrid',                                                   &
                                nx1=0, nx2=ni, ny1=0, ny2=nj, nz1=0, nz2=nk)
         error = vtk%xml_writer%write_fielddata(action='open')
         error = vtk%xml_writer%write_fielddata(data_name='Morton', x=self%field%code(b))
         error = vtk%xml_writer%write_fielddata(action='close')
         error = vtk%xml_writer%write_piece(nx1=0, nx2=ni, ny1=0, ny2=nj, nz1=0, nz2=nk)
         error = vtk%xml_writer%write_geo(x=self%field%compute_xyz(b, axis='x'), &
                                          y=self%field%compute_xyz(b, axis='y'), &
                                          z=self%field%compute_xyz(b, axis='z'))
         error = vtk%xml_writer%write_dataarray(location='cell', action='open')
         error = vtk%xml_writer%write_dataarray(data_name='u', x=[self%field%u(1:ni,1:nj,1:nk,b)])
         error = vtk%xml_writer%write_dataarray(data_name='myrn', x=[(((node%myrank_new, k=1,nk),j=1,nj),i=1,ni)])
         error = vtk%xml_writer%write_dataarray(location='cell', action='close')
         error = vtk%xml_writer%write_piece()
         error = vtk%finalize()
      enddo vtr_loop

      ! save VTM file (only master process does)
      if (self%myrank == 0_I4P) then
         error = vtm%initialize(filename=directory_//trim(basename)//'.vtm', scratch_units_number=max_level)
         vtm_group_loop : do l=1, max_level
            error = vtm%write_block(scratch=l, action='open', name='level-'//trim(str(l,.true.)))
         enddo vtm_group_loop
         vtm_filenames_loop : do while(self%tree%loop(node=node))
            b = node%block_index
            l = self%tree%level(code=node%code)
            error = vtm%write_block(scratch=l, action='write', filename=trim(basename)//'-block-'//trim(str(b,.true.))//&
                                                                        '-proc-'//trim(str(node%myrank,.true.))//'.vtr')
         enddo vtm_filenames_loop
         error = vtm%finalize()
      endif
   endassociate
   endsubroutine save_vtk

   ! operators
   ! =
   subroutine adam_assign_adam(lhs, rhs)
   !< Operator `=`.
   class(adam_object), intent(inout) :: lhs !< Left hand side.
   type(adam_object),  intent(in)    :: rhs !< Right hand side.

   lhs%grid         = rhs%grid
   lhs%tree         = rhs%tree
   lhs%field        = rhs%field
   lhs%error        = rhs%error
   lhs%procs_number = rhs%procs_number
   lhs%myrank       = rhs%myrank
   endsubroutine adam_assign_adam
endmodule adam_adam_object
