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
   type(grid_object)  :: grid  !< The grid.
   type(tree_object)  :: tree  !< The tree.
   type(field_object) :: field !< The field.
   ! MPI data
   integer(I4P) :: error=0_I4P        !< Error traping flag.
   integer(I4P) :: procs_number=1_I4P !< MPI Number of processes.
   integer(I4P) :: myrank=0_I4P       !< MPI rank process.
   contains
      ! public methods
      procedure, pass(self) :: adapt                        !< Adapt tree/field accordingly to refine/derefine necessity.
      procedure, pass(self) :: amr_update                   !< Update AMR status.
      procedure, pass(self) :: blocks_reorder               !< Reorder blocks (for asyncrhonous MPI)
      procedure, pass(self) :: destroy                      !< Destroy ADAM.
      procedure, pass(self) :: finalize                     !< Finalize ADAM.
      procedure, pass(self) :: initialize                   !< Initialize ADAM.
      procedure, pass(self) :: make_comm_local_maps_ghost   !< Make communication/local maps of ghost cells.
      procedure, pass(self) :: mpi_gather_refinement_needed !< Gather refinement needed.
      procedure, pass(self) :: mpi_redistribute             !< Redistribute nodes/blocks to processes, load balancing.
      procedure, pass(self) :: save_hdf5                    !< Save ADAM in HDF5 format.
      procedure, pass(self) :: save_vtk                     !< Save ADAM in VTK  format.
      procedure, pass(self) :: update_ghost                 !< Update ghost cells.
      ! operators
      generic :: assignment(=) => adam_assign_adam      !< Overload `=`.
      procedure, pass(lhs), private :: adam_assign_adam !< Operator `=`.
endtype adam_object

contains
   ! public methods
   subroutine adapt(self)
   !< Adapt tree/field accordingly to refine/derefine necessity.
   class(adam_object), intent(inout) :: self !< ADAM.

   call self%tree%adapt
   call self%field%adapt(ratio=self%tree%ratio,                                                            &
                         block_to_refine=self%tree%block_to_refine, block_refined=self%tree%block_refined, &
                         block_to_derefine=self%tree%block_to_derefine, block_derefined=self%tree%block_derefined)
   endsubroutine adapt

   subroutine amr_update(self, is_marked_by_field, is_marked_by_tree, do_mpi_redistribute, do_blocks_reorder, &
                         print_mpi_stats, is_grid_changed)
   !< Update AMR status.
   class(adam_object), intent(inout)         :: self                 !< ADAM.
   logical,            intent(in),  optional :: is_marked_by_field   !< Flag to check if marker is field.
   logical,            intent(in),  optional :: is_marked_by_tree    !< Flag to check if marker is tree.
   logical,            intent(in),  optional :: do_mpi_redistribute  !< Flag to activate MPI redistribute.
   logical,            intent(in),  optional :: do_blocks_reorder    !< Flag to activate blocks reorder.
   logical,            intent(in),  optional :: print_mpi_stats      !< Flag to activate MPI statistics print.
   logical,            intent(out), optional :: is_grid_changed      !< Flag to check if grid is changed.
   logical                                   :: do_mpi_redistribute_ !< Flag to activate MPI redistribute, local var.
   logical                                   :: do_blocks_reorder_   !< Flag to activate blocks reorder, local var.

   do_mpi_redistribute_ = .true.  ; if (present(do_mpi_redistribute )) do_mpi_redistribute_ = do_mpi_redistribute
   do_blocks_reorder_ = .true.  ; if (present(do_blocks_reorder)) do_blocks_reorder_ = do_blocks_reorder

   call self%mpi_gather_refinement_needed(is_marked_by_field=is_marked_by_field, is_marked_by_tree=is_marked_by_tree)

   call self%update_ghost

   call self%adapt

   if (present(is_grid_changed)) is_grid_changed = (size(self%tree%node_to_refine,   dim=1)>0_I4P).or.&
                                                   (size(self%tree%node_to_derefine, dim=1)>0_I4P)

   if (do_mpi_redistribute_) call self%mpi_redistribute(print_mpi_stats=print_mpi_stats)

   if (do_blocks_reorder_) call self%blocks_reorder

   call self%make_comm_local_maps_ghost
   endsubroutine amr_update

   subroutine blocks_reorder(self)
   !< Reorder blocks (for asyncrhonous MPI)
   class(adam_object), intent(inout) :: self !< ADAM.

   call self%tree%blocks_reorder
   call self%field%blocks_reorder(inner_outer_block_map=self%tree%inner_outer_block_map, &
                                  inner_blocks_number=self%tree%inner_blocks_number)
   endsubroutine blocks_reorder

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
   call self%amr_update
   endsubroutine initialize

   subroutine make_comm_local_maps_ghost(self)
   !< Make communication/local maps of ghost cells.
   class(adam_object), intent(inout) :: self !< ADAM.

   call self%tree%make_comm_local_maps_ghost
   call self%field%prepare_comm_local_ghost(local_map_ghost         = self%tree%local_map_ghost,         &
                                            comm_map_n_send_ghost   = self%tree%comm_map_n_send_ghost,   &
                                            comm_map_n_recv_ghost   = self%tree%comm_map_n_recv_ghost,   &
                                            comm_map_send_ptr_ghost = self%tree%comm_map_send_ptr_ghost, &
                                            comm_map_recv_ptr_ghost = self%tree%comm_map_recv_ptr_ghost, &
                                            comm_map_send_ghost     = self%tree%comm_map_send_ghost,     &
                                            comm_map_recv_ghost     = self%tree%comm_map_recv_ghost)
   endsubroutine make_comm_local_maps_ghost

   subroutine mpi_gather_refinement_needed(self, is_marked_by_field, is_marked_by_tree)
   !< Gather refinement needed.
   class(adam_object), intent(inout)         :: self                !< ADAM.
   logical,            intent(in),  optional :: is_marked_by_field  !< Flag to check if marker is field.
   logical,            intent(in),  optional :: is_marked_by_tree   !< Flag to check if marker is tree.
   logical                                   :: is_marked_by_field_ !< Flag to check if marker is field, local var.
   logical                                   :: is_marked_by_tree_  !< Flag to check if marker is tree, local var.

   is_marked_by_field_  = .false. ; if (present(is_marked_by_field  )) is_marked_by_field_  = is_marked_by_field
   is_marked_by_tree_   = .false. ; if (present(is_marked_by_tree   )) is_marked_by_tree_   = is_marked_by_tree

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

   subroutine save_hdf5(self, basename, directory, with_ghost, with_cell_morton, with_ls)
   !< Save ADAM in HDF5 format.
   class(adam_object), intent(inout)        :: self              !< ADAM.
   character(*),       intent(in)           :: basename          !< Base name of output files.
   character(*),       intent(in), optional :: directory         !< Directory name of output files.
   logical,            intent(in), optional :: with_ghost        !< Flag to save ghost cells.
   logical,            intent(in), optional :: with_cell_morton  !< Flag to save Morton code also in cells.
   logical,            intent(in), optional :: with_ls           !< Flag to save level set field.
   character(:), allocatable                :: directory_        !< Directory name of output files, local var.
   logical                                  :: with_ghost_       !< Flag to save ghost cells, local var.
   logical                                  :: with_cell_morton_ !< Flag to save Morton code also in cells, local var.
   logical                                  :: with_ls_          !< Flag to save level set field, local var.
   type(tree_node_object), pointer          :: node              !< Pointer to node.
   real(R8P)                                :: emin(3)           !< Minimum abscissa of current block.
   integer(I4P)                             :: b                 !< Counter.
   integer(I4P)                             :: xdmf              !< XDMF file handler.
   character(len=:), allocatable            :: h5_file_name      !< H5 Dataset name.
   character(len=:), allocatable            :: h5_dset_name      !< H5 Dataset name.
   integer(HID_T)                           :: h5_file_id        !< H5 File identifier.
   integer(HID_T)                           :: h5_dspace_id      !< H5 Dataspace identifier.
   integer(HID_T)                           :: h5_dset_id        !< H5 Dataset identifier.
   real(R8P)                                :: dx, dy, dz        !< Space steps.
   character(len=:), allocatable            :: grid_dims         !< Grid dimensions.
   integer(I4P)                             :: gci, gcj, gck     !< Ghost cells saved.
   integer(I4P)                             :: i, j, k, l        !< Counter.

   directory_ = '' ; if (present(directory)) directory_ = trim(directory)
   with_ghost_ = .false. ; if (present(with_ghost)) with_ghost_ = with_ghost
   with_cell_morton_ = .false. ; if (present(with_cell_morton)) with_cell_morton_ = with_cell_morton
   with_ls_ = .false. ; if (present(with_ls)) with_ls_ = with_ls.and.allocated(self%field%ls)
   if (with_ghost_) then
      gci = self%grid%gci
      gcj = self%grid%gcj
      gck = self%grid%gck
   else
      gci = 0_I4P
      gcj = 0_I4P
      gck = 0_I4P
   endif
   associate(grid=>self%grid, ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk)
   ! save H5 file (one for each process)
   ! open fortran interface
   call h5open_f(self%error)
   ! create a new file using default properties
   h5_file_name = directory_//trim(basename)//'-proc'//trim(strz(self%myrank,6))//'.h5'
   call h5fcreate_f(h5_file_name, H5F_ACC_TRUNC_F, h5_file_id, self%error)

   ! create the dataspace for 3D fields
   call h5screate_simple_f(3_I4P, [int(ni+2*gci,I8P),int(nj+2*gcj,I8P),int(nk+2*gck,I8P)], h5_dspace_id, self%error)
   ! save all blocks in process
   do b=1, self%field%blocks_number
      h5_dset_name = 'u-'//trim(str(self%myrank,.true.))//'-'//trim(str(b,.true.))
      call h5dcreate_f(h5_file_id, h5_dset_name, H5T_NATIVE_DOUBLE, h5_dspace_id, h5_dset_id, self%error)
      call h5dwrite_f(h5_dset_id, H5T_NATIVE_DOUBLE, self%field%u(1-gci:ni+gci,1-gcj:nj+gcj,1-gck:nk+gck,b), &
                      [int(ni+2*gci,I8P),int(nj+2*gcj,I8P),int(nk+2*gck,I8P)], self%error)
      call h5dclose_f(h5_dset_id, self%error)

      if (with_cell_morton_) then
         h5_dset_name = 'morton-'//trim(str(self%myrank,.true.))//'-'//trim(str(b,.true.))
         call h5dcreate_f(h5_file_id, h5_dset_name, H5T_NATIVE_DOUBLE, h5_dspace_id, h5_dset_id, self%error)
         call h5dwrite_f(h5_dset_id, H5T_NATIVE_DOUBLE, &
                         reshape([(real(self%field%code(b),R8P),i=1,(ni+2*gci)*(nj+2*gcj)*(nk+2*gck))], &
                                 [ni+2*gci,nj+2*gcj,nk+2*gck]),                              &
                         [int(ni+2*gci,I8P),int(nj+2*gcj,I8P),int(nk+2*gck,I8P)], self%error)
         call h5dclose_f(h5_dset_id, self%error)
      endif

      if (with_ls_) then
         h5_dset_name = 'level-set-'//trim(str(self%myrank,.true.))//'-'//trim(str(b,.true.))
         call h5dcreate_f(h5_file_id, h5_dset_name, H5T_NATIVE_DOUBLE, h5_dspace_id, h5_dset_id, self%error)
         call h5dwrite_f(h5_dset_id, H5T_NATIVE_DOUBLE, self%field%ls(1-gci:ni+gci,1-gcj:nj+gcj,1-gck:nk+gck,b), &
                         [int(ni+2*gci,I8P),int(nj+2*gcj,I8P),int(nk+2*gck,I8P)], self%error)
         call h5dclose_f(h5_dset_id, self%error)
      endif
   enddo
   ! terminate access to the data space
   call h5sclose_f(h5_dspace_id, self%error)

   ! close the file
   call h5fclose_f(h5_file_id, self%error)
   ! close FORTRAN interface
   call h5close_f(self%error)

   ! save XDMF file (only master process does)
   if (self%myrank == 0_I4P) then
      grid_dims = trim(str([ni+2*gci+1,nj+2*gcj+1,nk+2*gck+1],separator=' '))
      open(newunit=xdmf, file=directory_//trim(basename)//'.xdmf')
      write(xdmf, '(A)') '<?xml version="1.0" encoding="utf-8"?>'
      write(xdmf, '(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="3.0">'
      write(xdmf, '(A)') '  <Domain>'
      write(xdmf, '(A)') '    <Grid Name="ADAM" GridType="Collection">'
      do while(self%tree%loop(node=node))
         b = node%block_index
         h5_file_name = trim(basename)//'-proc'//trim(strz(node%myrank,6))//'.h5'
         call self%tree%morton_to_coordinates(code=node%code, i=i, j=j, k=k, l=l)
         call grid%compute_metrics(coordinates=[i,j,k,l], emin=emin, dx=dx, dy=dy, dz=dz)
         write(xdmf, '(A)') '      <Grid Name="'//trim(str(node%code))//'">'

         write(xdmf, '(A)') '        <Geometry Origin="" Type="ORIGIN_DXDYDZ">'
         write(xdmf, '(A)') '          <DataItem DataType="Float" Dimensions="3" Format="XML" Precision="8">'// &
                                        trim(str([emin(3)-gck*dz,emin(2)-gcj*dy,emin(1)-gci*dx],separator=' '))//'</DataItem>'
         write(xdmf, '(A)') '          <DataItem DataType="Float" Dimensions="3" Format="XML" Precision="8">'// &
                                        trim(str([dz,dy,dx],separator=' '))//'</DataItem>'
         write(xdmf, '(A)') '        </Geometry>'
         write(xdmf, '(A)') '        <Topology Dimensions="'//grid_dims//'" Type="3DCoRectMesh"/>'

         h5_dset_name = 'u-'//trim(str(node%myrank,.true.))//'-'//trim(str(b,.true.))
         write(xdmf, '(A)') '        <Attribute Name="u" Center="Cell" ElementDegree="0" Type="Scalar">'
         write(xdmf, '(A)') '          <DataItem DataType="Float" Dimensions="'//grid_dims//'" Format="HDF" Precision="8">'// &
                                       h5_file_name//':'//h5_dset_name//'</DataItem>'
         write(xdmf, '(A)') '        </Attribute>'

         if (with_cell_morton_) then
            h5_dset_name = 'morton-'//trim(str(node%myrank,.true.))//'-'//trim(str(b,.true.))
            write(xdmf, '(A)') '        <Attribute Name="morton" Center="Cell" ElementDegree="0" Type="Scalar">'
            write(xdmf, '(A)') '          <DataItem DataType="Float" Dimensions="'//grid_dims//'" Format="HDF" Precision="8">'// &
                                          h5_file_name//':'//h5_dset_name//'</DataItem>'
            write(xdmf, '(A)') '        </Attribute>'
         endif

         if (with_ls_) then
            h5_dset_name = 'level-set-'//trim(str(node%myrank,.true.))//'-'//trim(str(b,.true.))
            write(xdmf, '(A)') '        <Attribute Name="level-set" Center="Cell" ElementDegree="0" Type="Scalar">'
            write(xdmf, '(A)') '          <DataItem DataType="Float" Dimensions="'//grid_dims//'" Format="HDF" Precision="8">'// &
                                          h5_file_name//':'//h5_dset_name//'</DataItem>'
            write(xdmf, '(A)') '        </Attribute>'
         endif

         write(xdmf, '(A)') '        <Attribute Name="Morton" Center="Grid">'
         write(xdmf, '(A)') '          <DataItem Dimensions="1" Format="XML" DataType="Int">'//trim(str(node%code))//'</DataItem>'
         write(xdmf, '(A)') '        </Attribute>'

         write(xdmf, '(A)') '        <Attribute Name="block-index" Center="Grid">'
         write(xdmf, '(A)') '          <DataItem Dimensions="1" Format="XML" DataType="Int">'//trim(str(node%block_index))//&
                                       '</DataItem>'
         write(xdmf, '(A)') '        </Attribute>'

         write(xdmf, '(A)') '        <Attribute Name="myrank" Center="Grid">'
         write(xdmf, '(A)') '          <DataItem Dimensions="1" Format="XML" DataType="Int">'//trim(str(node%myrank))//'</DataItem>'
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

   subroutine save_vtk(self, basename, directory, with_ghost)
   !< Save ADAM in VTK files.
   class(adam_object), intent(inout)        :: self                                          !< ADAM.
   character(*),       intent(in)           :: basename                                      !< Base name of output files.
   character(*),       intent(in), optional :: directory                                     !< Output directory name.
   logical,            intent(in), optional :: with_ghost                                    !< Flag to save ghost cells.
   character(:), allocatable                :: directory_                                    !< Output directory name, local var.
   logical                                  :: with_ghost_                                   !< Flag to save ghost cells, local var.
   type(vtk_file)                           :: vtk                                           !< VTK file handler.
   type(vtm_file)                           :: vtm                                           !< VTM file handler.
   type(tree_node_object), pointer          :: node                                          !< Pointer to node.
   integer(I4P)                             :: b, l                                          !< Counter.
   integer(I4P)                             :: i, j, k                                       !< Counter.
   integer(I4P)                             :: max_level                                     !< Maximum level.
   integer(I4P)                             :: gci, gcj, gck                                 !< Ghost cells saved.
   real(R8P)                                :: x(0-self%grid%gci:self%grid%ni+self%grid%gci) !< X coordinates.
   real(R8P)                                :: y(0-self%grid%gcj:self%grid%nj+self%grid%gcj) !< Y coordinates.
   real(R8P)                                :: z(0-self%grid%gck:self%grid%nk+self%grid%gck) !< Z coordinates.

   directory_ = '' ; if (present(directory)) directory_ = trim(directory)
   with_ghost_ = .false. ; if (present(with_ghost)) with_ghost_ = with_ghost
   if (with_ghost_) then
      gci = self%grid%gci
      gcj = self%grid%gcj
      gck = self%grid%gck
   else
      gci = 0_I4P
      gcj = 0_I4P
      gck = 0_I4P
   endif
   associate(ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk)
      max_level = 0_I4P
      vtr_loop : do b=1, self%field%blocks_number
         call self%grid%compute_metrics(coordinates=self%field%coordinates(:,b), x_node=x, y_node=y, z_node=z)
         max_level = max(max_level, self%field%coordinates(4,b))
         self%error = vtk%initialize(format='raw', filename=directory_//trim(basename)//'-block-'//trim(str(b,.true.))//&
                                                            '-proc-'//trim(str(self%myrank,.true.))//'.vtr',            &
                                     mesh_topology='RectilinearGrid',                                                   &
                                     nx1=0-gci, nx2=ni+gci, ny1=0-gcj, ny2=nj+gcj, nz1=0-gck, nz2=nk+gck)
         self%error = vtk%xml_writer%write_fielddata(action='open')
         self%error = vtk%xml_writer%write_fielddata(data_name='Morton', x=self%field%code(b))
         self%error = vtk%xml_writer%write_fielddata(data_name='myrank', x=self%myrank)
         self%error = vtk%xml_writer%write_fielddata(action='close')
         self%error = vtk%xml_writer%write_piece(nx1=0-gci, nx2=ni+gci, ny1=0-gcj, ny2=nj+gcj, nz1=0-gck, nz2=nk+gck)
         self%error = vtk%xml_writer%write_geo(x=x(0-gci:ni+gci), y=y(0-gcj:nj+gcj), z=z(0-gck:nk+gck))
         self%error = vtk%xml_writer%write_dataarray(location='cell', action='open')
         self%error = vtk%xml_writer%write_dataarray(data_name='u', x=[self%field%u(1-gci:ni+gci,1-gcj:nj+gcj,1-gck:nk+gck,b)])
         self%error = vtk%xml_writer%write_dataarray(location='cell', action='close')
         self%error = vtk%xml_writer%write_piece()
         self%error = vtk%finalize()
      enddo vtr_loop

      ! save VTM file (only master process does)
      if (self%myrank == 0_I4P) then
         self%error = vtm%initialize(filename=directory_//trim(basename)//'.vtm', scratch_units_number=max_level)
         vtm_group_loop : do l=1, max_level
            self%error = vtm%write_block(scratch=l, action='open', name='level-'//trim(str(l,.true.)))
         enddo vtm_group_loop
         vtm_filenames_loop : do while(self%tree%loop(node=node))
            b = node%block_index
            l = self%tree%level(code=node%code)
            self%error = vtm%write_block(scratch=l, action='write', filename=trim(basename)//'-block-'//trim(str(b,.true.))//&
                                                                             '-proc-'//trim(str(node%myrank,.true.))//'.vtr')
         enddo vtm_filenames_loop
         self%error = vtm%finalize()
      endif
   endassociate
   endsubroutine save_vtk

   subroutine update_ghost(self)
   !< Update ghost cells.
   class(adam_object), intent(inout) :: self !< ADAM.

   if (self%tree%nodes_number > 1) call self%field%update_ghost(q=self%field%u(:,:,:,:))
   endsubroutine update_ghost

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
