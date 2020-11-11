!< ADAM, objects and procedurs.
module adam_objects
!< ADAM, objects and procedurs.

use adam_field_object
use adam_tree_node_object
use adam_tree_bucket_object
use adam_tree_object
use PENF
use stringifor
use vtk_fortran
use HDF5
use MPI

implicit none
private
public :: destroy_tree_node
public :: field_object
public :: tree_node_object, NODE_TO_BE_REFINED, NODE_TO_BE_DEREFINED, NODE_TO_NOT_TOUCH
public :: tree_bucket_object, len
public :: iterator_interface
public :: tree_object
public :: mark_sphere_nodes
public :: field_save_vtk
public :: save_hdf5

contains
   subroutine mark_sphere_nodes(tree, field, center, radius, threshold)
   !< Mark all nodes inside a sphere to be refined.
   type(tree_object),  intent(inout)        :: tree            !< The tree.
   type(field_object), intent(in)           :: field           !< The field.
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
   do while(tree%loop(node=node))
      if (tree%myrank /= node%myrank) cycle ! mark only my nodes
      block_center = (field%emax(:,node%block_index) + field%emin(:,node%block_index)) / 2._R8P
      block_diagonal = sqrt((field%emax(1,node%block_index) - field%emin(1,node%block_index))**2 + &
                            (field%emax(2,node%block_index) - field%emin(2,node%block_index))**2 + &
                            (field%emax(3,node%block_index) - field%emin(3,node%block_index))**2)

      associate (emin=>field%emin(:,node%block_index), emax=>field%emax(:,node%block_index), &
                 ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk)
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

      max_cell_delta = field%max_cell_delta(distance=distance(0))

      if (block_diagonal/min(ni,nj,nk) > max_cell_delta) then
         node%refinement_needed = NODE_TO_BE_REFINED
      elseif (block_diagonal/min(ni,nj,nk) * threshold_ < max_cell_delta) then
         node%refinement_needed = NODE_TO_BE_DEREFINED
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

   subroutine field_save_vtk(tree, field, basename, directory)
   !< Save field in VTK files.
   type(tree_object),  intent(in)           :: tree       !< The tree.
   type(field_object), intent(in)           :: field      !< The field.
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
   associate(emin=>field%emin, emax=>field%emax, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, &
             gc1=>field%grid%gc1, gc2=>field%grid%gc2, gc3=>field%grid%gc3, &
             gc4=>field%grid%gc4, gc5=>field%grid%gc5, gc6=>field%grid%gc6)

      max_level = 0_I4P
      vtr_loop : do while(tree%loop(node=node))
         if (field%myrank /= node%myrank) cycle ! only the process having node can save VTR file
         b = node%block_index
         max_level = max(max_level, tree%level(code=node%code))
         error = vtk%initialize(format='raw', filename=directory_//trim(basename)//'-block-'//trim(str(b,.true.))//&
                                                       '-proc-'//trim(str(node%myrank,.true.))//'.vtr',            &
                                mesh_topology='RectilinearGrid',                                                   &
                                nx1=0, nx2=ni, ny1=0, ny2=nj, nz1=0, nz2=nk)
         error = vtk%xml_writer%write_fielddata(action='open')
         error = vtk%xml_writer%write_fielddata(data_name='Morton', x=field%code(b))
         error = vtk%xml_writer%write_fielddata(action='close')
         error = vtk%xml_writer%write_piece(nx1=0, nx2=ni, ny1=0, ny2=nj, nz1=0, nz2=nk)
         error = vtk%xml_writer%write_geo(x=field%compute_xyz(b, axis='x'), &
                                          y=field%compute_xyz(b, axis='y'), &
                                          z=field%compute_xyz(b, axis='z'))
         error = vtk%xml_writer%write_dataarray(location='cell', action='open')
         error = vtk%xml_writer%write_dataarray(data_name='u', x=[field%u(1:ni,1:nj,1:nk,b)])
         error = vtk%xml_writer%write_dataarray(data_name='myrn', x=[(((node%myrank_new, k=1,nk),j=1,nj),i=1,ni)])
         error = vtk%xml_writer%write_dataarray(location='cell', action='close')
         error = vtk%xml_writer%write_piece()
         error = vtk%finalize()
      enddo vtr_loop

      ! only myrank == 0 save VTM file
      if (tree%myrank == 0_I4P) then
         error = vtm%initialize(filename=directory_//trim(basename)//'.vtm', scratch_units_number=max_level)
         vtm_group_loop : do l=1, max_level
            error = vtm%write_block(scratch=l, action='open', name='level-'//trim(str(l,.true.)))
         enddo vtm_group_loop
         vtm_filenames_loop : do while(tree%loop(node=node))
            b = node%block_index
            l = tree%level(code=node%code)
            error = vtm%write_block(scratch=l, action='write', filename=trim(basename)//'-block-'//trim(str(b,.true.))//&
                                                                        '-proc-'//trim(str(node%myrank,.true.))//'.vtr')
         enddo vtm_filenames_loop
         error = vtm%finalize()
      endif
   endassociate
   endsubroutine field_save_vtk

   subroutine save_hdf5(tree, field, basename, directory)
   !< Save ADAM data in HDF5 format.
   type(tree_object),  intent(in)           :: tree             !< The tree.
   type(field_object), intent(in)           :: field            !< The field.
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

   associate(grid=>field%grid, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk)
   ! save H5 file (one for each process)
   ! open fortran interface
   call h5open_f(error)
   ! create a new file using default properties
   h5_file_name = trim(basename)//trim(str(field%myrank,.true.))//'.h5'
   call h5fcreate_f(h5_file_name, H5F_ACC_TRUNC_F, h5_file_id, error)
   ! create the dataspace
   call h5screate_simple_f(3_I4P, [int(ni,I8P),int(nj,I8P),int(nk,I8P)], h5_dspace_id, error)
   ! save all blocks in process
   do b=1, field%blocks_number
      h5_dset_name = 'u-'//trim(str(field%myrank,.true.))//'-'//trim(str(b,.true.))
      call h5dcreate_f(h5_file_id, h5_dset_name, H5T_NATIVE_DOUBLE, h5_dspace_id, h5_dset_id, error)
      call h5dwrite_f(h5_dset_id, H5T_NATIVE_DOUBLE, field%u(1:ni,1:nj,1:nk,b), [int(ni,I8P),int(nj,I8P),int(nk,I8P)], error)
      call h5dclose_f(h5_dset_id, error)
   enddo
   ! terminate access to the data space
   call h5sclose_f(h5_dspace_id, error)
   ! close the file
   call h5fclose_f(h5_file_id, error)
   ! close FORTRAN interface
   call h5close_f(error)

   ! save XDMF file (only master process does)
   if (tree%myrank == 0_I4P) then
      open(newunit=xdmf, file=trim(basename)//'.xdmf')
      write(xdmf, '(A)') '<?xml version="1.0" encoding="utf-8"?>'
      write(xdmf, '(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="3.0">'
      write(xdmf, '(A)') '  <Domain>'
      write(xdmf, '(A)') '    <Grid Name="ADAM" GridType="Collection">'
      do while(tree%loop(node=node))
         b = node%block_index
         h5_file_name = trim(basename)//trim(str(node%myrank,.true.))//'.h5'
         h5_dset_name = 'u-'//trim(str(node%myrank,.true.))//'-'//trim(str(b,.true.))
         call tree%morton_to_coordinates(code=node%code, i=i, j=j, k=k, l=l)
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
endmodule adam_objects
