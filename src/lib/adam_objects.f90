!< ADAM, objects classes.
module adam_objects
!< ADAM, objects classes.

use adam_field_object
use adam_tree_node_object
use adam_tree_bucket_object
use adam_tree_object
use PENF, only : R8P, I4P, I8P, str
use stringifor, only : string
use vtk_fortran, only : vtK_file, vtm_file

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

   ! real(R8P)                          :: dx, dy, dz    !< Domain delta space.
   ! real(R8P)                          :: dxl, dyl, dzl !< Local delta space.
   ! real(R8P)                          :: emin(3), emax(3)
   integer(I4P)                       :: i, j, k, l, b !< Counter.

   ! dx = field%domain_emax(1) - field%domain_emin(1)
   ! dy = field%domain_emax(2) - field%domain_emin(2)
   ! dz = field%domain_emax(3) - field%domain_emin(3)

   threshold_ = 2200000000000_R8P ; if (present(threshold)) threshold_ = threshold
   do while(tree%loop(node=node))
      if (node%code==8_I8P) then
      print*, 'cazzo code:'//trim(str(node%code,.true.))//' ',trim(str(field%emin(:,node%block_index) ,.true.))
      print*, 'cazzo code:'//trim(str(node%code,.true.))//' ',trim(str(field%emax(:,node%block_index) ,.true.))
      call tree%morton_to_coordinates(code=node%code, i=i, j=j, k=k, l=l)
      print*, 'cazzo code:'//trim(str(node%code,.true.))//' ',trim(str([i,j,k,l] ,.true.))
      print*, 'cazzo code:'//trim(str(node%code,.true.))//' ',trim(str(node%block_index ,.true.))
      endif
      block_center = (field%emax(:,node%block_index) + field%emin(:,node%block_index)) / 2._R8P
      block_diagonal = sqrt((field%emax(1,node%block_index) - field%emin(1,node%block_index))**2 + &
                            (field%emax(2,node%block_index) - field%emin(2,node%block_index))**2 + &
                            (field%emax(3,node%block_index) - field%emin(3,node%block_index))**2)

      ! call tree%morton_to_coordinates(code=node%code, i=i, j=j, k=k, l=l)
      ! dxl = dx / 2**l
      ! dyl = dy / 2**l
      ! dzl = dz / 2**l
      ! emin(1) = i * dxl ; emax(1) = emin(1) + dxl
      ! emin(2) = j * dyl ; emax(2) = emin(2) + dyl
      ! emin(3) = k * dzl ; emax(3) = emin(3) + dzl
      ! block_center = (emax(:) + emin(:)) / 2._R8P
      ! block_diagonal = sqrt((emax(1) - emin(1))**2 + &
      !                       (emax(2) - emin(2))**2 + &
      !                       (emax(3) - emin(3))**2)

      ! associate (&
      associate (emin=>field%emin(:,node%block_index), emax=>field%emax(:,node%block_index), &
                 ni=>field%ni, nj=>field%nj, nk=>field%nk)
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
   integer(I4P)                             :: i,j,k      !< Counter.
   integer(I4P)                             :: max_level  !< Maximum level.

   directory_ = '' ; if (present(directory)) directory_ = trim(directory)
   associate(emin=>field%emin, emax=>field%emax, ni=>field%ni, nj=>field%nj, nk=>field%nk, &
             gc1=>field%gc1, gc2=>field%gc2, gc3=>field%gc3,  gc4=>field%gc4, gc5=>field%gc5, gc6=>field%gc6)

      max_level = 0_I4P
      vtk_loop : do while(tree%loop(node=node))
         b = node%block_index
         max_level = max(max_level, tree%level(code=node%code))
         ! only the process having node can save VTR file
         if (field%myrank == node%myrank) then
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
         endif
      enddo vtk_loop

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
endmodule adam_objects
