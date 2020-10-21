!< ADAM, objects classes.
module adam_objects
!< ADAM, objects classes.

use adam_field_object
use adam_tree_node_object
use adam_tree_bucket_object
use adam_tree_object
use PENF, only : I4P, str
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
public :: field_save_vtk

contains
   subroutine field_save_vtk(tree, field, basename)
   !< Save field in VTK files.
   type(tree_object),  intent(in)  :: tree         !< The tree.
   type(field_object), intent(in)  :: field        !< The field.
   character(*),       intent(in)  :: basename     !< Base name of output files.
   integer(I4P)                    :: error        !< Error trapping flag.
   ! character(:), allocatable       :: filenames !< File names list.
   type(string), allocatable       :: filenames(:) !< File names list, per level.
   type(vtk_file)                  :: vtk          !< VTK file handler.
   type(vtm_file)                  :: vtm          !< VTM file handler.
   type(tree_node_object), pointer :: node         !< Pointer to node.
   integer(I4P)                    :: b, l         !< Counter.

   allocate(filenames(tree%max_level))
   do l=1, tree%max_level
      filenames(l) = ''
   enddo
   associate(emin=>field%emin, emax=>field%emax, ni=>field%ni, nj=>field%nj, nk=>field%nk, &
             gc1=>field%gc1, gc2=>field%gc2, gc3=>field%gc3,  gc4=>field%gc4, gc5=>field%gc5, gc6=>field%gc6)
      do while(tree%loop(node=node))
         b = node%block_index
         l = tree%level(code=node%code)
         filenames(l) = filenames(l)//trim(basename)//'-block-'//trim(str(b,.true.))//'.vtr '
         error = vtk%initialize(format='raw', filename=trim(basename)//'-block-'//trim(str(b,.true.))//'.vtr', &
                                mesh_topology='RectilinearGrid',                                               &
                                nx1=0, nx2=ni, ny1=0, ny2=nj, nz1=0, nz2=nk)
         error = vtk%xml_writer%write_fielddata(action='open')
         error = vtk%xml_writer%write_fielddata(data_name='Morton', x=field%code(b))
         error = vtk%xml_writer%write_fielddata(action='close')
         error = vtk%xml_writer%write_piece(nx1=0, nx2=ni, ny1=0, ny2=nj, nz1=0, nz2=nk)
         error = vtk%xml_writer%write_geo(x=field%x(0:ni,b), &
                                          y=field%y(0:nj,b), &
                                          z=field%z(0:nk,b))
         error = vtk%xml_writer%write_dataarray(location='cell', action='open')
         error = vtk%xml_writer%write_dataarray(data_name='u', x=[field%u(1:ni,1:nj,1:nk,b)])
         error = vtk%xml_writer%write_dataarray(location='cell', action='close')
         error = vtk%xml_writer%write_piece()
         error = vtk%finalize()
      enddo
   endassociate

   error = vtm%initialize(filename=trim(basename)//'.vtm')
   do l=1, tree%max_level
      if (filenames(l)/='') &
      error = vtm%write_block(filenames=trim(filenames(l)%chars()), name='level-'//trim(str(l,.true.)))
   enddo
   error = vtm%finalize()
   endsubroutine field_save_vtk
endmodule adam_objects
