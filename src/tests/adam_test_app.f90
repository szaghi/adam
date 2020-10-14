program adam_test_app

! TODO shift dell'ancestor e find finest_code

use adam_objects
use PENF, only : R8P, I8P, I4P, str
use vtk_fortran, only : vtK_file, vtm_file

implicit none

type(tree_object)               :: tree               !< The tree.
type(tree_node_object), pointer :: node               !< Pointer to node.
type(field_object)              :: field              !< Field.
type(vtk_file)                  :: vtk                !< VTK file.
type(vtm_file)                  :: vtm                !< VTM file.
integer(I8P)                    :: code               !< Counter.
integer(I8P), allocatable       :: block_to_refine(:) !< List of field blocks to be refined.
integer(I8P), allocatable       :: block_refined(:)   !< List of field refined blocks.
real(R8P)                       :: emin(3,100), emax(3,100), dx, dy, dz
integer(I4P)                    :: b, ib, i, j, k, ic, ii
integer(I4P)                    :: ic1
integer(I4P)                    :: ic2
integer(I4P)                    :: ic3
integer(I4P)                    :: ic4
integer(I4P)                    :: ic5
integer(I4P)                    :: ic6
integer(I4P)                    :: ic7
integer(I4P)                    :: ic8
integer(I4P)                    :: error
character(:), allocatable       :: filenames

print '(A)', 'initialize'
call tree%initialize

emin = 0._R8P
emax = 1._R8P
call field%initialize(nb=100, emin=emin, emax=emax)

print '(A)', 'add ancestor node -1'
call tree%add_node(code=-1_I8P, content=-1_I8P)
print '(A)', 'loop in tree'
do while(tree%loop(node=node))
   call tree%print_code_topology(code=node%code, block_index=.true., coordinates=.true.)
enddo
print*, ''

print '(A)', 'refine level 0'
call tree%refine(force_all=.true., block_to_refine=block_to_refine, block_refined=block_refined)
call field_refine
print '(A)', 'loop in tree'
do while(tree%loop(node=node))
   call tree%print_code_topology(code=node%code, block_index=.true., coordinates=.true.)
enddo
print*, ''

print '(A)', 'blocks index to refine'
print '(A)', trim(str(block_to_refine))
print '(A)', 'blocks index refined'
print '(A)', trim(str(block_refined))
print*, ''
print '(A)', 'refine level 1'
node => tree%node(code=3_I8P)
node%refinement_needed = TO_BE_REFINED
call tree%sanitize
call tree%refine(block_to_refine=block_to_refine, block_refined=block_refined)
call field_refine
print '(A)', 'loop in tree'
do while(tree%loop(node=node))
   call tree%print_code_topology(code=node%code, block_index=.true., coordinates=.true.)
enddo
print*, ''
print '(A)', 'blocks index to refine'
print '(A)', trim(str(block_to_refine))
print '(A)', 'blocks refined'
print '(A)', trim(str(block_refined))
print*, ''
print '(A)', 'refine index level 2'
node => tree%node(code=37_I8P)
node%refinement_needed = TO_BE_REFINED
! call tree%sanitize
call tree%refine(block_to_refine=block_to_refine, block_refined=block_refined)
call field_refine
print '(A)', 'loop in tree'
do while(tree%loop(node=node))
   call tree%print_code_topology(code=node%code, block_index=.true., coordinates=.true.)
enddo
print*, ''
print '(A)', 'blocks index to refine'
print '(A)', trim(str(block_to_refine))
print '(A)', 'blocks index refined'
print '(A)', trim(str(block_refined))

filenames = ''
do b=1, tree%nodes_number
   filenames = filenames//'block-'//trim(str(b,.true.))//'.vtr '

   error = vtk%initialize(format='binary', filename='block-'//trim(str(b,.true.))//'.vtr', &
                          mesh_topology='RectilinearGrid', &
                          nx1=0, nx2=1, ny1=0, ny2=1, nz1=0, nz2=1)
   error = vtk%xml_writer%write_piece(nx1=0, nx2=1, ny1=0, ny2=1, nz1=0, nz2=1)
   error = vtk%xml_writer%write_geo(x=[field%emin(1,b), field%emax(1,b)], &
                                    y=[field%emin(2,b), field%emax(2,b)], &
                                    z=[field%emin(3,b), field%emax(3,b)])
   error = vtk%xml_writer%write_dataarray(location='cell', action='open')
   error = vtk%xml_writer%write_dataarray(data_name='morton', x=[field%u(1,1,1,b)])
   error = vtk%xml_writer%write_dataarray(location='cell', action='close')
   error = vtk%xml_writer%write_piece()
   error = vtk%finalize()
enddo

error = vtm%initialize(filename='adam.vtm')
error = vtm%write_block(filenames=trim(filenames), name='adam')
error = vtm%finalize()

contains
   subroutine field_refine

   do b=1, size(block_to_refine, dim=1)
      ib = block_to_refine(b)
      ic1 = block_refined((b-1)*tree%ratio+1)
      field%u(:,:,:,ic1) = field%u(:,:,:,ib)*tree%ratio
      ic2 = block_refined((b-1)*tree%ratio+2)
      field%u(:,:,:,ic2) = field%u(:,:,:,ib)*tree%ratio
      ic3 = block_refined((b-1)*tree%ratio+3)
      field%u(:,:,:,ic3) = field%u(:,:,:,ib)*tree%ratio
      ic4 = block_refined((b-1)*tree%ratio+4)
      field%u(:,:,:,ic4) = field%u(:,:,:,ib)*tree%ratio
      ic5 = block_refined((b-1)*tree%ratio+5)
      field%u(:,:,:,ic5) = field%u(:,:,:,ib)*tree%ratio
      ic6 = block_refined((b-1)*tree%ratio+6)
      field%u(:,:,:,ic6) = field%u(:,:,:,ib)*tree%ratio
      ic7 = block_refined((b-1)*tree%ratio+7)
      field%u(:,:,:,ic7) = field%u(:,:,:,ib)*tree%ratio
      ic8 = block_refined((b-1)*tree%ratio+8)
      field%u(:,:,:,ic8) = field%u(:,:,:,ib)*tree%ratio
   enddo

   do b=1, size(block_to_refine, dim=1)
      ib = block_to_refine(b)
      dx = field%emax(1,ib)-field%emin(1,ib)
      dy = field%emax(2,ib)-field%emin(2,ib)
      dz = field%emax(3,ib)-field%emin(3,ib)
      ii = 1
      do k=0, 1
         do j=0, 1
            do i=0, 1
               ic = block_refined((b-1)*tree%ratio+ii)
               field%emin(1,ic) = field%emin(1,ib) + i * dx/2
               field%emin(2,ic) = field%emin(2,ib) + j * dy/2
               field%emin(3,ic) = field%emin(3,ib) + k * dz/2

               field%emax(1,ic) = field%emin(1,ic) + dx/2
               field%emax(2,ic) = field%emin(2,ic) + dy/2
               field%emax(3,ic) = field%emin(3,ic) + dz/2

               ii = ii + 1
            enddo
         enddo
      enddo
   enddo

   endsubroutine field_refine
endprogram adam_test_app
