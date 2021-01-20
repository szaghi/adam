!< ADAM, test track a sphere.
program adam_test_track_sphere
!< ADAM, test track a sphere.

use adam_adam_object
use adam_parameters
use PENF

implicit none

type(adam_object)  :: adam            !< ADAM.
logical            :: is_grid_changed !< Flag to check grid changes.
integer(I4P)       :: t, st           !< Counter.

print '(A)', 'sphere tracking, initialize grid'
call adam%initialize(nb=190000, emin=[0._R8P,0._R8P,0._R8P], emax=[1._R8P,1._R8P,1._R8P])

do t=1,2
   call adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
   call adam%amr_update(do_mpi_redistribute=.false.)
enddo

print '(A)', 'catch first position'
do t=1,8
   print*, ''
   print '(A)', 'track iteration '//trim(str(t, .true.))
   call adam%tree%mark_sphere(center=[0.2_R8P,0.5_R8P,0.5_R8P], radius=0.1_R8P)
   call adam%amr_update(do_mpi_redistribute=.false.)
enddo

print*, ''
print '(A)', 'move sphere'
do t=1,10
   print*, ''
   print '(A)', 'track iteration '//trim(str(t, .true.))//' position x='//trim(str(0.2_R8P + t*0.05_R8P))
   sub_iteration_loop : do st=1, 10
      print '(A)', '  track su-iteration '//trim(str(st, .true.))
      call adam%tree%mark_sphere(center=[0.2_R8P+t*0.05_R8P,0.5_R8P,0.5_R8P], radius=0.1_R8P)
      call adam%amr_update(do_mpi_redistribute=.false., is_grid_changed=is_grid_changed)
      if (.not.is_grid_changed) exit sub_iteration_loop
   enddo sub_iteration_loop
enddo
call adam%save_vtk(basename='sphere-'//trim(strz(t,9)))
endprogram adam_test_track_sphere
