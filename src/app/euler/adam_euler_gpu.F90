!< ADAM, ADAM for Euler equation, GPU backend.
program adam_euler_gpu
!< ADAM, ADAM for Euler equation, GPU backend.

use adam_adam_object
use adam_equation_euler_gpu_object
use adam_parameters
use PENF
use CUDAFOR
use MPI

implicit none

type(adam_object)               :: adam              !< ADAM.
type(equation_euler_gpu_object) :: euler             !< Euler equations system.
integer(I4P)                    :: l, t, st          !< Counter.
logical                         :: is_grid_changed   !< Flag to check grid changes.
real(R8P)                       :: time              !< Time.
real(R8P)                       :: time_max=0.25_R8P !< Maximum time of integration.
integer(I4P)                    :: n_save=10         !< Frequency of saving output.
real(R8P)                       :: timing(1:2)       !< Tic toc timing.

print '(A)', 'Euler equation integration'
call adam%initialize(max_level=7,                                              &
                     emin=[0._R8P,0._R8P,0._R8P], emax=[1._R8P,1._R8P,1._R8P], &
                     ni=16_I4P, nj=16_I4P, nk=16_I4P, ngc=2_I4P,               &
                     bc_type=[BC_EXTRAPOLATION,BC_EXTRAPOLATION,               &
                              BC_EXTRAPOLATION,BC_EXTRAPOLATION,               &
                              BC_EXTRAPOLATION,BC_EXTRAPOLATION],              &
                     nb=1500, nv=5, nodes_number=50000_I8P)

call euler%initialize(field=adam%field, ns=1, CFL=0.3_R8P,              &
                      null_xyz=[.false.,.false.,.false.], weno_s=2_I4P, &
                      fields_gpu_number=14)

 print '(A)', 'create 3 levels of refinement'
 do l=1, 3
    call adam%tree%mark_all_nodes(mark=TO_BE_REFINED)
    call adam%amr_update(do_blocks_reorder=.false., do_mpi_redistribute=.true.)
    print '(A)', 'refine ADAM at level '//trim(str(l))
    print *, 'blocks_number: ',adam%tree%nodes_number
 enddo
print '(A)', 'set initial conditions'
call euler%set_initial_conditions(ic_type='kt-02')

call euler%copy_cpu_gpu
call euler%copy_gpu_cpu(compute_q_aux=.true.)
call adam%save_hdf5(basename='euler-sod-'//trim(strz(0,9)),           &
                    q=euler%field%q,                                  &
                    q_aux=euler%q_aux,                                &
                    q_name=['rho  ','rho-u','rho-v','rho-w','rho-E'], &
                    q_aux_name=['c1','r ','u ','v ','w ','g ','p '],  &
                    with_cell_morton=.true.)

print '(A)', 'track initial discontinuity'
track: do t=1, -1
   if (mod(t,1)==0.and.adam%myrank==0) print '(A)', 'track iteration '//trim(str(t, .true.))
   call euler%mark_by_grad_rho(grad_tol=0.05_R8P, delta_fine=0.006_R8P, delta_coarse=0.015_R8P)
   call euler%update_ghost_gpu(q_gpu=euler%q_gpu)
   call euler%copy_gpu_cpu()
   call adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false., is_grid_changed=is_grid_changed)
   call euler%copy_cpu_gpu
   if (.not.is_grid_changed) exit track
enddo track

time = 0._R8P
t = 0
call MPI_BARRIER(MPI_COMM_WORLD, adam%error) ; timing(1) = MPI_Wtime()
integration: do
   t = t + 1
!    ! adapt grids tracking discontinuities
!    if (mod(t,10)==0) then
!       sub_track: do st=1, 10
!          if (adam%myrank==0) print '(A)', '  track discontinuities sub-iteration '//trim(str(st, .true.))
!          call euler%mark_by_grad_rho(grad_tol=2.5_R8P, delta_fine=0.010_R8P, delta_coarse=0.1_R8P)
!          call euler%update_ghost(q=adam%field%q)
!          call adam%amr_update(is_marked_by_field=.true., do_blocks_reorder=.false., is_grid_changed=is_grid_changed)
!          if (.not.is_grid_changed) exit sub_track
!       enddo sub_track
!    endif
   ! integrate Euler equations
   call euler%compute_dt
   if (mod(t,1)==0.and.adam%myrank==0) then
      print '(A)',       'blocks number: '//trim(str(adam%tree%nodes_number, .true.))
      print '(A,F12.7)', 'time step:     ', euler%dt
      print '(A,F12.7)', 'time:          ', time
      print '(A)',       't:             '//trim(str(t,.true.))
   endif
   call euler%integrate(t=time)
   if (mod(t,n_save)==0) then
      call euler%copy_gpu_cpu(compute_q_aux=.true.)
      call adam%save_hdf5(basename='euler-sod-'//trim(strz(t,9)),           &
                          q=euler%field%q,                                  &
                          q_aux=euler%q_aux,                                &
                          q_name=['rho  ','rho-u','rho-v','rho-w','rho-E'], &
                          q_aux_name=['r1','r ','u ','v ','w ','g ','p '],  &
                          with_cell_morton=.true.)
   endif
   time = time + euler%dt
   if (time>=time_max .or. t > 100) exit integration
   !exit integration
enddo integration
call MPI_BARRIER(MPI_COMM_WORLD, adam%error) ; timing(2) = MPI_Wtime()
print '(A, F18.10)', 'timing: ', (timing(2) - timing(1))/t

!call euler%copy_gpu_cpu(compute_q_aux=.true.)
!call adam%save_hdf5(basename='euler-sod-'//trim(strz(t,9)),           &
!                    q=euler%field%q,                                  &
!                    q_aux=euler%q_aux,                                &
!                    q_name=['rho  ','rho-u','rho-v','rho-w','rho-E'], &
!                    q_aux_name=['r1','r ','u ','v ','w ','g ','p '],  &
!                    with_cell_morton=.true.)

call adam%finalize
endprogram adam_euler_gpu
