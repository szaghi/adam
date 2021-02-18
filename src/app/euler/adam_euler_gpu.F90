!< ADAM, ADAM for Euler equation, GPU backend.
program adam_euler_gpu
!< ADAM, ADAM for Euler equation, GPU backend.

use adam_adam_object, only : adam_object
use adam_equation_euler_gpu_object, only : equation_euler_gpu_object, BC_EXTRAPOLATION, BC_INFLOW
use adam_parameters
use FINER, only : file_ini
use PENF
use CUDAFOR
use MPI

implicit none

type(adam_object)               :: adam                         !< ADAM.
type(equation_euler_gpu_object) :: euler                        !< Euler equations system.
integer(I4P)                    :: t                            !< Counter.
real(R8P)                       :: time                         !< Time.
real(R8P)                       :: time_max=0.25_R8P            !< Maximum time of integration.
integer(I4P)                    :: t_max=huge(1_I4P)            !< Maximum time iteration.
integer(I4P)                    :: n_save=1                     !< Frequency of saving output.
character(999)                  :: output_basename              !< Output base name.
integer(I4P)                    :: nic_regions                  !< Number of initial conditions regions.
real(R8P), allocatable          :: emin_icr(:,:), emax_icr(:,:) !< Initial conditions regions extents.
real(R8P), allocatable          :: rho_icr(:)                   !< Initial conditions regions density.
real(R8P), allocatable          :: velocity_icr(:,:)            !< Initial conditions regions velocity.
real(R8P), allocatable          :: pressure_icr(:)              !< Initial conditions regions pressure.
real(R8P)                       :: timing(1:2)                  !< Tic toc timing.

call initialize(filename='src/app/euler/adam_euler_sod.ini')

! call euler%refine_uniform(refinement_levels=3)

call euler%set_initial_conditions(nic_regions=nic_regions,   &
                                  emin_icr=emin_icr,         &
                                  emax_icr=emax_icr,         &
                                  rho_icr=rho_icr,           &
                                  velocity_icr=velocity_icr, &
                                  pressure_icr=pressure_icr)

call euler%save_hdf5(output_basename=output_basename, t=0, time=0._R8P)

! call euler%amr_update(iterations=1)

time = 0._R8P
t = 0
call MPI_BARRIER(MPI_COMM_WORLD, adam%error) ; timing(1) = MPI_Wtime()
integration: do
   t = t + 1
   ! call euler%amr_update(iterations=1)
   call euler%compute_dt
   if (time + euler%dt > time_max) euler%dt = time_max - time
   if (mod(t,1)==0.and.adam%myrank==0) call euler%print_progress(t=t, time=time, time_max=time_max)
   call euler%integrate(t=time)
   if (mod(t,n_save)==0) call euler%save_hdf5(output_basename=output_basename, t=t, time=time)
   time = time + euler%dt
   if (time>=time_max.or.t>=t_max) exit integration
enddo integration
call MPI_BARRIER(MPI_COMM_WORLD, adam%error) ; timing(2) = MPI_Wtime()
print '(A, F18.10)', 'timing: ', (timing(2) - timing(1))/t

call euler%save_hdf5(output_basename=output_basename, t=t, time=time)

call adam%finalize
contains
   subroutine initialize(filename)
   !< Parse parameters file getting simulation input data.
   character(*), intent(in)  :: filename                     !< Parameters file name.
   type(file_ini)            :: file_parameters              !< Simulation parameters file handler.
   integer(I4P)              :: ni, nj, nk                   !< Cells number.
   integer(I4P)              :: ngc                          !< Ghost cells number.
   real(R8P)                 :: emin(3), emax(3)             !< Domain extents.
   integer(I4P)              :: bc_type(6)                   !< Type of boundary conditions in the 6 faces of grid.
   integer(I8P)              :: nodes_number                 !< Total number of nodes.
   integer(I4P)              :: max_level                    !< Maximum level of refinement.
   integer(I4P)              :: nv                           !< Number of variables.
   integer(I4P)              :: nb                           !< Number of blocks for each MPI process.
   integer(I4P)              :: ns                           !< Number of species.
   real(R8P)                 :: CFL                          !< Stability CFL coeffiecinet.
   integer(I4P)              :: nrk                          !< Number of Runge-Kutta stages.
   logical                   :: null_xyz(3)                  !< Flag to nullify fluxes computation for 1D or 2D sim.
   integer(I4P)              :: weno_stencils                !< WENO stencils number/dimension.
   integer(I4P)              :: fields_gpu_number            !< Number of fields allocated on GPU (for checkup).
   character(:), allocatable :: section_name                 !< Section name.
   integer(I4P)              :: icr                          !< Counter.

   call file_parameters%initialize(filename=trim(filename))
   call file_parameters%load

   call file_parameters%get(section_name='adam', option_name='ni'          , val=ni          )
   call file_parameters%get(section_name='adam', option_name='nj'          , val=nj          )
   call file_parameters%get(section_name='adam', option_name='nk'          , val=nk          )
   call file_parameters%get(section_name='adam', option_name='ngc'         , val=ngc         )
   call file_parameters%get(section_name='adam', option_name='emin_x'      , val=emin(1)     )
   call file_parameters%get(section_name='adam', option_name='emin_y'      , val=emin(2)     )
   call file_parameters%get(section_name='adam', option_name='emin_z'      , val=emin(3)     )
   call file_parameters%get(section_name='adam', option_name='emax_x'      , val=emax(1)     )
   call file_parameters%get(section_name='adam', option_name='emax_y'      , val=emax(2)     )
   call file_parameters%get(section_name='adam', option_name='emax_z'      , val=emax(3)     )
   call file_parameters%get(section_name='adam', option_name='nodes_number', val=nodes_number)
   call file_parameters%get(section_name='adam', option_name='max_level'   , val=max_level   )
   call file_parameters%get(section_name='adam', option_name='nv'          , val=nv          )
   call file_parameters%get(section_name='adam', option_name='nb'          , val=nb          )

   call adam%initialize(ni=ni, nj=nj, nk=nk, ngc=ngc, max_level=max_level, emin=emin, emax=emax,  &
                        bc_type=[BC_EXTRAPOLATION,BC_EXTRAPOLATION,                               &
                                 BC_EXTRAPOLATION,BC_EXTRAPOLATION,                               &
                                 BC_EXTRAPOLATION,BC_EXTRAPOLATION],                              &
                        nodes_number=nodes_number, nv=nv, nb=nb)

   call file_parameters%get(section_name='euler', option_name='time_max'         , val=time_max         )
   call file_parameters%get(section_name='euler', option_name='t_max'            , val=t_max            )
   call file_parameters%get(section_name='euler', option_name='output_basename'  , val=output_basename  )
   call file_parameters%get(section_name='euler', option_name='n_save'           , val=n_save           )
   call file_parameters%get(section_name='euler', option_name='ns'               , val=ns               )
   call file_parameters%get(section_name='euler', option_name='CFL'              , val=CFL              )
   call file_parameters%get(section_name='euler', option_name='nrk'              , val=nrk              )
   call file_parameters%get(section_name='euler', option_name='null_x'           , val=null_xyz(1)      )
   call file_parameters%get(section_name='euler', option_name='null_y'           , val=null_xyz(2)      )
   call file_parameters%get(section_name='euler', option_name='null_z'           , val=null_xyz(3)      )
   call file_parameters%get(section_name='euler', option_name='weno_stencils'    , val=weno_stencils    )
   call file_parameters%get(section_name='euler', option_name='fields_gpu_number', val=fields_gpu_number)
   call euler%initialize(adam=adam, ns=ns, CFL=CFL, nrk=nrk,             &
                         null_xyz=null_xyz, weno_stencils=weno_stencils, &
                         fields_gpu_number=fields_gpu_number)

   call file_parameters%get(section_name='initial_conditions', option_name='nic_regions', val=nic_regions)
   allocate(    emin_icr(1:3, 1:nic_regions))
   allocate(    emax_icr(1:3, 1:nic_regions))
   allocate(     rho_icr(     1:nic_regions))
   allocate(velocity_icr(1:3, 1:nic_regions))
   allocate(pressure_icr(     1:nic_regions))
   do icr=1, nic_regions
      section_name = 'initial_conditions_'//trim(str(icr,.true.))
      call file_parameters%get(section_name=section_name, option_name='emin_x'    , val=    emin_icr(1,icr))
      call file_parameters%get(section_name=section_name, option_name='emin_y'    , val=    emin_icr(2,icr))
      call file_parameters%get(section_name=section_name, option_name='emin_z'    , val=    emin_icr(3,icr))
      call file_parameters%get(section_name=section_name, option_name='emax_x'    , val=    emax_icr(1,icr))
      call file_parameters%get(section_name=section_name, option_name='emax_y'    , val=    emax_icr(2,icr))
      call file_parameters%get(section_name=section_name, option_name='emax_z'    , val=    emax_icr(3,icr))
      call file_parameters%get(section_name=section_name, option_name='rho'       , val=     rho_icr(  icr))
      call file_parameters%get(section_name=section_name, option_name='velocity_x', val=velocity_icr(1,icr))
      call file_parameters%get(section_name=section_name, option_name='velocity_y', val=velocity_icr(2,icr))
      call file_parameters%get(section_name=section_name, option_name='velocity_z', val=velocity_icr(3,icr))
      call file_parameters%get(section_name=section_name, option_name='pressure'  , val=pressure_icr(  icr))
   enddo
   endsubroutine initialize
endprogram adam_euler_gpu
