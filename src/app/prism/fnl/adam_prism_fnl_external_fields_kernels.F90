!< ADAM, PRISM external fields definition, FNL backend kernels.

#include "fundal.H"

module adam_prism_fnl_external_fields_kernels
!< ADAM, PRISM external fields definition, FNL backend kernels.

! ADAM modules
use :: adam_fnl_field_object
! PRISM modules
use :: adam_prism_external_fields_object
use :: adam_prism_parameters
! third party modules
use :: penf

implicit none
private
public :: external_fields_initialize_dev
public :: add_external_fields_dev
public :: sub_external_fields_dev
public :: add_external_fields_dev_interface
public :: sub_external_fields_dev_interface
public :: add_external_fields_rmf_dev
public :: sub_external_fields_rmf_dev
public :: add_external_fields_uniform_dev
public :: sub_external_fields_uniform_dev

! pointer (abstract) procedures
procedure(add_external_fields_dev_interface), pointer :: add_external_fields_dev=>null() !< Add external fields.
procedure(sub_external_fields_dev_interface), pointer :: sub_external_fields_dev=>null() !< Subtract external fields.

interface
   subroutine add_external_fields_dev_interface(external_fields, field_gpu, dt, time, q_gpu, gamm)
   import :: prism_external_fields_object, field_fnl_object, R8P
   type(prism_external_fields_object), intent(in)    :: external_fields !< External fields handler.
   type(field_fnl_object),             intent(in)    :: field_gpu       !< Field.
   real(R8P),                          intent(in)    :: dt              !< Time step.
   real(R8P),                          intent(in)    :: time            !< Current time.
   real(R8P),                          intent(inout) :: q_gpu(1:,              &
                                                              1-field_gpu%ngc:,&
                                                              1-field_gpu%ngc:,&
                                                              1-field_gpu%ngc:,&
                                                              1:)       !< Conservative variables.
   real(R8P), optional,                intent(in)    :: gamm            !< Gamma values of RK.
   endsubroutine add_external_fields_dev_interface

   subroutine sub_external_fields_dev_interface(external_fields, field_gpu, dt, time, q_gpu, gamm)
   import :: prism_external_fields_object, field_fnl_object, R8P
   type(prism_external_fields_object), intent(in)    :: external_fields !< External fields handler.
   type(field_fnl_object),             intent(in)    :: field_gpu       !< Field.
   real(R8P),                          intent(in)    :: dt              !< Time step.
   real(R8P),                          intent(in)    :: time            !< Current time.
   real(R8P),                          intent(inout) :: q_gpu(1:,              &
                                                              1-field_gpu%ngc:,&
                                                              1-field_gpu%ngc:,&
                                                              1-field_gpu%ngc:,&
                                                              1:)       !< Conservative variables.
   real(R8P), optional,                intent(in)    :: gamm            !< Gamma values of RK.
   endsubroutine sub_external_fields_dev_interface
endinterface

contains
   subroutine external_fields_initialize_dev(external_fields)
   !< Initialize external fields device kernels.
   type(prism_external_fields_object), intent(in) :: external_fields !< External fields handler.

   select case(external_fields%ef_type)
   case(EF_TYPE_RMF)
      add_external_fields_dev => add_external_fields_rmf_dev
      sub_external_fields_dev => sub_external_fields_rmf_dev
   case(EF_TYPE_UNIFORM_FIELD)
      add_external_fields_dev => add_external_fields_uniform_dev
      sub_external_fields_dev => sub_external_fields_uniform_dev
   !case(EF_TYPE_MAGNETIC_NOZZLE)
   !   add_external_fields => self%external_fields%add_external_fields_magnetic_nozzle
   !case(EF_TYPE_RMF_AND_MAGNETIC_NOZZLE)
   !   add_external_fields => self%external_fields%add_external_fields_rmf_and_magnetic_nozzle
   endselect
   endsubroutine external_fields_initialize_dev

   subroutine add_external_fields_rmf_dev(external_fields, field_gpu, dt, time, q_gpu, gamm)
   !< Add rotating magnetic field to the field.
   type(prism_external_fields_object), intent(in)    :: external_fields !< External fields handler.
   type(field_fnl_object),             intent(in)    :: field_gpu       !< Field.
   real(R8P),                          intent(in)    :: dt              !< Time step.
   real(R8P),                          intent(in)    :: time            !< Current time.
   real(R8P),                          intent(inout) :: q_gpu(1:,              &
                                                              1-field_gpu%ngc:,&
                                                              1-field_gpu%ngc:,&
                                                              1-field_gpu%ngc:,&
                                                              1:)       !< Conservative variables.
   real(R8P), optional,                intent(in)    :: gamm            !< Gamma values of RK.
   real(R8P)                                         :: time_next       !< Time at the next sub-step.
   real(R8P)                                         :: omega           !< Omega frequency.

   associate(blocks_number=>field_gpu%blocks_number, ni=>field_gpu%ni, nj=>field_gpu%nj, nk=>field_gpu%nk, ngc=>field_gpu%ngc, &
             RMF_frequency=>external_fields%RMF_frequency, RMF_B_amplitude=>external_fields%RMF_B_amplitude,                   &
             ef_alpha=>external_fields%alpha, ef_beta=>external_fields%beta, ef_gamma=>external_fields%gamm)
   time_next = time + dt ; if (present(gamm)) time_next = time + dt*gamm
   omega = 2.0_R8P*PI*RMF_frequency
   call add_external_fields_rmf_dev_kernel(ni              = ni                  ,&
                                           nj              = nj                  ,&
                                           nk              = nk                  ,&
                                           ngc             = ngc                 ,&
                                           blocks_number   = blocks_number       ,&
                                           time_next       = time_next           ,&
                                           ef_alpha        = ef_alpha            ,&
                                           ef_beta         = ef_beta             ,&
                                           ef_gamma        = ef_gamma            ,&
                                           omega           = omega               ,&
                                           RMF_B_amplitude = RMF_B_amplitude     ,&
                                           displacement_scale = EPS0,&
                                           x_cell_gpu      = field_gpu%x_cell_gpu,&
                                           y_cell_gpu      = field_gpu%y_cell_gpu,&
                                           z_cell_gpu      = field_gpu%z_cell_gpu,&
                                           q_gpu           = q_gpu)
   call correct_pec_external_ghosts(external_fields, field_gpu, time_next, 1._R8P, q_gpu)
   endassociate
   contains
      subroutine add_external_fields_rmf_dev_kernel(ni,nj,nk,ngc,blocks_number,                               &
                                                    time_next,ef_alpha,ef_beta,ef_gamma,omega,RMF_B_amplitude,displacement_scale,&
                                                    x_cell_gpu,y_cell_gpu,z_cell_gpu,q_gpu)
      !< Add rotating magnetic field to the field, device kernel.
      integer(I4P), intent(in)    :: ni,nj,nk,ngc,blocks_number        !< Grids dimensions.
      real(R8P),    intent(in)    :: time_next                         !< Time at the next sub-step.
      real(R8P),    intent(in)    :: RMF_B_amplitude                   !< Rotating magnetic field amplitude.
      real(R8P),    intent(in)    :: displacement_scale                 !< Displacement to magnetic field scale.
	   integer(I4P), intent(in)    :: ef_alpha                          !< RMF rotation axis coordinate 1
	   integer(I4P), intent(in)    :: ef_beta                           !< RMF rotation axis coordinate 2
	   integer(I4P), intent(in)    :: ef_gamma                          !< RMF rotation axis coordinate 3
      real(R8P),    intent(in)    :: omega                             !< Omega frequency.
      real(R8P),    intent(in)    :: x_cell_gpu(1:,1-ngc:)             !< Cells x coordinates on GPU [nb,1-ngc:ni+ngc].
      real(R8P),    intent(in)    :: y_cell_gpu(1:,1-ngc:)             !< Cells y coordinates on GPU [nb,1-ngc:nj+ngc].
      real(R8P),    intent(in)    :: z_cell_gpu(1:,1-ngc:)             !< Cells z coordinates on GPU [nb,1-ngc:nk+ngc].
      real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field cell centered variables.
      real(R8P)                   :: B_r, B_theta                      !< Radial/azimuthal components of the rotating B field.
      real(R8P)                   :: cell_coord(3)                     !< Cell coordinates.
      real(R8P)                   :: phase                             !< Phase.
      real(R8P)                   :: theta                             !< Angle in cylindrical coordinates.
      real(R8P)                   :: x,y,r,c,s                         !< Buffer.
      integer(I4P)                :: i,j,k,b                           !< Counter.

      !$acc parallel loop independent gang vector collapse(4)                                               &
      !$acc& DEVICEVAR(x_cell_gpu,y_cell_gpu,z_cell_gpu,q_gpu)                                              &
      !$acc& firstprivate(ni,nj,nk,blocks_number,time_next,ef_alpha,ef_beta,ef_gamma,omega,RMF_B_amplitude,displacement_scale) &
      !$acc& private(B_r,B_theta,cell_coord,phase,theta,x,y,r,c,s)
      !$omp OMPLOOP collapse(4) &
      !$omp& DEVICEPTR(x_cell_gpu,y_cell_gpu,z_cell_gpu,q_gpu) &
      !$omp& firstprivate(ni,nj,nk,blocks_number,time_next,ef_alpha,ef_beta,ef_gamma,omega,RMF_B_amplitude,displacement_scale) &
      !$omp& private(B_r,B_theta,cell_coord,phase,theta,x,y,r,c,s)
      do b = 1, blocks_number
      do k = 1 - ngc, nk + ngc
      do j = 1 - ngc, nj + ngc
      do i = 1 - ngc, ni + ngc
         cell_coord = [x_cell_gpu(b,i), y_cell_gpu(b,j), z_cell_gpu(b,k)]
         x = cell_coord(ef_alpha)
         y = cell_coord(ef_beta)
         r = sqrt(x*x + y*y)
         theta = atan2(y, x)
         phase = omega*time_next - theta
         B_r     = RMF_B_amplitude*cos(phase)
         B_theta = RMF_B_amplitude*sin(phase)
         c = cos(theta)
         s = sin(theta)
         q_gpu(b,i,j,k,ef_alpha+3) = q_gpu(b,i,j,k,ef_alpha+3) + B_r*c - B_theta*s
         q_gpu(b,i,j,k,ef_beta +3) = q_gpu(b,i,j,k,ef_beta +3) + B_r*s + B_theta*c
         q_gpu(b,i,j,k,ef_gamma  ) = q_gpu(b,i,j,k,ef_gamma  ) + r*omega*RMF_B_amplitude*cos(phase)*displacement_scale
      enddo
      enddo
      enddo
      enddo
      endsubroutine add_external_fields_rmf_dev_kernel
   endsubroutine add_external_fields_rmf_dev

   subroutine sub_external_fields_rmf_dev(external_fields, field_gpu, dt, time, q_gpu, gamm)
   !< Subtract rotating magnetic field to the field, device kernel.
   type(prism_external_fields_object), intent(in)    :: external_fields !< External fields handler.
   type(field_fnl_object),             intent(in)    :: field_gpu       !< Field.
   real(R8P),                          intent(in)    :: dt              !< Time step.
   real(R8P),                          intent(in)    :: time            !< Current time.
   real(R8P),                          intent(inout) :: q_gpu(1:,              &
                                                              1-field_gpu%ngc:,&
                                                              1-field_gpu%ngc:,&
                                                              1-field_gpu%ngc:,&
                                                              1:)       !< Conservative variables.
   real(R8P), optional,                intent(in)    :: gamm            !< Gamma values of RK.
   real(R8P)                                         :: time_next       !< Time at the next sub-step.
   real(R8P)                                         :: omega           !< Omega frequency.

   associate(blocks_number=>field_gpu%blocks_number, ni=>field_gpu%ni, nj=>field_gpu%nj, nk=>field_gpu%nk, ngc=>field_gpu%ngc, &
             RMF_frequency=>external_fields%RMF_frequency, RMF_B_amplitude=>external_fields%RMF_B_amplitude,                   &
             ef_alpha=>external_fields%alpha, ef_beta=>external_fields%beta, ef_gamma=>external_fields%gamm)
   time_next = time + dt ; if (present(gamm)) time_next = time + dt*gamm
   omega = 2.0_R8P*PI*RMF_frequency
   call sub_external_fields_rmf_dev_kernel(ni              = ni                  ,&
                                           nj              = nj                  ,&
                                           nk              = nk                  ,&
                                           ngc             = ngc                 ,&
                                           blocks_number   = blocks_number       ,&
                                           time_next       = time_next           ,&
                                           ef_alpha        = ef_alpha            ,&
                                           ef_beta         = ef_beta             ,&
                                           ef_gamma        = ef_gamma            ,&
                                           omega           = omega               ,&
                                           RMF_B_amplitude = RMF_B_amplitude     ,&
                                           displacement_scale = EPS0,&
                                           x_cell_gpu      = field_gpu%x_cell_gpu,&
                                           y_cell_gpu      = field_gpu%y_cell_gpu,&
                                           z_cell_gpu      = field_gpu%z_cell_gpu,&
                                           q_gpu           = q_gpu)
   call correct_pec_external_ghosts(external_fields, field_gpu, time_next, -1._R8P, q_gpu)
   endassociate
   contains
      subroutine sub_external_fields_rmf_dev_kernel(ni,nj,nk,ngc,blocks_number,                               &
                                                    time_next,ef_alpha,ef_beta,ef_gamma,omega,RMF_B_amplitude,displacement_scale,&
                                                    x_cell_gpu,y_cell_gpu,z_cell_gpu,q_gpu)
      !< Subtract rotating magnetic field to the field, device kernel.
      integer(I4P), intent(in)    :: ni,nj,nk,ngc,blocks_number        !< Grids dimensions.
      real(R8P),    intent(in)    :: time_next                         !< Time at the next sub-step.
      real(R8P),    intent(in)    :: RMF_B_amplitude                   !< Rotating magnetic field amplitude.
      real(R8P),    intent(in)    :: displacement_scale                 !< Displacement to magnetic field scale.
	   integer(I4P), intent(in)    :: ef_alpha                          !< RMF rotation axis coordinate 1
	   integer(I4P), intent(in)    :: ef_beta                           !< RMF rotation axis coordinate 2
	   integer(I4P), intent(in)    :: ef_gamma                          !< RMF rotation axis coordinate 3
      real(R8P),    intent(in)    :: omega                             !< Omega frequency.
      real(R8P),    intent(in)    :: x_cell_gpu(1:,1-ngc:)             !< Cells x coordinates on GPU [nb,1-ngc:ni+ngc].
      real(R8P),    intent(in)    :: y_cell_gpu(1:,1-ngc:)             !< Cells y coordinates on GPU [nb,1-ngc:nj+ngc].
      real(R8P),    intent(in)    :: z_cell_gpu(1:,1-ngc:)             !< Cells z coordinates on GPU [nb,1-ngc:nk+ngc].
      real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field cell centered variables.
      real(R8P)                   :: B_r, B_theta                      !< Radial/azimuthal components of the rotating B field.
      real(R8P)                   :: cell_coord(3)                     !< Cell coordinates.
      real(R8P)                   :: phase                             !< Phase.
      real(R8P)                   :: theta                             !< Angle in cylindrical coordinates.
      real(R8P)                   :: x,y,r,c,s                         !< Buffer.
      integer(I4P)                :: i,j,k,b                           !< Counter.

      !!$acc& private(cell_coord)
      !$acc parallel loop independent gang vector collapse(4)                                               &
      !$acc& DEVICEVAR(x_cell_gpu,y_cell_gpu,z_cell_gpu,q_gpu)                                              &
      !$acc& firstprivate(ni,nj,nk,blocks_number,time_next,ef_alpha,ef_beta,ef_gamma,omega,RMF_B_amplitude,displacement_scale) &
      !$acc& private(B_r,B_theta,cell_coord,phase,theta,x,y,r,c,s)
      !$omp OMPLOOP collapse(4) &
      !$omp& DEVICEPTR(x_cell_gpu,y_cell_gpu,z_cell_gpu,q_gpu) &
      !$omp& firstprivate(ni,nj,nk,blocks_number,time_next,ef_alpha,ef_beta,ef_gamma,omega,RMF_B_amplitude,displacement_scale) &
      !$omp& private(B_r,B_theta,cell_coord,phase,theta,x,y,r,c,s)
      do b = 1, blocks_number
      do k = 1 - ngc, nk + ngc
      do j = 1 - ngc, nj + ngc
      do i = 1 - ngc, ni + ngc
         cell_coord = [x_cell_gpu(b,i), y_cell_gpu(b,j), z_cell_gpu(b,k)]
         x = cell_coord(ef_alpha)
         y = cell_coord(ef_beta)
         r = sqrt(x*x + y*y)
         theta = atan2(y, x)
         phase = omega*time_next - theta
         B_r     = RMF_B_amplitude*cos(phase)
         B_theta = RMF_B_amplitude*sin(phase)
         c = cos(theta)
         s = sin(theta)
         q_gpu(b,i,j,k,ef_alpha+3) = q_gpu(b,i,j,k,ef_alpha+3) - (B_r*c - B_theta*s)
         q_gpu(b,i,j,k,ef_beta +3) = q_gpu(b,i,j,k,ef_beta +3) - (B_r*s + B_theta*c)
         q_gpu(b,i,j,k,ef_gamma  ) = q_gpu(b,i,j,k,ef_gamma  ) - r*omega*RMF_B_amplitude*cos(phase)*displacement_scale
      enddo
      enddo
      enddo
      enddo
      endsubroutine sub_external_fields_rmf_dev_kernel
   endsubroutine sub_external_fields_rmf_dev

   subroutine add_external_fields_uniform_dev(external_fields, field_gpu, dt, time, q_gpu, gamm)
   !< Add uniform external electric displacement and magnetic field to the field.
   type(prism_external_fields_object), intent(in)    :: external_fields !< External fields handler.
   type(field_fnl_object),             intent(in)    :: field_gpu       !< Field.
   real(R8P),                          intent(in)    :: dt              !< Time step.
   real(R8P),                          intent(in)    :: time            !< Current time.
   real(R8P),                          intent(inout) :: q_gpu(1:,              &
                                                              1-field_gpu%ngc:,&
                                                              1-field_gpu%ngc:,&
                                                              1-field_gpu%ngc:,&
                                                              1:)       !< Conservative variables.
   real(R8P), optional,                intent(in)    :: gamm            !< Gamma values of RK.

   associate(blocks_number=>field_gpu%blocks_number, ni=>field_gpu%ni, nj=>field_gpu%nj, nk=>field_gpu%nk, ngc=>field_gpu%ngc, &
             axis=>external_fields%uniform_axis, Uniform_D_amplitude=>external_fields%Uniform_D_amplitude,                    &
             Uniform_B_amplitude=>external_fields%Uniform_B_amplitude)
   if (present(gamm)) continue
   associate(dt_unused=>dt, time_unused=>time)
   endassociate
   call add_external_fields_uniform_dev_kernel(ni                  = ni,                  &
                                               nj                  = nj,                  &
                                               nk                  = nk,                  &
                                               ngc                 = ngc,                 &
                                               blocks_number       = blocks_number,       &
                                               axis                = axis,                &
                                               Uniform_D_amplitude = Uniform_D_amplitude, &
                                               Uniform_B_amplitude = Uniform_B_amplitude, &
                                               q_gpu               = q_gpu)
   call correct_pec_external_ghosts(external_fields, field_gpu, time, 1._R8P, q_gpu)
   endassociate
   contains
      subroutine add_external_fields_uniform_dev_kernel(ni,nj,nk,ngc,blocks_number,axis, &
                                                        Uniform_D_amplitude,Uniform_B_amplitude,q_gpu)
      !< Add uniform external field to the field, device kernel.
      integer(I4P), intent(in)    :: ni,nj,nk,ngc,blocks_number        !< Grids dimensions.
      integer(I4P), intent(in)    :: axis                              !< Uniform field direction index.
      real(R8P),    intent(in)    :: Uniform_D_amplitude               !< Uniform electric displacement amplitude.
      real(R8P),    intent(in)    :: Uniform_B_amplitude               !< Uniform magnetic field amplitude.
      real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field cell centered variables.
      integer(I4P)                :: i,j,k,b                           !< Counter.

      !$acc parallel loop independent gang vector collapse(4)                                      &
      !$acc& DEVICEVAR(q_gpu)                                                                       &
      !$acc& firstprivate(ni,nj,nk,blocks_number,axis,Uniform_D_amplitude,Uniform_B_amplitude)
      !$omp OMPLOOP collapse(4) &
      !$omp& DEVICEPTR(q_gpu) &
      !$omp& firstprivate(ni,nj,nk,blocks_number,axis,Uniform_D_amplitude,Uniform_B_amplitude)
      do b = 1, blocks_number
      do k = 1 - ngc, nk + ngc
      do j = 1 - ngc, nj + ngc
      do i = 1 - ngc, ni + ngc
         q_gpu(b,i,j,k,axis  ) = q_gpu(b,i,j,k,axis  ) + Uniform_D_amplitude
         q_gpu(b,i,j,k,axis+3) = q_gpu(b,i,j,k,axis+3) + Uniform_B_amplitude
      enddo
      enddo
      enddo
      enddo
      endsubroutine add_external_fields_uniform_dev_kernel
   endsubroutine add_external_fields_uniform_dev

   subroutine sub_external_fields_uniform_dev(external_fields, field_gpu, dt, time, q_gpu, gamm)
   !< Subtract uniform external electric displacement and magnetic field from the field.
   type(prism_external_fields_object), intent(in)    :: external_fields !< External fields handler.
   type(field_fnl_object),             intent(in)    :: field_gpu       !< Field.
   real(R8P),                          intent(in)    :: dt              !< Time step.
   real(R8P),                          intent(in)    :: time            !< Current time.
   real(R8P),                          intent(inout) :: q_gpu(1:,              &
                                                              1-field_gpu%ngc:,&
                                                              1-field_gpu%ngc:,&
                                                              1-field_gpu%ngc:,&
                                                              1:)       !< Conservative variables.
   real(R8P), optional,                intent(in)    :: gamm            !< Gamma values of RK.

   associate(blocks_number=>field_gpu%blocks_number, ni=>field_gpu%ni, nj=>field_gpu%nj, nk=>field_gpu%nk, ngc=>field_gpu%ngc, &
             axis=>external_fields%uniform_axis, Uniform_D_amplitude=>external_fields%Uniform_D_amplitude,                    &
             Uniform_B_amplitude=>external_fields%Uniform_B_amplitude)
   if (present(gamm)) continue
   associate(dt_unused=>dt, time_unused=>time)
   endassociate
   call sub_external_fields_uniform_dev_kernel(ni                  = ni,                  &
                                               nj                  = nj,                  &
                                               nk                  = nk,                  &
                                               ngc                 = ngc,                 &
                                               blocks_number       = blocks_number,       &
                                               axis                = axis,                &
                                               Uniform_D_amplitude = Uniform_D_amplitude, &
                                               Uniform_B_amplitude = Uniform_B_amplitude, &
                                               q_gpu               = q_gpu)
   call correct_pec_external_ghosts(external_fields, field_gpu, time, -1._R8P, q_gpu)
   endassociate
   contains
      subroutine sub_external_fields_uniform_dev_kernel(ni,nj,nk,ngc,blocks_number,axis, &
                                                        Uniform_D_amplitude,Uniform_B_amplitude,q_gpu)
      !< Subtract uniform external field from the field, device kernel.
      integer(I4P), intent(in)    :: ni,nj,nk,ngc,blocks_number        !< Grids dimensions.
      integer(I4P), intent(in)    :: axis                              !< Uniform field direction index.
      real(R8P),    intent(in)    :: Uniform_D_amplitude               !< Uniform electric displacement amplitude.
      real(R8P),    intent(in)    :: Uniform_B_amplitude               !< Uniform magnetic field amplitude.
      real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field cell centered variables.
      integer(I4P)                :: i,j,k,b                           !< Counter.

      !$acc parallel loop independent gang vector collapse(4)                                      &
      !$acc& DEVICEVAR(q_gpu)                                                                       &
      !$acc& firstprivate(ni,nj,nk,blocks_number,axis,Uniform_D_amplitude,Uniform_B_amplitude)
      !$omp OMPLOOP collapse(4) &
      !$omp& DEVICEPTR(q_gpu) &
      !$omp& firstprivate(ni,nj,nk,blocks_number,axis,Uniform_D_amplitude,Uniform_B_amplitude)
      do b = 1, blocks_number
      do k = 1 - ngc, nk + ngc
      do j = 1 - ngc, nj + ngc
      do i = 1 - ngc, ni + ngc
         q_gpu(b,i,j,k,axis  ) = q_gpu(b,i,j,k,axis  ) - Uniform_D_amplitude
         q_gpu(b,i,j,k,axis+3) = q_gpu(b,i,j,k,axis+3) - Uniform_B_amplitude
      enddo
      enddo
      enddo
      enddo
      endsubroutine sub_external_fields_uniform_dev_kernel
   endsubroutine sub_external_fields_uniform_dev

   subroutine correct_pec_external_ghosts(external_fields, field_gpu, time_stage, factor, q_gpu)
   !< Replace the raw external contribution on physical PEC ghosts by its mirrored value.
   type(prism_external_fields_object), intent(in)    :: external_fields
   type(field_fnl_object),             intent(in)    :: field_gpu
   real(R8P),                          intent(in)    :: time_stage, factor
   real(R8P),                          intent(inout) :: q_gpu(1:,1-field_gpu%ngc:,1-field_gpu%ngc:,1-field_gpu%ngc:,1:)
   integer(I4P)                                      :: b, faces(6), field_kind

   if (.not.allocated(external_fields%pec_faces)) return
   field_kind = 1_I4P
   if (external_fields%ef_type == EF_TYPE_UNIFORM_FIELD) field_kind = 2_I4P
   do b=1, field_gpu%blocks_number
      faces = external_fields%pec_faces(:,b)
      if (all(faces == 0_I4P)) cycle
      call correct_block(b, faces)
   enddo
   contains
      subroutine correct_block(b, faces)
      integer(I4P), intent(in) :: b, faces(6)
      integer(I4P)             :: i,j,k,axis,m,idx(3),sample(3)
      real(R8P)                :: sign_D(3),sign_B(3),value(6,2),coord(3),x,y,r,theta,phase,omega,br,bt,c,s
      real(R8P)                :: d_amp,b_amp
      real(R8P)                :: rmf_amp
      real(R8P), pointer       :: x_cell_gpu(:,:), y_cell_gpu(:,:), z_cell_gpu(:,:)
      integer(I4P)             :: alpha,beta,gamma,uniform_axis,ni,nj,nk,ngc,naxis(3)
      x_cell_gpu => field_gpu%x_cell_gpu
      y_cell_gpu => field_gpu%y_cell_gpu
      z_cell_gpu => field_gpu%z_cell_gpu
      ni=field_gpu%ni; nj=field_gpu%nj; nk=field_gpu%nk; ngc=field_gpu%ngc
      naxis=[ni,nj,nk]
      alpha = external_fields%alpha
      beta = external_fields%beta
      gamma = external_fields%gamm
      uniform_axis = external_fields%uniform_axis
      omega = 2._R8P*PI*external_fields%RMF_frequency
      d_amp = external_fields%Uniform_D_amplitude
      b_amp = external_fields%Uniform_B_amplitude
      rmf_amp = external_fields%RMF_B_amplitude
      !$acc parallel loop independent gang vector collapse(3) &
      !$acc& DEVICEVAR(q_gpu,x_cell_gpu,y_cell_gpu,z_cell_gpu) &
      !$acc& firstprivate(b,faces,field_kind,factor,time_stage,alpha,beta,gamma,uniform_axis,omega,d_amp,b_amp,rmf_amp,naxis) &
      !$acc& private(idx,sample,sign_D,sign_B,value,coord,x,y,r,theta,phase,br,bt,c,s,axis,m)
      !$omp OMPLOOP collapse(3) &
      !$omp& DEVICEPTR(q_gpu,x_cell_gpu,y_cell_gpu,z_cell_gpu) &
      !$omp& firstprivate(b,faces,field_kind,factor,time_stage,alpha,beta,gamma,uniform_axis,omega,d_amp,b_amp,rmf_amp,naxis) &
      !$omp& private(idx,sample,sign_D,sign_B,value,coord,x,y,r,theta,phase,br,bt,c,s,axis,m)
      do k=1-ngc,nk+ngc
      do j=1-ngc,nj+ngc
      do i=1-ngc,ni+ngc
         idx = [i,j,k]
         sample = idx
         sign_D = 1._R8P
         sign_B = 1._R8P
         do axis=1,3
            if (idx(axis)<1 .and. faces(2*axis-1)/=0) then
               sample(axis)=1-idx(axis)
               sign_D=-sign_D; sign_D(axis)=-sign_D(axis); sign_B(axis)=-sign_B(axis)
            elseif (idx(axis)>naxis(axis) .and. faces(2*axis)/=0) then
               sample(axis)=2*naxis(axis)+1-idx(axis)
               sign_D=-sign_D; sign_D(axis)=-sign_D(axis); sign_B(axis)=-sign_B(axis)
            endif
         enddo
         if (sample(1)==i .and. sample(2)==j .and. sample(3)==k) cycle
         value=0._R8P
         if (field_kind==2_I4P) then
            value(uniform_axis,1)=d_amp
            value(uniform_axis+3,1)=b_amp
            value(:,2)=value(:,1)
         else
            do m=1,2
               if (m==1) then
                  coord=[x_cell_gpu(b,i),y_cell_gpu(b,j),z_cell_gpu(b,k)]
               else
                  coord=[x_cell_gpu(b,sample(1)),y_cell_gpu(b,sample(2)),z_cell_gpu(b,sample(3))]
               endif
               x=coord(alpha); y=coord(beta); r=sqrt(x*x+y*y)
               theta=atan2(y,x); phase=omega*time_stage-theta
               br=rmf_amp*cos(phase)
               bt=rmf_amp*sin(phase)
               c=cos(theta); s=sin(theta)
               value(alpha+3,m)=br*c-bt*s
               value(beta+3,m)=br*s+bt*c
               value(gamma,m)=r*omega*rmf_amp*cos(phase)*EPS0
            enddo
         endif
         do m=1,3
            q_gpu(b,i,j,k,m)=q_gpu(b,i,j,k,m)+factor*(sign_D(m)*value(m,2)-value(m,1))
            q_gpu(b,i,j,k,m+3)=q_gpu(b,i,j,k,m+3)+factor*(sign_B(m)*value(m+3,2)-value(m+3,1))
         enddo
      enddo
      enddo
      enddo
      endsubroutine correct_block
   endsubroutine correct_pec_external_ghosts

endmodule adam_prism_fnl_external_fields_kernels
