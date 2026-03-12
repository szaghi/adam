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
                                           x_cell_gpu      = field_gpu%x_cell_gpu,&
                                           y_cell_gpu      = field_gpu%y_cell_gpu,&
                                           z_cell_gpu      = field_gpu%z_cell_gpu,&
                                           q_gpu           = q_gpu)
   endassociate
   contains
      subroutine add_external_fields_rmf_dev_kernel(ni,nj,nk,ngc,blocks_number,                               &
                                                    time_next,ef_alpha,ef_beta,ef_gamma,omega,RMF_B_amplitude,&
                                                    x_cell_gpu,y_cell_gpu,z_cell_gpu,q_gpu)
      !< Add rotating magnetic field to the field, device kernel.
      integer(I4P), intent(in)    :: ni,nj,nk,ngc,blocks_number        !< Grids dimensions.
      real(R8P),    intent(in)    :: time_next                         !< Time at the next sub-step.
      real(R8P),    intent(in)    :: RMF_B_amplitude                   !< Rotating magnetic field amplitude.
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

      !$acc parallel loop independent gang vector collapse(4)  &
      !$acc& DEVICEVAR(x_cell_gpu,y_cell_gpu,z_cell_gpu,q_gpu) &
      !$acc& firstprivate(time_next,ef_alpha,ef_beta,ef_gamma,omega,RMF_B_amplitude)
      do b = 1, blocks_number
      do k = 1, nk
      do j = 1, nj
      do i = 1, ni
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
         q_gpu(b,i,j,k,ef_gamma  ) = q_gpu(b,i,j,k,ef_gamma  ) + r*omega*RMF_B_amplitude*cos(phase)*EPS0
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
                                           x_cell_gpu      = field_gpu%x_cell_gpu,&
                                           y_cell_gpu      = field_gpu%y_cell_gpu,&
                                           z_cell_gpu      = field_gpu%z_cell_gpu,&
                                           q_gpu           = q_gpu)
   endassociate
   contains
      subroutine sub_external_fields_rmf_dev_kernel(ni,nj,nk,ngc,blocks_number,                               &
                                                    time_next,ef_alpha,ef_beta,ef_gamma,omega,RMF_B_amplitude,&
                                                    x_cell_gpu,y_cell_gpu,z_cell_gpu,q_gpu)
      !< Subtract rotating magnetic field to the field, device kernel.
      integer(I4P), intent(in)    :: ni,nj,nk,ngc,blocks_number        !< Grids dimensions.
      real(R8P),    intent(in)    :: time_next                         !< Time at the next sub-step.
      real(R8P),    intent(in)    :: RMF_B_amplitude                   !< Rotating magnetic field amplitude.
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

      !$acc parallel loop independent gang vector collapse(4)  &
      !$acc& DEVICEVAR(x_cell_gpu,y_cell_gpu,z_cell_gpu,q_gpu) &
      !$acc& firstprivate(time_next,ef_alpha,ef_beta,ef_gamma,omega,RMF_B_amplitude)
      do b = 1, blocks_number
      do k = 1, nk
      do j = 1, nj
      do i = 1, ni
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
         q_gpu(b,i,j,k,ef_alpha+3) = q_gpu(b,i,j,k,ef_alpha+3) - B_r*c - B_theta*s
         q_gpu(b,i,j,k,ef_beta +3) = q_gpu(b,i,j,k,ef_beta +3) - B_r*s + B_theta*c
         q_gpu(b,i,j,k,ef_gamma  ) = q_gpu(b,i,j,k,ef_gamma  ) - r*omega*RMF_B_amplitude*cos(phase)*EPS0
      enddo
      enddo
      enddo
      enddo
      endsubroutine sub_external_fields_rmf_dev_kernel
   endsubroutine sub_external_fields_rmf_dev
endmodule adam_prism_fnl_external_fields_kernels
