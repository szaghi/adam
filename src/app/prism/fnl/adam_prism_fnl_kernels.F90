!< ADAM, PRISM FNL application kernels.

#include "fundal.H"

module adam_prism_fnl_kernels
!< ADAM, PRISM FNL application kernels.

use adam_prism_parameters
use adam_weno_fnl_kernels
! use adam_prism_fnl_cns_kernels
use fundal
use penf, only : I4P, I8P, R8P

implicit none
private
public :: compute_coils_current_dev
public :: compute_div_d_b_dev
public :: compute_fluxes_convective_dev
public :: compute_fluxes_difference_dev
public :: compute_dxyz_min_dev
public :: set_bc_q_gpu_dev
public :: set_sir_dev

contains
   ! public procedures
   subroutine compute_coils_current_dev(ni,nj,nk,ngc,blocks_number,time,A_gpu,d_gpu,f_gpu,phase_gpu,coil_flag_gpu,&
                                        td,J_vec_gpu,dx_gpu,q_gpu)
   !< Compute current coils sources.
   integer(I4P), intent(in)    :: ni                                     !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                     !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                     !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                    !< Ghost cells number.
   integer(I4P), intent(in)    :: blocks_number                          !< Number of blocks.
   integer(I4P), intent(in)    :: coil_flag_gpu(1:,1-ngc:,1-ngc:,1-ngc:) !< Matrice contenente informazioni su quale spira passa.
   real(R8P),    intent(in)    :: time                                   !< Simulation time, to compute current value if AC
   real(R8P),    intent(in)    :: A_gpu(1:)                              !< Current amplitude (A)
   real(R8P),    intent(in)    :: f_gpu(1:)                              !< Current frequency, if AC (Hz)
   real(R8P),    intent(in)    :: phase_gpu(1:)                          !< Current initial phase, if AC
   real(R8P),    intent(in)    :: d_gpu(1:)                              !< Wire diameter
   real(R8P),    intent(in)    :: td                                     !< Delay di accensione della spira
   real(R8P),    intent(in)    :: J_vec_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Matrice versori di corrente delle spire nelle celle
   real(R8P),    intent(in)    :: dx_gpu                                 !< Space step in x direction (m)
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Field variables.
   real(R8P)                   :: current_density                        !< Current density
   real(R8P)                   :: g                                      !< Polinomio caratteristico transitorio accensione spira
   integer(I4P)                :: coil_id                                !< ID per identificare spira
   integer(I4P)                :: i,j,k,b,n                              !< Counter

   g = 10._R8P*(time/td)**3 - 15._R8P*(time/td)**4 + 6._R8P*(time/td)**5
   !$acc parallel loop independent gang vector collapse(4)                           &
   !$acc DEVICEVAR(coil_flag_gpu,A_gpu,f_gpu,phase_gpu,d_gpu,j_vec_gpu,q_gpu,dx_gpu) &
   !$acc firstprivate(g,time,td)
   !$omp OMPLOOP DEVICEVAR(coil_flag_gpu,A_gpu,f_gpu,phase_gpu,d_gpu,j_vec_gpu,q_gpu,dx_gpu) firstprivate(g,time,td)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               coil_id = coil_flag_gpu(b,i,j,k)
               if (coil_id /= 0_I4P) then
                  !Per DC frequenza e fase sono nulle, quindi se uso la funzione coseno
                  !mi rispramio anche il selectcase

                  !Densità di corrente al tempo t della spira n-esima identificata da (coil_id)
                  !current_density = 4*A(coil_id)/(pi*d(coil_id)**2)*cos(2*pi*f(coil_id)*time + phase(coil_id)*pi/180.0_R8P)

                  !Modifico calcolo densità di corrente considerando sezione quadrata, per coerenza con calcolo Filippo
                  !E aggiungo transitorio di corrente
                  if (time < td) then
                     current_density = g*A_gpu(coil_id)/((d_gpu(coil_id)-dx_gpu)**2)*cos(phase_gpu(coil_id)*PI/180.0_R8P)
                  else
                     current_density = A_gpu(coil_id)/((d_gpu(coil_id)-dx_gpu)**2)*cos(2*PI*f_gpu(coil_id)*(time-td) + &
                     phase_gpu(coil_id)*PI/180.0_R8P)
                  endif

                  do n=7, 9
                     q_gpu(b,i,j,k,n) = current_density* J_vec_gpu(b,i,j,k,n-6)
                  enddo

                  !if (sq_norm(q(7:9,i,j,k,b)) == 0._R8P) then
                     !Se la densità di corrente è nulla non faccio nulla
                  !   q(7:9,i,j,k,b) = current_density*q(7:9,i,j,k,b)
                  !else
                     !Se la densità di corrente è diversa da zero allora rinormalizzo il vettore corrente
                     !e lo moltiplico per la densità di corrente
                  !   q(7,i,j,k,b) = q(7,i,j,k,b)/(sq_norm(q(7:9,i,j,k,b)))**0.5 !lo devo rinormalizzare ogni volta
                  !   q(8,i,j,k,b) = q(8,i,j,k,b)/(sq_norm(q(7:9,i,j,k,b)))**0.5
                  !   q(9,i,j,k,b) = q(9,i,j,k,b)/(sq_norm(q(7:9,i,j,k,b)))**0.5
                  !   q(7:9,i,j,k,b) = current_density*q(7:9,i,j,k,b)
                  !endif
               endif
            enddo
         enddo
      enddo
   enddo
   endsubroutine compute_coils_current_dev

   subroutine compute_div_d_b_dev(ni, nj, nk, ngc, blocks_number, dx_gpu, q_gpu, div_gpu)
   !< Compute div(D), div(B).
   integer(I4P), intent(in)    :: ni                                  !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                  !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                  !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in)    :: blocks_number                       !< Number of blocks.
   real(R8P),    intent(in)    :: dx_gpu                              !< Space step in x direction (m)
   real(R8P),    intent(in)    ::   q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field variables.
   real(R8P),    intent(inout) :: div_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Divergence of D, B.
   integer(I4P)                :: i,j,k,b                             !< Counter

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(div_gpu,q_gpu,dx_gpu)
   !$omp OMPLOOP DEVICEVAR(div_gpu,q_gpu,dx_gpu)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               div_gpu(b,i,j,k,1) = 0.5_R8P*((q_gpu(b,i+1,j,k,1) - q_gpu(b,i-1,j,k,1))/dx_gpu + &
                                             (q_gpu(b,i,j+1,k,2) - q_gpu(b,i,j-1,k,2))/dx_gpu + &
                                             (q_gpu(b,i,j,k+1,3) - q_gpu(b,i,j,k-1,3))/dx_gpu)
               div_gpu(b,i,j,k,2) = 0.5_R8P*((q_gpu(b,i+1,j,k,4) - q_gpu(b,i-1,j,k,4))/dx_gpu + &
                                             (q_gpu(b,i,j+1,k,5) - q_gpu(b,i,j-1,k,5))/dx_gpu + &
                                             (q_gpu(b,i,j,k+1,6) - q_gpu(b,i,j,k-1,6))/dx_gpu)

            enddo
         enddo
      enddo
   enddo
   endsubroutine compute_div_d_b_dev

   subroutine compute_fluxes_convective_dev(dir,blocks_number,ni,nj,nk,ngc,nv,evmax,si_gpu,sir_gpu,q_gpu,fluxes_gpu)
   !< Compute convective fluxes along x direction.
   !< @NOTE Be carefull with `si` and `sir` variables: are not present on GPU, probably a mapping
   !< is done by compiler, check it!
   integer(I4P), intent(in)    :: dir                                    !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)    :: blocks_number                          !< Number of blocks.
   integer(I4P), intent(in)    :: ni                                     !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                     !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                     !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                    !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                                     !< Number of conservative varibales.
   real(R8P),    intent(in)    :: evmax                                  !< Maximum waves speeds estimation.
   integer(I4P), intent(in)    :: si_gpu(3)                              !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)    :: sir_gpu(3)                             !< Directional (1=x,2=y,3=z) increment, real cast.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Fields variables.
   real(R8P),    intent(inout) :: fluxes_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Fluxes.
   integer(I4P)                :: b, i, j, k                             !< Counter.

   select case(dir)
   case(1)
      !$acc parallel loop independent gang vector collapse(3) &
      !$acc default(none)                                     &
      !$acc DEVICEVAR(si_gpu,sir_gpu,q_gpu,fluxes_gpu)        &
      !$acc firstprivate(dir,blocks_number,ngc,nv,evmax)      &
      !$acc private(b,i,j,k)
      do b=1, blocks_number
      do j=1, nj
      do k=1, nk
      do i=0, ni
         call compute_fluxes_convective_ri_dev(dir=dir,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv, &
                                               evmax=evmax,si_gpu=si_gpu,sir_gpu=sir_gpu,q_gpu=q_gpu,fluxes_gpu=fluxes_gpu)
      enddo
      enddo
      enddo
      enddo
   case(2)
      !$acc parallel loop independent gang vector collapse(3) &
      !$acc default(none)                                     &
      !$acc DEVICEVAR(si_gpu,sir_gpu,q_gpu,fluxes_gpu)        &
      !$acc firstprivate(dir,blocks_number,ngc,nv,evmax)      &
      !$acc private(b,i,j,k)
      do b=1, blocks_number
      do i=1, ni
      do k=1, nk
      do j=0, nj
         call compute_fluxes_convective_ri_dev(dir=dir,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv, &
                                               evmax=evmax,si_gpu=si_gpu,sir_gpu=sir_gpu,q_gpu=q_gpu,fluxes_gpu=fluxes_gpu)
      enddo
      enddo
      enddo
      enddo
   case(3)
      !$acc parallel loop independent gang vector collapse(3) &
      !$acc default(none)                                     &
      !$acc DEVICEVAR(si_gpu,sir_gpu,q_gpu,fluxes_gpu)        &
      !$acc firstprivate(dir,blocks_number,ngc,nv,evmax)      &
      !$acc private(b,i,j,k)
      do b=1, blocks_number
      do i=1, ni
      do j=1, nj
      do k=0, nk
         call compute_fluxes_convective_ri_dev(dir=dir,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv, &
                                               evmax=evmax,si_gpu=si_gpu,sir_gpu=sir_gpu,q_gpu=q_gpu,fluxes_gpu=fluxes_gpu)
      enddo
      enddo
      enddo
      enddo
   endselect
   endsubroutine compute_fluxes_convective_dev

   subroutine compute_fluxes_difference_dev(null_x, null_y, null_z, blocks_number, ni, nj, nk, ngc, nv, ib_eps, &
                                            dx_gpu, dy_gpu, dz_gpu, flx_gpu, fly_gpu, flz_gpu, phi_gpu,         &
                                            q_gpu, eta, chi, d_divergence_cleaner, b_divergence_cleaner, dq_gpu)
   !< Compute fluxes difference.
   logical,      intent(in)           :: null_x, null_y, null_z              !< Nullified directions tags.
   integer(I4P), intent(in)           :: blocks_number                       !< Number of blocks.
   integer(I4P), intent(in)           :: ni                                  !< Grid cells number in I direction.
   integer(I4P), intent(in)           :: nj                                  !< Grid cells number in J direction.
   integer(I4P), intent(in)           :: nk                                  !< Grid cells number in K direction.
   integer(I4P), intent(in)           :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in)           :: nv                                  !< Number of conservative varibales.
   real(R8P),    intent(in)           :: ib_eps                              !< Tolerance IB delta ratio.
   real(R8P),    intent(in)           :: dx_gpu(1:), dy_gpu(1:), dz_gpu(1:)  !< Space steps.
   real(R8P),    intent(in)           :: flx_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< X direction fluxes.
   real(R8P),    intent(in)           :: fly_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Y direction fluxes.
   real(R8P),    intent(in)           :: flz_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Z direction fluxes.
   real(R8P),    intent(in), optional :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< IB distance function.
   real(R8P),    intent(in)           :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Fields.
   real(R8P),    intent(in)           :: eta, chi                            !< Coefficiente modello correzione divergenza campi.
   logical,      intent(in)           :: D_divergence_cleaner                !< Flag to perform electric field divergence cleaning.
   logical,      intent(in)           :: B_divergence_cleaner                !< Flag to perform magnetic field divergence cleaning.
   real(R8P),    intent(inout)        :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Fluxes differences.
   real(R8P)                          :: delta_x, delta_y, delta_z           !< Space steps.
   real(R8P)                          :: dx_locale, dy_locale, dz_locale     !< Local space steps.
   integer(I4P)                       :: b, i, j, k, v                       !< Counter.
   integer(I4P)                       :: all_solids                          !< Last phi index, all solids summary.
   real(R8P)                          :: qmx, qmy, qmz                       !< Momentum nullification scalar.

   qmx = 1._R8P ; if (null_x) qmx = 0._R8P
   qmy = 1._R8P ; if (null_y) qmy = 0._R8P
   qmz = 1._R8P ; if (null_z) qmz = 0._R8P
   if (present(phi_gpu)) then
      all_solids = ubound(phi_gpu, dim=5)
      !$acc parallel loop independent gang vector collapse(4) &
      !$acc DEVICEVAR(dx_gpu, dy_gpu, dz_gpu, flx_gpu, fly_gpu, flz_gpu, phi_gpu, dq_gpu, q_gpu)
      !$omp OMPLOOP DEVICEVAR(dx_gpu, dy_gpu, dz_gpu, flx_gpu, fly_gpu, flz_gpu, phi_gpu, dq_gpu, q_gpu)
      do b=1,blocks_number
      do k=1,nk
      do j=1,nj
      do i=1,ni
         dx_locale = dx_gpu(b)
         if (phi_gpu(b,i,j,k,all_solids)<0.) then
            if (phi_gpu(b,i+1,j,k,all_solids)*phi_gpu(b,i-1,j,k,all_solids)<0) then
               if (phi_gpu(b,i+1,j,k,all_solids)>0.) then
                  delta_x= -phi_gpu(b,i,j,k,all_solids)/(phi_gpu(b,i+1,j,k,all_solids)-phi_gpu(b,i,j,k,all_solids)+ib_eps)*dx_gpu(b)
                  dx_locale = dx_gpu(b)/2 + delta_x
               else
                  delta_x= -phi_gpu(b,i,j,k,all_solids)/(phi_gpu(b,i-1,j,k,all_solids)-phi_gpu(b,i,j,k,all_solids)+ib_eps)*dx_gpu(b)
                  dx_locale = dx_gpu(b)/2 + delta_x
               endif
            endif
         endif
         dy_locale = dy_gpu(b)
         if (phi_gpu(b,i,j,k,all_solids)<0.) then
            if (phi_gpu(b,i,j+1,k,all_solids)*phi_gpu(b,i,j-1,k,all_solids)<0) then
               if (phi_gpu(b,i,j+1,k,all_solids)>0.) then
                  delta_y= -phi_gpu(b,i,j,k,all_solids)/(phi_gpu(b,i,j+1,k,all_solids)-phi_gpu(b,i,j,k,all_solids)+ib_eps)*dy_gpu(b)
                  dy_locale = dy_gpu(b)/2 + delta_y
               else
                  delta_y= -phi_gpu(b,i,j,k,all_solids)/(phi_gpu(b,i,j-1,k,all_solids)-phi_gpu(b,i,j,k,all_solids)+ib_eps)*dy_gpu(b)
                  dy_locale = dy_gpu(b)/2 + delta_y
               endif
            endif
         endif
         dz_locale = dz_gpu(b)
         if (phi_gpu(b,i,j,k,all_solids)<0.) then
            if (phi_gpu(b,i,j,k+1,all_solids)*phi_gpu(b,i,j,k-1,all_solids)<0) then
               if (phi_gpu(b,i,j,k+1,all_solids)>0.) then
                  delta_z= -phi_gpu(b,i,j,k,all_solids)/(phi_gpu(b,i,j,k+1,all_solids)-phi_gpu(b,i,j,k,all_solids)+ib_eps)*dz_gpu(b)
                  dz_locale = dz_gpu(b)/2 + delta_z
               else
                  delta_z= -phi_gpu(b,i,j,k,all_solids)/(phi_gpu(b,i,j,k-1,all_solids)-phi_gpu(b,i,j,k,all_solids)+ib_eps)*dz_gpu(b)
                  dz_locale = dz_gpu(b)/2 + delta_z
               endif
            endif
         endif
         do v=1, nv
            dq_gpu(b,i,j,k,v) = - (flx_gpu(b,i,j,k,v)-flx_gpu(b,i-1,j,k,v))/dx_locale &
                                - (fly_gpu(b,i,j,k,v)-fly_gpu(b,i,j-1,k,v))/dy_locale &
                                - (flz_gpu(b,i,j,k,v)-flz_gpu(b,i,j,k-1,v))/dz_locale
         enddo
         dq_gpu(b,i,j,k,2) = dq_gpu(b,i,j,k,2) * qmx
         dq_gpu(b,i,j,k,3) = dq_gpu(b,i,j,k,3) * qmy
         dq_gpu(b,i,j,k,4) = dq_gpu(b,i,j,k,4) * qmz

         !Completo calcolo aggiungendo termini sorgenti legato alle correnti delle spire (per ora)
         !E ad una eventuale correzione parabolica della divergenza (parametro eta diverso da zero)

         dq_gpu(b,i,j,k,1) = dq_gpu(b,i,j,k,1) - q_gpu(b,i,j,k,7)
         dq_gpu(b,i,j,k,2) = dq_gpu(b,i,j,k,2) - q_gpu(b,i,j,k,8)
         dq_gpu(b,i,j,k,3) = dq_gpu(b,i,j,k,3) - q_gpu(b,i,j,k,9)

         if (d_divergence_cleaner .and. .not.b_divergence_cleaner .and. eta>0._R8P) then
            dq_gpu(b,i,j,k,10) = dq_gpu(b,i,j,k,10) - chi/eta*chi/eta*q_gpu(b,i,j,k,10)
         elseif (D_divergence_cleaner .and. B_divergence_cleaner .and. eta>0._R8P) then
            dq_gpu(b,i,j,k,10) = dq_gpu(b,i,j,k,10) - chi/eta*chi/eta*q_gpu(b,i,j,k,10)
            dq_gpu(b,i,j,k,11) = dq_gpu(b,i,j,k,11) - chi/eta*chi/eta*q_gpu(b,i,j,k,11)
         endif
      enddo
      enddo
      enddo
      enddo
   else
      !$acc parallel loop independent gang vector collapse(4) &
      !$acc DEVICEVAR(dx_gpu, dy_gpu, dz_gpu, flx_gpu, fly_gpu, flz_gpu, phi_gpu, dq_gpu, q_gpu)
      !$omp OMPLOOP DEVICEVAR(dx_gpu, dy_gpu, dz_gpu, flx_gpu, fly_gpu, flz_gpu, phi_gpu, dq_gpu, q_gpu)
      do b=1,blocks_number
      do k=1,nk
      do j=1,nj
      do i=1,ni
         do v=1, nv
            dq_gpu(b,i,j,k,v) = - (flx_gpu(b,i,j,k,v)-flx_gpu(b,i-1,j,k,v))/dx_gpu(b) &
                                - (fly_gpu(b,i,j,k,v)-fly_gpu(b,i,j-1,k,v))/dy_gpu(b) &
                                - (flz_gpu(b,i,j,k,v)-flz_gpu(b,i,j,k-1,v))/dz_gpu(b)
         enddo
         dq_gpu(b,i,j,k,2) = dq_gpu(b,i,j,k,2) * qmx
         dq_gpu(b,i,j,k,3) = dq_gpu(b,i,j,k,3) * qmy
         dq_gpu(b,i,j,k,4) = dq_gpu(b,i,j,k,4) * qmz

         !Completo calcolo aggiungendo termini sorgenti legato alle correnti delle spire (per ora)
         !E ad una eventuale correzione parabolica della divergenza (parametro eta diverso da zero)

         dq_gpu(b,i,j,k,1) = dq_gpu(b,i,j,k,1) - q_gpu(b,i,j,k,7)
         dq_gpu(b,i,j,k,2) = dq_gpu(b,i,j,k,2) - q_gpu(b,i,j,k,8)
         dq_gpu(b,i,j,k,3) = dq_gpu(b,i,j,k,3) - q_gpu(b,i,j,k,9)

         if (d_divergence_cleaner .and. .not.b_divergence_cleaner .and. eta>0._R8P) then
            dq_gpu(b,i,j,k,10) = dq_gpu(b,i,j,k,10) - chi/eta*chi/eta*q_gpu(b,i,j,k,10)
         elseif (D_divergence_cleaner .and. B_divergence_cleaner .and. eta>0._R8P) then
            dq_gpu(b,i,j,k,10) = dq_gpu(b,i,j,k,10) - chi/eta*chi/eta*q_gpu(b,i,j,k,10)
            dq_gpu(b,i,j,k,11) = dq_gpu(b,i,j,k,11) - chi/eta*chi/eta*q_gpu(b,i,j,k,11)
         endif
      enddo
      enddo
      enddo
      enddo
   endif
   endsubroutine compute_fluxes_difference_dev

   subroutine compute_dxyz_min_dev(blocks_number, dxyz_gpu, dxyz_min)
   !< Compute minimum dxyz space step.
   integer(I4P), intent(in)  :: blocks_number !< Number of blocks.
   real(R8P),    intent(in)  :: dxyz_gpu(:,:) !< XYZ space steps.
   real(R8P),    intent(out) :: dxyz_min      !< Minimum space step.
   integer(I4P)              :: b             !< Counter.

   dxyz_min = huge(0._R8P)
   !$acc parallel loop independent gang vector DEVICEVAR(dxyz_gpu) reduction(min: dxyz_min)
   !$omp OMPLOOP DEVICEVAR(dxyz_gpu) reduction(min: dxyz_min)
   do b=1, blocks_number
      dxyz_min = min(dxyz_min, dxyz_gpu(b,1), dxyz_gpu(b,2), dxyz_gpu(b,3))
   enddo
   dxyz_min = dxyz_min * 0.5_R8P
   endsubroutine compute_dxyz_min_dev

   subroutine set_bc_q_gpu_dev(BC_EXTRAPOLATION, BC_fWLAYER, nv, ni, nj, nk, ngc,    &
                               D_divergence_cleaner, B_divergence_cleaner, dxyz_gpu, &
                               l_map_bc_gpu, fec_1_6_array_gpu, q_gpu)
   !< Set BC over q.
   integer(I4P), intent(in)    :: BC_EXTRAPOLATION        !< Extrapolation BC parameter.
   integer(I4P), intent(in)    :: BC_fWLAYER              !< fWLAYER BC parameter.
   integer(I4P), intent(in)    :: nv                      !< Number of variables.
   integer(I4P), intent(in)    :: ni, nj, nk              !< Cells number.
   integer(I4P), intent(in)    :: ngc                     !< Ghost cells number.
   logical                     :: D_divergence_cleaner    !< Enable electric field divergence cleaning.
   logical                     :: B_divergence_cleaner    !< Enable magnetic field divergence cleaning.
   real(R8P),    intent(in)    :: dxyz_gpu(:,:)           !< Delta cells GPU.
   integer(I8P), intent(in)    :: l_map_bc_gpu(:,:,:)     !< Local map for BC ghost cells.
   integer(I4P), intent(in)    :: fec_1_6_array_gpu(:)    !< Local map for BC ghost cells.
   real(R8P),    intent(inout) :: q_gpu(1:,    &
                                        1-ngc:,&
                                        1-ngc:,&
                                        1-ngc:,1:)        !< Conservative variables.
   integer(I4P)                :: b                       !< Counter.
   integer(I4P)                :: c, i, j, k, v           !< Counter.
   integer(I4P)                :: idelta                  !< IJK delta step for extrapolation.
   integer(I4P)                :: jdelta                  !< IJK delta step for extrapolation.
   integer(I4P)                :: kdelta                  !< IJK delta step for extrapolation.
   integer(I4P)                :: bc_type                 !< Boundary condition type.
   integer(I4P)                :: crown                   !< Crown counter.
   integer(I4P)                :: fec                     !< Boundary fec (1 to 26).
   integer(I4P)                :: fec_1_6                 !< Boundary fec (1 to 6).
   integer(I4P)                :: alfa_D, beta_D, gamma_D !< Indici alfa beta gamma come in Barbas.
   integer(I4P)                :: alfa_B, beta_B, gamma_B !< Indici alfa beta gamma come in Barbas.
   real(R8P)                   :: s1                      !< Coefficiente pari a +-1.
   real(R8P)                   :: ds                      !< Distanza tra le celle in x, y o z.
   real(R8P)                   :: ngc_r, crown_r          !< Numero di gc totale, reale
   real(R8P)                   :: ref(1:9)                !< Vettore di stato di riferimento per assegnazione gc.
   real(R8P)                   :: fi, f                   !< Variabili phi e f fWL.

   do crown=1, ngc
      !$acc parallel loop independent gang vector DEVICEVAR(l_map_bc_gpu, fec_1_6_array_gpu, q_gpu)
      !$omp OMPLOOP DEVICEVAR(l_map_bc_gpu, fec_1_6_array_gpu, q_gpu)
      do c=1, size(l_map_bc_gpu, dim=1)
         b = l_map_bc_gpu(c, 1 ,crown)
         if (b>0) then
            i       = l_map_bc_gpu(c, 2 ,crown)
            j       = l_map_bc_gpu(c, 3 ,crown)
            k       = l_map_bc_gpu(c, 4 ,crown)
            idelta  = l_map_bc_gpu(c, 5 ,crown)
            jdelta  = l_map_bc_gpu(c, 6 ,crown)
            kdelta  = l_map_bc_gpu(c, 7 ,crown)
            bc_type = l_map_bc_gpu(c, 8 ,crown)
            fec     = l_map_bc_gpu(c, 9 ,crown)
            fec_1_6 = fec_1_6_array_gpu(fec)
            if (bc_type == BC_EXTRAPOLATION) then
               do v=1, 9
                  q_gpu(b,i,j,k,v) = q_gpu(b,i-idelta,j-jdelta,k-kdelta,v) !ni,j,k coordinate della cella da cui prendo i valori
               enddo
               if (D_divergence_cleaner) then
                  q_gpu(b,i,j,k,10) = 0._R8P
               endif
               if (B_divergence_cleaner) then
                  q_gpu(b,i,j,k,11) = 0._R8P
               endif
            elseif (bc_type == BC_fWLayer) then
               !if (fec <= 6) then
               !   select case(fec)
               !   !Identifico gli alfa beta gamma come nel paper di Barbas, distinguendo tra alfa_D e alfa_B ecc
               !   case(1) ! x- face alfa = 2, beta = 3, gamma = 1
               !      s1 = 1.0_R8P
               !      alfa_D = 2_I4P
               !      beta_D = 3_I4P
               !      gamma_D = 1_I4P
               !      alfa_B = 5_I4P
               !      beta_B = 6_I4P
               !      gamma_B = 4_I4P
               !      ds = dxyz_gpu(b,1) !distanza tra le celle in x
               !      ref = q_gpu(b,1,j,k,:) !vettore di stato di riferimento per assegnazione gc
               !   case(2) ! x+ face
               !      s1 = -1.0_R8P
               !      alfa_D = 2_I4P
               !      beta_D = 3_I4P
               !      gamma_D = 1_I4P
               !      alfa_B = 5_I4P
               !      beta_B = 6_I4P
               !      gamma_B = 4_I4P
               !      ds = dxyz_gpu(b,1) !distanza tra le celle in x
               !      ref = q_gpu(b,ni,j,k,:)
               !   case(3) ! y- face
               !      s1 = 1.0_R8P
               !      alfa_D = 3_I4P
               !      beta_D = 1_I4P
               !      gamma_D = 2_I4P
               !      alfa_B = 6_I4P
               !      beta_B = 4_I4P
               !      gamma_B = 5_I4P
               !      ds = dxyz_gpu(b,2) !distanza tra le celle in y
               !      ref = q_gpu(b,i,1,k,:)
               !   case(4) ! y+ face
               !      s1 = -1.0_R8P
               !      alfa_D = 3_I4P
               !      beta_D = 1_I4P
               !      gamma_D = 2_I4P
               !      alfa_B = 6_I4P
               !      beta_B = 4_I4P
               !      gamma_B = 5_I4P
               !      ds = dxyz_gpu(b,2) !distanza tra le celle in y
               !      ref = q_gpu(b,i,nj,k,:)
               !   case(5) ! z- face
               !      s1 = 1.0_R8P
               !      alfa_D = 1_I4P
               !      beta_D = 2_I4P
               !      gamma_D = 3_I4P
               !      alfa_B = 4_I4P
               !      beta_B = 5_I4P
               !      gamma_B = 6_I4P
               !      ds = dxyz_gpu(b,3) !distanza tra le celle in z
               !      ref = q_gpu(b,i,j,1,:)
               !   case(6) ! z+ face
               !      s1 = -1.0_R8P
               !      alfa_D = 1_I4P
               !      beta_D = 2_I4P
               !      gamma_D = 3_I4P
               !      alfa_B = 4_I4P
               !      beta_B = 5_I4P
               !      gamma_B = 6_I4P
               !      ds = dxyz_gpu(b,3) !distanza tra le celle in z
               !      ref = q_gpu(b,i,j,nk,:)
               !   endselect
               !   ngc_r = real(ngc,R8P)
               !   crown_r = real(crown, R8P)
               !   if (ngc < 40_I4P) then
               !      fi = 1/150._R8P*(-7.0_R8P*ngc_r**2 + 255._R8P*ngc_r + 250._R8P) !polinomio di Barbas
               !   else
               !      fi = 25.0_R8P
               !   endif
               !   !x - xa è la distanza tra il centro della gc considerata e il bordo del dominio (fatto col centro cella, vedr
               !   !è pari quindi a (ngc_r - crown_r) * ds

               !   !xb - xa è la distanza tra il centro della gc più esterna considerata e il bordo del dominio (fatto col centr
               !   !è pari quindi a C * ds

               !   f = 1._R8P/fi*LOG10(((ngc_r-crown_r)*ds)/(ngc_r*ds)*(10._R8P**fi-1._R8P)+1._R8P) !funzione f

               !   q_gpu(b,i,j,k,alfa_D) = 1/(2*MU0**0.5_R8P)*(s1*(f-1._R8P)*ref(beta_B)*EPS0**0.5_R8P + &
               !                           (f+1._R8P)*ref(alfa_D)*MU0**0.5_R8P)

               !   q_gpu(b,i,j,k,beta_B) = 1/(2*EPS0**0.5_R8P)*((f+1._R8P)*ref(beta_B)*EPS0**0.5_R8P + &
               !                           s1*(f-1._R8P)*ref(alfa_D)*MU0**0.5_R8P)

               !   q_gpu(b,i,j,k,beta_D) = 1/(2*MU0**0.5_R8P)*(-s1*(f-1._R8P)*ref(alfa_B)*EPS0**0.5_R8P + &
               !                           (f+1._R8P)*ref(beta_D)*MU0**0.5_R8P)

               !   q_gpu(b,i,j,k,alfa_B) = 1/(2*EPS0**0.5_R8P)*((f+1._R8P)*ref(alfa_B)*EPS0**0.5_R8P - &
               !                           s1*(f-1._R8P)*ref(beta_D)*MU0**0.5_R8P)

               !   q_gpu(b,i,j,k,gamma_D) = ref(gamma_D)

               !   q_gpu(b,i,j,k,gamma_B) = ref(gamma_B)

               !   q_gpu(b,i,j,k,7:9) = 0._R8P
               !else
               !   do v=1, nv
               !     q_gpu(b,i,j,k,v) = 0.0_R8P
               !   enddo
               !endif
            endif
         endif
      enddo
   enddo
   endsubroutine set_bc_q_gpu_dev

   subroutine set_sir_dev(si_gpu, sir_gpu)
   !< Set directional increment si and sir on device.
   integer(I4P), intent(inout) :: si_gpu(3,3)  !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(inout) :: sir_gpu(3,3) !< Directional (1=x,2=y,3=z) increment, real cast.

   !$acc serial DEVICEVAR(si_gpu, sir_gpu)
   si_gpu (1:3,1) = [1,     0,     0     ]
   sir_gpu(1:3,1) = [1._R8P,0._R8P,0._R8P]
   si_gpu (1:3,2) = [0,     1,     0     ]
   sir_gpu(1:3,2) = [0._R8P,1._R8P,0._R8P]
   si_gpu (1:3,3) = [0,     0     ,1     ]
   sir_gpu(1:3,3) = [0._R8P,0._R8P,1._R8P]
   !$acc end serial
   endsubroutine set_sir_dev

   ! private procedures
   subroutine compute_fluxes_convective_ri_dev(dir,b,i,j,k,ngc,nv,evmax,si_gpu,sir_gpu,q_gpu,fluxes_gpu)
   !< Compute convective fluxes at right interface of b,i,j,k.
   integer(I4P), intent(in)    :: dir                                       !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)    :: b, i, j, k                                !< Counter.
   integer(I4P), intent(in)    :: ngc                                       !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                                        !< Number of conservative varibales.
   real(R8P),    intent(in)    :: evmax                                     !< Maximum waves speeds estimation.
   integer(I4P), intent(in)    :: si_gpu(3)                                 !< Stencil increment.
   real(R8P),    intent(in)    :: sir_gpu(3)                                !< Stencil increment, real cast.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)         !< Fields variables.
   real(R8P),    intent(inout) :: fluxes_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Fluxes.
   real(R8P)                   :: fpmr(1:2,1:11)                            !< Fluxes +- reconstructed.
   real(R8P)                   :: fpmr_(1:2,1:11)                           !< Fluxes +- reconstructed (to be removed).
   integer(I4P)                :: v                                         !< Counter.
   !$acc routine(compute_fluxes_convective_ri_dev)
   !$omp declare target(compute_fluxes_convective_ri_dev)

   call decompose_fluxes_convective_dev(sir   = sir_gpu, &
                                        b     = b,       &
                                        i     = i,       &
                                        j     = j,       &
                                        k     = k,       &
                                        ngc   = ngc,     &
                                        nv    = nv,      &
                                        evmax = evmax,   &
                                        q_gpu = q_gpu,   &
                                        fmp   = fpmr)
   call decompose_fluxes_convective_dev(sir   = sir_gpu,     &
                                        b     = b,           &
                                        i     = i+si_gpu(1), &
                                        j     = j+si_gpu(2), &
                                        k     = k+si_gpu(3), &
                                        ngc   = ngc,         &
                                        nv    = nv,          &
                                        evmax = evmax,       &
                                        q_gpu = q_gpu,       &
                                        fmp   = fpmr_)
   do v=1,nv
      fluxes_gpu(b,i,j,k,v) = fpmr(2,v) + fpmr_(1,v)
   enddo
   endsubroutine compute_fluxes_convective_ri_dev

   subroutine decompose_fluxes_convective_dev(sir,b,i,j,k,ngc,nv,evmax,q_gpu,fmp)
   !< Decompose convective fluxes.
   real(R8P),    intent(in)    :: sir(3)                            !< Stencil increment, real cast.
   integer(I4P), intent(in)    :: b, i, j, k                        !< Counter.
   integer(I4P), intent(in)    :: ngc                               !< Ghost cells number.
   real(R8P),    intent(in)    :: evmax                             !< Maximum waves speeds estimation.
   integer(I4P), intent(in)    :: nv                                !< Number of conservative varibales.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: fmp(1:,1:)                        !< Fluxes -+ decomposition.
   real(R8P)                   :: f(11)                             !< Conservative fluxes.
   integer(I4P)                :: v                                 !< Counter.
   !$acc routine(decompose_fluxes_convective_dev)
   !$omp declare target(decompose_fluxes_convective_dev)

   call compute_convective_fluxes_maxwell_dev(sir=sir,q_gpu=q_gpu(b,i,j,k,:),f_gpu=f)
   do v=1, nv
      fmp(2,v) = 0.5_R8P * (f(v) + evmax * q_gpu(b,i,j,k,v))
      fmp(1,v) = f(v) - fmp(2,v)
   enddo
   do v=1, 2
      fmp(v,7) = 0._R8P
      fmp(v,8) = 0._R8P
      fmp(v,9) = 0._R8P
   enddo
   endsubroutine decompose_fluxes_convective_dev

   subroutine compute_convective_fluxes_maxwell_dev(sir, q_gpu, f_gpu)
   !< Compute convective fluxes.
   !<
   !<```
   !< Fx(1) = 0        Fy(1) = -Bz/muz    Fz(1) = By/muy
   !< Fx(2) = Bz/muz   Fy(2) = 0          Fz(2) = -Bx/mux
   !< Fx(3) = -By/muy  Fy(3) = Bx/mux     Fz(3) = 0
   !< Fx(4) = 0        Fy(4) = Dz/epsz    Fz(4) = -Dy/epsy
   !< Fx(5) = -Dz/epsz Fy(5) = 0          Fz(5) = Dx/epsx
   !< Fx(6) = Dy/epsy  Fy(6) = -Dx/epsx   Fz(6) = 0
   !< Fx(7) = 0        Fy(7) = 0          Fz(7) = 0
   !< Fx(8) = 0        Fy(8) = 0          Fz(8) = 0
   !< Fx(9) = 0        Fy(9) = 0          Fz(9) = 0
   !<```
   real(R8P), intent(in)    :: sir(3)    !< Directional (1=x,2=y,3=z) increment.
   real(R8P), intent(in)    :: q_gpu(1:) !< Auxiliary variables.
   real(R8P), intent(inout) :: f_gpu(1:) !< Conservative fluxes.
   !$acc routine(compute_conservative_fluxes_maxwell_dev)
   !$omp declare target(compute_conservative_fluxes_maxwell_dev)

   if (sir(1)==1._R8P) then
      f_gpu(1) =  0.0_R8P
      f_gpu(2) =  q_gpu(6)/MU0
      f_gpu(3) = -q_gpu(5)/MU0
      f_gpu(4) =  0.0_R8P
      f_gpu(5) = -q_gpu(3)/EPS0
      f_gpu(6) =  q_gpu(2)/EPS0
      f_gpu(7) = 0._R8P
      f_gpu(8) = 0._R8P
      f_gpu(9) = 0._R8P
   elseif(sir(2)==1._R8P) then
      f_gpu(1) = -q_gpu(6)/MU0
      f_gpu(2) =  0.0_R8P
      f_gpu(3) =  q_gpu(4)/MU0
      f_gpu(4) =  q_gpu(3)/EPS0
      f_gpu(5) =  0.0_R8P
      f_gpu(6) = -q_gpu(1)/EPS0
      f_gpu(7) = 0._R8P
      f_gpu(8) = 0._R8P
      f_gpu(9) = 0._R8P
   elseif(sir(3)==1._R8P) then
      f_gpu(1) =  q_gpu(5)/MU0
      f_gpu(2) = -q_gpu(4)/MU0
      f_gpu(3) =  0.0_R8P
      f_gpu(4) = -q_gpu(2)/EPS0
      f_gpu(5) =  q_gpu(1)/EPS0
      f_gpu(6) =  0.0_R8P
      f_gpu(7) = 0._R8P
      f_gpu(8) = 0._R8P
      f_gpu(9) = 0._R8P
   endif
   endsubroutine compute_convective_fluxes_maxwell_dev

   subroutine compute_convective_fluxes_maxwell_div_d_dev(sir,q_gpu,chi,f_gpu)
   !< Compute convective fluxes with div(D) correction.
   !<```
   !< Fx(1) = phi      Fy(1) = -Bz/muz    Fz(1) = By/muy
   !< Fx(2) = Bz/muz   Fy(2) = phi        Fz(2) = -Bx/mux
   !< Fx(3) = -By/muy  Fy(3) = Bx/mux     Fz(3) = phi
   !< Fx(4) = 0        Fy(4) = Dz/epsz    Fz(4) = -Dy/epsy
   !< Fx(5) = -Dz/epsz Fy(5) = 0          Fz(5) = Dx/epsx
   !< Fx(6) = Dy/epsy  Fy(6) = -Dx/epsx   Fz(6) = 0
   !< Fx(7) = 0        Fy(7) = 0          Fz(7) = 0
   !< Fx(8) = 0        Fy(8) = 0          Fz(8) = 0
   !< Fx(9) = 0        Fy(9) = 0          Fz(9) = 0
   !< Fx(10) = ch^2*Dx Fy(10) = ch^2*Dy   Fz(10) = ch^2*Dz
   !<```
   real(R8P),    intent(in)    :: sir(3)    !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)    :: q_gpu(1:) !< Auxiliary variables.
   real(R8P),    intent(in)    :: chi       !< Coefficiente velocità trasporto errori divergenza campi
   real(R8P),    intent(inout) :: f_gpu(1:) !< Conservative fluxes.
   real(R8P)                   :: ch        !< Velocità trasporto errori divergenza campi modello iperbolico
   !$acc routine(compute_conservative_fluxes_maxwell_div_d_dev)
   !$omp declare target(compute_conservative_fluxes_maxwell_div_d_dev)

   ch = chi/sqrt(EPS0*MU0)
   if (sir(1)==1._R8P) then
      f_gpu(1)  =  q_gpu(10)
      f_gpu(2)  =  q_gpu(6)/MU0
      f_gpu(3)  = -q_gpu(5)/MU0
      f_gpu(4)  =  0.0_R8P
      f_gpu(5)  = -q_gpu(3)/EPS0
      f_gpu(6)  =  q_gpu(2)/EPS0
      f_gpu(7)  = 0._R8P
      f_gpu(8)  = 0._R8P
      f_gpu(9)  = 0._R8P
      f_gpu(10) = ch*ch*q_gpu(1)
   elseif(sir(2)==1._R8P) then
      f_gpu(1)  = -q_gpu(6)/MU0
      f_gpu(2)  =  q_gpu(10)
      f_gpu(3)  =  q_gpu(4)/MU0
      f_gpu(4)  =  q_gpu(3)/EPS0
      f_gpu(5)  =  0.0_R8P
      f_gpu(6)  = -q_gpu(1)/EPS0
      f_gpu(7)  = 0._R8P
      f_gpu(8)  = 0._R8P
      f_gpu(9)  = 0._R8P
      f_gpu(10) = ch*ch*q_gpu(2)
   elseif(sir(3)==1._R8P) then
      f_gpu(1)  =  q_gpu(5)/MU0
      f_gpu(2)  = -q_gpu(4)/MU0
      f_gpu(3)  =  q_gpu(10)
      f_gpu(4)  = -q_gpu(2)/EPS0
      f_gpu(5)  =  q_gpu(1)/EPS0
      f_gpu(6)  =  0.0_R8P
      f_gpu(7)  = 0._R8P
      f_gpu(8)  = 0._R8P
      f_gpu(9)  = 0._R8P
      f_gpu(10) = ch*ch*q_gpu(3)
   endif
   endsubroutine compute_convective_fluxes_maxwell_div_d_dev

   subroutine compute_convective_fluxes_maxwell_div_d_b_dev(sir,q_gpu,chi,f_gpu)
   !< Compute convective fluxes with div(D) and div(B) correction.
   ! Fx(1) = phi      Fy(1) = -Bz/muz    Fz(1) = By/muy
   ! Fx(2) = Bz/muz   Fy(2) = phi        Fz(2) = -Bx/mux
   ! Fx(3) = -By/muy  Fy(3) = Bx/mux     Fz(3) = phi
   ! Fx(4) = csi      Fy(4) = Dz/epsz    Fz(4) = -Dy/epsy
   ! Fx(5) = -Dz/epsz Fy(5) = csi        Fz(5) = Dx/epsx
   ! Fx(6) = Dy/epsy  Fy(6) = -Dx/epsx   Fz(6) = csi
   ! Fx(7) = 0        Fy(7) = 0          Fz(7) = 0
   ! Fx(8) = 0        Fy(8) = 0          Fz(8) = 0
   ! Fx(9) = 0        Fy(9) = 0          Fz(9) = 0
   ! Fx(10) = ch^2*Dx Fy(10) = ch^2*Dy   Fz(10) = ch^2*Dz
   ! Fx(11) = ch^2*Bx Fy(10) = ch^2*By   Fz(10) = ch^2*Bz
   real(R8P),    intent(in)    :: sir(3)    !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)    :: q_gpu(1:) !< Auxiliary variables.
   real(R8P),    intent(in)    :: chi       !< Coefficiente velocità trasporto errori divergenza campi
   real(R8P),    intent(inout) :: f_gpu(1:) !< Conservative fluxes.
   real(R8P)                   :: ch        !< Velocità trasporto errori divergenza campi modello iperbolico
   !$acc routine(compute_conservative_fluxes_maxwell_div_d_b_dev)
   !$omp declare target(compute_conservative_fluxes_maxwell_div_d_b_dev)

   ch = chi/sqrt(EPS0*MU0)
   if (sir(1)==1._R8P) then
      f_gpu(1)  =  q_gpu(10)
      f_gpu(2)  =  q_gpu(6)/MU0
      f_gpu(3)  = -q_gpu(5)/MU0
      f_gpu(4)  =  q_gpu(11)
      f_gpu(5)  = -q_gpu(3)/EPS0
      f_gpu(6)  =  q_gpu(2)/EPS0
      f_gpu(7)  = 0._R8P
      f_gpu(8)  = 0._R8P
      f_gpu(9)  = 0._R8P
      f_gpu(10) = ch*ch*q_gpu(1)
      f_gpu(11) = ch*ch*q_gpu(4)
   elseif(sir(2)==1._R8P) then
      f_gpu(1)  = -q_gpu(6)/MU0
      f_gpu(2)  =  q_gpu(10)
      f_gpu(3)  =  q_gpu(4)/MU0
      f_gpu(4)  =  q_gpu(3)/EPS0
      f_gpu(5)  =  q_gpu(11)
      f_gpu(6)  = -q_gpu(1)/EPS0
      f_gpu(7)  = 0._R8P
      f_gpu(8)  = 0._R8P
      f_gpu(9)  = 0._R8P
      f_gpu(10) = ch*ch*q_gpu(2)
      f_gpu(11) = ch*ch*q_gpu(5)
   elseif(sir(3)==1._R8P) then
      f_gpu(1)  =  q_gpu(5)/MU0
      f_gpu(2)  = -q_gpu(4)/MU0
      f_gpu(3)  =  q_gpu(10)
      f_gpu(4)  = -q_gpu(2)/EPS0
      f_gpu(5)  =  q_gpu(1)/EPS0
      f_gpu(6)  =  q_gpu(11)
      f_gpu(7)  = 0._R8P
      f_gpu(8)  = 0._R8P
      f_gpu(9)  = 0._R8P
      f_gpu(10) = ch*ch*q_gpu(3)
      f_gpu(11) = ch*ch*q_gpu(6)
   endif
   endsubroutine compute_convective_fluxes_maxwell_div_d_b_dev
endmodule adam_prism_fnl_kernels
