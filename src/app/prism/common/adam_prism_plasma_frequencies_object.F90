!< PRISM plasma characteristic frequencies and the time-step limit they set for the PIC models.
!<
!< Reference: `frequenze_caratteristiche_plasma.pdf` (rev. 2), sections 1-3 (frequencies) and 9 (Runge-Kutta criterion).
!< The maxima are evaluated at every step (and at t0), after the step, for the next time step:
!< - omega_p^2 per cell with NGP: sum over the particles of each cell of q_p^2/(m_p eps0 V), electrons (q < 0) and ions
!<   (q > 0) apart; q_p and m_p are the macro-particle charge and mass (both scaled by the weight), so q_p^2/m_p = w q^2/m.
!<   NGP needs no ghost fold and no communication; on a noisy density it overestimates the maximum (conservative).
!< - Omega_c per particle: |q_p| |B(x_p)| / m_p (non-relativistic pusher), with B the field gathered at the particle, the
!<   external fields included; electrons and ions apart.
!< - transit: max over the particles of max_d |v_d| / dx_d (cells per unit time).
!< From the global maxima: omega_UH, omega_LH, the cutoffs omega_R, omega_L, omega_O; omega_R >= omega_UH >=
!< max(omega_pe, Omega_ce) is the highest, omega_max. Each frequency grows with omega_pe and Omega_ce, so taking the two
!< maxima where they are (possibly in different cells) bounds omega_max from above.
!< Time step: with dt_CFL the stability-limited step of the fields, dt = dt_CFL if omega_max dt_CFL < omega_dt_max
!< (`[PIC] omega_dt_max`, default 0.05, accuracy), otherwise dt = min(omega_dt_max/omega_max,
!< (1/dt_CFL^2 + omega_max^2/beta_I^2)^(-1/2)), beta_I the imaginary-axis stability bound of the RK scheme (section 9a);
!< the second term never prevails while omega_dt_max < beta_I, it is kept as the stability statement. The step source,
!< 'CFL' or 'plasma_frequencies', is printed at every step.
module adam_prism_plasma_frequencies_object

use :: adam_common_library
use :: adam_mpih_global, only : mpih
use :: adam_prism_parameters, only : EPS0
use :: adam_prism_physics_object, only : prism_physics_object, is_adim_model, is_pic_model
use :: adam_prism_pic_object, only : prism_pic_object, NUM_SCHEME_TIME_PIC_RUNGE_KUTTA
use :: adam_prism_numerics_object, only : prism_numerics_object, NUM_SCHEME_TIME_RUNGE_KUTTA
use :: mpi
use :: penf

implicit none
private
public :: prism_plasma_frequencies_object
public :: DT_SOURCE_CFL
public :: DT_SOURCE_PLASMA

character(len=3),  parameter :: DT_SOURCE_CFL    = 'CFL'                !< dt set by the fields (CFL, GLM, PML).
character(len=18), parameter :: DT_SOURCE_PLASMA = 'plasma_frequencies' !< dt set by the plasma frequencies.

type :: prism_plasma_frequencies_object
   logical      :: active = .false.         !< PIC model: frequencies evaluated and the limit applied.
   real(R8P)    :: omega_dt_max = 0.05_R8P  !< Accuracy threshold on omega_max dt.
   real(R8P)    :: beta_i = 0._R8P          !< Imaginary-axis stability bound of the RK scheme (0: none).
   logical      :: single_step = .true.     !< Field and PIC time schemes are single-step RK (dt may vary).
   real(R8P)    :: inv_eps = 1._R8P/EPS0    !< 1/eps0 (1 for the non-dimensional models).
   real(R8P)    :: wpe = 0._R8P             !< Max electron plasma frequency [rad/s].
   real(R8P)    :: wpi = 0._R8P             !< Max ion plasma frequency [rad/s].
   real(R8P)    :: wce = 0._R8P             !< Max electron cyclotron frequency [rad/s].
   real(R8P)    :: wci = 0._R8P             !< Max ion cyclotron frequency [rad/s].
   real(R8P)    :: wuh = 0._R8P             !< Upper hybrid frequency [rad/s].
   real(R8P)    :: wlh = 0._R8P             !< Lower hybrid frequency [rad/s].
   real(R8P)    :: wr  = 0._R8P             !< R cutoff [rad/s].
   real(R8P)    :: wl  = 0._R8P             !< L cutoff [rad/s].
   real(R8P)    :: wo  = 0._R8P             !< O cutoff (= omega_pe) [rad/s].
   real(R8P)    :: wmax = 0._R8P            !< Highest characteristic frequency (omega_R) [rad/s].
   real(R8P)    :: transit = 0._R8P         !< Max over the particles of max_d |v_d|/dx_d [1/s].
   real(R8P), allocatable :: wp2(:,:,:,:,:) !< CPU work: omega_p^2 per interior cell [2,ni,nj,nk,nb].
   contains
      procedure, pass(self) :: initialize
      procedure, pass(self) :: description
      procedure, pass(self) :: compute_cpu
      procedure, pass(self) :: set_maxima
      procedure, pass(self) :: dt_limit
      procedure, pass(self) :: dt_source
      procedure, pass(self) :: check_step
endtype prism_plasma_frequencies_object

contains
   subroutine initialize(self, pic, physics, numerics, rk_scheme)
   !< Activate for the PIC models, read the threshold from pic, set beta_I from the RK scheme.
   class(prism_plasma_frequencies_object), intent(inout) :: self      !< Plasma frequencies.
   type(prism_pic_object),                 intent(in)    :: pic       !< PIC (omega_dt_max, time scheme).
   type(prism_physics_object),             intent(in)    :: physics   !< Physics.
   type(prism_numerics_object),            intent(in)    :: numerics  !< Numerics (field time scheme).
   character(*),                           intent(in)    :: rk_scheme !< RK scheme name.

   self%active = is_pic_model(physics%physical_model)
   if (.not.self%active) return
   self%omega_dt_max = pic%omega_dt_max
   self%inv_eps = 1._R8P/EPS0
   if (is_adim_model(physics%physical_model)) self%inv_eps = 1._R8P
   self%single_step = trim(numerics%scheme_time) == NUM_SCHEME_TIME_RUNGE_KUTTA .and. &
                      trim(pic%scheme_time) == NUM_SCHEME_TIME_PIC_RUNGE_KUTTA
   ! end of the stability region on the imaginary axis: |R(i y)| <= 1 for |y| <= beta_I (0: none, R(i y) > 1 at once)
   select case(trim(rk_scheme))
   case(RK_SSP_54)
      self%beta_i = 3.278_R8P            ! computed from the Shu-Osher coefficients of ssp_tableau
   case(RK_3, RK_SSP_33)
      self%beta_i = sqrt(3._R8P)
   case default
      self%beta_i = 0._R8P
   endselect
   endsubroutine initialize

   function description(self) result(desc)
   !< Frequencies, omega dt and the transit, for the log.
   class(prism_plasma_frequencies_object), intent(in) :: self !< Plasma frequencies.
   character(len=:), allocatable                      :: desc !< Description.

   desc = 'plasma frequencies [rad/s]: w_pe='//trim(str(self%wpe))//' w_pi='//trim(str(self%wpi))// &
          ' W_ce='//trim(str(self%wce))//' W_ci='//trim(str(self%wci))//' w_UH='//trim(str(self%wuh))// &
          ' w_LH='//trim(str(self%wlh))//' w_R='//trim(str(self%wr))//' w_L='//trim(str(self%wl))// &
          ' w_O='//trim(str(self%wo))//'; max |v|/dx='//trim(str(self%transit))//' 1/s'
   endfunction description

   subroutine compute_cpu(self, pic, field, q_pic, pic_fields, blocks_number, ni, nj, nk)
   !< Local maxima on the CPU, then the global ones (set_maxima, collective).
   class(prism_plasma_frequencies_object), intent(inout) :: self              !< Plasma frequencies.
   type(prism_pic_object),                 intent(in)    :: pic               !< PIC (particles, neighbour list).
   type(field_object),                     intent(in)    :: field             !< Field (cell sizes).
   real(R8P),                              intent(in)    :: q_pic(1:,1:)      !< Particles [8,np].
   real(R8P),                              intent(in)    :: pic_fields(1:,1:) !< Fields at the particles [6,np].
   integer(I4P),                           intent(in)    :: blocks_number     !< Blocks of this rank.
   integer(I4P),                           intent(in)    :: ni, nj, nk        !< Block cells.
   real(R8P)                                             :: wpe2, wpi2, wce, wci, tr, qp, w, val
   integer(I4P)                                          :: n, b, i, j, k, s

   if (.not.self%active) return
   wce = 0._R8P ; wci = 0._R8P ; tr = 0._R8P ; wpe2 = 0._R8P ; wpi2 = 0._R8P
   if (allocated(self%wp2)) then
      if (size(self%wp2, dim=5) < blocks_number) deallocate(self%wp2)
   endif
   if (.not.allocated(self%wp2)) allocate(self%wp2(2,ni,nj,nk,max(1_I4P,blocks_number)))
   self%wp2 = 0._R8P
   if (pic%particle_number > 0_I4P .and. allocated(pic%neighbour_list)) then
      !$omp parallel do private(qp, w, b, i, j, k, s, val) reduction(max: wce, wci, tr)
      do n = 1, pic%particle_number
         qp = q_pic(7,n)
         if (qp == 0._R8P) cycle
         w = abs(qp)*sqrt(pic_fields(4,n)**2 + pic_fields(5,n)**2 + pic_fields(6,n)**2)/q_pic(8,n)
         if (qp < 0._R8P) then
            wce = max(wce, w)
         else
            wci = max(wci, w)
         endif
         b = pic%neighbour_list(1,n)
         if (b <= 0_I4P) cycle
         tr = max(tr, abs(q_pic(4,n))/field%dxyz(1,b), abs(q_pic(5,n))/field%dxyz(2,b), abs(q_pic(6,n))/field%dxyz(3,b))
         i = pic%neighbour_list(2,n) ; j = pic%neighbour_list(3,n) ; k = pic%neighbour_list(4,n)
         ! a particle in the ghost halo (wrapped across a periodic face at the end of the step, one rank) is skipped
         if (i < 1 .or. i > ni .or. j < 1 .or. j > nj .or. k < 1 .or. k > nk) cycle
         s = 1_I4P ; if (qp > 0._R8P) s = 2_I4P
         val = qp*qp*self%inv_eps/(q_pic(8,n)*field%dxyz(1,b)*field%dxyz(2,b)*field%dxyz(3,b))
         !$omp atomic update
         self%wp2(s,i,j,k,b) = self%wp2(s,i,j,k,b) + val
      enddo
      if (blocks_number > 0_I4P) then
         wpe2 = maxval(self%wp2(1,:,:,:,1:blocks_number))
         wpi2 = maxval(self%wp2(2,:,:,:,1:blocks_number))
      endif
   endif
   call self%set_maxima(wpe2=wpe2, wpi2=wpi2, wce=wce, wci=wci, transit=tr)
   endsubroutine compute_cpu

   subroutine set_maxima(self, wpe2, wpi2, wce, wci, transit)
   !< Reduce the local maxima over the ranks (collective) and derive the composite frequencies (section 3).
   class(prism_plasma_frequencies_object), intent(inout) :: self    !< Plasma frequencies.
   real(R8P),                              intent(in)    :: wpe2    !< Local max omega_pe^2.
   real(R8P),                              intent(in)    :: wpi2    !< Local max omega_pi^2.
   real(R8P),                              intent(in)    :: wce     !< Local max Omega_ce.
   real(R8P),                              intent(in)    :: wci     !< Local max Omega_ci.
   real(R8P),                              intent(in)    :: transit !< Local max |v|/dx.
   real(R8P)                                             :: m(5)    !< Maxima.
   real(R8P)                                             :: d       !< sqrt(Omega_ce^2 + 4 omega_pe^2).
   real(R8P)                                             :: inv     !< 1/omega_LH^2.

   if (.not.self%active) return
   m = [wpe2, wpi2, wce, wci, transit]
   call MPI_ALLREDUCE(MPI_IN_PLACE, m, 5, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, mpih%error)
   self%wpe = sqrt(m(1)) ; self%wpi = sqrt(m(2)) ; self%wce = m(3) ; self%wci = m(4) ; self%transit = m(5)
   self%wuh = sqrt(self%wpe**2 + self%wce**2)
   inv = 0._R8P
   if (self%wce*self%wci > 0._R8P) inv = inv + 1._R8P/(self%wce*self%wci)
   if (self%wpi**2 + self%wci**2 > 0._R8P) inv = inv + 1._R8P/(self%wpi**2 + self%wci**2)
   self%wlh = 0._R8P ; if (inv > 0._R8P) self%wlh = 1._R8P/sqrt(inv)
   d = sqrt(self%wce**2 + 4._R8P*self%wpe**2)
   self%wr = 0.5_R8P*(d + self%wce)
   self%wl = 0.5_R8P*(d - self%wce)
   self%wo = self%wpe
   self%wmax = max(self%wr, self%wuh, self%wpe, self%wce, self%wpi, self%wci, self%wlh)
   endsubroutine set_maxima

   pure function dt_limit(self, dt_cfl) result(dt)
   !< The step: dt_CFL while omega_max dt_CFL < omega_dt_max, otherwise the accuracy limit (and the RK stability bound).
   class(prism_plasma_frequencies_object), intent(in) :: self   !< Plasma frequencies.
   real(R8P),                              intent(in) :: dt_cfl !< Step limited by the fields.
   real(R8P)                                          :: dt     !< Step.

   dt = dt_cfl
   if (.not.self%active .or. self%wmax <= 0._R8P) return
   if (self%wmax*dt_cfl < self%omega_dt_max) return
   dt = self%omega_dt_max/self%wmax
   if (self%beta_i > 0._R8P) dt = min(dt, 1._R8P/sqrt(1._R8P/dt_cfl**2 + (self%wmax/self%beta_i)**2))
   endfunction dt_limit

   pure function dt_source(self, dt) result(source)
   !< 'plasma_frequencies' when the global step is the plasma limit, 'CFL' otherwise.
   class(prism_plasma_frequencies_object), intent(in) :: self   !< Plasma frequencies.
   real(R8P),                              intent(in) :: dt     !< Global step (forest reduction).
   character(len=:), allocatable                      :: source !< Step source.

   source = DT_SOURCE_CFL
   if (.not.self%active .or. self%wmax <= 0._R8P) return
   ! the CFL step satisfies omega_max dt < omega_dt_max strictly; the plasma step has omega_max dt = omega_dt_max
   if (self%wmax*dt >= self%omega_dt_max*(1._R8P - 1.e-12_R8P)) source = DT_SOURCE_PLASMA
   endfunction dt_source

   subroutine check_step(self, dt)
   !< Stop when the particles cross more than one cell per step, or when the plasma limit acts on a constant-dt scheme.
   class(prism_plasma_frequencies_object), intent(in) :: self !< Plasma frequencies.
   real(R8P),                              intent(in) :: dt   !< Global step.

   if (.not.self%active) return
   if (self%transit*dt > 1._R8P) &
      call mpih%error_stop(msg=': PIC particles cross '//trim(str(self%transit*dt))//' cells in a step (max |v| dt/dx > 1, '// &
                               'dt='//trim(str(dt))//'): reduce [time].(CFL)')
   if (.not.self%single_step .and. self%dt_source(dt) == DT_SOURCE_PLASMA) &
      call mpih%error_stop(msg=': the plasma frequencies limit dt (omega_max dt_CFL >= [PIC].(omega_dt_max)), but the '// &
                               'field or PIC time scheme is not a single-step Runge-Kutta and needs a constant dt: reduce '// &
                               '[time].(CFL) (leapfrog + Boris is not supported yet)')
   endsubroutine check_step
endmodule adam_prism_plasma_frequencies_object
