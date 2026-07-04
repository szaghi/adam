!< ADAM, RK class definition.
module adam_rk_object
!< ADAM, RK class definition.

! ADAM singleton objects
use :: adam_field_object, only : field_object
use :: adam_mpih_global,  only : mpih
use :: adam_grid_object,  only : grid_object
! third party modules
use :: finer
use :: penf

implicit none
private
public :: rk_object
public :: RK_1
public :: RK_2
public :: RK_3
public :: RK_SSP_11
public :: RK_SSP_22
public :: RK_SSP_33
public :: RK_SSP_54
public :: RK_YOSHIDA

character(len=13), parameter :: RK_1       ="runge-kutta-1"       !< Parameter of time scheme, Runge-Kutta 1.
character(len=13), parameter :: RK_2       ="runge-kutta-2"       !< Parameter of time scheme, Runge-Kutta 2.
character(len=13), parameter :: RK_3       ="runge-kutta-3"       !< Parameter of time scheme, Runge-Kutta 3.
character(len=18), parameter :: RK_SSP_11  ="runge-kutta-ssp-11"  !< Parameter of time scheme, Runge-Kutta SSP 11 (forward Euler).
character(len=18), parameter :: RK_SSP_22  ="runge-kutta-ssp-22"  !< Parameter of time scheme, Runge-Kutta SSP 22.
character(len=18), parameter :: RK_SSP_33  ="runge-kutta-ssp-33"  !< Parameter of time scheme, Runge-Kutta SSP 33.
character(len=18), parameter :: RK_SSP_54  ="runge-kutta-ssp-54"  !< Parameter of time scheme, Runge-Kutta SSP 54.
character(len=19), parameter :: RK_YOSHIDA ="runge-kutta-yoshida" !< Parameter of time scheme, Runge-Kutta Yoshida, symplectic 4.

character(len=11), parameter :: INI_SECTION_NAME="runge_kutta" !< INI (config) file section name containing configs.

type :: rk_object
   !< RK class definition.
   character(:), allocatable :: scheme    !< RK scheme.
   integer(I4P)              :: nrk=3_I4P !< Runge-Kutta stages number.
   ! classic, Butcher schemes
   real(R8P), allocatable    :: ark(:)    !< Runge-Kutta low storage alpha coefficients.
   real(R8P), allocatable    :: brk(:)    !< Runge-Kutta low storage beta coefficients.
   real(R8P), allocatable    :: crk(:)    !< Runge-Kutta low storage beta coefficients.
   real(R8P), allocatable    :: alph(:,:) !< Runge-Kutta SSP alpha coefficients.
   real(R8P), allocatable    :: beta(:)   !< Runge-Kutta SSP beta coefficients.
   real(R8P), allocatable    :: gamm(:)   !< Runge-Kutta SSP gamma coefficients.
   ! symplectic (splitting) schemes
   real(R8P), allocatable :: ssa(:) !< Runge-Kutta sympletic-splitting part A coefficients.
   real(R8P), allocatable :: ssb(:) !< Runge-Kutta sympletic-splitting part B coefficients.
   ! RK data
   real(R8P), allocatable    :: q_rk(:,:,:,:,:,:) !< Field cell centered variables, RK stages.
   ! grid data replica for easy handling
   integer(I4P),       pointer :: ngc=>null()           !< Number of ghost cells.
   integer(I4P),       pointer :: ni=>null()            !< Number of cells in i direction.
   integer(I4P),       pointer :: nj=>null()            !< Number of cells in j direction.
   integer(I4P),       pointer :: nk=>null()            !< Number of cells in k direction.
   contains
      ! public methods
      procedure, pass(self) :: assign_stage      !< Assign q to RK stage.
      procedure, pass(self) :: compute_stage     !< Compute RK stage.
      procedure, pass(self) :: compute_stage_ls  !< Compute RK stage, low storage scheme.
      procedure, pass(self) :: description       !< Return pretty-printed object description.
      procedure, pass(self) :: initialize        !< Initialize class.
      procedure, pass(self) :: initialize_stages !< Initialize RK stages.
      procedure, pass(self) :: load_from_file    !< Load config from file.
      procedure, pass(self) :: update_q          !< Update RK q.
endtype rk_object
contains
   ! public methods
   subroutine assign_stage(self, field, s, q, phi)
   !< Assign q to RK stage.
   class(rk_object), intent(inout)        :: self          !< RK object.
   type(field_object), intent(in)         :: field         !< Field (sibling realm component, threaded in).
   integer(I4P),     intent(in)           :: s             !< Current stage number.
   real(R8P),        intent(in)           :: q(1:     ,     &
                                               1-self%ngc:, &
                                               1-self%ngc:, &
                                               1-self%ngc:, &
                                               1:)         !< Conservative variables.
   real(R8P),        intent(in), optional :: phi(1:,          &
                                                 1-self%ngc:, &
                                                 1-self%ngc:, &
                                                 1-self%ngc:, &
                                                 1:)       !< IB distance.
   integer(I4P)                           :: all_solids    !< Last phi index, all solids summary.
   integer(I4P)                           :: i, j, k, b, v !< Counter.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>field%nv, blocks_number=>field%blocks_number)
   if (present(phi)) then
      all_solids = ubound(phi, dim=1)
      !$omp parallel do collapse(5) default(firstprivate) shared(phi,q,self)
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do v=1, nv
                     if (phi(all_solids,i,j,k,b) < 0._R8P) then
                        self%q_rk(v,i,j,k,b,s) = q(v,i,j,k,b)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   else
      !$omp parallel do collapse(5) default(firstprivate) shared(q,self)
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do v=1, nv
                     self%q_rk(v,i,j,k,b,s) = q(v,i,j,k,b)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   endif
   endassociate
   endsubroutine assign_stage

   subroutine compute_stage(self, field, s, dt, phi)
   !< Compute RK stage.
   class(rk_object), intent(inout)        :: self    !< RK object.
   type(field_object), intent(in)         :: field   !< Field (sibling realm component, threaded in).
   integer(I4P),     intent(in)           :: s       !< Current stage number.
   real(R8P),        intent(in)           :: dt      !< Current time step.
   real(R8P),        intent(in), optional :: phi(1:,          &
                                                 1-self%ngc:, &
                                                 1-self%ngc:, &
                                                 1-self%ngc:, &
                                                 1:) !< IB distance.
   integer(I4P)                           :: all_solids        !< Last phi index, all solids summary.
   integer(I4P)                           :: i, j, k, b, v, ss !< Counter.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>field%nv, blocks_number=>field%blocks_number)
   if (present(phi)) then
      all_solids = ubound(phi, dim=1)
      !$omp parallel do collapse(6) default(firstprivate) shared(phi,self)
      do ss=1, s-1
         do b=1, blocks_number
            do k=1, nk
               do j=1, nj
                  do i=1, ni
                     do v=1, nv
                        if (phi(all_solids,i,j,k,b) < 0._R8P) then
                           self%q_rk(v,i,j,k,b,s) = self%q_rk(v,i,j,k,b,s) + dt * self%alph(s,ss) * self%q_rk(v,i,j,k,b,ss)
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   else
      !$omp parallel do collapse(6) default(firstprivate) shared(self)
      do ss=1, s-1
         do b=1, blocks_number
            do k=1, nk
               do j=1, nj
                  do i=1, ni
                     do v=1, nv
                        self%q_rk(v,i,j,k,b,s) = self%q_rk(v,i,j,k,b,s) + dt * self%alph(s, ss) * self%q_rk(v,i,j,k,b,ss)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   endif
   endassociate
   endsubroutine compute_stage

   subroutine compute_stage_ls(self, field, s, dt, phi, dq, q)
   !< Compute RK stage, low storage scheme.
   !< The first (only) stage is assumed to be the previous time step q solution.
   class(rk_object), intent(in)           :: self          !< RK object.
   type(field_object), intent(in)         :: field         !< Field (sibling realm component, threaded in).
   integer(I4P),     intent(in)           :: s             !< Current RK stage.
   real(R8P),        intent(in)           :: dt            !< Current time step.
   real(R8P),        intent(in), optional :: phi(1:,          &
                                                 1-self%ngc:, &
                                                 1-self%ngc:, &
                                                 1-self%ngc:, &
                                                 1:)       !< IB distance.
   real(R8P),        intent(in)           :: dq(1:,          &
                                                1-self%ngc:, &
                                                1-self%ngc:, &
                                                1-self%ngc:, &
                                                1:)        !< Conservative variables residuals.
   real(R8P),        intent(inout)        :: q(1:,          &
                                               1-self%ngc:, &
                                               1-self%ngc:, &
                                               1-self%ngc:, &
                                               1:)         !< Conservative variables stage.
   integer(I4P)                           :: all_solids    !< Last phi index, all solids summary.
   integer(I4P)                           :: i, j, k, b, v !< Counter.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>field%nv, blocks_number=>field%blocks_number)

   !print *, maxval(dq(7,:,:,:,:)), 'Stampa il valore del residuo massimo della densità di corrente lungo x, RK input'
   !print *, maxval(dq(8,:,:,:,:)), 'Stampa il valore del residuo massimo della densità di corrente lungo y, RK input'
   !print *, maxval(dq(9,:,:,:,:)), 'Stampa il valore del residuo massimo della densità di corrente lungo z, RK input'
   !print *, minval(dq(7,:,:,:,:)), 'Stampa il valore del residuo minimo della densità di corrente lungo x, RK input'
   !print *, minval(dq(8,:,:,:,:)), 'Stampa il valore del residuo minimo della densità di corrente lungo y, RK input'
   !print *, minval(dq(9,:,:,:,:)), 'Stampa il valore del residuo minimo della densità di corrente lungo z, RK input'
!
!
   !print *, maxval(q(7,:,:,:,:)), 'Stampa il valore massimo della densità di corrente lungo x RK input'
   !print *, maxval(q(8,:,:,:,:)), 'Stampa il valore massimo della densità di corrente lungo y RK input'
   !print *, maxval(q(9,:,:,:,:)), 'Stampa il valore massimo della densità di corrente lungo z RK input'
   !print *, minval(q(7,:,:,:,:)), 'Stampa il valore minimo della densità di corrente lungo x  RK input'
   !print *, minval(q(8,:,:,:,:)), 'Stampa il valore minimo della densità di corrente lungo y  RK input'
   !print *, minval(q(9,:,:,:,:)), 'Stampa il valore minimo della densità di corrente lungo z  RK input'

   !print *, maxval(self%q_rk(7,:,:,:,:,1)), 'Stampa il valore massimo della densità di corrente lungo x RK input2'
   !print *, maxval(self%q_rk(8,:,:,:,:,1)), 'Stampa il valore massimo della densità di corrente lungo y RK input2'
   !print *, maxval(self%q_rk(9,:,:,:,:,1)), 'Stampa il valore massimo della densità di corrente lungo z RK input2'
   !print *, minval(self%q_rk(7,:,:,:,:,1)), 'Stampa il valore minimo della densità di corrente lungo x  RK input2'
   !print *, minval(self%q_rk(8,:,:,:,:,1)), 'Stampa il valore minimo della densità di corrente lungo y  RK input2'
   !print *, minval(self%q_rk(9,:,:,:,:,1)), 'Stampa il valore minimo della densità di corrente lungo z  RK input2'

   if (present(phi)) then
      all_solids = ubound(phi, dim=1)
      !$omp parallel do collapse(5) default(firstprivate) shared(phi,q,dq,self)
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do v=1, nv
                     if (phi(all_solids,i,j,k,b) < 0._R8P) then
                        q(v,i,j,k,b) = self%ark(s) * self%q_rk(v,i,j,k,b,1) + &
                                       self%brk(s) * q(v,i,j,k,b) + dt * self%crk(s) * dq(v,i,j,k,b)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   else
      !$omp parallel do collapse(5) default(firstprivate) shared(q,dq,self)
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do v=1, nv
                     q(v,i,j,k,b) = self%ark(s) * self%q_rk(v,i,j,k,b,1) + &
                                    self%brk(s) * q(v,i,j,k,b) + dt * self%crk(s) * dq(v,i,j,k,b)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   endif

!   print *, maxval(dq(7,:,:,:,:)), 'Stampa il valore del residuo massimo della densità di corrente lungo x, RK output'
!   print *, maxval(dq(8,:,:,:,:)), 'Stampa il valore del residuo massimo della densità di corrente lungo y, RK output'
!   print *, maxval(dq(9,:,:,:,:)), 'Stampa il valore del residuo massimo della densità di corrente lungo z, RK output'
!   print *, minval(dq(7,:,:,:,:)), 'Stampa il valore del residuo minimo della densità di corrente lungo x, RK output'
!   print *, minval(dq(8,:,:,:,:)), 'Stampa il valore del residuo minimo della densità di corrente lungo y, RK output'
!   print *, minval(dq(9,:,:,:,:)), 'Stampa il valore del residuo minimo della densità di corrente lungo z, RK output'
!
!   print *, maxval(q(7,:,:,:,:)), 'Stampa il valore massimo della densità di corrente lungo x RK output'
!   print *, maxval(q(8,:,:,:,:)), 'Stampa il valore massimo della densità di corrente lungo y RK output'
!   print *, maxval(q(9,:,:,:,:)), 'Stampa il valore massimo della densità di corrente lungo z RK output'
!   print *, minval(q(7,:,:,:,:)), 'Stampa il valore minimo della densità di corrente lungo x  RK output'
!   print *, minval(q(8,:,:,:,:)), 'Stampa il valore minimo della densità di corrente lungo y  RK output'
!   print *, minval(q(9,:,:,:,:)), 'Stampa il valore minimo della densità di corrente lungo z  RK output'

   endassociate
   endsubroutine compute_stage_ls

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(rk_object), intent(in)  :: self             !< RK object.
   character(len=:), allocatable :: desc             !< Description.
   character(len=1), parameter   :: NL=new_line('a') !< New line character.
   integer(I4P)                  :: s                !< Counter.

   desc =       mpih%myrankstr//'Runge-Kutta scheme main data'//NL
   if (allocated(self%ark)) &
   desc = desc//mpih%myrankstr//'  ark:                             '//trim(str(self%ark                ))//NL
   if (allocated(self%brk)) &
   desc = desc//mpih%myrankstr//'  brk:                             '//trim(str(self%brk                ))//NL
   if (allocated(self%crk)) &
   desc = desc//mpih%myrankstr//'  crk:                             '//trim(str(self%crk                ))//NL
   if (allocated(self%alph)) then
   do s=1, self%nrk
   desc = desc//mpih%myrankstr//'  alph('//trim(str(s,.true.))//'): '//trim(str(self%alph(:,s)          ))//NL
   enddo
   endif
   if (allocated(self%beta)) &
   desc = desc//mpih%myrankstr//'  beta:                            '//trim(str(self%beta               ))//NL
   if (allocated(self%gamm)) &
   desc = desc//mpih%myrankstr//'  gamm:                            '//trim(str(self%gamm               ))//NL
   desc = desc//mpih%myrankstr//'  nrk:                             '//trim(str(self%nrk                ))
   endfunction description

   subroutine initialize(self, field, grid, file_parameters, scheme)
   !< Initialize class.
   class(rk_object),   intent(inout)        :: self            !< RK object.
   type(field_object), intent(in)           :: field           !< Field (sibling realm component, threaded in).
   type(grid_object),  intent(in), target   :: grid            !< Grid (sibling realm component, threaded in).
   type(file_ini),     intent(in), optional :: file_parameters !< Simulation parameters ini file handler.
   character(*),       intent(in), optional :: scheme          !< Runge-Kutta scheme.
   real(R8P)                                :: w0, w1          !< Sympletic RK coefficients.

   call mpih%print_message('rk_object%initialize start')
   call associate_adam_data
   if (present(file_parameters)) then
      call self%load_from_file(file_parameters=file_parameters)
   elseif (present(scheme)) then
      self%scheme = trim(adjustl(scheme))
   else
      call mpih%error_stop(msg=': failed to initialize rk object, one between file parameters and scheme must be passed')
   endif
   select case(self%scheme)
   case(RK_1) ! 1 stage, 1st order, Euler
      self%nrk = 1
      allocate(self%ark(self%nrk), self%brk(self%nrk), self%crk(self%nrk))
      self%ark(1) = 1._R8P ; self%brk(1) = 0._R8P ; self%crk(1) = 1._R8P
   case(RK_2) ! 2 stages, low storage, 2nd order TVD
      self%nrk = 2
      allocate(self%ark(self%nrk), self%brk(self%nrk), self%crk(self%nrk))
      self%ark(1) = 1._R8P  ; self%brk(1) = 0._R8P  ; self%crk(1) = 1._R8P
      self%ark(2) = 0.5_R8P ; self%brk(2) = 0.5_R8P ; self%crk(2) = 0.5_R8P
   case(RK_3) ! 3 stages, low storage, 3rd order TVD
      self%nrk = 3
      allocate(self%ark(self%nrk), self%brk(self%nrk), self%crk(self%nrk))
      self%ark(1) = 1._R8P        ; self%brk(1) = 0._R8P        ; self%crk(1) = 1._R8P
      self%ark(2) = 0.75_R8P      ; self%brk(2) = 0.25_R8P      ; self%crk(2) = 0.25_R8P
      self%ark(3) = 1._R8P/3._R8P ; self%brk(3) = 2._R8P/3._R8P ; self%crk(3) = 2._R8P/3._R8P
   case(RK_SSP_11) ! 1 stage, 1st order SSP (forward Euler): the 1-stage debugging
                   ! instrument of the staged forest path (issue #25); alph = Butcher A = 0,
                   ! beta = Butcher b = 1, gamm = Butcher c = 0.
      self%nrk = 1
      allocate(self%alph(self%nrk,self%nrk), self%beta(self%nrk), self%gamm(self%nrk))
      self%alph = 0._R8P
      self%beta = 0._R8P
      self%gamm = 0._R8P

      self%beta(1) = 1._R8P
   case(RK_SSP_22) ! 2 stages, 2nd order SSP
      self%nrk = 2
      allocate(self%alph(self%nrk,self%nrk), self%beta(self%nrk), self%gamm(self%nrk))
      self%alph = 0._R8P
      self%beta = 0._R8P
      self%gamm = 0._R8P

      self%beta(1) = 0.5_R8P
      self%beta(2) = 0.5_R8P

      self%alph(2,1) = 1._R8P

      self%gamm(2) = 1._R8P
   case(RK_SSP_33) ! 3 stages, 3rd order SSP
      self%nrk = 3 ! was 2 (issue #25 finding): beta(3)/alph(3,:)/gamm(3) below wrote OUT OF
                   ! BOUNDS of nrk=2 allocations -- silent heap corruption for any SSP-33 user
      allocate(self%alph(self%nrk,self%nrk), self%beta(self%nrk), self%gamm(self%nrk))
      self%alph = 0._R8P
      self%beta = 0._R8P
      self%gamm = 0._R8P

      self%beta(1) = 1._R8P/6._R8P
      self%beta(2) = 1._R8P/6._R8P
      self%beta(3) = 2._R8P/3._R8P

      self%alph(2,1) = 1._R8P
      self%alph(3,1) = 0.25_R8P ; self%alph(3,2) = 0.25_R8P

      self%gamm(2) = 1._R8P
      self%gamm(3) = 0.5_R8P
   case(RK_SSP_54) ! 5 stages, 4th order SSP
      self%nrk = 5
      allocate(self%alph(self%nrk,self%nrk), self%beta(self%nrk), self%gamm(self%nrk))
      self%alph = 0._R8P
      self%beta = 0._R8P
      self%gamm = 0._R8P

      self%beta(1) = 0.14681187618661_R8P
      self%beta(2) = 0.24848290924556_R8P
      self%beta(3) = 0.10425883036650_R8P
      self%beta(4) = 0.27443890091960_R8P
      self%beta(5) = 0.22600748319395_R8P

      self%alph(2,1) = 0.39175222700392_R8P
      self%alph(3,1) = 0.21766909633821_R8P ; self%alph(3,2) = 0.36841059262959_R8P
      self%alph(4,1) = 0.08269208670950_R8P ; self%alph(4,2) = 0.13995850206999_R8P ; self%alph(4,3) = 0.25189177424738_R8P
      self%alph(5,1) = 0.06796628370320_R8P ; self%alph(5,2) = 0.11503469844438_R8P ; self%alph(5,3) = 0.20703489864929_R8P
      self%alph(5,4) = 0.54497475021237_R8P

      self%gamm(2) = 0.39175222700392_R8P
      self%gamm(3) = 0.58607968896780_R8P
      self%gamm(4) = 0.47454236302687_R8P
      self%gamm(5) = 0.93501063100924_R8P
   case(RK_YOSHIDA)
      self%nrk = 4
      w0 = -1.702414383919315_R8P
      w1 =  1.351207191959658_R8P
      allocate(self%ssa(self%nrk), self%ssb(self%nrk-1))
      self%ssa = [w1/2.0_R8P,(w0+w1)/2.0_R8P,(w0+w1)/2.0_R8P,w1/2.0_R8P]
      self%ssb = [w1,w0,w1]
   case default
      ! issue #25: an unknown scheme used to fall through SILENTLY (default nrk, no
      ! coefficient arrays) and crash far from the cause; fail fast at initialization.
      call mpih%error_stop(msg=': unknown Runge-Kutta scheme "'//trim(adjustl(self%scheme))//'"')
   endselect

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>field%nv, nb=>field%nb, nrk=>self%nrk)
   select case(self%scheme)
   case(RK_1, RK_2, RK_3) ! low storage, only stage 1 is necessary
      call allocate_variable(var=self%q_rk,         &
                             ulb=reshape([1,nv,         &
                                          1-ngc,ni+ngc, &
                                          1-ngc,nj+ngc, &
                                          1-ngc,nk+ngc, &
                                          1,nb,         &
                                          1,1],[2,6]),  &
                             msg=mpih%myrankstr//'rk_object%initialize allocate q_rk')
   case(RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54)
      call allocate_variable(var=self%q_rk,              &
                             ulb=reshape([1,nv,          &
                                          1-ngc,ni+ngc,  &
                                          1-ngc,nj+ngc,  &
                                          1-ngc,nk+ngc,  &
                                          1,nb,          &
                                          1,nrk],[2,6]), &
                             msg=mpih%myrankstr//'rk_object%initialize allocate q_rk')
   endselect
   endassociate
   print '(A)', self%description()
   call mpih%print_message('rk_object%initialize finish')
   contains
      subroutine associate_adam_data
      !< Associate grid data pointers for easy handling.

      self%ni  => grid%ni
      self%nj  => grid%nj
      self%nk  => grid%nk
      self%ngc => grid%ngc
      endsubroutine associate_adam_data
   endsubroutine initialize

   subroutine initialize_stages(self, field, q)
   !< Initialize RK stages.
   class(rk_object), intent(inout) :: self             !< RK object.
   type(field_object), intent(in)  :: field            !< Field (sibling realm component, threaded in).
   real(R8P),        intent(in)    :: q(1:,          &
                                        1-self%ngc:, &
                                        1-self%ngc:, &
                                        1-self%ngc:, &
                                        1:)            !< Conservative variables.
   integer(I4P)                    :: i, j, k, b, v, s !< Counter.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>field%nv, blocks_number=>field%blocks_number)
   !$omp parallel do collapse(6) default(firstprivate) shared(q,self)
   do s=lbound(self%q_rk,dim=6),ubound(self%q_rk,dim=6)
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do v=1, nv
                     self%q_rk(v,i,j,k,b,s) = q(v,i,j,k,b)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endassociate
   endsubroutine initialize_stages

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(rk_object), intent(inout)        :: self            !< RK object.
   type(file_ini),   intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,          intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                :: go_on_fail_     !< Go on if load fails.
   character(99)                          :: buff_c          !< Character buffer.
   integer(I4P)                           :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='scheme', val=buff_c, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(scheme)')
   self%scheme = trim(adjustl(buff_c))
   endsubroutine load_from_file

   subroutine update_q(self, field, dt, phi, q, dq)
   !< Update RK q.
   class(rk_object), intent(in)              :: self                 !< RK object.
   type(field_object), intent(in)            :: field                !< Field (sibling realm component, threaded in).
   real(R8P),        intent(in)              :: dt                   !< Current time step.
   real(R8P),        intent(in), optional    :: phi(1:,          &
                                                    1-self%ngc:, &
                                                    1-self%ngc:, &
                                                    1-self%ngc:, &
                                                    1:)              !< IB distance.
   real(R8P),        intent(inout)           :: q(1:,          &
                                               1-self%ngc:, &
                                               1-self%ngc:, &
                                               1-self%ngc:, &
                                               1:)                   !< Conservative variables.
   real(R8P),        intent(inout), optional :: dq(1:,          &
                                                   1-self%ngc:, &
                                                   1-self%ngc:, &
                                                   1-self%ngc:, &
                                                   1:)               !< Conservative variables residuals.
   integer(I4P)                              :: all_solids           !< Last phi index, all solids summary.
   integer(I4P)                              :: i, j, k, b, v, s     !< Counter.

   associate(nrk=>self%nrk, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>field%nv, blocks_number=>field%blocks_number)
   if (present(phi)) then
      all_solids = ubound(phi, dim=1)
      !$omp parallel do collapse(6) default(firstprivate) shared(phi,q,self)
      do s=1, nrk
         do b=1, blocks_number
            do k=1, nk
               do j=1, nj
                  do i=1, ni
                     do v=1, nv
                        if (phi(all_solids,i,j,k,b) < 0._R8P) then
                           q(v,i,j,k,b) = q(v,i,j,k,b) + dt * self%beta(s) * self%q_rk(v,i,j,k,b,s)
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   else
      !$omp parallel do collapse(6) default(firstprivate) shared(q,self)
      if (present(dq)) dq = 0._R8P
      do s=1, nrk
         do b=1, blocks_number
            do k=1, nk
               do j=1, nj
                  do i=1, ni
                     do v=1, nv
                        if (present(dq)) dq(v,i,j,k,b) = dq(v,i,j,k,b) + self%beta(s)*self%q_rk(v,i,j,k,b,s)
                        q(v,i,j,k,b) = q(v,i,j,k,b) + dt * self%beta(s) * self%q_rk(v,i,j,k,b,s)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   endif
   endassociate
   endsubroutine update_q
endmodule adam_rk_object
