!< ADAM, Runge Kutta class definition, CPU backend.
module adam_rk_cpu_object
!< ADAM, Runge Kutta class definition, CPU backend.

use adam_field_object, only : field_object
use adam_grid_object,  only : grid_object
use adam_mpih_object,  only : mpih_object
use finer
use penf
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none
private
public :: rk_cpu_object

character(len=7), parameter :: INI_SECTION_NAME='schemes' !< INI (config) file section name containing RK configs.

character(len=11), parameter :: RK_SUPPORTED_SCHEMES(7)=['rk-ls-qn-11', &
                                                         'rk-ls-qn-22', &
                                                         'rk-ls-qn-33', &
                                                         'rk-ls-qs-54', &
                                                         'rk-ls-qs-64', &
                                                         'rk-ls-qs-74', &
                                                         'rk-ns-qs-54'] !< List of supported schemes.
integer(I4P),      parameter :: RK_NRK_SUPPORTED_SCHEMES(7)=[1, &
                                                             2, &
                                                             3, &
                                                             5, &
                                                             6, &
                                                             7, &
                                                             5]         !< List RK stages number of supported schemes.

type :: rk_cpu_object
   !< Runge Kutta class definition, CPU backend.
   type(mpih_object)         :: mpih      !< MPI handler.
   character(:), allocatable :: scheme    !< RK scheme, ['rk-11', 'rk-33', 'rk-43', 'rk-54']
   integer(I4P)              :: nrk=1_I4P !< Runge-Kutta stages number.
   ! low storage coefficients
   real(R8P), allocatable    :: ark(:) !< RK a coefficients.
   real(R8P), allocatable    :: brk(:) !< RK b coefficients.
   real(R8P), allocatable    :: crk(:) !< RK c coefficients.
   ! normal storage coefficients
   real(R8P), allocatable    :: nark(:,:) !< RK a coefficients.
   real(R8P), allocatable    :: nbrk(:)   !< RK b coefficients.
   real(R8P), allocatable    :: ncrk(:)   !< RK c coefficients.
   ! large arrays
   real(R8P), allocatable    :: q_n(:,:,:,:,:)   !< Field at t(ns), previous stage, low storage algorithms.
   real(R8P), allocatable    :: q_s(:,:,:,:,:,:) !< Field at t(ns), previous stage, normal storage algorithms.
   contains
      ! public methods
      procedure, pass(self) :: compute_stage    !< Compute RK stage.
      procedure, pass(self) :: description      !< Return pretty-printed object description.
      procedure, pass(self) :: initialize       !< Initialize RK.
      procedure, pass(self) :: initialize_stage !< Initialize RK stage.
      procedure, pass(self) :: is_supported     !< Return true if input scheme is supported.
      procedure, pass(self) :: load_from_file   !< Load config from file.
      procedure, pass(self) :: sum_stages       !< Sum stages to q for computing q(tn+1). Used only by normal storage schemes.
      ! private methods
      procedure, pass(self), private :: compute_stage_ls_phi_q_n !< Compute RK (low storage) stage with IB, q_n=q(tn).
      procedure, pass(self), private :: compute_stage_ls_q_n     !< Compute RK (low storage) stage without IB, q_n=q(tn).
      procedure, pass(self), private :: compute_stage_ls_phi_q_s !< Compute RK (low storage) stage with IB, q_n=q(s).
      procedure, pass(self), private :: compute_stage_ls_q_s     !< Compute RK (low storage) stage without IB, q_n=q(s).
      ! procedure, pass(self), private :: compute_stage_ns_phi_q_s !< Compute RK (normal storage) stage with IB, q_n=q(s).
      procedure, pass(self), private :: compute_stage_ns_q_s     !< Compute RK (normal storage) stage without IB, q_n=q(s).
endtype rk_cpu_object

contains
   ! public methods
   subroutine compute_stage(self, ni, nj, nk, ngc, nv, blocks_number, dt, s, phi, rq, q)
   !< Compute RK stage.
   class(rk_cpu_object), intent(inout)           :: self                            !< RK.
   integer(I4P),         intent(in)              :: ni                              !< Grid cells number in I direction.
   integer(I4P),         intent(in)              :: nj                              !< Grid cells number in J direction.
   integer(I4P),         intent(in)              :: nk                              !< Grid cells number in K direction.
   integer(I4P),         intent(in)              :: ngc                             !< Ghost cells number.
   integer(I4P),         intent(in)              :: nv                              !< Number of conservative variables.
   integer(I4P),         intent(in)              :: blocks_number                   !< Number of blocks.
   real(R8P),            intent(in)              :: dt                              !< Time step.
   integer(I4P),         intent(in)              :: s                               !< RK stage.
   real(R8P),            intent(in), allocatable :: phi( :,     :,     :,     :, :) !< IB distance function.
   real(R8P),            intent(in)              ::  rq(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Residuals of previous stage.
   real(R8P),            intent(inout)           ::   q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.

   select case(self%scheme(1:5))
   case('rk-ls') ! low storage schemes
      select case(self%scheme(7:8))
      case('qn') ! scheme with q_n=q(tn) fixed
         if (allocated(phi)) then ! with IB
            call self%compute_stage_ls_phi_q_n(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                               dt=dt, s=s, phi=phi, rq=rq, q=q)
         else                     ! without IB
            call self%compute_stage_ls_q_n(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                           dt=dt, s=s, rq=rq, q=q)
         endif
      case('qs') ! scheme with q_n=q(s) evolved with q
         if (allocated(phi)) then ! with IB
            call self%compute_stage_ls_phi_q_s(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                               dt=dt, s=s, phi=phi, rq=rq, q=q)
         else                     ! without IB
            call self%compute_stage_ls_q_s(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                           dt=dt, s=s, rq=rq, q=q)
         endif
      endselect
   case('rk-ns') ! normal storage schemes
      if (allocated(phi)) then ! with IB
         ! call self%compute_stage_phi_ns_q_s(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                            ! dt=dt, s=s, phi=phi, rq=rq, q=q)
      else                     ! without IB
         call self%compute_stage_ns_q_s(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, dt=dt, s=s, q=q)
      endif
   endselect
   endsubroutine compute_stage

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(rk_cpu_object), intent(in) :: self             !< RK.
   character(len=:), allocatable    :: desc             !< Description.
   character(len=1), parameter      :: NL=new_line('a') !< New line character.

   desc =       self%mpih%myrankstr//'RK main data'//NL
   desc = desc//self%mpih%myrankstr//'  scheme:            '//trim(adjustl(self%scheme))//NL
   desc = desc//self%mpih%myrankstr//'  RK stages number:  '//trim(str(self%nrk)       )//NL
   select case(self%scheme(1:5))
   case('rk-ls') ! low storage schemes
      desc = desc//self%mpih%myrankstr//'  RK A coefficients: '//trim(str(self%ark))//NL
      desc = desc//self%mpih%myrankstr//'  RK B coefficients: '//trim(str(self%brk))//NL
      desc = desc//self%mpih%myrankstr//'  RK C coefficients: '//trim(str(self%crk))
   case('rk-ns') ! normal storage schemes
      desc = desc//self%mpih%myrankstr//'  RK A coefficients: '//trim(str(reshape(self%nark,[self%nrk*self%nrk])))//NL
      desc = desc//self%mpih%myrankstr//'  RK B coefficients: '//trim(str(self%nbrk))//NL
      desc = desc//self%mpih%myrankstr//'  RK C coefficients: '//trim(str(self%ncrk))
   endselect
   endfunction description

   subroutine initialize(self, file_parameters, grid, field)
   !< Initialize RK.
   class(rk_cpu_object), intent(inout) :: self            !< RK.
   type(file_ini),       intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   type(grid_object),    intent(in)    :: grid            !< The grid.
   type(field_object),   intent(in)    :: field           !< The field.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'rk_cpu_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   select case(trim(adjustl(self%scheme)))
   case('rk-ls-qn-11') ! low storage euler 1st order
      allocate(self%ark(self%nrk), self%brk(self%nrk), self%crk(self%nrk))
      self%ark(1) = 1._R8P ; self%brk(1) = 0._R8P ; self%crk(1) = 1._R8P
      allocate(self%q_n(1:field%nv,1-grid%ngc:grid%ni+grid%ngc,1-grid%ngc:grid%nj+grid%ngc,1-grid%ngc:grid%nk+grid%ngc,1:field%nb))
   case('rk-ls-qn-22') ! low storage 2nd order TVD
      allocate(self%ark(self%nrk), self%brk(self%nrk), self%crk(self%nrk))
      self%ark(1) = 1._R8P  ; self%brk(1) = 0._R8P  ; self%crk(1) = 1._R8P
      self%ark(2) = 0.5_R8P ; self%brk(2) = 0.5_R8P ; self%crk(2) = 0.5_R8P
      allocate(self%q_n(1:field%nv,1-grid%ngc:grid%ni+grid%ngc,1-grid%ngc:grid%nj+grid%ngc,1-grid%ngc:grid%nk+grid%ngc,1:field%nb))
   case('rk-ls-qn-33') ! low storage 3rd order TVD
      allocate(self%ark(self%nrk), self%brk(self%nrk), self%crk(self%nrk))
      self%ark(1) = 1._R8P        ; self%brk(1) = 0._R8P        ; self%crk(1) = 1._R8P
      self%ark(2) = 0.75_R8P      ; self%brk(2) = 0.25_R8P      ; self%crk(2) = 0.25_R8P
      self%ark(3) = 1._R8P/3._R8P ; self%brk(3) = 2._R8P/3._R8P ; self%crk(3) = 2._R8P/3._R8P
      allocate(self%q_n(1:field%nv,1-grid%ngc:grid%ni+grid%ngc,1-grid%ngc:grid%nj+grid%ngc,1-grid%ngc:grid%nk+grid%ngc,1:field%nb))
   case('rk-ls-qs-54') ! low storage 4th order SSP 5 stages
      allocate(self%ark(self%nrk), self%brk(self%nrk), self%crk(self%nrk))
      self%ark(1) =  0._R8P
      self%ark(2) = -real(567301805773_I8P,  kind=R8P) / real(1357537059087_I8P, kind=R8P)
      self%ark(3) = -real(2404267990393_I8P, kind=R8P) / real(2016746695238_I8P, kind=R8P)
      self%ark(4) = -real(3550918686646_I8P, kind=R8P) / real(2091501179385_I8P, kind=R8P)
      self%ark(5) = -real(1275806237668_I8P, kind=R8P) / real(842570457699_I8P,  kind=R8P)

      self%brk(1) =  real(1432997174477_I8P, kind=R8P) / real(9575080441755_I8P,  kind=R8P)
      self%brk(2) =  real(5161836677717_I8P, kind=R8P) / real(13612068292357_I8P, kind=R8P)
      self%brk(3) =  real(1720146321549_I8P, kind=R8P) / real(2090206949498_I8P,  kind=R8P)
      self%brk(4) =  real(3134564353537_I8P, kind=R8P) / real(4481467310338_I8P,  kind=R8P)
      self%brk(5) =  real(2277821191437_I8P, kind=R8P) / real(14882151754819_I8P, kind=R8P)

      self%crk(1) =  0._R8P
      self%crk(2) =  real(1432997174477_I8P, kind=R8P) / real(9575080441755_I8P, kind=R8P)
      self%crk(3) =  real(2526269341429_I8P, kind=R8P) / real(6820363962896_I8P, kind=R8P)
      self%crk(4) =  real(2006345519317_I8P, kind=R8P) / real(3224310063776_I8P, kind=R8P)
      self%crk(5) =  real(2802321613138_I8P, kind=R8P) / real(2924317926251_I8P, kind=R8P)
      allocate(self%q_n(1:field%nv,1-grid%ngc:grid%ni+grid%ngc,1-grid%ngc:grid%nj+grid%ngc,1-grid%ngc:grid%nk+grid%ngc,1:field%nb))
   case('rk-ls-qs-64') ! low storage 4th order 6 stages
      allocate(self%ark(self%nrk), self%brk(self%nrk), self%crk(self%nrk))
      self%ark(1) =  0._R8P             ; self%brk(1) = 0.122000000000_R8P ; self%crk(1) = 0._R8P
      self%ark(2) = -0.691750960670_R8P ; self%brk(2) = 0.477263056358_R8P ; self%crk(2) = 0.122000000000_R8P
      self%ark(3) = -1.727127405211_R8P ; self%brk(3) = 0.381941220320_R8P ; self%crk(3) = 0.269115878630_R8P
      self%ark(4) = -0.694890150986_R8P ; self%brk(4) = 0.447757195744_R8P ; self%crk(4) = 0.447717183551_R8P
      self%ark(5) = -1.039942756197_R8P ; self%brk(5) = 0.498614246822_R8P ; self%crk(5) = 0.749979795490_R8P
      self%ark(6) = -1.531977447611_R8P ; self%brk(6) = 0.186648570846_R8P ; self%crk(6) = 0.898555413085_R8P
      allocate(self%q_n(1:field%nv,1-grid%ngc:grid%ni+grid%ngc,1-grid%ngc:grid%nj+grid%ngc,1-grid%ngc:grid%nk+grid%ngc,1:field%nb))
   case('rk-ls-qs-74') ! low storage 4th order 7 stages
      allocate(self%ark(self%nrk), self%brk(self%nrk), self%crk(self%nrk))
      self%ark(1) =  0._R8P             ; self%brk(1) = 0.117322146869_R8P ; self%crk(1) = 0._R8P
      self%ark(2) = -0.647900745934_R8P ; self%brk(2) = 0.503270262127_R8P ; self%crk(2) = 0.117322146869_R8P
      self%ark(3) = -2.704760863204_R8P ; self%brk(3) = 0.233663281658_R8P ; self%crk(3) = 0.294523230758_R8P
      self%ark(4) = -0.460080550118_R8P ; self%brk(4) = 0.283419634625_R8P ; self%crk(4) = 0.305658622131_R8P
      self%ark(5) = -0.500581787785_R8P ; self%brk(5) = 0.540367414023_R8P ; self%crk(5) = 0.582864148403_R8P
      self%ark(6) = -1.906532255913_R8P ; self%brk(6) = 0.371499414620_R8P ; self%crk(6) = 0.858664273599_R8P
      self%ark(7) = -1.450000000000_R8P ; self%brk(7) = 0.136670099385_R8P ; self%crk(7) = 0.868664273599_R8P
      allocate(self%q_n(1:field%nv,1-grid%ngc:grid%ni+grid%ngc,1-grid%ngc:grid%nj+grid%ngc,1-grid%ngc:grid%nk+grid%ngc,1:field%nb))
   case('rk-ns-qs-54') ! normal storage 4th order SSP 5 stages
      allocate(self%nbrk(1:self%nrk            )) ; self%nbrk = 0._R_P
      allocate(self%nark(1:self%nrk, 1:self%nrk)) ; self%nark = 0._R_P
      allocate(self%ncrk(            1:self%nrk)) ; self%ncrk = 0._R_P

      self%nbrk(1) = 0.14681187618661_R8P
      self%nbrk(2) = 0.24848290924556_R8P
      self%nbrk(3) = 0.10425883036650_R8P
      self%nbrk(4) = 0.27443890091960_R8P
      self%nbrk(5) = 0.22600748319395_R8P

      self%nark(2,1) = 0.39175222700392_R8P
      self%nark(3,1) = 0.21766909633821_R8P ; self%nark(3,2) = 0.36841059262959_R8P
      self%nark(4,1) = 0.08269208670950_R8P ; self%nark(4,2) = 0.13995850206999_R8P ; self%nark(4,3) = 0.25189177424738_R8P
      self%nark(5,1) = 0.06796628370320_R8P ; self%nark(5,2) = 0.11503469844438_R8P ; self%nark(5,3) = 0.20703489864929_R8P
      self%nark(5,4) = 0.54497475021237_R8P

      self%ncrk(2) = 0.39175222700392_R8P
      self%ncrk(3) = 0.58607968896780_R8P
      self%ncrk(4) = 0.47454236302687_R8P
      self%ncrk(5) = 0.93501063100924_R8P
      allocate(self%q_s(1:field%nv,&
                        1-grid%ngc:grid%ni+grid%ngc,1-grid%ngc:grid%nj+grid%ngc,1-grid%ngc:grid%nk+grid%ngc,1:field%nb,1:self%nrk))
   endselect
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'rk_cpu_object%initialize finish'
   endsubroutine initialize

   subroutine initialize_stage(self, ngc, q)
   !< Initialize RK stage.
   class(rk_cpu_object), intent(inout) :: self                          !< RK.
   integer(I4P),         intent(in)    :: ngc                           !< Number of ghost cells.
   real(R8P),            intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field at previoust time step.
   integer(I4P)                        :: s                             !< Counter.

   select case(self%scheme(1:5))
   case('rk-ls') ! low storage schemes
      select case(self%scheme(7:8))
      case('qn')
         self%q_n = q ! store solution at time t(n)
      case('qs')
         self%q_n = 0._R8P
      endselect
   case('rk-ns') ! normal storage schemes
      do s=1, self%nrk
         self%q_s(:,:,:,:,:,s) = q ! initialize with solution at time t(n)
      enddo
   endselect
   endsubroutine initialize_stage

   function is_supported(self, scheme, nrk)
   !< Return true if input scheme is supported.
   class(rk_cpu_object), intent(in)            :: self         !< RK.
   character(*),         intent(in)            :: scheme       !< Selected scheme.
   integer(I4P),         intent(out), optional :: nrk          !< Number of RK stages of the input scheme.
   logical                                     :: is_supported !< Inquire result.
   integer(I4P)                                :: s            !< Counter.

   is_supported = .false.
   do s=lbound(RK_SUPPORTED_SCHEMES, dim=1), ubound(RK_SUPPORTED_SCHEMES, dim=1)
     if (trim(adjustl(scheme)) == trim(adjustl(RK_SUPPORTED_SCHEMES(s)))) then
       is_supported = .true.
       if (present(nrk)) nrk = RK_NRK_SUPPORTED_SCHEMES(s)
       return
     endif
   enddo
   endfunction is_supported

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(rk_cpu_object), intent(inout)        :: self            !< RK.
   type(file_ini),       intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,              intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                    :: go_on_fail_     !< Go on if load fails.
   character(99)                              :: char_buff       !< Character buffer.
   integer(I4P)                               :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='runge_kutta', val=char_buff, error=error)
   ! if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//INI_SECTION_NAME//'].(runge_kutta)'
   self%scheme = trim(adjustl(char_buff))
   if (.not.self%is_supported(scheme=trim(adjustl(self%scheme)), nrk=self%nrk)) then ! check if scheme is supported and assign nrk
      write(stderr, '(A)') self%mpih%myrankstr//' error: RK schemes "'//trim(adjustl(self%scheme)) //'" is unknown'
      call self%mpih%abort
      stop
   endif
   endsubroutine load_from_file

   subroutine sum_stages(self, ni, nj, nk, ngc, blocks_number, dt, q)
   !< Sum stages to q for computing q(tn+1). Used only by normal storage schemes.
   class(rk_cpu_object), intent(in)    :: self                          !< RK.
   integer(I4P),         intent(in)    :: ni                            !< Grid cells number in I direction.
   integer(I4P),         intent(in)    :: nj                            !< Grid cells number in J direction.
   integer(I4P),         intent(in)    :: nk                            !< Grid cells number in K direction.
   integer(I4P),         intent(in)    :: ngc                           !< Ghost cells number.
   integer(I4P),         intent(in)    :: blocks_number                 !< Number of blocks.
   real(R8P),            intent(in)    :: dt                            !< Time step.
   real(R8P),            intent(inout) :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field at previoust time step.
   integer(I4P)                        :: i, j, k, b, s                 !< Counter.

   do s=1, self%nrk
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  q(:,i,j,k,b) = q(:,i,j,k,b) + self%q_s(:,i,j,k,b,s) * dt * self%nbrk(s)
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine sum_stages

   ! private methods
   subroutine compute_stage_ls_phi_q_n(self, ni, nj, nk, ngc, nv, blocks_number, dt, s, phi, rq, q)
   !< Compute RK (low storage) stage with IB, q_n being the q(tn) time step not evolved within q.
   class(rk_cpu_object), intent(inout) :: self                            !< RK.
   integer(I4P),         intent(in)    :: ni                              !< Grid cells number in I direction.
   integer(I4P),         intent(in)    :: nj                              !< Grid cells number in J direction.
   integer(I4P),         intent(in)    :: nk                              !< Grid cells number in K direction.
   integer(I4P),         intent(in)    :: ngc                             !< Ghost cells number.
   integer(I4P),         intent(in)    :: nv                              !< Number of conservative variables.
   integer(I4P),         intent(in)    :: blocks_number                   !< Number of blocks.
   real(R8P),            intent(in)    :: dt                              !< Time step.
   integer(I4P),         intent(in)    :: s                               !< RK stage.
   real(R8P),            intent(in)    :: phi(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< IB distance function.
   real(R8P),            intent(in)    ::  rq(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Residuals of previous stage.
   real(R8P),            intent(inout) ::   q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.
   real(R8P)                           :: dtcrk                           !< Dt times crk.
   integer(I4P)                        :: i, j, k, b                      !< Counter.

   associate(ark=>self%ark(s), brk=>self%brk(s), crk=>self%crk(s))
   dtcrk = dt * crk
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               if (all(phi(:,i,j,k,b) < 0._R8P)) then
                  q(:,i,j,k,b) = ark * self%q_n(:,i,j,k,b) + brk * q(:,i,j,k,b) + dtcrk * rq(:,i,j,k,b)
               endif
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine compute_stage_ls_phi_q_n

   subroutine compute_stage_ls_q_n(self, ni, nj, nk, ngc, nv, blocks_number, dt, s, rq, q)
   !< Compute RK (low storage) stage without IB, q_n being the q(tn) time step not evolved within q.
   class(rk_cpu_object), intent(inout) :: self                            !< RK.
   integer(I4P),         intent(in)    :: ni                              !< Grid cells number in I direction.
   integer(I4P),         intent(in)    :: nj                              !< Grid cells number in J direction.
   integer(I4P),         intent(in)    :: nk                              !< Grid cells number in K direction.
   integer(I4P),         intent(in)    :: ngc                             !< Ghost cells number.
   integer(I4P),         intent(in)    :: nv                              !< Number of conservative variables.
   integer(I4P),         intent(in)    :: blocks_number                   !< Number of blocks.
   real(R8P),            intent(in)    :: dt                              !< Time step.
   integer(I4P),         intent(in)    :: s                               !< RK stage.
   real(R8P),            intent(in)    ::  rq(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Residuals of previous stage.
   real(R8P),            intent(inout) ::   q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.
   real(R8P)                           :: dtcrk                           !< Dt times crk.
   integer(I4P)                        :: i, j, k, b                      !< Counter.

   associate(ark=>self%ark(s), brk=>self%brk(s), crk=>self%crk(s))
   dtcrk = dt * crk
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               q(:,i,j,k,b) = ark * self%q_n(:,i,j,k,b) + brk * q(:,i,j,k,b) + dtcrk * rq(:,i,j,k,b)
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine compute_stage_ls_q_n

   subroutine compute_stage_ls_phi_q_s(self, ni, nj, nk, ngc, nv, blocks_number, dt, s, phi, rq, q)
   !< Compute RK (low storage) stage with IB, q_n being the s-th stage evolved within q.
   class(rk_cpu_object), intent(inout) :: self                            !< RK.
   integer(I4P),         intent(in)    :: ni                              !< Grid cells number in I direction.
   integer(I4P),         intent(in)    :: nj                              !< Grid cells number in J direction.
   integer(I4P),         intent(in)    :: nk                              !< Grid cells number in K direction.
   integer(I4P),         intent(in)    :: ngc                             !< Ghost cells number.
   integer(I4P),         intent(in)    :: nv                              !< Number of conservative variables.
   integer(I4P),         intent(in)    :: blocks_number                   !< Number of blocks.
   real(R8P),            intent(in)    :: dt                              !< Time step.
   integer(I4P),         intent(in)    :: s                               !< RK stage.
   real(R8P),            intent(in)    :: phi(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< IB distance function.
   real(R8P),            intent(in)    ::  rq(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Residuals of previous stage.
   real(R8P),            intent(inout) ::   q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.
   integer(I4P)                        :: i, j, k, b                      !< Counter.

   associate(ark=>self%ark(s), brk=>self%brk(s), crk=>self%crk(s))
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               if (all(phi(:,i,j,k,b) < 0._R8P)) then
                  self%q_n(:,i,j,k,b) = self%q_n(:,i,j,k,b) * ark + rq(:,i,j,k,b) * dt
                         q(:,i,j,k,b) =        q(:,i,j,k,b) + self%q_n(:,i,j,k,b) * brk
               endif
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine compute_stage_ls_phi_q_s

   subroutine compute_stage_ls_q_s(self, ni, nj, nk, ngc, nv, blocks_number, dt, s, rq, q)
   !< Compute RK (low storage) stage without IB, q_n being the s-th stage evolved within q.
   class(rk_cpu_object), intent(inout) :: self                            !< RK.
   integer(I4P),         intent(in)    :: ni                              !< Grid cells number in I direction.
   integer(I4P),         intent(in)    :: nj                              !< Grid cells number in J direction.
   integer(I4P),         intent(in)    :: nk                              !< Grid cells number in K direction.
   integer(I4P),         intent(in)    :: ngc                             !< Ghost cells number.
   integer(I4P),         intent(in)    :: nv                              !< Number of conservative variables.
   integer(I4P),         intent(in)    :: blocks_number                   !< Number of blocks.
   real(R8P),            intent(in)    :: dt                              !< Time step.
   integer(I4P),         intent(in)    :: s                               !< RK stage.
   real(R8P),            intent(in)    ::  rq(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Residuals of previous stage.
   real(R8P),            intent(inout) ::   q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.
   integer(I4P)                        :: i, j, k, b                      !< Counter.

   associate(ark=>self%ark(s), brk=>self%brk(s), crk=>self%crk(s))
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               self%q_n(:,i,j,k,b) = self%q_n(:,i,j,k,b) * ark + rq(:,i,j,k,b) * dt
                      q(:,i,j,k,b) =        q(:,i,j,k,b) + self%q_n(:,i,j,k,b) * brk
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine compute_stage_ls_q_s

   subroutine compute_stage_ns_q_s(self, ni, nj, nk, ngc, nv, blocks_number, dt, s, q)
   !< Compute RK (normal storage) stage without IB, q_n being the s-th stage evolved within q.
   class(rk_cpu_object), intent(inout) :: self                          !< RK.
   integer(I4P),         intent(in)    :: ni                            !< Grid cells number in I direction.
   integer(I4P),         intent(in)    :: nj                            !< Grid cells number in J direction.
   integer(I4P),         intent(in)    :: nk                            !< Grid cells number in K direction.
   integer(I4P),         intent(in)    :: ngc                           !< Ghost cells number.
   integer(I4P),         intent(in)    :: nv                            !< Number of conservative variables.
   integer(I4P),         intent(in)    :: blocks_number                 !< Number of blocks.
   real(R8P),            intent(in)    :: dt                            !< Time step.
   integer(I4P),         intent(in)    :: s                             !< RK stage.
   real(R8P),            intent(inout) :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.
   integer(I4P)                        :: i, j, k, b, ss                !< Counter.

   associate(nark=>self%nark)
   do ss=1, s - 1
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  self%q_s(:,i,j,k,b,s) = self%q_s(:,i,j,k,b,s) + self%q_s(:,i,j,k,b,ss) * dt * nark(s, ss)
               enddo
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine compute_stage_ns_q_s
endmodule adam_rk_cpu_object
