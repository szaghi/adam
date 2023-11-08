!< ADAM, RK class definition.
module adam_rk_object
!< ADAM, RK class definition.

use adam_field_object
use adam_grid_object
use adam_mpih_object
use finer
use penf

implicit none
save
private
public :: rk_object
public :: RK_1
public :: RK_2
public :: RK_3
public :: RK_SSP_22
public :: RK_SSP_33
public :: RK_SSP_54

character(len=13), parameter :: RK_1     ="runge-kutta-1"      !< Parameter of time scheme, Runge-Kutta 1.
character(len=13), parameter :: RK_2     ="runge-kutta-2"      !< Parameter of time scheme, Runge-Kutta 2.
character(len=13), parameter :: RK_3     ="runge-kutta-3"      !< Parameter of time scheme, Runge-Kutta 3.
character(len=18), parameter :: RK_SSP_22="runge-kutta-ssp-22" !< Parameter of time scheme, Runge-Kutta SSP 22.
character(len=18), parameter :: RK_SSP_33="runge-kutta-ssp-33" !< Parameter of time scheme, Runge-Kutta SSP 33.
character(len=18), parameter :: RK_SSP_54="runge-kutta-ssp-54" !< Parameter of time scheme, Runge-Kutta SSP 54.

character(len=11), parameter :: INI_SECTION_NAME="runge_kutta" !< INI (config) file section name containing time configs.

type :: rk_object
   !< RK class definition.
   type(mpih_object)         :: mpih          !< MPI handler.
   character(:), allocatable :: scheme        !< RK scheme.
   integer(I4P)              :: nrk=3_I4P     !< Runge-Kutta stages number.
   real(R8P), allocatable    :: ark(:)        !< Runge-Kutta low storage alpha coefficients.
   real(R8P), allocatable    :: brk(:)        !< Runge-Kutta low storage beta coefficients.
   real(R8P), allocatable    :: crk(:)        !< Runge-Kutta low storage beta coefficients.
   real(R8P), allocatable    :: alph(:,:)     !< Runge-Kutta SSP alpha coefficients.
   real(R8P), allocatable    :: beta(:)       !< Runge-Kutta SSP beta coefficients.
   real(R8P), allocatable    :: gamm(:)       !< Runge-Kutta SSP gamma coefficients.
   ! grid/field data replica for easy handling
   type(field_object), pointer :: field=>null()         !< The field.
   type(grid_object),  pointer :: grid=>null()          !< The grid.
   integer(I4P),       pointer :: ngc=>null()           !< Number of ghost cells.
   integer(I4P),       pointer :: ni=>null()            !< Number of cells in i direction.
   integer(I4P),       pointer :: nj=>null()            !< Number of cells in j direction.
   integer(I4P),       pointer :: nk=>null()            !< Number of cells in k direction.
   integer(I4P),       pointer :: nb=>null()            !< Total blocks number for MPI.
   integer(I4P),       pointer :: blocks_number=>null() !< Actual blocks number.
   integer(I4P),       pointer :: ns=>null()            !< Number of fluids specie.
   integer(I4P),       pointer :: nv=>null()            !< Number of conservative variables.
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize class.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype rk_object
contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(rk_object), intent(in)  :: self             !< RK object.
   character(len=:), allocatable :: desc             !< Description.
   character(len=1), parameter   :: NL=new_line('a') !< New line character.
   integer(I4P)                  :: s                !< Counter.

   desc =       self%mpih%myrankstr//'Runge-Kutta scheme main data'//NL
   if (allocated(self%ark)) &
   desc = desc//self%mpih%myrankstr//'  ark:                             '//trim(str(self%ark                ))//NL
   if (allocated(self%brk)) &
   desc = desc//self%mpih%myrankstr//'  brk:                             '//trim(str(self%brk                ))//NL
   if (allocated(self%crk)) &
   desc = desc//self%mpih%myrankstr//'  crk:                             '//trim(str(self%crk                ))//NL
   if (allocated(self%alph)) then
   do s=1, self%nrk
   desc = desc//self%mpih%myrankstr//'  alph('//trim(str(s,.true.))//'): '//trim(str(self%alph(:,s)          ))//NL
   enddo
   endif
   if (allocated(self%beta)) &
   desc = desc//self%mpih%myrankstr//'  beta:                            '//trim(str(self%beta               ))//NL
   if (allocated(self%gamm)) &
   desc = desc//self%mpih%myrankstr//'  gamm:                            '//trim(str(self%gamm               ))//NL
   desc = desc//self%mpih%myrankstr//'  nrk:                             '//trim(str(self%nrk                ))
   endfunction description

   subroutine initialize(self, file_parameters, scheme, grid, field)
   !< Initialize class.
   class(rk_object),   intent(inout)        :: self            !< RK object.
   type(file_ini),     intent(in), optional :: file_parameters !< Simulation parameters ini file handler.
   character(*),       intent(in), optional :: scheme          !< Runge-Kutta scheme.
   type(grid_object),  intent(in), target   :: grid            !< The grid.
   type(field_object), intent(in), target   :: field           !< The field.

   call self%mpih%initialize(do_mpi_init=.false.)
   call self%mpih%print_message('rk_object%initialize start')
   call associate_adam_data(grid=grid, field=field)
   if (present(file_parameters)) then
      call self%load_from_file(file_parameters=file_parameters)
   elseif (present(scheme)) then
      self%scheme = trim(adjustl(scheme))
   else
      call self%mpih%error_stop(msg=': failed to initialize rk object, one between file parameters and scheme must be passed')
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
      self%nrk = 2
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
   endselect
   call self%mpih%print_message('rk_object%initialize finish')
   contains
      subroutine associate_adam_data(grid, field)
      !< Associate objects data to equation for easy handling.
      type(grid_object),          intent(in), target :: grid    !< The grid.
      type(field_object),         intent(in), target :: field   !< The field.

      self%grid          => grid
      self%field         => field
      self%blocks_number => field%blocks_number
      self%ni            => field%grid%ni
      self%nj            => field%grid%nj
      self%nk            => field%grid%nk
      self%ngc           => field%grid%ngc
      self%nb            => field%nb
      self%nv            => field%nv
      endsubroutine associate_adam_data
   endsubroutine initialize

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
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(scheme)')
   self%scheme = trim(adjustl(buff_c))
   endsubroutine load_from_file
endmodule adam_rk_object
