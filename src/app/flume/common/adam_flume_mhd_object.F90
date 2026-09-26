!< ADAM, FLUME MHD class definition: the `[mhd]` section.
module adam_flume_mhd_object
!< ADAM, FLUME MHD class definition: the `[mhd]` section.
!<
!< Loaded only when `[physics].(physical_model) = mhd-ideal` (issue #41, section 6.2). `divergence_control` selects
!< the model variant, hence the state width: `glm` (mixed GLM cleaning, `nv = 9`) or `none` (`nv = 8`). The GLM keys
!< are required with `glm` only: `glm_ch` is the constant, uniform cleaning speed; the damping is
!< `c_h^2 / c_p^2 = glm_alpha c_h / L`, `L = glm_damping_length`, a positive length or `min-cell` (the minimum cell
!< spacing of the realm over the active directions, Mignone & Tzeferacos 2010; set by `set_glm_damping` once the grid
!< exists, issue #41, section 3.5); `glm_ch_check` decides what happens when the fastest wave
!< outruns `c_h` (`warning` or `error`). `divb_tol` arms the div(B) monitor (0 disables it) and `divb_error` makes an
!< exceeded tolerance fatal. `rho_floor` and `p_floor` are the positivity floors (0 disables them).

! ADAM singleton objects
use :: adam_mpih_global,      only : mpih
! FLUME modules
use :: adam_flume_parameters, only : DIVERGENCE_CONTROL_GLM, DIVERGENCE_CONTROL_NONE, GLM_CH_CHECK_ERROR, &
                                     GLM_CH_CHECK_WARNING, GLM_DAMPING_LENGTH_MIN_CELL, strip_control
! third party modules
use :: finer,                 only : file_ini
use :: penf,                  only : I4P, R8P, str

implicit none
private
public :: flume_mhd_object

character(len=3), parameter :: INI_SECTION_NAME="mhd" !< INI (config) file section name containing MHD configs.

type :: flume_mhd_object
   !< FLUME MHD class definition.
   character(:), allocatable :: divergence_control        !< Divergence control: glm or none.
   logical                   :: has_glm=.false.           !< GLM cleaning active (host-side flag only).
   real(R8P)                 :: glm_ch=0._R8P             !< GLM cleaning speed (constant, uniform).
   real(R8P)                 :: glm_alpha=0._R8P          !< GLM damping parameter.
   real(R8P)                 :: glm_damping_length=0._R8P !< GLM damping length (min-cell: set by set_glm_damping).
   logical                   :: glm_damping_min_cell=.false. !< The damping length is the minimum cell spacing.
   real(R8P)                 :: glm_damping=0._R8P        !< Damping rate c_h^2 / c_p^2 (set by set_glm_damping).
   character(:), allocatable :: glm_ch_check              !< Action when the fastest wave outruns c_h.
   real(R8P)                 :: divb_tol=0._R8P           !< div(B) monitor tolerance (0: disabled).
   logical                   :: divb_error=.false.        !< Stop when the div(B) tolerance is exceeded.
   real(R8P)                 :: rho_floor=0._R8P          !< Density floor (0: disabled).
   real(R8P)                 :: p_floor=0._R8P            !< Pressure floor (0: disabled).
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize MHD configs.
      procedure, pass(self) :: load_from_file !< Load config from file.
      procedure, pass(self) :: set_glm_damping !< Set the damping length (min-cell) and the damping rate.
endtype flume_mhd_object

contains
   ! public methods
   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(flume_mhd_object), intent(in) :: self             !< MHD configs.
   character(len=:), allocatable       :: desc             !< Description.
   character(len=1), parameter         :: NL=new_line('a') !< New line character.

   desc =       mpih%myrankstr//'MHD main data'//NL
   desc = desc//mpih%myrankstr//'  divergence_control: '//self%divergence_control//NL
   if (self%has_glm) then
      desc = desc//mpih%myrankstr//'  glm_ch:             '//trim(str(self%glm_ch))//NL
      desc = desc//mpih%myrankstr//'  glm_alpha:          '//trim(str(self%glm_alpha))//NL
      if (self%glm_damping_min_cell) then
         desc = desc//mpih%myrankstr//'  glm_damping_length: '//GLM_DAMPING_LENGTH_MIN_CELL//NL
      else
         desc = desc//mpih%myrankstr//'  glm_damping_length: '//trim(str(self%glm_damping_length))//NL
      endif
      desc = desc//mpih%myrankstr//'  glm_ch_check:       '//self%glm_ch_check//NL
   endif
   desc = desc//mpih%myrankstr//'  divb_tol:           '//trim(str(self%divb_tol))//NL
   desc = desc//mpih%myrankstr//'  divb_error:         '//trim(merge('.true. ', '.false.', self%divb_error))//NL
   desc = desc//mpih%myrankstr//'  rho_floor:          '//trim(str(self%rho_floor))//NL
   desc = desc//mpih%myrankstr//'  p_floor:            '//trim(str(self%p_floor))
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize MHD configs.
   class(flume_mhd_object), intent(inout) :: self            !< MHD configs.
   type(file_ini),          intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   print '(A)', mpih%myrankstr//'flume_mhd_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', mpih%myrankstr//'flume_mhd_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters)
   !< Load config from file; every key in use is required and an unknown or out-of-range value is fatal.
   class(flume_mhd_object), intent(inout) :: self            !< MHD configs.
   type(file_ini),          intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   character(999)                         :: buff            !< Option value buffer.
   integer(I4P)                           :: error           !< Error status.

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='divergence_control', val=buff, error=error)
   if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(divergence_control)')
   self%divergence_control = trim(adjustl(strip_control(buff)))
   select case(self%divergence_control)
   case(DIVERGENCE_CONTROL_GLM)
      self%has_glm = .true.
   case(DIVERGENCE_CONTROL_NONE)
      self%has_glm = .false.
   case default
      call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(divergence_control) "'//self%divergence_control// &
                               '"; expected one of '//DIVERGENCE_CONTROL_GLM//', '//DIVERGENCE_CONTROL_NONE)
   endselect

   if (self%has_glm) then
      call load_real(key='glm_ch', val=self%glm_ch)
      if (.not.(self%glm_ch > 0._R8P)) call range_error(key='glm_ch', val=self%glm_ch, expected='> 0')
      call load_real(key='glm_alpha', val=self%glm_alpha)
      if (.not.(self%glm_alpha >= 0._R8P)) call range_error(key='glm_alpha', val=self%glm_alpha, expected='>= 0')
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='glm_damping_length', val=buff, error=error)
      if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(glm_damping_length)')
      buff = trim(adjustl(strip_control(buff)))
      self%glm_damping_min_cell = trim(buff) == GLM_DAMPING_LENGTH_MIN_CELL
      if (.not.self%glm_damping_min_cell) then
         read(buff, *, iostat=error) self%glm_damping_length
         if (error /= 0) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(glm_damping_length) = "'//trim(buff)// &
                                                  '"; expected '//GLM_DAMPING_LENGTH_MIN_CELL//' or a positive length')
         if (.not.(self%glm_damping_length > 0._R8P)) &
            call range_error(key='glm_damping_length', val=self%glm_damping_length, expected='> 0')
      endif
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='glm_ch_check', val=buff, error=error)
      if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(glm_ch_check)')
      self%glm_ch_check = trim(adjustl(strip_control(buff)))
      select case(self%glm_ch_check)
      case(GLM_CH_CHECK_WARNING, GLM_CH_CHECK_ERROR)
      case default
         call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(glm_ch_check) "'//self%glm_ch_check// &
                                  '"; expected one of '//GLM_CH_CHECK_WARNING//', '//GLM_CH_CHECK_ERROR)
      endselect
   else
      self%glm_ch_check = GLM_CH_CHECK_WARNING
   endif

   call load_real(key='divb_tol', val=self%divb_tol)
   if (.not.(self%divb_tol >= 0._R8P)) call range_error(key='divb_tol', val=self%divb_tol, expected='>= 0')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='divb_error', val=self%divb_error, error=error)
   if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(divb_error)')
   call load_real(key='rho_floor', val=self%rho_floor)
   if (.not.(self%rho_floor >= 0._R8P)) call range_error(key='rho_floor', val=self%rho_floor, expected='>= 0')
   call load_real(key='p_floor', val=self%p_floor)
   if (.not.(self%p_floor >= 0._R8P)) call range_error(key='p_floor', val=self%p_floor, expected='>= 0')
   contains
      subroutine load_real(key, val)
      !< Load one required real key.
      character(*), intent(in)  :: key !< Option name.
      real(R8P),    intent(out) :: val !< Option value.

      call file_parameters%get(section_name=INI_SECTION_NAME, option_name=key, val=val, error=error)
      if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].('//key//')')
      endsubroutine load_real

      subroutine range_error(key, val, expected)
      !< Stop on an out-of-range value.
      character(*), intent(in) :: key      !< Option name.
      real(R8P),    intent(in) :: val      !< Option value.
      character(*), intent(in) :: expected !< Accepted range.

      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].('//key//') = '//trim(str(val))//'; expected '//expected)
      endsubroutine range_error
   endsubroutine load_from_file

   subroutine set_glm_damping(self, min_cell)
   !< Set the damping length (with `min-cell`, the minimum cell spacing of the realm) and the damping rate
   !< `c_h^2 / c_p^2 = glm_alpha c_h / L` (issue #41, section 3.5); a no-op without GLM.
   class(flume_mhd_object), intent(inout) :: self     !< MHD configs.
   real(R8P),               intent(in)    :: min_cell !< Minimum cell spacing of the realm (active directions).

   if (.not.self%has_glm) return
   if (self%glm_damping_min_cell) self%glm_damping_length = min_cell
   self%glm_damping = self%glm_alpha * self%glm_ch / self%glm_damping_length
   print '(A)', mpih%myrankstr//'MHD GLM damping: length '//trim(str(self%glm_damping_length))//', rate c_h^2/c_p^2 '// &
                trim(str(self%glm_damping))
   endsubroutine set_glm_damping
endmodule adam_flume_mhd_object
