!< ADAM, the `[turbulence]` section: algebraic subgrid-scale model of a realm (issue #119, T1).
module adam_turbulence_object
!< ADAM, the `[turbulence]` section: algebraic subgrid-scale model of a realm (issue #119, T1).
!<
!< `sgs_model = none | vreman | wale | sigma` (default `none`; any other value is fatal), `sgs_constant` (default 0.07 Vreman,
!< 0.5 WALE, 1.35 Sigma; `> 0`), `prandtl_t` (default 0.9; `> 0`). The app owns the object (a realm component), the field and the
!< selection of the kernel of the model: `adam_turbulence_kernels` (CPU), `adam_fnl_turbulence_kernels` (FNL).

! ADAM modules
use :: adam_mpih_global, only : mpih
! third party modules
use :: finer,            only : file_ini
use :: penf,             only : I4P, R8P, str

implicit none
private
public :: SGS_NONE
public :: SGS_SIGMA_MODEL
public :: SGS_VREMAN_MODEL
public :: SGS_WALE_MODEL
public :: turbulence_object

integer(I4P), parameter :: SGS_NONE=0_I4P         !< No model.
integer(I4P), parameter :: SGS_VREMAN_MODEL=1_I4P !< Vreman (2004).
integer(I4P), parameter :: SGS_WALE_MODEL=2_I4P   !< WALE (Nicoud & Ducros 1999).
integer(I4P), parameter :: SGS_SIGMA_MODEL=3_I4P  !< Sigma (Nicoud et al. 2011).

character(len=10), parameter :: INI_SECTION_NAME='turbulence' !< Input section.

type :: turbulence_object
   !< Subgrid-scale model parameters.
   integer(I4P) :: model=SGS_NONE      !< Model id.
   real(R8P)    :: constant=0._R8P     !< Model constant.
   real(R8P)    :: prandtl_t=0.9_R8P   !< Turbulent Prandtl number.
   logical      :: is_active=.false.   !< A model is on.
   contains
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: load_from_file !< Load the section.
endtype turbulence_object

contains
   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(turbulence_object), intent(in) :: self !< Turbulence.
   character(len=:), allocatable        :: desc !< Description.

   desc = mpih%myrankstr//'  sgs model id:    '//trim(str(self%model))//', constant '//trim(str(self%constant))// &
          ', Pr_t '//trim(str(self%prandtl_t))
   endfunction description

   subroutine load_from_file(self, file_parameters)
   !< Load the section: without it, or with `sgs_model = none`, no model.
   class(turbulence_object), intent(inout) :: self            !< Turbulence.
   type(file_ini),           intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   character(99)                           :: model           !< Model name.
   integer(I4P)                            :: error           !< Error status.

   model = 'none'
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='sgs_model', val=model, error=error)
   if (error > 0) model = 'none'
   select case(trim(adjustl(model)))
   case('none')
      self%model = SGS_NONE ; self%constant = 0._R8P
   case('vreman')
      self%model = SGS_VREMAN_MODEL ; self%constant = 0.07_R8P
   case('wale')
      self%model = SGS_WALE_MODEL ; self%constant = 0.5_R8P
   case('sigma')
      self%model = SGS_SIGMA_MODEL ; self%constant = 1.35_R8P
   case default
      call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(sgs_model) "'//trim(adjustl(model))// &
                               '": none, vreman, wale or sigma')
   endselect
   self%is_active = self%model /= SGS_NONE
   if (.not.self%is_active) return
   call read_positive('sgs_constant', self%constant)
   self%prandtl_t = 0.9_R8P
   call read_positive('prandtl_t', self%prandtl_t)

   contains
      subroutine read_positive(key, val)
      !< Read an optional positive value: absent keeps the default, not positive is fatal.
      character(*), intent(in)    :: key   !< Option name.
      real(R8P),    intent(inout) :: val   !< Value.
      real(R8P)                   :: given !< Value read.
      integer(I4P)                :: err   !< Error status.

      given = val
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name=key, val=given, error=err)
      if (err > 0) return
      if (.not.(given > 0._R8P)) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].('//key//') must be positive')
      val = given
      endsubroutine read_positive
   endsubroutine load_from_file
endmodule adam_turbulence_object
