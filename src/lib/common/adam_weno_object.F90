!< ADAM, WENO configuration and the run-time dispatchers over SETH's WENO kernels.
module adam_weno_object
!< ADAM, WENO configuration and the run-time dispatchers over SETH's WENO kernels.
!<
!< [SETH](https://github.com/szaghi/SETH) is the only source of WENO, TENO and hybrid operators in ADAM (issue #139):
!< every coefficient table, smoothness indicator and nonlinear weight is SETH's, in kernels with compile-time
!< coefficients (`pure`, `acc routine seq`, `omp declare target`: nothing to mirror on a device). ADAM keeps the policy:
!< `weno_object` holds the `[weno]` configuration (the scheme, its half-width `S`, the weights `js` or `si` and their
!< descaler switch `sigma`, the regulariser `eps`), and the module holds three run-time dispatchers that select SETH's
!< kernel of the width `S` for the callers whose `S` is a run-time value per call (RIPPLE, the tutorial application).
!< FLUME and PRISM call SETH's kernel of each width directly, from kernels compiled with `S` a named constant.
!<
!< SETH's stencil convention: cells by offset from the face i+1/2 (cell 0 left of it), `vp(1-S:S)` the values
!< reconstructed with a left bias, `vm(1-S:S)` those reconstructed with a right bias, both unshifted (`vp(S)` and
!< `vm(1-S)` are not read). The upwind schemes `weno-u-1` to `weno-u-9` (S = 1 ... 5, order 2S-1) are accepted; the
!< centred names `weno-c-*` are refused until an application calls SETH's centred kernels (#139, D6).

! ADAM singleton objects
use :: adam_ini_library, only : ini_check, ini_stop
use :: adam_mpih_global, only : mpih
! third party modules
use :: finer
use :: penf
use :: seth, only : SETH_WENO_JS_EPS,                                                                                      &
                    seth_weno1_js, seth_weno3_js, seth_weno5_js, seth_weno7_js, seth_weno9_js,                             &
                    seth_weno1_js_wr, seth_weno3_js_wr, seth_weno5_js_wr, seth_weno7_js_wr, seth_weno9_js_wr,              &
                    seth_weno1_js_interp_wr, seth_weno3_js_interp_wr, seth_weno5_js_interp_wr, seth_weno7_js_interp_wr,    &
                    seth_weno9_js_interp_wr

implicit none
save
private
public :: weno_object
public :: weno_reconstruct
public :: weno_reconstruct_wr
public :: weno_reconstruct_interp_wr
public :: WENO_WEIGHTS_JS
public :: WENO_WEIGHTS_SI
public :: S_MAX
public :: S_MAX_M1
public :: WENO_U_1
public :: WENO_U_3
public :: WENO_U_5
public :: WENO_U_7
public :: WENO_U_9
public :: WENO_C_2
public :: WENO_C_4
public :: WENO_C_6
public :: WENO_C_8

integer(I4P), parameter :: S_MAX=5    !< Maximum number/dimensions of stencils.
integer(I4P), parameter :: S_MAX_M1=4 !< Maximum number/dimensions of stencils minus 1.

character(len=8), parameter :: WENO_U_1="weno-u-1" !< Parameter of WENO scheme, upwind   S=1, order 1.
character(len=8), parameter :: WENO_U_3="weno-u-3" !< Parameter of WENO scheme, upwind   S=2, order 3.
character(len=8), parameter :: WENO_U_5="weno-u-5" !< Parameter of WENO scheme, upwind   S=3, order 5.
character(len=8), parameter :: WENO_U_7="weno-u-7" !< Parameter of WENO scheme, upwind   S=4, order 7.
character(len=8), parameter :: WENO_U_9="weno-u-9" !< Parameter of WENO scheme, upwind   S=5, order 9.
character(len=8), parameter :: WENO_C_2="weno-c-2" !< Centred WENO name, S=1, order 2: refused (#139 D6).
character(len=8), parameter :: WENO_C_4="weno-c-4" !< Centred WENO name, S=2, order 4: refused (#139 D6).
character(len=8), parameter :: WENO_C_6="weno-c-6" !< Centred WENO name, S=3, order 6: refused (#139 D6).
character(len=8), parameter :: WENO_C_8="weno-c-8" !< Centred WENO name, S=4, order 8: refused (#139 D6).

character(len=2), parameter :: WENO_WEIGHTS_JS="js" !< Nonlinear weights of Jiang & Shu, absolute `eps` (default).
character(len=2), parameter :: WENO_WEIGHTS_SI="si" !< Scale-invariant weights, `eps` relative to the data (#49).

character(len=4), parameter :: INI_SECTION_NAME="weno" !< INI (config) file section name containing time configs.

type :: weno_object
   !< WENO configuration: the scheme of `[weno]`, the kernels are SETH's.
   character(:), allocatable :: scheme                    !< WENO scheme.
   integer(I4P)              :: S=0_I4P                   !< Stencils number/dimensions, 2S-1 order of accuracy.
   character(:), allocatable :: weights                   !< Nonlinear weights, WENO_WEIGHTS_JS or WENO_WEIGHTS_SI.
   real(R8P)                 :: sigma=0._R8P              !< Descaler switch of the weights: 0 Jiang-Shu, 1 scale-invariant.
   real(R8P)                 :: eps=SETH_WENO_JS_EPS      !< Regulariser of the smoothness indicators of SETH's kernels.
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize class.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype weno_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(weno_object), intent(in) :: self             !< WENO object.
   character(len=:), allocatable  :: desc             !< Description.
   character(len=1), parameter    :: NL=new_line('a') !< New line character.

   desc =       mpih%myrankstr//'WENO scheme main data'//NL
   desc = desc//mpih%myrankstr//'  scheme:                         '//trim(    self%scheme              )//NL
   desc = desc//mpih%myrankstr//'  S:                              '//trim(str(self%S                  ))//NL
   desc = desc//mpih%myrankstr//'  eps:                            '//trim(str(self%eps                ))//NL
   if (allocated(self%weights)) &
   desc = desc//mpih%myrankstr//'  weights:                        '//trim(    self%weights             )//NL
   desc = desc//mpih%myrankstr//'  kernels:                        SETH (seth_weno'//trim(str(2*self%S-1, .true.))//'_js*)'
   endfunction description

   subroutine initialize(self, file_parameters, scheme)
   !< Initialize class: the configuration of `[weno]` (or the scheme alone, the Jiang-Shu weights).
   class(weno_object), intent(inout)        :: self            !< WENO object.
   type(file_ini),     intent(in), optional :: file_parameters !< Simulation parameters ini file handler.
   character(*),       intent(in), optional :: scheme          !< WENO scheme.

   call mpih%print_message('weno_object%initialize start')
   if (present(file_parameters)) then
      call self%load_from_file(file_parameters=file_parameters)
   elseif (present(scheme)) then
      self%scheme = trim(adjustl(scheme))
      self%weights = WENO_WEIGHTS_JS
   else
      call mpih%error_stop(msg=': failed to initialize weno object, file parameters or scheme type must be passed')
   endif
   select case(self%scheme)
   case(WENO_U_1)
      self%S = 1
   case(WENO_U_3)
      self%S = 2
   case(WENO_U_5)
      self%S = 3
   case(WENO_U_7)
      self%S = 4
   case(WENO_U_9)
      self%S = 5
   case(WENO_C_2, WENO_C_4, WENO_C_6, WENO_C_8)
      ! refused until an application calls SETH's centred kernels seth_weno{4,6,8}_c (issue #139, D6): PRISM used to run
      ! the upwind primitive at that S under a centred name (F3), and no input names a centred scheme
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(scheme)="'//self%scheme//'": the centred WENO schemes are '// &
                               'refused (issue #139); accepted: '//WENO_U_1//', '//WENO_U_3//', '//WENO_U_5//', '//      &
                               WENO_U_7//', '//WENO_U_9)
   case default
      call mpih%error_stop(msg=': failed to initialize weno object, WENO scheme "'//self%scheme//'" unknown; accepted: '// &
                               WENO_U_1//', '//WENO_U_3//', '//WENO_U_5//', '//WENO_U_7//', '//WENO_U_9)
   endselect
   self%eps = SETH_WENO_JS_EPS
   select case(self%weights)
   case(WENO_WEIGHTS_JS)
      self%sigma = 0._R8P
   case(WENO_WEIGHTS_SI)
      self%sigma = 1._R8P
   case default
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(weights)="'//self%weights//'" unknown: accepted "'// &
                               WENO_WEIGHTS_JS//'" (Jiang-Shu, default) and "'//WENO_WEIGHTS_SI//'" (scale-invariant)')
   endselect
   print '(A)', self%description()
   call mpih%print_message('weno_object%initialize finish')
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters)
   !< Load config from file.
   class(weno_object), intent(inout)        :: self            !< WENO object.
   type(file_ini),     intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   character(99)                            :: buff_c          !< Character buffer.
   character(:), allocatable                :: sname           !< Section name.
   integer(I4P)                             :: error           !< Error status.

   sname = INI_SECTION_NAME
   call file_parameters%get(section_name=sname, option_name='scheme', val=buff_c, error=error)
   if (error>0) call ini_stop(msg='failed to load ['//INI_SECTION_NAME//'].(scheme)')
   self%scheme = trim(adjustl(buff_c))
   ! optional: the Jiang-Shu weights unless the scale-invariant ones are asked for
   buff_c = WENO_WEIGHTS_JS
   call file_parameters%get(section_name=sname, option_name='weights', val=buff_c, error=error)
   self%weights = trim(adjustl(buff_c))
   ! The order reduction (ROR) and the IB order reduction keys are no longer read (issue #80; their state left with
   ! ADAM's WENO kernels, #139 S6): an input that sets them reports them unread.
   call ini_check()
   endsubroutine load_from_file

   pure subroutine weno_reconstruct(S, vp, vm, eps, rmu, rp, rm)
   !< WENO-JS reconstruction of one face by SETH's kernel of width `S` (issue #139): the run-time dispatcher of the callers
   !< whose `S` is a run-time value per call (RIPPLE, the tutorial application); FLUME and PRISM call SETH's kernel of each
   !< width directly from a kernel compiled with `S` a named constant. SETH's stencil convention: cells by offset from the
   !< face (cell 0 left of it), `vp` the left-biased values, `vm` the right-biased ones, both unshifted.
   integer(I4P), intent(in)  :: S         !< Candidates of the scheme, WENO-(2S-1).
   real(R8P),    intent(in)  :: vp(1-S:S) !< Values reconstructed with a left bias.
   real(R8P),    intent(in)  :: vm(1-S:S) !< Values reconstructed with a right bias.
   real(R8P),    intent(in)  :: eps       !< Regulariser of the smoothness indicators.
   real(R8P),    intent(in)  :: rmu       !< Inverse descaler (1: Jiang-Shu weights).
   real(R8P),    intent(out) :: rp        !< Left state at the face.
   real(R8P),    intent(out) :: rm        !< Right state at the face.
   !$acc routine seq
   !$omp declare target

   select case(S)
   case(1) ; call seth_weno1_js(vp=vp, vm=vm, eps=eps, rmu=rmu, rp=rp, rm=rm)
   case(2) ; call seth_weno3_js(vp=vp, vm=vm, eps=eps, rmu=rmu, rp=rp, rm=rm)
   case(3) ; call seth_weno5_js(vp=vp, vm=vm, eps=eps, rmu=rmu, rp=rp, rm=rm)
   case(4) ; call seth_weno7_js(vp=vp, vm=vm, eps=eps, rmu=rmu, rp=rp, rm=rm)
   case(5) ; call seth_weno9_js(vp=vp, vm=vm, eps=eps, rmu=rmu, rp=rp, rm=rm)
   endselect
   endsubroutine weno_reconstruct

   pure subroutine weno_reconstruct_wr(S, vp, vm, eps, rmu, rp, rm, wrp, wrm)
   !< `weno_reconstruct` with the weights ratios `min_k w_k/d_k` of the two states.
   integer(I4P), intent(in)  :: S         !< Candidates of the scheme, WENO-(2S-1).
   real(R8P),    intent(in)  :: vp(1-S:S) !< Values reconstructed with a left bias.
   real(R8P),    intent(in)  :: vm(1-S:S) !< Values reconstructed with a right bias.
   real(R8P),    intent(in)  :: eps       !< Regulariser of the smoothness indicators.
   real(R8P),    intent(in)  :: rmu       !< Inverse descaler (1: Jiang-Shu weights).
   real(R8P),    intent(out) :: rp        !< Left state at the face.
   real(R8P),    intent(out) :: rm        !< Right state at the face.
   real(R8P),    intent(out) :: wrp       !< Weights ratio of the left state.
   real(R8P),    intent(out) :: wrm       !< Weights ratio of the right state.
   !$acc routine seq
   !$omp declare target

   select case(S)
   case(1) ; call seth_weno1_js_wr(vp=vp, vm=vm, eps=eps, rmu=rmu, rp=rp, rm=rm, wrp=wrp, wrm=wrm)
   case(2) ; call seth_weno3_js_wr(vp=vp, vm=vm, eps=eps, rmu=rmu, rp=rp, rm=rm, wrp=wrp, wrm=wrm)
   case(3) ; call seth_weno5_js_wr(vp=vp, vm=vm, eps=eps, rmu=rmu, rp=rp, rm=rm, wrp=wrp, wrm=wrm)
   case(4) ; call seth_weno7_js_wr(vp=vp, vm=vm, eps=eps, rmu=rmu, rp=rp, rm=rm, wrp=wrp, wrm=wrm)
   case(5) ; call seth_weno9_js_wr(vp=vp, vm=vm, eps=eps, rmu=rmu, rp=rp, rm=rm, wrp=wrp, wrm=wrm)
   endselect
   endsubroutine weno_reconstruct_wr

   pure subroutine weno_reconstruct_interp_wr(S, vp, vm, eps, rmu, rp, rm, wrp, wrm)
   !< `weno_reconstruct_wr` from point values (the interpolation coefficients of SETH's `_interp` kernels).
   integer(I4P), intent(in)  :: S         !< Candidates of the scheme, WENO-(2S-1).
   real(R8P),    intent(in)  :: vp(1-S:S) !< Values interpolated with a left bias.
   real(R8P),    intent(in)  :: vm(1-S:S) !< Values interpolated with a right bias.
   real(R8P),    intent(in)  :: eps       !< Regulariser of the smoothness indicators.
   real(R8P),    intent(in)  :: rmu       !< Inverse descaler (1: Jiang-Shu weights).
   real(R8P),    intent(out) :: rp        !< Left state at the face.
   real(R8P),    intent(out) :: rm        !< Right state at the face.
   real(R8P),    intent(out) :: wrp       !< Weights ratio of the left state.
   real(R8P),    intent(out) :: wrm       !< Weights ratio of the right state.
   !$acc routine seq
   !$omp declare target

   select case(S)
   case(1) ; call seth_weno1_js_interp_wr(vp=vp, vm=vm, eps=eps, rmu=rmu, rp=rp, rm=rm, wrp=wrp, wrm=wrm)
   case(2) ; call seth_weno3_js_interp_wr(vp=vp, vm=vm, eps=eps, rmu=rmu, rp=rp, rm=rm, wrp=wrp, wrm=wrm)
   case(3) ; call seth_weno5_js_interp_wr(vp=vp, vm=vm, eps=eps, rmu=rmu, rp=rp, rm=rm, wrp=wrp, wrm=wrm)
   case(4) ; call seth_weno7_js_interp_wr(vp=vp, vm=vm, eps=eps, rmu=rmu, rp=rp, rm=rm, wrp=wrp, wrm=wrm)
   case(5) ; call seth_weno9_js_interp_wr(vp=vp, vm=vm, eps=eps, rmu=rmu, rp=rp, rm=rm, wrp=wrp, wrm=wrm)
   endselect
   endsubroutine weno_reconstruct_interp_wr
endmodule adam_weno_object
