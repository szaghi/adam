!< ADAM, the time integrator adapter: ADAM's `[runge_kutta]` schemes on a SELAH stepper (SELAH M7, issue #132).

#include "fundal.H"

module adam_stepper_object
!< ADAM, the time integrator adapter: ADAM's `[runge_kutta]` schemes on a SELAH stepper (SELAH M7, issue #132).
!<
!< [SELAH](https://github.com/szaghi/SELAH) owns the time operators (the coefficients, the registers and the stage
!< kernels); this object is the thin ADAM side of it:
!<
!< - the input: `[runge_kutta] scheme`, read as `rk_object` reads it (no other key), and the mapping of ADAM's scheme
!<   names onto a SELAH scheme and storage form (`stepper_map`);
!< - the field: ADAM's `q(nv, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, nb)` (`LAYOUT_V_FIRST`, the CPU backends) or
!<   `q_gpu(nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, nv)` (`LAYOUT_V_LAST`, the FNL backends), the box of each call
!<   built from the active blocks (`field%blocks_number`, never stored: it changes at every regrid);
!< - the immersed boundary: `ib%phi` is SELAH's mask (a cell is active where the last component, the all-solids
!<   distance, is negative, as in `rk_object`);
!< - the registers: allocated here through FUNDAL (device memory in a device build) and passed to SELAH (`registers=`);
!<   on the host zeroed block by block the first time a block is active (issue #45: zeroing the whole fill-the-node
!<   pool at allocation would commit all of it), on the device zeroed at allocation (as the FNL stage buffers always
!<   were); reachable between steps as scratch (`scratch(n)`, the D_flux diagnostic and the equilibrium residual, #132
!<   rule 3).
!<
!< **The scheme mapping** (#132): `runge-kutta-1/2/3` are the two-register Shu-Osher forms of Euler, SSP-22 and SSP-33
!< (ADAM's low-storage coefficients are these, bitwise); `runge-kutta-ssp-11/22/33/54` are the Butcher forms of Euler,
!< SSP-22, SSP-33 and SSP-54 (bitwise ADAM's tableaux), every stage on its own register (`stages_on_q` off), as ADAM's
!< stages are; `runge-kutta-yoshida` is SELAH's Yoshida-4 splitting (the PRISM integrator of phase 4: it needs the
!< components of its two partitions).
!<
!< **The stage window.** Between `begin_stage(k)` and `end_stage(k)` the state of stage k is `q` itself when
!< `state_is_q(k)` (every stage of the two-register forms), else `stage_state()`; the caller fills its ghost cells and
!< writes the residual into `rhs_register()`. The masked (solid) cells of every stage state hold `q`, as every stage
!< buffer of `rk_object` did: SELAH's mask policy (`masked_q=.true.`, szaghi/SELAH#5), written by its stage kernel in
!< the same pass; the immersed-boundary extrapolation of the residual starts from them.
!<
!< The same source serves the CPU backends (host build) and the FNL backends (OpenACC build, `-DDEV_OAC`): the arrays
!< passed (`q`, `phi`, `dq`) are host arrays in the first, device arrays (FUNDAL pointers) in the second.
use :: adam_field_object, only : field_object
use :: adam_grid_object,  only : grid_object
use :: adam_ini_library,  only : ini_check, ini_get
use :: adam_mpih_global,  only : mpih
use :: fundal,            only : dev_alloc, dev_free
use :: selah,             only : stepper, box, box_interior, registers_needed, R_P, LAYOUT_V_FIRST, LAYOUT_V_LAST,  &
                                 selah_copy, butcher, butcher_tableau,                                                &
                                 FORM_BUTCHER, FORM_SHU_OSHER_2R, FORM_SPLITTING,                                     &
                                 SCHEME_EULER, SCHEME_SSP_22, SCHEME_SSP_33, SCHEME_SSP_54, SCHEME_YOSHIDA_4
! third party modules
use :: finer
use :: penf

implicit none
private
public :: stepper_object
public :: stepper_map
public :: stepper_stored_fields
public :: stepper_butcher_stages
public :: LAYOUT_V_FIRST
public :: LAYOUT_V_LAST
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

type :: stepper_object
   !< ADAM's `[runge_kutta]` scheme on a SELAH stepper.
   character(:), allocatable          :: scheme                       !< ADAM scheme name (`[runge_kutta] scheme`).
   character(:), allocatable          :: selah_scheme                 !< SELAH scheme.
   character(:), allocatable          :: form                         !< SELAH storage form.
   type(stepper)                      :: st                           !< SELAH stepper.
   real(R8P), pointer, contiguous     :: regs(:,:,:,:,:,:) => null()  !< Registers (v, i, j, k, b, slot).
   integer(I4P)                       :: nregs=0_I4P                  !< Registers number.
   integer(I4P)                       :: blocks_zeroed=0_I4P          !< Block slots of `regs` zeroed so far (issue #45).
   integer(I4P)                       :: nv=0_I4P                     !< Number of variables.
   integer(I4P)                       :: ni=0_I4P, nj=0_I4P, nk=0_I4P !< Cells per block, per direction.
   integer(I4P)                       :: ngc=0_I4P                    !< Ghost cells.
   integer(I4P)                       :: nb=0_I4P                     !< Block slots of the pool.
   integer(I4P)                       :: layout=LAYOUT_V_FIRST        !< Field layout: (v, i, j, k, b) or (b, i, j, k, v).
   contains
      ! public methods
      procedure, pass(self) :: begin_step      !< Open a step.
      procedure, pass(self) :: begin_stage     !< Open stage k: build its state.
      procedure, pass(self) :: description     !< Return pretty-printed object description.
      procedure, pass(self) :: destroy         !< Free the registers.
      procedure, pass(self) :: end_stage       !< Close stage k.
      procedure, pass(self) :: finish_step     !< Close the step.
      procedure, pass(self) :: initialize      !< Initialize the stepper and its registers.
      procedure, pass(self) :: imag_axis_bound !< Stability bound on the imaginary axis (SELAH).
      procedure, pass(self) :: is_butcher      !< The form is Butcher's (the stages on their own registers).
      procedure, pass(self) :: is_splitting    !< The form is a splitting of two partitions.
      procedure, pass(self) :: is_two_register !< The form is the two-register (low-storage) Shu-Osher one.
      procedure, pass(self) :: load_from_file  !< Load config from file.
      procedure, pass(self) :: rhs_register    !< Residual register of the active stage.
      procedure, pass(self) :: scratch         !< A register as scratch, between steps.
      procedure, pass(self) :: set_rhs         !< Copy a residual into the residual register of the active stage.
      procedure, pass(self) :: ssp_coefficient !< SSP coefficient of the scheme (SELAH; 0 if not SSP).
      procedure, pass(self) :: stage_euler_fraction !< Forward-Euler step of stage k, a fraction of the step (SELAH).
      procedure, pass(self) :: stage_state     !< State register of the active stage (null when it is q).
      procedure, pass(self) :: stage_time      !< Fraction of the step at which stage k is evaluated.
      procedure, pass(self) :: stage_weight    !< Effective Butcher weight of stage k.
      procedure, pass(self) :: stages          !< Stages per step.
      procedure, pass(self) :: state_is_q      !< The state of stage k is q itself.
      ! private methods
      procedure, pass(self), private :: interior_box !< Box of the interior cells of the active blocks.
endtype stepper_object

contains
   ! public procedures
   pure subroutine stepper_map(scheme, selah_scheme, form)
   !< Map an ADAM `[runge_kutta]` scheme name onto a SELAH scheme and storage form; an unknown name gives empty strings.
   character(*),              intent(in)  :: scheme       !< ADAM scheme name.
   character(:), allocatable, intent(out) :: selah_scheme !< SELAH scheme.
   character(:), allocatable, intent(out) :: form         !< SELAH storage form.

   selah_scheme = '' ; form = ''
   select case(trim(adjustl(scheme)))
   case(RK_1)       ; selah_scheme = SCHEME_EULER     ; form = FORM_SHU_OSHER_2R
   case(RK_2)       ; selah_scheme = SCHEME_SSP_22    ; form = FORM_SHU_OSHER_2R
   case(RK_3)       ; selah_scheme = SCHEME_SSP_33    ; form = FORM_SHU_OSHER_2R
   case(RK_SSP_11)  ; selah_scheme = SCHEME_EULER     ; form = FORM_BUTCHER
   case(RK_SSP_22)  ; selah_scheme = SCHEME_SSP_22    ; form = FORM_BUTCHER
   case(RK_SSP_33)  ; selah_scheme = SCHEME_SSP_33    ; form = FORM_BUTCHER
   case(RK_SSP_54)  ; selah_scheme = SCHEME_SSP_54    ; form = FORM_BUTCHER
   case(RK_YOSHIDA) ; selah_scheme = SCHEME_YOSHIDA_4 ; form = FORM_SPLITTING
   endselect
   endsubroutine stepper_map

   pure function stepper_stored_fields(scheme) result(fields)
   !< Return the block-sized fields (of nv variables) the stepper of a scheme allocates: SELAH's `registers_needed`, the
   !< residual register included. It lets an app size its memory budget before the stepper exists; an unknown scheme
   !< returns -1.
   character(*), intent(in)  :: scheme       !< ADAM scheme name.
   integer(I4P)              :: fields       !< Block-sized fields.
   character(:), allocatable :: selah_scheme !< SELAH scheme.
   character(:), allocatable :: form         !< SELAH storage form.

   call stepper_map(scheme=scheme, selah_scheme=selah_scheme, form=form)
   fields = -1_I4P
   if (len(form) > 0) fields = registers_needed(selah_scheme, form)
   endfunction stepper_stored_fields

   pure function stepper_butcher_stages(scheme) result(stages)
   !< Return the stages of a scheme of the Butcher form (the SSP schemes), from SELAH's tableau; -1 for any other scheme.
   !< It lets an app size what follows the field's stages (PRISM's particle stages) before the stepper exists.
   character(*), intent(in)  :: scheme       !< ADAM scheme name.
   integer(I4P)              :: stages       !< Stages per step.
   character(:), allocatable :: selah_scheme !< SELAH scheme.
   character(:), allocatable :: form         !< SELAH storage form.
   type(butcher_tableau)     :: bt           !< SELAH tableau.

   call stepper_map(scheme=scheme, selah_scheme=selah_scheme, form=form)
   stages = -1_I4P
   if (form /= FORM_BUTCHER) return
   bt = butcher(selah_scheme)
   stages = bt%stages
   endfunction stepper_butcher_stages

   ! public methods
   subroutine begin_step(self, field, q, phi)
   !< Open a step: zero the block slots active for the first time, then let SELAH keep what its form needs of `q`.
   class(stepper_object), intent(inout)        :: self    !< Stepper.
   type(field_object),    intent(in)           :: field   !< Field (sibling realm component, threaded in).
   real(R8P),             intent(in)           :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)   !< State.
   real(R8P),             intent(in), optional :: phi(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< IB distance.

   ! the host registers are zeroed block by block (issue #45); the device ones at allocation (blocks_zeroed = nb)
   if (field%blocks_number > self%blocks_zeroed) then
      if (self%layout == LAYOUT_V_FIRST) then
         self%regs(:,:,:,:,self%blocks_zeroed+1:field%blocks_number,:) = 0._R8P
      else
         self%regs(self%blocks_zeroed+1:field%blocks_number,:,:,:,:,:) = 0._R8P
      endif
      self%blocks_zeroed = field%blocks_number
   endif
   call self%st%begin_step(self%interior_box(field), q, mask=phi)
   endsubroutine begin_step

   subroutine begin_stage(self, k, field, dt, q, phi)
   !< Open stage `k`: SELAH builds its state (its masked cells `q`, the mask policy).
   class(stepper_object), intent(inout)        :: self    !< Stepper.
   integer(I4P),          intent(in)           :: k       !< Stage.
   type(field_object),    intent(in)           :: field   !< Field (sibling realm component, threaded in).
   real(R8P),             intent(in)           :: dt      !< Time step.
   real(R8P),             intent(in)           :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)   !< State.
   real(R8P),             intent(in), optional :: phi(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< IB distance.

   call self%st%begin_stage(k, self%interior_box(field), dt, q, mask=phi)
   endsubroutine begin_stage

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(stepper_object), intent(in) :: self             !< Stepper.
   character(len=:), allocatable     :: desc             !< Description.
   character(len=1), parameter       :: NL=new_line('a') !< New line character.
   integer(I4P)                      :: k                !< Counter.

   desc = mpih%myrankstr//'Time integrator (SELAH stepper)'//NL
   desc = desc//mpih%myrankstr//'  scheme:    '//self%scheme//' -> SELAH '//self%selah_scheme//', form '//self%form//NL
   desc = desc//mpih%myrankstr//'  stages:    '//trim(str(self%st%stages(), .true.))//NL
   desc = desc//mpih%myrankstr//'  registers: '//trim(str(self%nregs, .true.))
   if (self%form /= FORM_SPLITTING) then
      desc = desc//NL//mpih%myrankstr//'  b:        '
      do k=1, self%st%stages()
         desc = desc//' '//trim(str(self%st%stage_weight(k)))
      enddo
   endif
   endfunction description

   subroutine destroy(self)
   !< Free the stepper and its registers.
   class(stepper_object), intent(inout) :: self !< Stepper.

   call self%st%destroy
   if (associated(self%regs)) call dev_free(self%regs)
   nullify(self%regs)
   self%nregs = 0_I4P
   self%blocks_zeroed = 0_I4P
   endsubroutine destroy

   subroutine end_stage(self, k, field, dt, q, phi)
   !< Close stage `k`, its residual written into `rhs_register()`: the two-register forms advance `q` here.
   class(stepper_object), intent(inout)        :: self    !< Stepper.
   integer(I4P),          intent(in)           :: k       !< Stage.
   type(field_object),    intent(in)           :: field   !< Field (sibling realm component, threaded in).
   real(R8P),             intent(in)           :: dt      !< Time step.
   real(R8P),             intent(inout)        :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)   !< State.
   real(R8P),             intent(in), optional :: phi(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< IB distance.

   call self%st%end_stage(k, self%interior_box(field), dt, q, mask=phi)
   endsubroutine end_stage

   subroutine finish_step(self, field, dt, q, dq, phi)
   !< Close the step: the Butcher form advances `q` now and gives the step residual `dq = sum_k b_k R_k` (0 in the masked
   !< cells, #128 D4); the two-register forms already advanced `q` and leave `dq` as the caller wrote it.
   class(stepper_object), intent(inout)           :: self    !< Stepper.
   type(field_object),    intent(in)              :: field   !< Field (sibling realm component, threaded in).
   real(R8P),             intent(in)              :: dt      !< Time step.
   real(R8P),             intent(inout)           :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)   !< State.
   real(R8P),             intent(inout), optional :: dq(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)  !< Step residual.
   real(R8P),             intent(in),    optional :: phi(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< IB distance.

   if (self%is_butcher()) then
      call self%st%finish_step(self%interior_box(field), dt, q, dq=dq, mask=phi)
   else
      call self%st%finish_step(self%interior_box(field), dt, q, mask=phi)
   endif
   endsubroutine finish_step

   subroutine initialize(self, field, grid, file_parameters, scheme, partitions, layout)
   !< Initialize the stepper of `[runge_kutta] scheme` (or of `scheme`) and its registers, the field's pool shape.
   class(stepper_object), intent(inout)        :: self            !< Stepper.
   type(field_object),    intent(in)           :: field           !< Field (sibling realm component, threaded in).
   type(grid_object),     intent(in)           :: grid            !< Grid (sibling realm component, threaded in).
   type(file_ini),        intent(in), optional :: file_parameters !< Simulation parameters ini file handler.
   character(*),          intent(in), optional :: scheme          !< Runge-Kutta scheme.
   integer(I4P),          intent(in), optional :: partitions(2,2) !< Splitting: the variables of its two partitions.
   integer(I4P),          intent(in), optional :: layout          !< Field layout (default LAYOUT_V_FIRST).
   integer(I4P)                                :: ierr            !< Error status.
   integer(I4P)                                :: lb(5), ub(5)    !< Bounds of the field.

   call mpih%print_message('stepper_object%initialize start')
   call self%destroy
   if (R_P /= R8P) call mpih%error_stop(msg=': SELAH must be built with R8P reals in ADAM (no -D_R_P_IS_R4P)')
   if (present(file_parameters)) then
      call self%load_from_file(file_parameters=file_parameters)
   elseif (present(scheme)) then
      self%scheme = trim(adjustl(scheme))
   else
      call mpih%error_stop(msg=': failed to initialize the stepper, one between file parameters and scheme must be passed')
   endif
   call stepper_map(scheme=self%scheme, selah_scheme=self%selah_scheme, form=self%form)
   if (len(self%form) == 0) call mpih%error_stop(msg=': unknown Runge-Kutta scheme "'//self%scheme//'"')
   if (self%form == FORM_SPLITTING .and. .not.present(partitions)) &
      call mpih%error_stop(msg=': the Runge-Kutta scheme "'//self%scheme//'" is a splitting: the stepper needs the '// &
                               'variables of its two partitions')
   self%nv = field%nv ; self%nb = field%nb
   self%ni = grid%ni ; self%nj = grid%nj ; self%nk = grid%nk ; self%ngc = grid%ngc
   self%layout = LAYOUT_V_FIRST ; if (present(layout)) self%layout = layout
   self%nregs = registers_needed(self%selah_scheme, self%form)
   lb = [1_I4P, 1_I4P-self%ngc, 1_I4P-self%ngc, 1_I4P-self%ngc, 1_I4P]
   if (self%layout == LAYOUT_V_FIRST) then
      ub = [self%nv, self%ni+self%ngc, self%nj+self%ngc, self%nk+self%ngc, self%nb]
   else
      ub = [self%nb, self%ni+self%ngc, self%nj+self%ngc, self%nk+self%ngc, self%nv]
   endif
#if defined DEV_OAC || defined DEV_OMP
   ! device memory: zeroed once at allocation, as the FNL stage buffers always were
   call dev_alloc(fptr_dev=self%regs, lbounds=[lb, 1_I4P], ubounds=[ub, self%nregs], init_value=0._R8P, ierr=ierr)
   self%blocks_zeroed = self%nb
#else
   ! host memory, not zeroed here: begin_step zeroes each block slot the first time it is active (issue #45)
   call dev_alloc(fptr_dev=self%regs, lbounds=[lb, 1_I4P], ubounds=[ub, self%nregs], ierr=ierr)
   self%blocks_zeroed = 0_I4P
#endif
   if (ierr /= 0) call mpih%error_stop(msg=': stepper_object%initialize: the '//trim(str(self%nregs, .true.))// &
                                           ' registers could not be allocated')
   call self%st%initialize(scheme=self%selah_scheme, form=self%form, lb=lb, ub=ub, partitions=partitions, &
                           registers=self%regs, stages_on_q=.false., masked_q=.true.)
   print '(A)', self%description()
   call mpih%print_message('stepper_object%initialize finish')
   endsubroutine initialize

   pure function imag_axis_bound(self) result(beta)
   !< Return the end of the scheme's stability region on the imaginary axis, |R(iy)| <= 1 for |y| <= beta (SELAH, from
   !< its stability polynomial; 0 when R(iy) leaves the unit disc at once, as for Euler and SSP-22).
   class(stepper_object), intent(in) :: self !< Stepper.
   real(R8P)                         :: beta !< Imaginary-axis bound.

   beta = self%st%imag_axis_bound()
   endfunction imag_axis_bound

   elemental function is_butcher(self) result(yes)
   !< Return true if the form is Butcher's: every stage on its own state register, the step closed by `finish_step`.
   class(stepper_object), intent(in) :: self !< Stepper.
   logical                           :: yes  !< True for the Butcher form.

   yes = .false.
   if (allocated(self%form)) yes = self%form == FORM_BUTCHER
   endfunction is_butcher

   elemental function is_splitting(self) result(yes)
   !< Return true if the form is a splitting of two partitions (each stage a flow of one partition).
   class(stepper_object), intent(in) :: self !< Stepper.
   logical                           :: yes  !< True for a splitting.

   yes = .false.
   if (allocated(self%form)) yes = self%form == FORM_SPLITTING
   endfunction is_splitting

   elemental function is_two_register(self) result(yes)
   !< Return true if the form is the two-register Shu-Osher one (every stage on q, low storage).
   class(stepper_object), intent(in) :: self !< Stepper.
   logical                           :: yes  !< True for the two-register form.

   yes = .false.
   if (allocated(self%form)) yes = self%form == FORM_SHU_OSHER_2R
   endfunction is_two_register

   subroutine load_from_file(self, file_parameters)
   !< Load config from file: `[runge_kutta] scheme`, the only key (as `rk_object` reads it).
   class(stepper_object), intent(inout) :: self            !< Stepper.
   type(file_ini),        intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   character(99)                        :: buff_c          !< Character buffer.

   call ini_get(file_parameters, INI_SECTION_NAME, 'scheme', buff_c)
   self%scheme = trim(adjustl(buff_c))
   call ini_check()
   endsubroutine load_from_file

   function rhs_register(self) result(r)
   !< Return the residual register of the active stage, with the field's bounds (null outside a stage window).
   class(stepper_object), intent(in) :: self         !< Stepper.
   real(R8P), pointer                :: r(:,:,:,:,:) !< Residual register.

   r => self%st%rhs_register()
   if (.not.associated(r)) call mpih%error_stop(msg=': stepper_object%rhs_register: no active stage')
   endfunction rhs_register

   function scratch(self, n) result(r)
   !< Return register `n` as scratch, between steps (SELAH's registers are scratch after `finish_step`): every form has
   !< at least two registers (#132 rule 3, #128 D1). Writing them leaves the next step unchanged: each step rebuilds
   !< what it reads of them.
   class(stepper_object), intent(in) :: self         !< Stepper.
   integer(I4P),          intent(in) :: n            !< Register.
   real(R8P), pointer                :: r(:,:,:,:,:) !< The register, with the field's bounds.

   if (self%st%k_active /= 0_I4P) &
      call mpih%error_stop(msg=': stepper_object%scratch: register '//trim(str(n, .true.))//' asked inside a stage')
   if (n < 1_I4P .or. n > self%nregs) &
      call mpih%error_stop(msg=': stepper_object%scratch: register '//trim(str(n, .true.))//' of '// &
                               trim(str(self%nregs, .true.)))
   r(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) => self%regs(:,:,:,:,:,n)
   endfunction scratch

   subroutine set_rhs(self, field, dq)
   !< Copy the residual `dq` into the residual register of the active stage, the interior cells of the active blocks
   !< (SELAH's copy kernel: on the device in a device build). For a caller that evaluates the residual into its own
   !< array (the residual log reads it, a correction acts on it) before the stage update reads the register.
   class(stepper_object), intent(in) :: self    !< Stepper.
   type(field_object),    intent(in) :: field   !< Field (sibling realm component, threaded in).
   real(R8P),             intent(in) :: dq(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Residual.
   real(R8P), pointer                :: r(:,:,:,:,:) !< Residual register.

   r => self%rhs_register()
   call selah_copy(self%interior_box(field), dq, r)
   endsubroutine set_rhs

   pure function ssp_coefficient(self) result(c)
   !< Return the SSP coefficient C of the scheme (SELAH, from its Shu-Osher form): every stage a convex combination of
   !< forward-Euler steps when C > 0; 0 for a scheme that is not SSP and for a splitting.
   class(stepper_object), intent(in) :: self !< Stepper.
   real(R8P)                         :: c    !< SSP coefficient.

   c = self%st%ssp_coefficient()
   endfunction ssp_coefficient

   pure function stage_euler_fraction(self, k) result(f)
   !< Return the forward-Euler step the residual of stage k enters, as a fraction of the step (SELAH, f_k = max_i
   !< beta_ik / alpha_ik of the Shu-Osher form): a limiter that makes the state of stage k plus f_k dt times its residual
   !< admissible keeps every later stage and the step admissible; 0 when the scheme has no such decomposition.
   class(stepper_object), intent(in) :: self !< Stepper.
   integer(I4P),          intent(in) :: k    !< Stage.
   real(R8P)                         :: f    !< Fraction of the step.

   f = self%st%stage_euler_fraction(k)
   endfunction stage_euler_fraction

   function stage_state(self) result(u)
   !< Return the state register of the active stage, with the field's bounds; null when the stage's state is q itself.
   class(stepper_object), intent(in) :: self         !< Stepper.
   real(R8P), pointer                :: u(:,:,:,:,:) !< Stage state.

   u => self%st%stage_state()
   endfunction stage_state

   pure function stage_time(self, k) result(c)
   !< Return the fraction of the step at which the residual of stage `k` is evaluated (the abscissa c_k).
   class(stepper_object), intent(in) :: self !< Stepper.
   integer(I4P),          intent(in) :: k    !< Stage.
   real(R8P)                         :: c    !< Stage time, fraction of the step.

   c = self%st%stage_time(k)
   endfunction stage_time

   pure function stage_weight(self, k) result(w)
   !< Return the effective Butcher weight b_k of stage `k`: the weight of its fluxes in the step (the flux register's).
   class(stepper_object), intent(in) :: self !< Stepper.
   integer(I4P),          intent(in) :: k    !< Stage.
   real(R8P)                         :: w    !< Stage weight.

   w = self%st%stage_weight(k)
   endfunction stage_weight

   pure function stages(self) result(n)
   !< Return the stages per step.
   class(stepper_object), intent(in) :: self !< Stepper.
   integer(I4P)                      :: n    !< Stages.

   n = self%st%stages()
   endfunction stages

   pure function state_is_q(self, k) result(yes)
   !< Return true if the state of stage `k` is `q` itself (every stage of the two-register forms).
   class(stepper_object), intent(in) :: self !< Stepper.
   integer(I4P),          intent(in) :: k    !< Stage.
   logical                           :: yes  !< True if the stage state is q.

   yes = self%st%stage_state_is_q(k)
   endfunction state_is_q

   ! private methods
   pure function interior_box(self, field) result(bx)
   !< Return the box of the interior cells of the active blocks (`field%blocks_number`, at the time of the call).
   class(stepper_object), intent(in) :: self  !< Stepper.
   type(field_object),    intent(in) :: field !< Field.
   type(box)                         :: bx    !< Box.

   bx = box_interior(self%layout, self%nv, self%ni, self%nj, self%nk, self%ngc, field%blocks_number)
   endfunction interior_box

endmodule adam_stepper_object
