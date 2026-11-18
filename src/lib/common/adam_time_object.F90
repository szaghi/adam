!< ADAM, time handler class definition: the step counter, the time, the time step and the `[time]` limits.
module adam_time_object
!< ADAM, time handler class definition: the step counter, the time, the time step and the `[time]` limits.
!<
!< The state every application integrates in time and the keys of the `[time]` section, `it_max`, `time_max` and `CFL`,
!< all required (issue #80 P3: FLUME, PRISM and RIPPLE each carried a copy). An application extends it with its own
!< state and overrides the predicates whose semantics it needs (PRISM's `is_to_save` also saves at `it == it_max`).
!< `time_unit` divides `time_max` for a model that integrates in a scaled time (PRISM's adimensional models).

use :: adam_ini_library, only : ini_check, ini_get
use :: adam_mpih_global, only : mpih
use :: finer,            only : file_ini
use :: penf,             only : I4P, R8P, str

implicit none
private
public :: time_object

character(len=4), parameter :: INI_SECTION_NAME="time" !< INI (config) file section name containing time configs.

type :: time_object
   !< Time handler class definition.
   integer(I4P) :: it_max=-1_I4P   !< Maximum number of integration time steps (<= 0: time-driven run).
   real(R8P)    :: time_max=1._R8P !< Maximum integration time.
   real(R8P)    :: CFL=0.3_R8P     !< CFL number.
   integer(I4P) :: it=0_I4P        !< Time steps counter.
   real(R8P)    :: time=0._R8P     !< Time.
   real(R8P)    :: dt=0._R8P       !< Time step.
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize time handler.
      procedure, pass(self) :: is_done        !< Return true if the run has reached its end.
      procedure, pass(self) :: is_to_save     !< Return true if the current step is a save step.
      procedure, pass(self) :: load_from_file !< Load config from file.
      procedure, pass(self) :: print_progress !< Print simulation progress.
endtype time_object

contains
   ! public methods
   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(time_object), intent(in) :: self             !< Time handler.
   character(len=:), allocatable  :: desc             !< Description.
   character(len=1), parameter    :: NL=new_line('a') !< New line character.

   desc =       mpih%myrankstr//'Time main data'//NL
   desc = desc//mpih%myrankstr//'  it_max:   '//trim(str(self%it_max  ))//NL
   desc = desc//mpih%myrankstr//'  time_max: '//trim(str(self%time_max))//NL
   desc = desc//mpih%myrankstr//'  CFL:      '//trim(str(self%CFL     ))
   endfunction description

   subroutine initialize(self, file_parameters, time_unit, verbose)
   !< Initialize time handler: load the `[time]` keys; `time_unit` (if present) divides `time_max`.
   class(time_object), intent(inout)        :: self            !< Time handler.
   type(file_ini),     intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   real(R8P),          intent(in), optional :: time_unit       !< Unit of a scaled time: time_max/time_unit is used.
   logical,            intent(in), optional :: verbose         !< Trigger verbose output (default .true.).
   logical                                  :: verbose_        !< Trigger verbose output, local variable.

   verbose_ = .true. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) print '(A)', mpih%myrankstr//'time_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   if (present(time_unit)) self%time_max = self%time_max / time_unit
   if (verbose_) then
      print '(A)', self%description()
      print '(A)', mpih%myrankstr//'time_object%initialize finish'
   endif
   endsubroutine initialize

   function is_done(self) result(done)
   !< Return true if the run has reached its end: `it_max` steps (iteration-driven) or `time_max` (time-driven).
   class(time_object), intent(in) :: self !< Time handler.
   logical                        :: done !< Check result.

   done = ((self%it_max <= 0_I4P) .and. (self%time >= self%time_max)) .or. &
          ((self%it_max >  0_I4P) .and. (self%it >= self%it_max))
   endfunction is_done

   function is_to_save(self, cadence) result(is_to)
   !< Return true if the current step is a save step: every `cadence` steps and at the end of the run.
   !<
   !< A non-positive cadence disables saving (it never reaches `mod(it, 0)`). Every save site uses this predicate
   !< only, so the gate and the write cannot disagree.
   class(time_object), intent(in) :: self    !< Time handler.
   integer(I4P),       intent(in) :: cadence !< Save cadence in time steps.
   logical                        :: is_to   !< Check result.

   is_to = .false.
   if (cadence > 0_I4P) is_to = (mod(self%it, cadence) == 0_I4P) .or. self%is_done()
   endfunction is_to_save

   subroutine load_from_file(self, file_parameters)
   !< Load config from file; every key is required.
   class(time_object), intent(inout) :: self            !< Time handler.
   type(file_ini),     intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call ini_get(file_parameters, INI_SECTION_NAME, 'it_max', self%it_max)
   call ini_get(file_parameters, INI_SECTION_NAME, 'time_max', self%time_max)
   call ini_get(file_parameters, INI_SECTION_NAME, 'CFL', self%CFL)
   call ini_check()
   endsubroutine load_from_file

   subroutine print_progress(self, nodes_number)
   !< Print simulation progress.
   class(time_object), intent(in) :: self         !< Time handler.
   integer(I4P),       intent(in) :: nodes_number !< Nodes number, global blocks number.

   associate(r=>mpih%myrankstr, it=>self%it, time=>self%time, dt=>self%dt, it_max=>self%it_max, time_max=>self%time_max)
   print '(A)', r//''
   print '(A)', r//'t:            '//trim(str(it, .true.))
   print '(A)', r//'nodes number: '//trim(str(nodes_number, .true.))
   print '(A)', r//'time step:    '//trim(str(dt, .true.))
   print '(A)', r//'time:         '//trim(str(time, .true.))
   if (it_max <= 0_I4P) then
      print '(A)', r//'progress:     '//trim(str(int(time/time_max * 100), .true.))//'%'
   else
      print '(A)', r//'progress:     '//trim(str(int((it*1._R8P)/it_max * 100), .true.))//'%'
   endif
   print '(A)', r//''
   endassociate
   endsubroutine print_progress
endmodule adam_time_object
