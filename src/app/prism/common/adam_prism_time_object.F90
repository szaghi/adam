!< ADAM, PRISM time handler class definition.
module adam_prism_time_object
!< ADAM, PRISM time handler class definition.
!<
!< The library time handler (`adam_time_object`, issue #80 P3) plus what PRISM needs: the origin of the time step and
!< the plasma-frequency bound of a PIC run (printed with the progress), the initial `dt` of 1e-4, the cadence predicate
!< and a save predicate that also saves at `it == it_max`. The adimensional models pass `time_unit = T0`.

use :: adam_mpih_global, only : mpih
use :: adam_time_object, only : time_object
use :: finer,            only : file_ini
use :: penf,             only : I4P, R8P, str

implicit none
private
public :: prism_time_object

type, extends(time_object) :: prism_time_object
   !< PRISM time handler class definition.
   character(len=18) :: dt_source='CFL' !< Origin of dt: CFL (fields) or plasma_frequencies (PIC).
   real(R8P)         :: omega_dt=0._R8P  !< PIC: highest plasma frequency times dt (0: not a PIC run).
   contains
      ! public methods
      procedure, pass(self) :: initialize     !< Initialize time handler (dt = 1e-4 until the first compute_dt).
      procedure, pass(self) :: is_on_cadence  !< Return true if the current step is a multiple of a cadence.
      procedure, pass(self) :: is_to_save     !< Return true if data are to save.
      procedure, pass(self) :: print_progress !< Print simulation progress.
endtype prism_time_object

contains
   ! public methods
   subroutine initialize(self, file_parameters, time_unit, verbose)
   !< Initialize time handler: the library `[time]` keys, and `dt = 1e-4` until the first `compute_dt`.
   class(prism_time_object), intent(inout)        :: self            !< Time handler.
   type(file_ini),           intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   real(R8P),                intent(in), optional :: time_unit       !< Unit of a scaled time (T0 of the adim models).
   logical,                  intent(in), optional :: verbose         !< Trigger verbose output (default .true.).

   self%dt = 0.0001_R8P
   call self%time_object%initialize(file_parameters=file_parameters, time_unit=time_unit, verbose=verbose)
   endsubroutine initialize

   pure function is_on_cadence(self, cadence)
   !< Return true if the current step is a multiple of `cadence`.
   !<
   !< A non-positive cadence is never on cadence: it disables the output it drives and never reaches `mod(it, 0)`
   !< (issue #97). The test is nested, not `cadence > 0 .and. mod(...)`, because Fortran does not short-circuit.
   class(prism_time_object), intent(in) :: self          !< Time handler.
   integer(I4P),             intent(in) :: cadence       !< Cadence in time steps.
   logical                              :: is_on_cadence !< Check result.

   is_on_cadence = .false.
   if (cadence > 0_I4P) is_on_cadence = (mod(self%it, cadence) == 0_I4P)
   endfunction is_on_cadence

   function is_to_save(self, cadence) result(is_to)
   !< Return true if data are to save: every `cadence` steps, at `it == it_max` and at the end of the run.
   !<
   !< A non-positive cadence disables saving, end of run included (as FLUME, issue #97).
   class(prism_time_object), intent(in) :: self    !< Time handler.
   integer(I4P),             intent(in) :: cadence !< Save iterations frequency.
   logical                              :: is_to   !< Check result.

   is_to = .false.
   if (cadence <= 0_I4P) return
   if (self%is_on_cadence(cadence).or.self%it==self%it_max.or. &
      (((self%it_max<=0).and.(self%time>=self%time_max)).or. &
       ((self%it_max> 0).and.(self%it>=self%it_max)))) is_to = .true.
   endfunction is_to_save

   subroutine print_progress(self, nodes_number)
   !< Print simulation progress.
   class(prism_time_object), intent(in) :: self         !< Time handler.
   integer(I4P),             intent(in) :: nodes_number !< Nodes number, global blocks number.

   associate(r=>mpih%myrankstr, it=>self%it, time=>self%time, dt=>self%dt, it_max=>self%it_max, time_max=>self%time_max)
      print '(A)', r//''
      print '(A)', r//'t:            '//trim(str(it,.true.))
      print '(A)', r//'nodes number: '//trim(str(nodes_number, .true.))
      if (self%omega_dt > 0._R8P) then
         print '(A)', r//'time step:    '//trim(str(dt, .true.))//' ('//trim(self%dt_source)//', omega_max dt = '// &
                      trim(str(self%omega_dt, .true.))//')'
      else
         print '(A)', r//'time step:    '//trim(str(dt, .true.))//' ('//trim(self%dt_source)//')'
      endif
      print '(A)', r//'time:         '//trim(str(time, .true.))
   if (it_max <= 0) then
      print '(A)', r//'progress:     '//trim(str(int(time/time_max * 100), .true.))//'%'
   else
      print '(A)', r//'progress:     '//trim(str(int((it*1._R8P)/it_max * 100), .true.))//'%'
   endif
      print '(A)', r//''
   endassociate
   endsubroutine print_progress
endmodule adam_prism_time_object
