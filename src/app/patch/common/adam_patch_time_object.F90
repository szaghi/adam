!< ADAM, PATCH time handler class definition, common CPU backend.
module adam_patch_time_object
!< ADAM, PATCH time handler class definition, common CPU backend.

use adam_mpih_object
use finer
use penf

implicit none
private
public :: patch_time_object

character(len=4), parameter :: INI_SECTION_NAME="time" !< INI (config) file section name containing time configs.

character(9),  parameter, public :: SMOOTHING_MULTIGRID   ='MULTIGRID'    !< Smoothing multigrid parameter.
character(12), parameter, public :: SMOOTHING_GAUSS_SEIDEL='GAUSS-SEIDEL' !< Smoothing Gauss-Seidel parameter.
character(3),  parameter, public :: SMOOTHING_SOR         ='SOR'          !< Smoothing SOR parameter.
character(7),  parameter, public :: SMOOTHING_SOR_OMP     ='SOR-OMP'      !< Smoothing SOR-OpenMP parameter.

type :: patch_time_object
   !< PATCH time handler class definition, CPU backend.
   type(mpih_object) :: mpih              !< MPI handler.
   integer(I4P)      :: it_max=-1_I4P     !< Maximum number of integration time steps.
   real(R8P)         :: time_max=1._R8P   !< Maximum integration time.
   real(R8P)         :: CFL=0.3_R8P       !< CFL time limit.
   integer(I4P)      :: it=0_I4P          !< Time steps counter.
   real(R8P)         :: time=0._R8P       !< Time.
   real(R8P)         :: dt=0.0001_R8P     !< Maximum time step accordingly to CFL criterion.
   character(:), allocatable :: smoothing !< Iterative smoothing method.
   contains
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize time handler.
      procedure, pass(self) :: is_to_save     !< Return true if data are to save.
      procedure, pass(self) :: load_from_file !< Load config from file.
      procedure, pass(self) :: print_progress !< Print simulation progress.
endtype patch_time_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(patch_time_object), intent(in) :: self             !< Time handler.
   character(len=:), allocatable        :: desc             !< Description.
   character(len=1), parameter          :: NL=new_line('a') !< New line character.

   desc =       self%mpih%myrankstr//'Time main data'//NL
   desc = desc//self%mpih%myrankstr//'  it_max:    '//trim(str(self%it_max  ))//NL
   desc = desc//self%mpih%myrankstr//'  time_max:  '//trim(str(self%time_max))//NL
   desc = desc//self%mpih%myrankstr//'  CFL:       '//trim(str(self%CFL     ))//NL
   desc = desc//self%mpih%myrankstr//'  smoothing: '//self%smoothing
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize time handler.
   class(patch_time_object), intent(inout) :: self            !< Time handler.
   type(file_ini),           intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'patch_time_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'patch_time_object%initialize finish'
   endsubroutine initialize

   function is_to_save(self, it_save)
   !< Return true if slices are to save.
   class(patch_time_object), intent(inout) :: self       !< Time handler.
   integer(I4P),             intent(in)    :: it_save    !< Save iterations frequency.
   logical                                 :: is_to_save !< Check result.

   is_to_save = .false.
   if (mod(self%it,it_save)==0.or.self%it==self%it_max.or.   &
      (((self%it_max<=0).and.(self%time>=self%time_max)).or. &
       ((self%it_max> 0).and.(self%it>=self%it_max)))) is_to_save = .true.
   endfunction is_to_save

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(patch_time_object), intent(inout)        :: self            !< Time handler.
   type(file_ini),           intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                  intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                        :: go_on_fail_     !< Go on if load fails.
   character(99)                                  :: buffer          !< Buffer.
   integer(I4P)                                   :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='it_max', val=self%it_max, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(it_max)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='time_max', val=self%time_max, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(time_max)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='CFL', val=self%CFL, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(CFL)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='smoothing', val=buffer, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(smoothing)')
   select case(trim(adjustl(buffer)))
   case('MULTIGRID','multigrid','Multigrid')
      self%smoothing = SMOOTHING_MULTIGRID
   case('GAUSS-SEIDEL','gauss-seidel','Gauss-Seidel')
      self%smoothing = SMOOTHING_GAUSS_SEIDEL
   case('SOR','sor','Sor')
      self%smoothing = SMOOTHING_SOR
   case('SOR-OMP','sor-omp','Sor-omp')
      self%smoothing = SMOOTHING_SOR_OMP
   endselect
   endsubroutine load_from_file

   subroutine print_progress(self, nodes_number)
   !< Print simulation progress.
   class(patch_time_object), intent(in) :: self         !< Time handler.
   integer(I4P),             intent(in) :: nodes_number !< Nodes number, global blocks number.

   associate(r=>self%mpih%myrankstr, it=>self%it, time=>self%time, dt=>self%dt, it_max=>self%it_max, time_max=>self%time_max)
      print '(A)', r//''
      print '(A)', r//'t:            '//trim(str(it,.true.))
      print '(A)', r//'nodes number: '//trim(str(nodes_number, .true.))
      print '(A)', r//'time step:    '//trim(str(dt, .true.))
      print '(A)', r//'time:         '//trim(str(time, .true.))
   if (it_max <= 0) then
      print '(A)', r//'progress:     '//trim(str(int(time/time_max * 100), .true.))//'%'
   else
      print '(A)', r//'progress:     '//trim(str(int((it*1._R8P)/it_max * 100), .true.))//'%'
   endif
      print '(A)', r//''
   endassociate
   endsubroutine print_progress
endmodule adam_patch_time_object
