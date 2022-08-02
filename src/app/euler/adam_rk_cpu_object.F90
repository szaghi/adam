!< ADAM, Runge Kutta class definition, CPU backend.
module adam_rk_cpu_object
!< ADAM, Runge Kutta class definition, CPU backend.

use adam_mpih_object, only : mpih_object
use finer
use penf

implicit none
private
public :: rk_cpu_object

character(len=7), parameter :: INI_SECTION_NAME='schemes' !< INI (config) file section name containing RK configs.

type :: rk_cpu_object
   !< Runge Kutta class definition, CPU backend.
   type(mpih_object)         :: mpih      !< MPI handler.
   character(:), allocatable :: scheme    !< RK scheme, ['rk-11', 'rk-33', 'rk-43', 'rk-54']
   integer(I4P)              :: nrk=1_I4P !< Runge-Kutta stages number.
   real(R8P), allocatable    :: ark(:)    !< RK a coefficients.
   real(R8P), allocatable    :: brk(:)    !< RK b coefficients.
   real(R8P), allocatable    :: crk(:)    !< RK c coefficients.
   contains
      ! public methods
      procedure, pass(self) :: compute_stage  !< Compute RK stage.
      procedure, pass(self) :: initialize     !< Initialize RK.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype rk_cpu_object

contains
   ! public methods
   subroutine compute_stage(self, ni, nj, nk, ngc, nv, blocks_number, dt, s, phi, q_n, rq, q)
   !< Compute RK stage.
   class(rk_cpu_object), intent(in)    :: self                            !< RK.
   integer(I4P),         intent(in)    :: ni                              !< Grid cells number in I direction.
   integer(I4P),         intent(in)    :: nj                              !< Grid cells number in J direction.
   integer(I4P),         intent(in)    :: nk                              !< Grid cells number in K direction.
   integer(I4P),         intent(in)    :: ngc                             !< Ghost cells number.
   integer(I4P),         intent(in)    :: nv                              !< Number of conservative variables.
   integer(I4P),         intent(in)    :: blocks_number                   !< Number of blocks.
   real(R8P),            intent(in)    :: dt                              !< Time step.
   integer(I4P),         intent(in)    :: s                               !< RK stage.
   real(R8P),            intent(in)    :: phi(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< IB distance function.
   real(R8P),            intent(in)    :: q_n(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field at time n.
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
               if (phi(1,i,j,k,b) < 0._R8P) then
                  q(:,i,j,k,b) = ark * q_n(:,i,j,k,b) + brk * q(:,i,j,k,b) + dtcrk * rq(:,i,j,k,b)
               endif
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine compute_stage

   subroutine initialize(self, file_parameters)
   !< Initialize the equation.
   class(rk_cpu_object), intent(inout) :: self            !< RK.
   type(file_ini),       intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'rk_cpu_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   allocate(self%ark(self%nrk), self%brk(self%nrk), self%crk(self%nrk))
   select case(self%nrk)
   case (1_I4P) ! euler 1st order
       self%ark(1) = 1._R8P ; self%brk(1) = 0._R8P ; self%crk(1) = 1._R8P
   case (2_I4P) ! 2nd order TVD
       self%ark(1) = 1._R8P  ; self%brk(1) = 0._R8P  ; self%crk(1) = 1._R8P
       self%ark(2) = 0.5_R8P ; self%brk(2) = 0.5_R8P ; self%crk(2) = 0.5_R8P
   case (3_I4P) ! 3rd order TVD
       self%ark(1) = 1._R8P        ; self%brk(1) = 0._R8P        ; self%crk(1) = 1._R8P
       self%ark(2) = 0.75_R8P      ; self%brk(2) = 0.25_R8P      ; self%crk(2) = 0.25_R8P
       self%ark(3) = 1._R8P/3._R8P ; self%brk(3) = 2._R8P/3._R8P ; self%crk(3) = 2._R8P/3._R8P
   endselect
   print '(A)', self%mpih%myrankstr//'rk_cpu_object%initialize finish'
   endsubroutine initialize

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
   if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//INI_SECTION_NAME//'].(runge_kutta)'
   self%scheme = trim(adjustl(char_buff))
   select case(trim(adjustl(self%scheme)))
   case('rk-11')
      self%nrk = 1
   case('rk-33')
      self%nrk = 3
   case('rk-43')
      self%nrk = 4
   case('rk-54')
      self%nrk = 5
   endselect
   endsubroutine load_from_file
endmodule adam_rk_cpu_object
