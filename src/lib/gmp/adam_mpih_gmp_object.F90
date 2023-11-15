!< ADAM, MPI handler class definition, GMP backend.
module adam_mpih_gmp_object
!< ADAM, MPI handler class definition, GMP backend.
!<
!< Extend common mpih class adding GMP features.

use adam_mpih_object
use penf
use mpi
use omp_lib, only : omp_set_default_device, omp_get_initial_device

implicit none
private
public :: mpih_gmp_object

type, extends(mpih_object) :: mpih_gmp_object
   !< MPI handler class, GMP backend.
   integer(I4P) :: mydev=0_I4P      !< My GPU rank.
   integer(I4P) :: myhos=0_I4P      !< Host ID.
   integer(I4P) :: local_comm=0_I4P !< Local communicator.
   contains
      ! override methods
      procedure, pass(self) :: initialize !< Initialize MPI handler data.
      ! public methods
      procedure, pass(self) :: description !< Return pretty-printed object description.
endtype mpih_gmp_object

contains
    ! override methods
   subroutine initialize(self, do_mpi_init, do_device_init)
   !< Initialize MPI handler data.
   class(mpih_gmp_object) , intent(inout)        :: self              !< MPI handler.
   logical,                 intent(in), optional :: do_mpi_init       !< Flag to activate MPI init call.
   logical,                 intent(in), optional :: do_device_init    !< Flag to activate device init call.

   call self%mpih_object%initialize(do_mpi_init=do_mpi_init)
   call self%print_message('mpih_gmp_object%initialize start')
   if (present(do_device_init)) then
      if (do_device_init) then
         call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, self%local_comm, self%error)
         call MPI_COMM_RANK(self%local_comm, self%mydev, self%error)
         call omp_set_default_device(self%mydev)
         self%myhos = omp_get_initial_device()
         self%memory_avail = 8._R8P ! to be generalized!
      endif
   endif
   call self%print_message('mpih_gmp_object%initialize finish')
   endsubroutine initialize

   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(mpih_gmp_object) , intent(in) :: self             !< MPI handler.
   character(len=:), allocatable       :: desc             !< Description.
   character(len=1), parameter         :: NL=new_line('a') !< New line character.

   desc = self%mpih_object%description()//NL
   desc = desc//self%myrankstr//'MPIH GMP main data'//NL
   desc = desc//self%myrankstr//'  mydev:      '//trim(str(self%mydev     ))//NL
   desc = desc//self%myrankstr//'  myhos:      '//trim(str(self%myhos     ))//NL
   desc = desc//self%myrankstr//'  local_comm: '//trim(str(self%local_comm))
   endfunction description
endmodule adam_mpih_gmp_object
