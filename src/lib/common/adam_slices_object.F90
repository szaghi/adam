!< ADAM, slices (of domain) class definition, CPU backend.
module adam_slices_object
   !< ADAM, slices (of domain) class definition, CPU backend.

use adam_adam_object, only : adam_object
use adam_mpih_object, only : mpih_object
use finer
use penf

implicit none
private
public :: slices_object

character(len=6), parameter :: INI_SECTION_NAME='slices'!< INI (config) file section name containing slices configs.

type :: slice_object
   !< Single slice class definition, CPU backend.
   character(99)          :: itype           !< Slice interpolation type.
   integer(I4P)           :: n_save          !< Iteration interval between subsequent data-slice saves.
   integer(I4P)           :: nijk(3)         !< Slice number of points.
   real(R8P)              :: emin(3)         !< Slice minimum extents.
   real(R8P)              :: emax(3)         !< Slice maximum extents.
   real(R8P), allocatable :: points(:,:,:,:) !< Slice points coordinates [3,ni,nj,nk].
endtype slice_object

type :: slices_object
   !< Slices (of domain) class definition, CPU backend.
   type(mpih_object)               :: mpih            !< MPI handler.
   integer(I4P)                    :: slices_number=0 !< Number of slices to be saved.
   type(slice_object), allocatable :: slice(:)        !< Slices data.
   contains
      ! public methods
      procedure, pass(self) :: initialize     !< Initialize RK.
      procedure, pass(self) :: is_to_save     !< Return true if slices are to save.
      procedure, pass(self) :: load_from_file !< Load config from file.
      procedure, pass(self) :: save_mat       !< Save simulation data slices in mat format.
endtype slices_object

contains
   ! public methods
   subroutine initialize(self, file_parameters)
   !< Initialize the equation.
   class(slices_object), intent(inout) :: self            !< Slices.
   type(file_ini),       intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   call self%mpih%print_message('slices_object%initialize start')
   call self%load_from_file(file_parameters=file_parameters)
   call self%mpih%print_message('slices_object%initialize finish')
   endsubroutine initialize

   pure function is_to_save(self, it, it_max, time, time_max)
   !< Return true if slices are to save.
   class(slices_object), intent(in) :: self       !< Slices.
   integer(I4P),         intent(in) :: it         !< Time step iteration.
   integer(I4P),         intent(in) :: it_max     !< Time step iteration max.
   real(R8P),            intent(in) :: time       !< Time iteration.
   real(R8P),            intent(in) :: time_max   !< Time iteration max.
   logical                          :: is_to_save !< Check result.
   integer(I4P)                     :: s          !< Slices counter.

   is_to_save = .false.
   if (self%slices_number>0) then
      slices_loop : do s=1, self%slices_number
         if (it>0) then
            if (mod(it,self%slice(s)%n_save)==0.or.it==it_max.or.&
               (((it_max <= 0).and.(time >= time_max)).or.((it>=it_max).and.(it_max > 0)))) then
               is_to_save = .true.
               exit slices_loop
            endif
         endif
      enddo slices_loop
   endif
   endfunction is_to_save

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(slices_object), intent(inout)        :: self            !< Slices.
   type(file_ini),       intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,              intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                    :: go_on_fail_     !< Go on if load fails.
   character(999)                             :: sname           !< Section name.
   real(R8P)                                  :: dxyz(3)         !< Space steps.
   integer(I4P)                               :: i, j, k, s      !< Counter.
   integer(I4P)                               :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='slices_number', val=self%slices_number, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(slices_number)')
   if (self%slices_number > 0) then
      allocate(self%slice(self%slices_number))
      do s=1, self%slices_number
         sname = 'slice_'//trim(str(s,.true.))
         call file_parameters%get(section_name=sname, option_name='itype',  val=self%slice(s)%itype)
         call file_parameters%get(section_name=sname, option_name='n_save', val=self%slice(s)%n_save)
         call file_parameters%get(section_name=sname, option_name='ni',     val=self%slice(s)%nijk(1))
         call file_parameters%get(section_name=sname, option_name='nj',     val=self%slice(s)%nijk(2))
         call file_parameters%get(section_name=sname, option_name='nk',     val=self%slice(s)%nijk(3))
         call file_parameters%get(section_name=sname, option_name='emin_x', val=self%slice(s)%emin(1))
         call file_parameters%get(section_name=sname, option_name='emin_y', val=self%slice(s)%emin(2))
         call file_parameters%get(section_name=sname, option_name='emin_z', val=self%slice(s)%emin(3))
         call file_parameters%get(section_name=sname, option_name='emax_x', val=self%slice(s)%emax(1))
         call file_parameters%get(section_name=sname, option_name='emax_y', val=self%slice(s)%emax(2))
         call file_parameters%get(section_name=sname, option_name='emax_z', val=self%slice(s)%emax(3))
         allocate(self%slice(s)%points(3,self%slice(s)%nijk(1),self%slice(s)%nijk(2),self%slice(s)%nijk(3)))
         dxyz(1) = (self%slice(s)%emax(1) - self%slice(s)%emin(1)) / self%slice(s)%nijk(1)
         dxyz(2) = (self%slice(s)%emax(2) - self%slice(s)%emin(2)) / self%slice(s)%nijk(2)
         dxyz(3) = (self%slice(s)%emax(3) - self%slice(s)%emin(3)) / self%slice(s)%nijk(3)
         do k=1, self%slice(s)%nijk(3)
            do j=1, self%slice(s)%nijk(2)
               do i=1, self%slice(s)%nijk(1)
                  self%slice(s)%points(1,i,j,k) = self%slice(s)%emin(1) + (i - 0.5_R8P) * dxyz(1)
                  self%slice(s)%points(2,i,j,k) = self%slice(s)%emin(2) + (j - 0.5_R8P) * dxyz(2)
                  self%slice(s)%points(3,i,j,k) = self%slice(s)%emin(3) + (k - 0.5_R8P) * dxyz(3)
               enddo
            enddo
         enddo
      enddo
   endif
   endsubroutine load_from_file

   subroutine save_mat(self, basename, it, it_max, time, time_max, adam, q, q_name)
   !< Save simulation data slices in mat format.
   class(slices_object), intent(inout)        :: self      !< Slices.
   character(*),         intent(in)           :: basename  !< Output file basename.
   integer(I4P),         intent(in)           :: it        !< Time step iteration.
   integer(I4P),         intent(in)           :: it_max    !< Time step iteration max.
   real(R8P),            intent(in)           :: time      !< Time iteration.
   real(R8P),            intent(in)           :: time_max  !< Time iteration max.
   type(adam_object),    intent(inout)        :: adam      !< Adam object.
   real(R8P),            intent(in)           :: q(1:,               &
                                                   1-adam%grid%ngc:, &
                                                   1-adam%grid%ngc:, &
                                                   1-adam%grid%ngc:, &
                                                   1:)     !< Field variables.
   character(*),         intent(in), optional :: q_name(:) !< Variables names.
   integer(I4P)                               :: s         !< Slices counter.

   if (self%slices_number>0) then
      do s=1, self%slices_number
         if (it>0) then
            if (mod(it,self%slice(s)%n_save)==0.or.it==it_max.or.&
               (((it_max <= 0).and.(time >= time_max)).or.((it>=it_max).and.(it_max > 0)))) then
               call self%mpih%barrier(tictoc=.true.)
               call self%mpih%print_message('save slice n: '//trim(str(s,.true.))//&
                                            ', t: '//trim(str(it,.true.))//', time: '//trim(str(time,.true.)))
               call adam%save_slice(points=self%slice(s)%points,                                                &
                                    itype=trim(self%slice(s)%itype),                                            &
                                    basename=trim(basename)//'-slice_'//trim(strz(s,2))//'-'//trim(strz(it,9)), &
                                    q=q,                                                                        &
                                    q_name=q_name)
               call self%mpih%barrier(tictoc=.true.)
            endif
         endif
      enddo
   endif
   endsubroutine save_mat
endmodule adam_slices_object
