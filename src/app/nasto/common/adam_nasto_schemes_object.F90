!< ADAM, Navier-Stokes schemes handler class definition, CPU backend.
module adam_nasto_schemes_object
!< ADAM, Navier-Stokes schemes handler class definition, CPU backend.

use adam_mpih_object
use adam_nasto_parameters
use finer
use penf

implicit none
private
public :: nasto_schemes_object
public :: SCHEME_TIME_RK_1
public :: SCHEME_TIME_RK_2
public :: SCHEME_TIME_RK_3
public :: SCHEME_FCONV_WENO_UPWIND
public :: SCHEME_FCONV_WENO_CENTRAL_2
public :: SCHEME_FCONV_WENO_CENTRAL_4
public :: SCHEME_FCONV_WENO_CENTRAL_6
public :: SCHEME_FDIFF_CENTRAL_2
public :: SCHEME_FDIFF_CENTRAL_4
public :: SCHEME_FDIFF_CENTRAL_6

character(len=7), parameter :: INI_SECTION_NAME="schemes" !< INI (config) file section name containing time configs.

character(len=13), parameter :: SCHEME_TIME_RK_1           ="runge-kutta-1" !< Parameter of time scheme, Runge-Kutta 1.
character(len=13), parameter :: SCHEME_TIME_RK_2           ="runge-kutta-2" !< Parameter of time scheme, Runge-Kutta 2.
character(len=13), parameter :: SCHEME_TIME_RK_3           ="runge-kutta-3" !< Parameter of time scheme, Runge-Kutta 3.
character(len=11), parameter :: SCHEME_FCONV_WENO_UPWIND   ="weno-upwind"   !< Parameter of WENO upwind fluxes convective scheme.
character(len=14), parameter :: SCHEME_FCONV_WENO_CENTRAL_2="weno-central-2"!< Parameter of WENO central 2 fluxes convective scheme.
character(len=14), parameter :: SCHEME_FCONV_WENO_CENTRAL_4="weno-central-4"!< Parameter of WENO central 4 fluxes convective scheme.
character(len=14), parameter :: SCHEME_FCONV_WENO_CENTRAL_6="weno-central-6"!< Parameter of WENO central 6 fluxes convective scheme.
character(len=9),  parameter :: SCHEME_FDIFF_CENTRAL_2     ="central-2"     !< Parameter of central 2 fluxes diffusive scheme.
character(len=9),  parameter :: SCHEME_FDIFF_CENTRAL_4     ="central-4"     !< Parameter of central 4 fluxes diffusive scheme.
character(len=9),  parameter :: SCHEME_FDIFF_CENTRAL_6     ="central-6"     !< Parameter of central 6 fluxes diffusive scheme.

type :: nasto_schemes_object
   !< NASTO schemes handler class definition.
   type(mpih_object) :: mpih                              !< MPI handler.
   character(:), allocatable :: time                      !< Scheme for time integration.
   character(:), allocatable :: fluxes_convective         !< Scheme for computing conv fluxes (weno-upwind/weno-central-2/4/6).
   character(:), allocatable :: fluxes_diffusive          !< Scheme for computing diff fluxes.
   integer(I4P)              :: ror_number=0_I4P          !< Number of ROR iterations
   integer(I4P), allocatable :: ror_schemes(:)            !< Scheme for each ROR iteration (4=weno-7/3=weno-5/2=weno-3/1=weno-1).
   real(R8P)                 :: ror_threshold=0.9_R8P     !< ROR threshold triggering
   integer(I4P)              :: ror_vars_number=2         !< Number of variables to check in ROR iterations.
   integer(I4P), allocatable :: ror_ivar(:)               !< Index of each variable to check in ROR iterations.
   logical                   :: enable_ror_stats=.false.  !< Enable ror statistic saving.
   integer(I4P)              :: ib_reduction_extent=0_I4P !< Extent of order reduction close to IB solids.
   integer(I4P)              :: ib_reduced_order=1        !< Reduced order close to IB solids (4=weno-7/3=weno-5/2=weno-3/1=weno-1).
   integer(I4P)              :: lmax=2_I4P                !< Central convective half stencil.
   integer(I4P)              :: iweno=2_I4P               !< WENO (half) stencil lenght.
   real(R8P), allocatable    :: fc_coeff(:,:)             !< Convective fluxes integration coefficients.
   real(R8P), allocatable    :: fd_coeff1(:)              !< Diffusive fluxes integration coefficients, first order.
   real(R8P), allocatable    :: fd_coeff2(:)              !< Diffusive fluxes integration coefficients, second order.
   integer(I4P)              :: nrk=4_I4P                 !< Runge-Kutta stages number.
   real(R8P), allocatable    :: ark(:)                    !< Runge-Kutta alpha coefficients.
   real(R8P), allocatable    :: brk(:)                    !< Runge-Kutta beta coefficients.
   real(R8P), allocatable    :: crk(:)                    !< Runge-Kutta beta coefficients.
   ! cell-centered arrays
   integer(I4P), allocatable :: ror_stats(:,:,:,:,:)   !< ROR statistics.
   integer(I4P), allocatable :: cell_scheme(:,:,:,:,:) !< Local-cell WENO scheme: iweno everywhere, but modified close to solids.
   contains
      procedure, pass(self) :: allocate_cellc_arrays   !< Allocate cell-centered arrays.
      procedure, pass(self) :: description             !< Return pretty-printed object description.
      procedure, pass(self) :: initialize              !< Initialize schemes handler.
      procedure, pass(self) :: initialize_coefficients !< Initialize fluxes integration coefficients.
      procedure, pass(self) :: initialize_runge_kutta  !< Initialize Runge-Kutta data.
      procedure, pass(self) :: load_from_file          !< Load config from file.
endtype nasto_schemes_object

contains
   ! public methods
   subroutine allocate_cellc_arrays(self, nb, ngc, ni, nj, nk)
   !< Allocate cell-centered arrays.
   class(nasto_schemes_object), intent(inout) :: self !< Schemes handler.
   integer(I4P),                intent(in)    :: nb   !< Total blocks number for MPI.
   integer(I4P),                intent(in)    :: ngc  !< Number of ghost cells.
   integer(I4P),                intent(in)    :: ni   !< Number of cells in i direction.
   integer(I4P),                intent(in)    :: nj   !< Number of cells in j direction.
   integer(I4P),                intent(in)    :: nk   !< Number of cells in k direction.

   allocate(self%cell_scheme(1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:3))
   self%cell_scheme = self%iweno
   if (self%enable_ror_stats) allocate(self%ror_stats(1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:3))
   endsubroutine allocate_cellc_arrays

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(nasto_schemes_object), intent(in) :: self             !< Schemes handler.
   character(len=:), allocatable           :: desc             !< Description.
   character(len=1), parameter             :: NL=new_line('a') !< New line character.

   desc =       self%mpih%myrankstr//'Schemes main data'//NL
   if (allocated(self%time)) &
   desc = desc//self%mpih%myrankstr//'  time:                '//         self%time                 //NL
   if (allocated(self%fluxes_convective)) &
   desc = desc//self%mpih%myrankstr//'  fluxes convective:   '//         self%fluxes_convective    //NL
   if (allocated(self%fluxes_diffusive)) &
   desc = desc//self%mpih%myrankstr//'  fluxes diffusive:    '//         self%fluxes_diffusive     //NL
   desc = desc//self%mpih%myrankstr//'  ror number:          '//trim(str(self%ror_number         ))//NL
   if (allocated(self%ror_schemes)) &
   desc = desc//self%mpih%myrankstr//'  ror schemes:         '//trim(str(self%ror_schemes        ))//NL
   desc = desc//self%mpih%myrankstr//'  ror threshold:       '//trim(str(self%ror_threshold      ))//NL
   desc = desc//self%mpih%myrankstr//'  ror vars number      '//trim(str(self%ror_vars_number    ))//NL
   if (allocated(self%ror_ivar)) &
   desc = desc//self%mpih%myrankstr//'  ror ivar:            '//trim(str(self%ror_ivar           ))//NL
   desc = desc//self%mpih%myrankstr//'  enable ror stats:    '//trim(str(self%enable_ror_stats   ))//NL
   desc = desc//self%mpih%myrankstr//'  ib reduction extent: '//trim(str(self%ib_reduction_extent))//NL
   desc = desc//self%mpih%myrankstr//'  ib reduced order:    '//trim(str(self%ib_reduced_order   ))//NL
   desc = desc//self%mpih%myrankstr//'  lmax:                '//trim(str(self%lmax               ))//NL
   desc = desc//self%mpih%myrankstr//'  iweno:               '//trim(str(self%iweno              ))//NL
   ! if (allocated(self%fc_coeff)) &
   ! desc = desc//self%mpih%myrankstr//'  fc coeff:            '//trim(str(self%fc_coeff           ))//NL
   if (allocated(self%fd_coeff1)) &
   desc = desc//self%mpih%myrankstr//'  fd coeff1:           '//trim(str(self%fd_coeff1          ))//NL
   if (allocated(self%fd_coeff2)) &
   desc = desc//self%mpih%myrankstr//'  fd coeff2:           '//trim(str(self%fd_coeff2          ))
   endfunction description

   subroutine initialize(self, file_parameters, nb, ngc, ni, nj, nk)
   !< Initialize time handler.
   class(nasto_schemes_object), intent(inout) :: self            !< Schemes handler.
   type(file_ini),              intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   integer(I4P),                intent(in)    :: nb              !< Total blocks number for MPI.
   integer(I4P),                intent(in)    :: ngc             !< Number of ghost cells.
   integer(I4P),                intent(in)    :: ni              !< Number of cells in i direction.
   integer(I4P),                intent(in)    :: nj              !< Number of cells in j direction.
   integer(I4P),                intent(in)    :: nk              !< Number of cells in k direction.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'nasto_schemes_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   call self%allocate_cellc_arrays(nb=nb, ngc=ngc, ni=ni, nj=nj, nk=nk)
   call self%initialize_coefficients
   call self%initialize_runge_kutta
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'nasto_schemes_object%initialize finish'
   endsubroutine initialize

   subroutine initialize_coefficients(self)
   !< Initialize fluxes integration coefficients.
   class(nasto_schemes_object), intent(inout) :: self !< Schemes handler.

   ! convective fluxes coefficients
   select case(self%fluxes_convective)
   case(SCHEME_FCONV_WENO_CENTRAL_2,SCHEME_FCONV_WENO_CENTRAL_4,SCHEME_FCONV_WENO_CENTRAL_6)
      allocate(self%fc_coeff(4,4))
      self%fc_coeff(1,1) = 1._R8P/2._R8P

      self%fc_coeff(1,2) =  2._R8P/3._R8P
      self%fc_coeff(2,2) = -1._R8P/12._R8P

      self%fc_coeff(1,3) =  3._R8P/4._R8P
      self%fc_coeff(2,3) = -3._R8P/20._R8P
      self%fc_coeff(3,3) =  1._R8P/60._R8P

      self%fc_coeff(1,4) =  4._R8P/5._R8P
      self%fc_coeff(2,4) = -1._R8P/5._R8P
      self%fc_coeff(3,4) =  4._R8P/105._R8P
      self%fc_coeff(4,4) = -1._R8P/280._R8P
   endselect

   ! diffusive fluxes coefficients
   allocate(self%fd_coeff1(3), self%fd_coeff2(0:3))
   select case(self%fluxes_diffusive)
   case(SCHEME_FDIFF_CENTRAL_2)
                                            self%fd_coeff2(0) = -2._R8P
      self%fd_coeff1(1) = 0.5_R8P         ; self%fd_coeff2(1) =  1._R8P
   case(SCHEME_FDIFF_CENTRAL_4)
                                            self%fd_coeff2(0) = -2.5_R8P
      self%fd_coeff1(1) = 2._R8P/3._R8P   ; self%fd_coeff2(1) = 4._R8P/3._R8P
      self%fd_coeff1(2) = -1._R8P/12._R8P ; self%fd_coeff2(2) = -1._R8P/12._R8P
   case(SCHEME_FDIFF_CENTRAL_6)
                                            self%fd_coeff2(0) = -245._R8P/90._R8P
      self%fd_coeff1(1) = 0.75_R8P        ; self%fd_coeff2(1) = 1.5_R8P
      self%fd_coeff1(2) = -0.15_R8P       ; self%fd_coeff2(2) = -0.15_R8P
      self%fd_coeff1(3) = 1._R8P/60._R8P  ; self%fd_coeff2(3) = 1._R8P/90._R8P
   endselect
   endsubroutine initialize_coefficients

   subroutine initialize_runge_kutta(self)
   !< Initialize Runge-Kutta data.
   class(nasto_schemes_object), intent(inout) :: self !< Schemes handler.

   select case(self%time)
   case(SCHEME_TIME_RK_1)
      self%nrk = 1
   case(SCHEME_TIME_RK_2)
      self%nrk = 2
   case(SCHEME_TIME_RK_3)
      self%nrk = 3
   endselect
   allocate(self%ark(self%nrk), self%brk(self%nrk), self%crk(self%nrk))
   select case(self%nrk)
      case(1_I4P) ! first order, Euler
         self%ark(1) = 1d0  ; self%brk(1) = 0d0; self%crk(1) = 1d0
      case(2_I4P) ! second order TVD
         self%ark(1) = 1d0    ; self%brk(1) = 0d0  ; self%crk(1) = 1d0
         self%ark(2) = 0.5d0  ; self%brk(2) = 0.5d0; self%crk(2) = 0.5d0
      case(3_I4P) ! third order TVD
         self%ark(1) = 1d0     ; self%brk(1) = 0d0     ; self%crk(1) = 1d0
         self%ark(2) = 0.75d0  ; self%brk(2) = 0.25d0  ; self%crk(2) = 0.25d0
         self%ark(3) = 1d0/3d0 ; self%brk(3) = 2d0/3d0 ; self%crk(3) = 2d0/3d0
   endselect
   endsubroutine initialize_runge_kutta

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(nasto_schemes_object), intent(inout)        :: self            !< Schemes handler.
   type(file_ini),              intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                     intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                           :: go_on_fail_     !< Go on if load fails.
   character(99)                                     :: buff_c          !< Character buffer.
   character(:), allocatable                         :: sname           !< Section name.
   character(:), allocatable                         :: oname           !< Option name.
   integer(I4P)                                      :: error           !< Error status.
   integer(I4P)                                      :: r               !< Counter.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='time', val=buff_c, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(time)')
   self%time = trim(adjustl(buff_c))
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='fluxes_convective', val=buff_c, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(fluxes_convective)')
   self%fluxes_convective = trim(adjustl(buff_c))
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='fluxes_diffusive', val=buff_c, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(fluxes_diffusive)')
   self%fluxes_diffusive = trim(adjustl(buff_c))

   select case(self%fluxes_convective)
   case(SCHEME_FCONV_WENO_CENTRAL_2)
      self%lmax = 1_I4P
   case(SCHEME_FCONV_WENO_CENTRAL_4)
      self%lmax = 2_I4P
   case(SCHEME_FCONV_WENO_CENTRAL_6)
      self%lmax = 3_I4P
   case(SCHEME_FCONV_WENO_UPWIND)
      sname = INI_SECTION_NAME//'_weno_upwind'
      call file_parameters%get(section_name=sname, option_name='ror_number', val=self%ror_number, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(ror_number)')
      call file_parameters%get(section_name=sname, option_name='iweno', val=self%iweno, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(iweno)')
      if (self%ror_number>0) then
         allocate(self%ror_schemes(self%ror_number))
         do r=1, self%ror_number
            oname = 'ror_scheme_'//trim(str(r,.true.))
            call file_parameters%get(section_name=sname, option_name=oname, val=self%ror_schemes(r), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].('//oname//')')
         enddo
         self%iweno = IWENO_FROM_SCHEME(self%ror_schemes(1))
      endif
      call file_parameters%get(section_name=sname, option_name='ror_threshold', val=self%ror_threshold, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(ror_threshold)')
      call file_parameters%get(section_name=sname, option_name='ror_vars_number', val=self%ror_vars_number, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(ror_vars_number)')
      if (self%ror_vars_number>0) then
         allocate(self%ror_ivar(self%ror_vars_number))
         do r=1, self%ror_vars_number
            oname = 'ror_ivar_'//trim(str(r,.true.))
            call file_parameters%get(section_name=sname, option_name=oname, val=self%ror_ivar(r), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].('//oname//')')
         enddo
      endif
      call file_parameters%get(section_name=sname, option_name='enable_ror_stats', val=self%enable_ror_stats, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(enable_ror_stats)')
      call file_parameters%get(section_name=sname, option_name='ib_reduction_extent', val=self%ib_reduction_extent, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(ib_reduction_extent)')
      call file_parameters%get(section_name=sname, option_name='ib_reduced_order', val=self%ib_reduced_order, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(ib_reduced_order)')
   endselect
   endsubroutine load_from_file
endmodule adam_nasto_schemes_object
