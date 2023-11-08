!< ADAM, WENO class definition.
module adam_weno_object
!< ADAM, WENO class definition.
!< Detailed reference of this implementation in contained into NASA/CR-97-206253 ICASE Report No. 97-65 "Essentially
!< Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes for Hyperbolic Conservation Laws" of Chi-Wang Shu (1997).
!< The stencil is arragnged as the following where for the sake of simplicity the case S=2 is sketched:
!<
!<```
!<                                                ___________________
!<                                               | right stencil 1   |
!<                              ___________________________          |
!<                             |      right stencil 0      |         |
!<                             |___________________________|         |
!<                             |       left stencil 1      |         |
!<                    _________|_________________          |         |
!<                   |    left stencil 0         |         |         |
!<...|--0--|--1--|...|--(i-1)--|-------(i)-------|--(i+1)--|--(i+2)--|...|--ni--|--ni+1--|...
!<                             ^                 ^
!<                             vl(i)             vr(i)
!<                             v_{i-1/2}^+       v_{i+1/2}^-
!<```
!< The non TBP procedures [[weno_reconstruct_optimal]] and [[weno_reconstruct_upwind]] are also provided for use without class
!< object in device backend kernels.

use adam_mpih_object
use finer
use penf

implicit none
save
private
public :: weno_object
public :: weno_reconstruct_optimal
public :: weno_reconstruct_upwind
public :: S_max
public :: S_max_m1

integer(I4P), parameter :: S_max=5    !< Maximum number/dimensions of stencils.
integer(I4P), parameter :: S_max_m1=4 !< Maximum number/dimensions of stencils minus 1.

character(len=4), parameter :: INI_SECTION_NAME="weno" !< INI (config) file section name containing time configs.

type :: weno_object
   !< WENO class definition.
   type(mpih_object)         :: mpih                      !< MPI handler.
   integer(I4P)              :: S                         !< Stencils number/dimensions, 2S-1 order of accuracy.
   real(R8P), allocatable    :: a(:,:,:)                  !< Optimal weights                    [1:2,0:S-1,1:S].
   real(R8P), allocatable    :: p(:,:,:,:)                !< Polinomials coefficients           [1:2,0:S-1,0:S-1,1:S].
   real(R8P), allocatable    :: d(:,:,:,:)                !< Smoothness indicators coefficients [0:S-1,0:S-1,0:S-1,1:S].
   real(R8P)                 :: zeps                      !< Parameter for avoiding division by zero in computing IS.
   integer(I4P)              :: sodd                      !< Branching between odd and even number of stencils (mod(S,2)).
   integer(I4P)              :: wexp                      !< Exponent for growing the diffusive part of weights.
   integer(I4P)              :: ror_number=0_I4P          !< Number of ROR iterations.
   integer(I4P), allocatable :: ror_schemes(:)            !< Scheme (S value) for each ROR iteration.
   real(R8P)                 :: ror_threshold=0.9_R8P     !< ROR threshold triggering.
   integer(I4P)              :: ror_vars_number=2         !< Number of variables to check in ROR iterations.
   integer(I4P), allocatable :: ror_ivar(:)               !< Index of each variable to check in ROR iterations.
   logical                   :: enable_ror_stats=.false.  !< Enable ror statistic saving.
   integer(I4P)              :: ib_reduction_extent=0_I4P !< Extent of order reduction close to IB solids.
   integer(I4P)              :: ib_reduced_order=1        !< Reduced order (S value) close to IB solids.
   integer(I4P), allocatable :: ror_stats(:,:,:,:,:)      !< ROR statistics.
   integer(I4P), allocatable :: cell_scheme(:,:,:,:,:)    !< Local-cell WENO scheme: S everywhere, but modified close to solids.
   contains
      ! public methods
      procedure, pass(self) :: description         !< Return pretty-printed object description.
      procedure, pass(self) :: initialize          !< Initialize class.
      procedure, pass(self) :: load_from_file      !< Load config from file.
      procedure, pass(self) :: reconstruct_optimal !< Return WENO optimal reconstruction of 2S-1 order.
      procedure, pass(self) :: reconstruct_upwind  !< Return WENO upwind  reconstruction of 2S-1 order.
      ! private methods
      procedure, pass(self) :: initialize_S1       !< Initialize coefficients for S=1.
      procedure, pass(self) :: initialize_S2       !< Initialize coefficients for S=2.
      procedure, pass(self) :: initialize_S3       !< Initialize coefficients for S=3.
      procedure, pass(self) :: initialize_S4       !< Initialize coefficients for S=4.
      procedure, pass(self) :: initialize_S5       !< Initialize coefficients for S=5.
      procedure, pass(self) :: initialize_S6       !< Initialize coefficients for S=6.
      procedure, pass(self) :: compute_polynomials !< Compute WENO polynomials.
      procedure, pass(self) :: compute_weights     !< Compute WENO weights.
endtype weno_object
contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(weno_object), intent(in) :: self             !< WENO object.
   character(len=:), allocatable  :: desc             !< Description.
   character(len=1), parameter    :: NL=new_line('a') !< New line character.
   integer(I4P)                   :: s1, s2           !< Counter.

   desc =       self%mpih%myrankstr//'WENO upwind scheme main data'//NL
   desc = desc//self%mpih%myrankstr//'  S:                              '//trim(str(self%S                  ))//NL
   desc = desc//self%mpih%myrankstr//'  zeps:                           '//trim(str(self%zeps               ))//NL
   desc = desc//self%mpih%myrankstr//'  sodd:                           '//trim(str(self%sodd               ))//NL
   desc = desc//self%mpih%myrankstr//'  wexp:                           '//trim(str(self%wexp               ))//NL
   if (allocated(self%a)) then
   do s2=0, self%S -1
   desc = desc//self%mpih%myrankstr//'  a(:,'//trim(str(s2,.true.))//'):'//trim(str(self%a(:,s2,self%S)     ))//NL
   enddo
   endif
   if (allocated(self%p)) then
   do s2=0, self%S -1
   do s1=0, self%S -1
   desc = desc//self%mpih%myrankstr//'  p(:,'//trim(str(s1,.true.))//','// &
                                               trim(str(s2,.true.))//'):'//trim(str(self%p(:,s1,s2,self%S)  ))//NL
   enddo
   enddo
   endif
   if (allocated(self%d)) then
   do s2=0, self%S -1
   do s1=0, self%S -1
   desc = desc//self%mpih%myrankstr//'  d(:,'//trim(str(s1,.true.))//','// &
                                               trim(str(s2,.true.))//'):'//trim(str(self%d(:,s1,s2,self%S)  ))//NL
   enddo
   enddo
   endif
   desc = desc//self%mpih%myrankstr//'  ror number:                     '//trim(str(self%ror_number         ))//NL
   if (allocated(self%ror_schemes)) &
   desc = desc//self%mpih%myrankstr//'  ror schemes:                    '//trim(str(self%ror_schemes        ))//NL
   desc = desc//self%mpih%myrankstr//'  ror threshold:                  '//trim(str(self%ror_threshold      ))//NL
   desc = desc//self%mpih%myrankstr//'  ror vars number                 '//trim(str(self%ror_vars_number    ))//NL
   if (allocated(self%ror_ivar)) &
   desc = desc//self%mpih%myrankstr//'  ror ivar:                       '//trim(str(self%ror_ivar           ))//NL
   desc = desc//self%mpih%myrankstr//'  enable ror stats:               '//trim(str(self%enable_ror_stats   ))//NL
   desc = desc//self%mpih%myrankstr//'  ib reduction extent:            '//trim(str(self%ib_reduction_extent))//NL
   desc = desc//self%mpih%myrankstr//'  ib reduced order:               '//trim(str(self%ib_reduced_order   ))
   endfunction description

   subroutine initialize(self, file_parameters, S, nb, ngc, ni, nj, nk)
   !< Initialize class.
   class(weno_object), intent(inout)        :: self            !< WENO object.
   type(file_ini),     intent(in), optional :: file_parameters !< Simulation parameters ini file handler.
   integer(I4P),       intent(in), optional :: S               !< Number of stencils used.
   integer(I4P),       intent(in)           :: nb              !< Total blocks number for MPI.
   integer(I4P),       intent(in)           :: ngc             !< Number of ghost cells.
   integer(I4P),       intent(in)           :: ni              !< Number of cells in i direction.
   integer(I4P),       intent(in)           :: nj              !< Number of cells in j direction.
   integer(I4P),       intent(in)           :: nk              !< Number of cells in k direction.

   call self%mpih%initialize(do_mpi_init=.false.)
   call self%mpih%print_message('weno_object%initialize start')
   if (present(file_parameters)) then
      call self%load_from_file(file_parameters=file_parameters)
   elseif (present(S)) then
      self%S = S
   else
      call self%mpih%error_stop(msg=': failed to initialize weno object, one between file parameters and S number must be passed')
   endif
   ! initialize coefficients for all schemes up to S
   if (allocated(self%a)) deallocate(self%a) ; allocate(self%a(1:2,0:self%S-1,1:self%S))
   if (allocated(self%p)) deallocate(self%p) ; allocate(self%p(1:2,0:self%S-1,0:self%S-1,1:self%S))
   if (allocated(self%d)) deallocate(self%d) ; allocate(self%d(0:self%S-1,0:self%S-1,0:self%S-1,1:self%S))
   select case(self%S)
   case(1) ! 1st order, godunov
      call self%initialize_S1
   case(2) ! 3rd order
      call self%initialize_S1
      call self%initialize_S2
   case(3) ! 5th order
      call self%initialize_S1
      call self%initialize_S2
      call self%initialize_S3
   case(4) ! 7th order
      call self%initialize_S1
      call self%initialize_S2
      call self%initialize_S3
      call self%initialize_S4
   case(5) ! 9th order
      call self%initialize_S1
      call self%initialize_S2
      call self%initialize_S3
      call self%initialize_S4
      call self%initialize_S5
   ! case(6) ! 11th order
   !    call self%initialize_S1
   !    call self%initialize_S2
   !    call self%initialize_S3
   !    call self%initialize_S4
   !    call self%initialize_S5
   !    call self%initialize_S6
   case default
      call self%mpih%error_stop(msg=': failed to initialize weno object, S must be in [1,'//trim(str(S_max,.true.))//']')
   endselect
   ! set constants
   self%wexp = self%S
   if (self%S>4) self%wexp = self%S - 1
   self%sodd = mod(self%S,2)
   self%zeps = 1.0e-6_R8P
   ! allocate cell-centered arrays
   allocate(self%cell_scheme(1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:3))
   self%cell_scheme = self%S
   if (self%enable_ror_stats) allocate(self%ror_stats(1:nb, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:3))
   print '(A)', self%description()
   call self%mpih%print_message('weno_object%initialize finish')
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(weno_object), intent(inout)        :: self            !< WENO object.
   type(file_ini),     intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,            intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                  :: go_on_fail_     !< Go on if load fails.
   character(99)                            :: buff_c          !< Character buffer.
   character(:), allocatable                :: sname           !< Section name.
   character(:), allocatable                :: oname           !< Option name.
   integer(I4P)                             :: error           !< Error status.
   integer(I4P)                             :: r               !< Counter.

   ! call self%weno%initialize(S=self%iweno, dxyz_min=1.0e-8_R8P)
   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   sname = INI_SECTION_NAME
   call file_parameters%get(section_name=sname, option_name='S', val=self%S, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(S)')
   call file_parameters%get(section_name=sname, option_name='ror_number', val=self%ror_number, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(ror_number)')
   if (self%ror_number>0) then
      allocate(self%ror_schemes(self%ror_number))
      do r=1, self%ror_number
         oname = 'ror_scheme_'//trim(str(r,.true.))
         call file_parameters%get(section_name=sname, option_name=oname, val=self%ror_schemes(r), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].('//oname//')')
      enddo
      self%S = self%ror_schemes(1) ! first ROR schemes must be the highest order used
   else
      ! allocate anyway a single element array for backends compatibility
      allocate(self%ror_schemes(1))
      self%ror_schemes = self%S
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
   else
      ! allocate anyway a single element array for backends compatibility
      allocate(self%ror_ivar(1))
      ! self%ror_ivar(1) = 1
   endif
   call file_parameters%get(section_name=sname, option_name='enable_ror_stats', val=self%enable_ror_stats, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(enable_ror_stats)')
   call file_parameters%get(section_name=sname, option_name='ib_reduction_extent', val=self%ib_reduction_extent, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(ib_reduction_extent)')
   call file_parameters%get(section_name=sname, option_name='ib_reduced_order', val=self%ib_reduced_order, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(ib_reduced_order)')
   endsubroutine load_from_file

   pure subroutine reconstruct_optimal(self, S, V, VR)
   !< Reconstruct by WENO method with optimal weigths (without smoothness indicators computations) of 2S-1 order.
   class(weno_object), intent(in)  :: self             !< WENO object.
   integer(I4P),       intent(in)  :: S                !< Number of stencils used.
   real(R8P),          intent(in)  :: V (1:2,1-S:-1+S) !< Variable to be reconstructed.
   real(R8P),          intent(out) :: VR(1:2         ) !< Left and right (1,2) interface value of reconstructed V.

   call weno_reconstruct_optimal(S=S, weno_a=self%a, weno_p=self%p, V=V, VR=VR)
   endsubroutine reconstruct_optimal

   pure subroutine reconstruct_upwind(self, S, V, VR)
   !< Reconstruct by WENO upwind method of 2S-1 order.
   class(weno_object), intent(in)  :: self             !< WENO object.
   integer(I4P),       intent(in)  :: S                !< Number of stencils used.
   real(R8P),          intent(in)  :: V (1:2,1-S:-1+S) !< Variables to be reconstructed.
   real(R8P),          intent(out) :: VR(1:2         ) !< Left and right (1,2) interface value of reconstructed V.

   call weno_reconstruct_upwind(S=S, weno_a=self%a, weno_p=self%p, weno_d=self%d, weno_zeps=self%zeps, V=V, VR=VR)
   endsubroutine reconstruct_upwind

   ! private methods
   subroutine initialize_S1(self)
   !< Initialize coefficients for S=1.
   class(weno_object), intent(inout) :: self !< WENO object.

   ! optimal weights
   self%a = 1._R8P

   ! polinomials coefficients
   self%p = 1._R8P

   ! smoothness indicators coefficients
   self%d = 1._R8P
   endsubroutine initialize_S1

   subroutine initialize_S2(self)
   !< Initialize coefficients for S=2, 2S-1=3rd order.
   class(weno_object), intent(inout) :: self !< WENO object.
   integer(I4P), parameter           :: S=2  !< Number of stencils used.

   associate(a=>self%a, p=>self%p, d=>self%d)
     ! optimal weights
     ! 1 => left interface (i-1/2)
     a(1,0,S) = 2._R8P/3._R8P ! stencil 0
     a(1,1,S) = 1._R8P/3._R8P ! stencil 1
     ! 2 => right interface (i+1/2)
     a(2,0,S) = 1._R8P/3._R8P ! stencil 0
     a(2,1,S) = 2._R8P/3._R8P ! stencil 1

     ! polinomials coefficients
     ! 1 => left interface (i-1/2)
     ! cell  0             ;   cell  1
     p(1,0,0,S) =  0.5_R8P ; p(1,1,0,S) =  0.5_R8P ! stencil 0
     p(1,0,1,S) = -0.5_R8P ; p(1,1,1,S) =  1.5_R8P ! stencil 1
     ! 2 => right interface (i+1/2)
     ! cell  0             ;   cell  1
     p(2,0,0,S) =  1.5_R8P ; p(2,1,0,S) = -0.5_R8P ! stencil 0
     p(2,0,1,S) =  0.5_R8P ; p(2,1,1,S) =  0.5_R8P ! stencil 1

     ! smoothness indicators coefficients
     ! stencil 0
     ! i*i               ;   (i-1)*i
     d(0,0,0,S) = 1._R8P ; d(1,0,0,S) =-2._R8P
     ! /                 ;   (i-1)*(i-1)
     d(0,1,0,S) = 0._R8P ; d(1,1,0,S) = 1._R8P
     ! stencil 1
     ! (i+1)*(i+1)       ;   (i+1)*i
     d(0,0,1,S) = 1._R8P ; d(1,0,1,S) =-2._R8P
     ! /                 ;   i*i
     d(0,1,1,S) = 0._R8P ; d(1,1,1,S) = 1._R8P
   endassociate
   endsubroutine initialize_S2

   subroutine initialize_S3(self)
   !< Initialize coefficients for S=3, 2S-1=5th order.
   class(weno_object), intent(inout) :: self !< WENO object.
   integer(I4P), parameter           :: S=3  !< Number of stencils used.

   associate(a=>self%a, p=>self%p, d=>self%d)
     ! optimal weights
     ! 1 => left interface (i-1/2)
     a(1,0,S) = 0.3_R8P ! stencil 0
     a(1,1,S) = 0.6_R8P ! stencil 1
     a(1,2,S) = 0.1_R8P ! stencil 2
     ! 2 => right interface (i+1/2)
     a(2,0,S) = 0.1_R8P ! stencil 0
     a(2,1,S) = 0.6_R8P ! stencil 1
     a(2,2,S) = 0.3_R8P ! stencil 2

     ! polinomials coefficients
     ! 1 => left interface (i-1/2)
     ! cell  0                   ;   cell  1                   ;   cell  2
     p(1,0,0,S) =  1._R8P/3._R8P ; p(1,1,0,S) =  5._R8P/6._R8P ; p(1,2,0,S) = -1._R8P/6._R8P ! stencil 0
     p(1,0,1,S) = -1._R8P/6._R8P ; p(1,1,1,S) =  5._R8P/6._R8P ; p(1,2,1,S) =  1._R8P/3._R8P ! stencil 1
     p(1,0,2,S) =  1._R8P/3._R8P ; p(1,1,2,S) = -7._R8P/6._R8P ; p(1,2,2,S) = 11._R8P/6._R8P ! stencil 2
     ! 2 => right interface (i+1/2)
     ! cell  0                   ;   cell  1                   ;   cell  2
     p(2,0,0,S) = 11._R8P/6._R8P ; p(2,1,0,S) = -7._R8P/6._R8P ; p(2,2,0,S) =  1._R8P/3._R8P ! stencil 0
     p(2,0,1,S) =  1._R8P/3._R8P ; p(2,1,1,S) =  5._R8P/6._R8P ; p(2,2,1,S) = -1._R8P/6._R8P ! stencil 1
     p(2,0,2,S) = -1._R8P/6._R8P ; p(2,1,2,S) =  5._R8P/6._R8P ; p(2,2,2,S) =  1._R8P/3._R8P ! stencil 2

     ! smoothness indicators coefficients
     ! stencil 0
     !      i*i                   ;   (i-1)*i                    ;   (i-2)*i
     d(0,0,0,S) =  10._R8P/3._R8P ; d(1,0,0,S) = -31._R8P/3._R8P ; d(2,0,0,S) =  11._R8P/3._R8P
     ! /                          ;   (i-1)*(i-1)                ;   (i-2)*(i-1)
     d(0,1,0,S) =   0._R8P        ; d(1,1,0,S) =  25._R8P/3._R8P ; d(2,1,0,S) = -19._R8P/3._R8P
     ! /                          ;   /                          ;   (i-2)*(i-2)
     d(0,2,0,S) =   0._R8P        ; d(1,2,0,S) =   0._R8P        ; d(2,2,0,S) =   4._R8P/3._R8P
     ! stencil 1
     ! (i+1)*(i+1)                ;   i*(i+1)                    ;   (i-1)*(i+1)
     d(0,0,1,S) =   4._R8P/3._R8P ; d(1,0,1,S) = -13._R8P/3._R8P ; d(2,0,1,S) =   5._R8P/3._R8P
     ! /                          ;   i*i                        ;   (i-1)*i
     d(0,1,1,S) =   0._R8P        ; d(1,1,1,S) =  13._R8P/3._R8P ; d(2,1,1,S) = -13._R8P/3._R8P
     ! /                          ;   /                          ;   (i-1)*(i-1)
     d(0,2,1,S) =   0._R8P        ; d(1,2,1,S) =   0._R8P        ; d(2,2,1,S) =   4._R8P/3._R8P
     ! stencil 2
     ! (i+2)*(i+2)                ;   (i+1)*(i+2)                ;   i*(i+2)
     d(0,0,2,S) =   4._R8P/3._R8P ; d(1,0,2,S) = -19._R8P/3._R8P ; d(2,0,2,S) =  11._R8P/3._R8P
     ! /                          ;   (i+1)*(i+1)                ;   i*(i+1)
     d(0,1,2,S) =   0._R8P        ; d(1,1,2,S) =  25._R8P/3._R8P ; d(2,1,2,S) = -31._R8P/3._R8P
     ! /                          ;   /                          ;   i*i
     d(0,2,2,S) =   0._R8P        ; d(1,2,2,S) =   0._R8P        ; d(2,2,2,S) =  10._R8P/3._R8P
   endassociate
   endsubroutine initialize_S3

   subroutine initialize_S4(self)
   !< Initialize coefficients for S=4, 2S-1=7th order..
   class(weno_object), intent(inout) :: self !< WENO object.
   integer(I4P), parameter           :: S=4  !< Number of stencils used.

   associate(a=>self%a, p=>self%p, d=>self%d)
      ! optimal weights
      ! 1 => left interface (i-1/2)
      a(1,0,S) =  4._R8P/35._R8P ! stencil 0
      a(1,1,S) = 18._R8P/35._R8P ! stencil 1
      a(1,2,S) = 12._R8P/35._R8P ! stencil 2
      a(1,3,S) =  1._R8P/35._R8P ! stencil 3
      ! 2 => right interface (i+1/2)
      a(2,0,S) =  1._R8P/35._R8P ! stencil 0
      a(2,1,S) = 12._R8P/35._R8P ! stencil 1
      a(2,2,S) = 18._R8P/35._R8P ! stencil 2
      a(2,3,S) =  4._R8P/35._R8P ! stencil 3

      ! polinomials coefficients
      ! 1 => left interface (i-1/2)
      ! cell  0                  ;   cell  1                  ;   cell  2                   ;   cell  3
      p(1,0,0,S)= 1._R8P/4._R8P  ; p(1,1,0,S)=13._R8P/12._R8P ; p(1,2,0,S)=-5._R8P /12._R8P ; p(1,3,0,S)= 1._R8P /12._R8P! sten 0
      p(1,0,1,S)=-1._R8P/12._R8P ; p(1,1,1,S)= 7._R8P/12._R8P ; p(1,2,1,S)=  7._R8P/12._R8P ; p(1,3,1,S)=-1._R8P /12._R8P! sten 1
      p(1,0,2,S)= 1._R8P/12._R8P ; p(1,1,2,S)=-5._R8P/12._R8P ; p(1,2,2,S)= 13._R8P/12._R8P ; p(1,3,2,S)= 1._R8P /4._R8P ! sten 2
      p(1,0,3,S)=-1._R8P/4._R8P  ; p(1,1,3,S)=13._R8P/12._R8P ; p(1,2,3,S)=-23._R8P/12._R8P ; p(1,3,3,S)= 25._R8P/12._R8P! sten 3
      ! 2 => right interface (i+1/2)
      ! cell  0                  ;   cell  1                   ;   cell  2                  ;   cell  3
      p(2,0,0,S)=25._R8P/12._R8P ; p(2,1,0,S)=-23._R8P/12._R8P ; p(2,2,0,S)=13._R8P/12._R8P ; p(2,3,0,S)=-1._R8P/4._R8P ! sten 0
      p(2,0,1,S)= 1._R8P/4._R8P  ; p(2,1,1,S)= 13._R8P/12._R8P ; p(2,2,1,S)=-5._R8P/12._R8P ; p(2,3,1,S)= 1._R8P/12._R8P! sten 1
      p(2,0,2,S)=-1._R8P/12._R8P ; p(2,1,2,S)= 7._R8P /12._R8P ; p(2,2,2,S)= 7._R8P/12._R8P ; p(2,3,2,S)=-1._R8P/12._R8P! sten 2
      p(2,0,3,S)= 1._R8P/12._R8P ; p(2,1,3,S)=-5._R8P /12._R8P ; p(2,2,3,S)=13._R8P/12._R8P ; p(2,3,3,S)= 1._R8P/4._R8P ! sten 3

      ! smoothness indicators coefficients
      ! stencil 0
      ! i*i                 ;  (i-1)*i               ;  (i-2)*i                ;  (i-3)*i
      d(0,0,0,S)= 2107._R8P ; d(1,0,0,S)=-9402._R8P  ; d(2,0,0,S)= 7042._R8P   ; d(3,0,0,S)=-1854._R8P
      !/                    ;  (i-1)*(i-1)           ;  (i-2)*(i-1)            ;  (i-3)*(i-1)
      d(0,1,0,S)= 0._R8P    ; d(1,1,0,S)= 11003._R8P ; d(2,1,0,S)=-17246._R8P  ; d(3,1,0,S)= 4642._R8P
      !/                    ;  /                     ;  (i-2)*(i-2)            ;  (i-3)*(i-2)
      d(0,2,0,S)= 0._R8P    ; d(1,2,0,S)= 0._R8P     ; d(2,2,0,S)= 7043._R8P   ; d(3,2,0,S)=-3882._R8P
      !/                    ;  /                     ;  /                      ;  (i-3)*(i-3)
      d(0,3,0,S)= 0._R8P    ; d(1,3,0,S)= 0._R8P     ; d(2,3,0,S)= 0._R8P      ; d(3,3,0,S)= 547._R8P
      ! stencil 1
      !(i+1)*(i+1)          ;   i*(i+1)              ;  (i-1)*(i+1)            ;  (i-2)*(i+1)
      d(0,0,1,S)= 547._R8P  ; d(1,0,1,S)=-2522._R8P  ; d(2,0,1,S)= 1922._R8P   ; d(3,0,1,S)=-494._R8P
      !/                    ;   i*i                  ;  (i-1)*i                ;  (i-2)*i
      d(0,1,1,S)= 0._R8P    ; d(1,1,1,S)= 3443._R8P  ; d(2,1,1,S)=-5966._R8P   ; d(3,1,1,S)= 1602._R8P
      !/                    ;  /                     ;  (i-1)*(i-1)            ;  (i-2)*(i-1)
      d(0,2,1,S)= 0._R8P    ; d(1,2,1,S)= 0._R8P     ; d(2,2,1,S)= 2843._R8P   ; d(3,2,1,S)=-1642._R8P
      !/                    ;  /                     ;  /                      ;  (i-2)*(i-2)
      d(0,3,1,S)= 0._R8P    ; d(1,3,1,S)= 0._R8P     ; d(2,3,1,S)= 0._R8P      ; d(3,3,1,S)= 267._R8P
      ! stencil 2
      !(i+2)*(i+2)          ;  (i+1)*(i+2)           ;   i*(i+2)               ;  (i-1)*(i+2)
      d(0,0,2,S)= 267._R8P  ; d(1,0,2,S)=-1642._R8P  ; d(2,0,2,S)= 1602._R8P   ; d(3,0,2,S) =-494._R8P
      !/                    ;  (i+1)*(i+1)           ;   i*(i+1)               ;  (i-1)*(i+1)
      d(0,1,2,S)= 0._R8P    ; d(1,1,2,S)= 2843._R8P  ; d(2,1,2,S)=-5966._R8P   ; d(3,1,2,S) = 1922._R8P
      !/                    ;  /                     ;   i*i                   ;  (i-1)*i
      d(0,2,2,S)= 0._R8P    ; d(1,2,2,S)= 0._R8P     ; d(2,2,2,S)= 3443._R8P   ; d(3,2,2,S) =-2522._R8P
      !/                    ;  /                     ;  /                      ;  (i-1)*(i-1)
      d(0,3,2,S)= 0._R8P    ; d(1,3,2,S)= 0._R8P     ; d(2,3,2,S)= 0._R8P      ; d(3,3,2,S) = 547._R8P
      ! stencil 3
      !(i+3)*(i+3)          ;  (i+2)*(i+3)           ;  (i+1)*(i+3)            ;  i*(i+3)
      d(0,0,3,S)= 547._R8P  ; d(1,0,3,S)=-3882._R8P  ; d(2,0,3,S)= 4642._R8P   ; d(3,0,3,S)=-1854._R8P
      !/                    ;  (i+2)*(i+2)           ;  (i+1)*(i+2)            ;  i*(i+2)
      d(0,1,3,S)= 0._R8P    ; d(1,1,3,S)= 7043._R8P  ; d(2,1,3,S)= -17246._R8P ; d(3,1,3,S)= 7042._R8P
      !/                    ;  /                     ;  (i+1)*(i+1)            ;  i*(i+1)
      d(0,2,3,S)= 0._R8P    ; d(1,2,3,S)= 0._R8P     ; d(2,2,3,S)= 11003._R8P  ; d(3,2,3,S)=-9402._R8P
      !/                    ;  /                     ;  /                      ;  i*i
      d(0,3,3,S)= 0._R8P    ; d(1,3,3,S)= 0._R8P     ; d(2,3,3,S)= 0._R8P      ; d(3,3,3,S)= 2107._R8P
   endassociate
   endsubroutine initialize_S4

   subroutine initialize_S5(self)
   !< Initialize coefficients for S=5, 2S-1=9th order..
   class(weno_object), intent(inout) :: self !< WENO object.
   integer(I4P), parameter           :: S=5  !< Number of stencils used.

   associate(a=>self%a, p=>self%p, d=>self%d)
      ! optimal weights
      ! 1 => left interface (i-1/2)
      a(1,0,S) =  5._R8P/126._R8P ! stencil 0
      a(1,1,S) = 20._R8P/63._R8P  ! stencil 1
      a(1,2,S) = 10._R8P/21._R8P  ! stencil 2
      a(1,3,S) = 10._R8P/63._R8P  ! stencil 3
      a(1,4,S) =  1._R8P/126._R8P ! stencil 4
      ! 2 => right interface (i+1/2)
      a(2,0,S) =  1._R8P/126._R8P ! stencil 0
      a(2,1,S) = 10._R8P/63._R8P  ! stencil 1
      a(2,2,S) = 10._R8P/21._R8P  ! stencil 2
      a(2,3,S) = 20._R8P/63._R8P  ! stencil 3
      a(2,4,S) =  5._R8P/126._R8P ! stencil 4

      ! polinomials coefficients
      ! 1 => left interface (i-1/2)
      !  cell  0                ;    cell  1                ;    cell  2                ;    cell  3
      p(1,0,0,S)=   1._R8P/5._R8P ;p(1,1,0,S)=  77._R8P/60._R8P;p(1,2,0,S)= -43._R8P/60._R8P;p(1,3,0,S)=  17._R8P/60._R8P! stencil 0
      p(1,0,1,S)=  -1._R8P/20._R8P;p(1,1,1,S)=   9._R8P/20._R8P;p(1,2,1,S)=  47._R8P/60._R8P;p(1,3,1,S)= -13._R8P/60._R8P! stencil 1
      p(1,0,2,S)=   1._R8P/30._R8P;p(1,1,2,S)= -13._R8P/60._R8P;p(1,2,2,S)=  47._R8P/60._R8P;p(1,3,2,S)=   9._R8P/20._R8P! stencil 2
      p(1,0,3,S)=  -1._R8P/20._R8P;p(1,1,3,S)=  17._R8P/60._R8P;p(1,2,3,S)= -43._R8P/60._R8P;p(1,3,3,S)=  77._R8P/60._R8P! stencil 3
      p(1,0,4,S)=   1._R8P/5._R8P ;p(1,1,4,S)= -21._R8P/20._R8P;p(1,2,4,S)= 137._R8P/60._R8P;p(1,3,4,S)=-163._R8P/60._R8P! stencil 4
      !  cell  4
      p(1,4,0,S)=  -1._R8P/20._R8P  ! stencil 0
      p(1,4,1,S)=   1._R8P/30._R8P  ! stencil 1
      p(1,4,2,S)=  -1._R8P/20._R8P  ! stencil 2
      p(1,4,3,S)=   1._R8P/5._R8P   ! stencil 3
      p(1,4,4,S)= 137._R8P/60._R8P  ! stencil 4
      ! 2 => right interface (i+1/2)
      !  cell  0               ;    cell  1               ;   cell  2                ;    cell  3
      p(2,0,0,S)= 137._R8P/60._R8P;p(2,1,0,S)=-163._R8P/60._R8P;p(2,2,0,S)= 137._R8P/60._R8P;p(2,3,0,S)= -21._R8P/20._R8P! stencil 0
      p(2,0,1,S)=   1._R8P/5._R8P ;p(2,1,1,S)=  77._R8P/60._R8P;p(2,2,1,S)= -43._R8P/60._R8P;p(2,3,1,S)=  17._R8P/60._R8P! stencil 1
      p(2,0,2,S)=  -1._R8P/20._R8P;p(2,1,2,S)=   9._R8P/20._R8P;p(2,2,2,S)=  47._R8P/60._R8P;p(2,3,2,S)= -13._R8P/60._R8P! stencil 2
      p(2,0,3,S)=   1._R8P/30._R8P;p(2,1,3,S)= -13._R8P/60._R8P;p(2,2,3,S)=  47._R8P/60._R8P;p(2,3,3,S)=   9._R8P/20._R8P! stencil 3
      p(2,0,4,S)=  -1._R8P/20._R8P;p(2,1,4,S)=  17._R8P/60._R8P;p(2,2,4,S)= -43._R8P/60._R8P;p(2,3,4,S)=  77._R8P/60._R8P! stencil 4
      !  cell  4
      p(2,4,0,S)=   1._R8P/5._R8P  ! stencil 0
      p(2,4,1,S)=  -1._R8P/20._R8P ! stencil 1
      p(2,4,2,S)=   1._R8P/30._R8P ! stencil 2
      p(2,4,3,S)=  -1._R8P/20._R8P ! stencil 3
      p(2,4,4,S)=   1._R8P/5._R8P  ! stencil 4

      ! smoothness indicators coefficients
      ! stencil 0
      ! i*i                  ; (i-1)*i                ; (i-2)*i                ; (i-3)*i                ; (i-4)*i
      d(0,0,0,S)= 107918._R8P;d(1,0,0,S)=-649501._R8P ;d(2,0,0,S)= 758823._R8P ;d(3,0,0,S)=-411487._R8P ;d(4,0,0,S)= 86329._R8P
      !/                     ; (i-1)*(i-1)            ; (i-2)*(i-1)            ; (i-3)*(i-1)            ; (i-4)*(i-1)
      d(0,1,0,S)= 0._R8P     ;d(1,1,0,S)= 1020563._R8P;d(2,1,0,S)=-2462076._R8P;d(3,1,0,S)= 1358458._R8P;d(4,1,0,S)=-288007._R8P
      !/                     ; /                      ; (i-2)*(i-2)            ; (i-3)*(i-2)            ; (i-4)*(i-2)
      d(0,2,0,S)= 0._R8P     ;d(1,2,0,S)= 0._R8P      ;d(2,2,0,S)= 1521393._R8P;d(3,2,0,S)=-1704396._R8P;d(4,2,0,S)= 364863._R8P
      !/                     ; /                      ; /                      ; (i-3)*(i-3)            ; (i-4)*(i-3)
      d(0,3,0,S)= 0._R8P     ;d(1,3,0,S)= 0._R8P      ;d(2,3,0,S)= 0._R8P      ;d(3,3,0,S)= 482963._R8P ;d(4,3,0,S)=-208501._R8P
      !/                     ; /                      ; /                      ; /                      ; (i-4)*(i-4)
      d(0,4,0,S)= 0._R8P     ;d(1,4,0,S)= 0._R8P      ;d(2,4,0,S)= 0._R8P      ;d(3,4,0,S)= 0._R8P      ;d(4,4,0,S)= 22658._R8P
      ! stencil 1
      !(i+1)*(i+1)          ;  i*(i+1)              ; (i-1)*(i+1)           ; (i-2)*(i+1)           ; (i-3)*(i+1)
      d(0,0,1,S)= 22658._R8P;d(1,0,1,S)=-140251._R8P;d(2,0,1,S)= 165153._R8P;d(3,0,1,S)=-88297._R8P ;d(4,0,1,S)= 18079._R8P
      !/                    ;  i*i                  ; (i-1)*i               ; (i-2)*i               ; (i-3)*i
      d(0,1,1,S)= 0._R8P    ;d(1,1,1,S)= 242723._R8P;d(2,1,1,S)=-611976._R8P;d(3,1,1,S)= 337018._R8P;d(4,1,1,S)=-70237._R8P
      !/                    ; /                     ; (i-1)*(i-1)           ; (i-2)*(i-1)           ; (i-3)*(i-1)
      d(0,2,1,S)= 0._R8P    ;d(1,2,1,S)= 0._R8P     ;d(2,2,1,S)= 406293._R8P;d(3,2,1,S)=-464976._R8P;d(4,2,1,S)= 99213._R8P
      !/                    ; /                     ; /                     ; (i-2)*(i-2)           ; (i-3)*(i-2)
      d(0,3,1,S)= 0._R8P    ;d(1,3,1,S)= 0._R8P     ;d(2,3,1,S)= 0._R8P     ;d(3,3,1,S)= 138563._R8P;d(4,3,1,S)=-60871._R8P
      !/                    ; /                     ; /                     ; /                     ; (i-3)*(i-3)
      d(0,4,1,S)= 0._R8P    ;d(1,4,1,S)= 0._R8P     ;d(2,4,1,S)= 0._R8P     ;d(3,4,1,S)= 0._R8P     ;d(4,4,1,S)= 6908._R8P
      ! stencil 2
      !(i+2)*(i+2)         ; (i+1)*(i+2)           ; i*(i+2)               ; (i-1)*(i+2)           ; (i-2)*(i+2)
      d(0,0,2,S)= 6908._R8P;d(1,0,2,S)=-51001._R8P;d(2,0,2,S)= 67923._R8P;d(3,0,2,S)=-38947._R8P;d(4,0,2,S)= 8209._R8P
      !/                   ; (i+1)*(i+1)           ; i*(i+1)               ; (i-1)*(i+1)           ; (i-2)*(i+1)
      d(0,1,2,S)= 0._R8P   ;d(1,1,2,S)= 104963._R8P;d(2,1,2,S)=-299076._R8P;d(3,1,2,S)= 179098._R8P;d(4,1,2,S)=-38947._R8P
      !/                   ; /                     ; i*i                   ; (i-1)*i               ; (i-2)*i
      d(0,2,2,S)= 0._R8P   ;d(1,2,2,S)= 0._R8P     ;d(2,2,2,S)= 231153._R8P;d(3,2,2,S)=-299076._R8P;d(4,2,2,S)= 67923._R8P
      !/                   ; /                     ; /                     ; (i-1)*(i-1)           ; (i-2)*(i-1)
      d(0,3,2,S)= 0._R8P   ;d(1,3,2,S)= 0._R8P     ;d(2,3,2,S)= 0._R8P     ;d(3,3,2,S)= 104963._R8P;d(4,3,2,S)=-51001._R8P
      !/                   ; /                     ; /                     ; /                     ; (i-2)*(i-2)
      d(0,4,2,S)= 0._R8P   ;d(1,4,2,S)= 0._R8P     ;d(2,4,2,S)= 0._R8P     ;d(3,4,2,S)= 0._R8P     ;d(4,4,2,S)= 6908._R8P
      ! stencil 3
      !(i+3)*(i+3)         ; (i+2)*(i+3)           ; (i+1)*(i+3)           ; i*(i+3)               ; (i-1)*(i+3)
      d(0,0,3,S)= 6908._R8P;d(1,0,3,S)=-60871._R8P ;d(2,0,3,S)= 99213._R8P ;d(3,0,3,S)=-70237._R8P ;d(4,0,3,S)= 18079._R8P
      !/                   ; (i+2)*(i+2)           ; (i+1)*(i+2)           ; i*(i+2)               ; (i-1)*(i+2)
      d(0,1,3,S)= 0._R8P   ;d(1,1,3,S)= 138563._R8P;d(2,1,3,S)=-464976._R8P;d(3,1,3,S)= 337018._R8P;d(4,1,3,S)=-88297._R8P
      !/                   ; /                     ; (i+1)*(i+1)           ; i*(i+1)               ; (i-1)*(i+1)
      d(0,2,3,S)= 0._R8P   ;d(1,2,3,S)= 0._R8P     ;d(2,2,3,S)= 406293._R8P;d(3,2,3,S)=-611976._R8P;d(4,2,3,S)= 165153._R8P
      !/                   ; /                     ; /                     ; i*i                   ; (i-1)*i
      d(0,3,3,S)= 0._R8P   ;d(1,3,3,S)= 0._R8P     ;d(2,3,3,S)= 0._R8P     ;d(3,3,3,S)= 242723._R8P;d(4,3,3,S)=-140251._R8P
      !/                   ; /                     ; /                     ; /                     ; (i-1)*(i-1)
      d(0,4,3,S)= 0._R8P   ;d(1,4,3,S)= 0._R8P     ;d(2,4,3,S)= 0._R8P     ;d(3,4,3,S)= 0._R8P     ;d(4,4,3,S)= 22658._R8P
      ! stencil 4
      !(i+4)*(i+4)          ; (i+3)*(i+4)           ; (i+2)*(i+4)            ; (i+1)*(i+4)            ; i*(i+4)
      d(0,0,4,S)= 22658._R8P;d(1,0,4,S)=-208501._R8P;d(2,0,4,S)= 364863._R8P ;d(3,0,4,S)=-288007._R8P ;d(4,0,4,S)= 86329._R8P
      !/                    ; (i+3)*(i+3)           ; (i+2)*(i+3)            ; (i+1)*(i+3)            ; i*(i+3)
      d(0,1,4,S)= 0._R8P    ;d(1,1,4,S)= 482963._R8P;d(2,1,4,S)=-1704396._R8P;d(3,1,4,S)= 1358458._R8P;d(4,1,4,S)=-411487._R8P
      !/                    ; /                     ; (i+2)*(i+2)            ; (i+1)*(i+2)            ; i*(i+2)
      d(0,2,4,S)= 0._R8P    ;d(1,2,4,S)= 0._R8P     ;d(2,2,4,S)= 1521393._R8P;d(3,2,4,S)=-2462076._R8P;d(4,2,4,S)= 758823._R8P
      !/                    ; /                     ; /                      ; (i+1)*(i+1)            ; i*(i+1)
      d(0,3,4,S)= 0._R8P    ;d(1,3,4,S)= 0._R8P     ;d(2,3,4,S)= 0._R8P      ;d(3,3,4,S)= 1020563._R8P;d(4,3,4,S)=-649501._R8P
      !/                    ; /                     ; /                      ; /                      ; i*i
      d(0,4,4,S)= 0._R8P    ;d(1,4,4,S)= 0._R8P     ;d(2,4,4,S)= 0._R8P      ;d(3,4,4,S)= 0._R8P      ;d(4,4,4,S)= 107918._R8P
   endassociate
   endsubroutine initialize_S5

   subroutine initialize_S6(self)
   !< Initialize coefficients for S=6, 2S-1=11th order..
   class(weno_object), intent(inout) :: self !< WENO object.
   integer(I4P), parameter           :: S=6  !< Number of stencils used.

   associate(a=>self%a, p=>self%p, d=>self%d)
      ! optimal weights
      ! 1 => left interface (i-1/2)
      a(1,0,S) =   1._R8P/77._R8P  ! stencil 0
      a(1,1,S) =  25._R8P/154._R8P ! stencil 1
      a(1,2,S) = 100._R8P/231._R8P ! stencil 2
      a(1,3,S) =  25._R8P/77._R8P  ! stencil 3
      a(1,4,S) =   5._R8P/77._R8P  ! stencil 4
      a(1,5,S) =   1._R8P/462._R8P ! stencil 5
      ! 2 => right interface (i+1/2)
      a(2,0,S) =   1._R8P/462._R8P ! stencil 0
      a(2,1,S) =   5._R8P/77._R8P  ! stencil 1
      a(2,2,S) =  25._R8P/77._R8P  ! stencil 2
      a(2,3,S) = 100._R8P/231._R8P ! stencil 3
      a(2,4,S) =  25._R8P/154._R8P ! stencil 4
      a(2,5,S) =   1._R8P/77._R8P  ! stencil 5

      ! polinomials coefficients
      ! 1 => left interface (i-1/2)
      !  cell  0                ;    cell  1                ;    cell  2                ;    cell  3
      p(1,0,0,S)=   1._R8P/6._R8P ;p(1,1,0,S)=  29._R8P/20._R8P;p(1,2,0,S)= -21._R8P/20._R8P;p(1,3,0,S)=  37._R8P/60._R8P! stencil 0
      p(1,0,1,S)=  -1._R8P/30._R8P;p(1,1,1,S)=  11._R8P/30._R8P;p(1,2,1,S)=  19._R8P/20._R8P;p(1,3,1,S)= -23._R8P/60._R8P! stencil 1
      p(1,0,2,S)=   1._R8P/60._R8P;p(1,1,2,S)=  -2._R8P/15._R8P;p(1,2,2,S)=  37._R8P/60._R8P;p(1,3,2,S)=  37._R8P/60._R8P! stencil 2
      p(1,0,3,S)=  -1._R8P/60._R8P;p(1,1,3,S)=   7._R8P/60._R8P;p(1,2,3,S)= -23._R8P/60._R8P;p(1,3,3,S)=  19._R8P/20._R8P! stencil 3
      p(1,0,4,S)=   1._R8P/30._R8P;p(1,1,4,S)= -13._R8P/60._R8P;p(1,2,4,S)=  37._R8P/60._R8P;p(1,3,4,S)= -21._R8P/20._R8P! stencil 4
      p(1,0,5,S)=  -1._R8P/6._R8P ;p(1,1,5,S)=  31._R8P/30._R8P;p(1,2,5,S)=-163._R8P/60._R8P;p(1,3,5,S)=  79._R8P/20._R8P! stencil 5
      !  cell  4                ;    cell  5
      p(1,4,0,S)= -13._R8P/60._R8P; p(1,5,0,S)=   1._R8P/30._R8P  ! stencil 0
      p(1,4,1,S)=   7._R8P/60._R8P; p(1,5,1,S)=  -1._R8P/60._R8P  ! stencil 1
      p(1,4,2,S)=  -2._R8P/15._R8P; p(1,5,2,S)=   1._R8P/60._R8P  ! stencil 2
      p(1,4,3,S)=  11._R8P/30._R8P; p(1,5,3,S)=  -1._R8P/30._R8P  ! stencil 3
      p(1,4,4,S)=  29._R8P/20._R8P; p(1,5,4,S)=   1._R8P/6._R8P   ! stencil 4
      p(1,4,5,S)= -71._R8P/20._R8P; p(1,5,5,S)=  49._R8P/20._R8P  ! stencil 5
      ! 2 => right interface (i+1/2)
      !  cell  0                ;    cell  1                ;   cell  2                 ;    cell  3
      p(2,0,0,S)=  49._R8P/20._R8P;p(2,1,0,S)= -71._R8P/20._R8P;p(2,2,0,S)=  79._R8P/20._R8P;p(2,3,0,S)=-163._R8P/60._R8P! stencil 0
      p(2,0,1,S)=   1._R8P/6._R8P ;p(2,1,1,S)=  29._R8P/20._R8P;p(2,2,1,S)= -21._R8P/20._R8P;p(2,3,1,S)=  37._R8P/60._R8P! stencil 1
      p(2,0,2,S)=  -1._R8P/30._R8P;p(2,1,2,S)=  11._R8P/30._R8P;p(2,2,2,S)=  19._R8P/20._R8P;p(2,3,2,S)= -23._R8P/60._R8P! stencil 2
      p(2,0,3,S)=   1._R8P/60._R8P;p(2,1,3,S)=  -2._R8P/15._R8P;p(2,2,3,S)=  37._R8P/60._R8P;p(2,3,3,S)=  37._R8P/60._R8P! stencil 3
      p(2,0,4,S)=  -1._R8P/60._R8P;p(2,1,4,S)=   7._R8P/60._R8P;p(2,2,4,S)= -23._R8P/60._R8P;p(2,3,4,S)=  19._R8P/20._R8P! stencil 4
      p(2,0,5,S)=   1._R8P/30._R8P;p(2,1,5,S)= -13._R8P/60._R8P;p(2,2,5,S)=  37._R8P/60._R8P;p(2,3,5,S)= -21._R8P/20._R8P! stencil 5
      !  cell  4                ;    cell  5
      p(2,4,0,S)=  31._R8P/30._R8P; p(2,5,0,S)=  -1._R8P/6._R8P   ! stencil 0
      p(2,4,1,S)= -13._R8P/60._R8P; p(2,5,1,S)=   1._R8P/30._R8P  ! stencil 1
      p(2,4,2,S)=   7._R8P/60._R8P; p(2,5,2,S)=  -1._R8P/60._R8P  ! stencil 2
      p(2,4,3,S)=  -2._R8P/15._R8P; p(2,5,3,S)=   1._R8P/60._R8P  ! stencil 3
      p(2,4,4,S)=  11._R8P/30._R8P; p(2,5,4,S)=  -1._R8P/30._R8P  ! stencil 4
      p(2,4,5,S)=  29._R8P/20._R8P; p(2,5,5,S)=   1._R8P/6._R8P   ! stencil 5

      ! smoothness indicators coefficients
      ! stencil 0
      !                 i*i                  ;                (i-1)*i                 ;                 (i-2)*i
      d(0,0,0,S) =  6150211._R8P / 120960._R8P; d(1,0,0,S) =  -2966279._R8P /   7560._R8P; d(2,0,0,S) =   4762921._R8P /   7560._R8P
      !             (i-3)*i                  ;                (i-4)*i                 ;                 (i-5)*i
      d(3,0,0,S) =-15848531._R8P /  30240._R8P; d(4,0,0,S) =   2706017._R8P /  12096._R8P; d(5,0,0,S) =   -235637._R8P /   6048._R8P

      !                  /                   ;                (i-1)*(i-1)             ;                 (i-2)*(i-1)
      d(0,1,0,S) =        0._R8P              ; d(1,1,0,S) =  31617079._R8P /  40320._R8P; d(2,1,0,S) = -25980937._R8P /  10080._R8P
      !             (i-3)*(i-1)              ;                (i-4)*(i-1)             ;                 (i-5)*(i-1)
      d(3,1,0,S) = 32862709._R8P /  15120._R8P; d(4,1,0,S) =  -1048211._R8P /   1120._R8P; d(5,1,0,S) =    661145._R8P /   4032._R8P

      !                  /                   ;                     /                  ;                 (i-2)*(i-2)
      d(0,2,0,S) =        0._R8P              ; d(1,2,0,S) =         0._R8P              ; d(2,2,0,S) =  21703781._R8P /  10080._R8P
      !              (i-3)*(i-2)             ;                (i-4)*(i-2)             ;                 (i-5)*(i-2)
      d(3,2,0,S) = -6937561._R8P /   1890._R8P; d(4,2,0,S) =   2674951._R8P /   1680._R8P; d(5,2,0,S) =   -314063._R8P /   1120._R8P

      !                  /                   ;                     /                  ;                      /
      d(0,3,0,S) =        0._R8P              ; d(1,3,0,S) =         0._R8P              ; d(2,3,0,S) =         0._R8P
      !             (i-3)*(i-3)              ;                (i-4)*(i-3)             ;                 (i-5)*(i-3)
      d(3,3,0,S) = 47689393._R8P /  30240._R8P; d(4,3,0,S) = -41615261._R8P /  30240._R8P; d(5,3,0,S) =   1840141._R8P /   7560._R8P

      !                  /                   ;                     /                  ;                      /
      d(0,4,0,S) =        0._R8P              ; d(1,4,0,S) =         0._R8P              ; d(2,4,0,S) =         0._R8P
      !                  /                   ;                (i-4)*(i-4)             ;                 (i-5)*(i-4)
      d(3,4,0,S) =        0._R8P              ; d(4,4,0,S) =  12160229._R8P /  40320._R8P; d(5,4,0,S) =   -539591._R8P /   5040._R8P

      !                  /                   ;                     /                  ;                      /
      d(0,5,0,S) =        0._R8P              ; d(1,5,0,S) =         0._R8P              ; d(2,5,0,S) =         0._R8P
      !                  /                   ;                     /                  ;                 (i-5)*(i-5)
      d(3,5,0,S) =        0._R8P              ; d(4,5,0,S) =         0._R8P              ; d(5,5,0,S) =    384187._R8P /  40320._R8P
      ! stencil 1
      !             (i+1)*(i+1)              ;                    i*(i+1)             ;                 (i-1)*(i+1)
      d(0,0,1,S) =   384187._R8P /  40320._R8P; d(1,0,1,S) =  -1139749._R8P /  15120._R8P; d(2,0,1,S) =     61427._R8P /   504._R8P
      !             (i-2)*(i+1)              ;                (i-3)*(i+1)             ;                 (i-4)*(i+1)
      d(3,0,1,S) = -1015303._R8P /  10080._R8P; d(4,0,1,S) =   2567287._R8P /  60480._R8P; d(5,0,1,S) =    -73379._R8P / 10080._R8P

      !                  /                   ;                    i*i                 ;                 (i-1)*i
      d(0,1,1,S) =        0._R8P              ; d(1,1,1,S) =  19365967._R8P / 120960._R8P; d(2,1,1,S) = -16306061._R8P /  30240._R8P
      !             (i-2)*i                  ;                (i-3)*i                 ;                 (i-4)*i
      d(3,1,1,S) =  6881719._R8P /  15120._R8P; d(4,1,1,S) =  -5877617._R8P /  30240._R8P; d(5,1,1,S) =   2033509._R8P /  60480._R8P

      !                  /                   ;                     /                  ;                 (i-1)*(i-1)
      d(0,2,1,S) =        0._R8P              ; d(1,2,1,S) =         0._R8P              ; d(2,2,1,S) =   4721851._R8P /  10080._R8P
      !             (i-2)*(i-1)              ;                (i-3)*(i-1)             ;                 (i-4)*(i-1)
      d(3,2,1,S) =  -169859._R8P /    210._R8P; d(4,2,1,S) =   5300629._R8P /  15120._R8P; d(5,2,1,S) =    -68601._R8P /   1120._R8P

      !                  /                   ;                     /                  ;                      /
      d(0,3,1,S) =        0._R8P              ; d(1,3,1,S) =         0._R8P              ; d(2,3,1,S) =         0._R8P
      !             (i-2)*(i-2)              ;                (i-3)*(i-2)             ;                 (i-4)*(i-2)
      d(3,3,1,S) =  1197047._R8P /   3360._R8P; d(4,3,1,S) =  -9478331._R8P /  30240._R8P; d(5,3,1,S) =    139471._R8P /   2520._R8P

      !                  /                   ;                     /                  ;                      /
      d(0,4,1,S) =        0._R8P              ; d(1,4,1,S) =         0._R8P              ; d(2,4,1,S) =         0._R8P
      !                  /                   ;                (i-3)*(i-3)             ;                 (i-4)*(i-3)
      d(3,4,1,S) =        0._R8P              ; d(4,4,1,S) =   8449957._R8P / 120960._R8P; d(5,4,1,S) =   -188483._R8P /   7560._R8P

      !                  /                   ;                     /                  ;                      /
      d(0,5,1,S) =        0._R8P              ; d(1,5,1,S) =         0._R8P              ; d(2,5,1,S) =         0._R8P
      !                  /                   ;                     /                  ;                 (i-4)*(i-4)
      d(3,5,1,S) =        0._R8P              ; d(4,5,1,S) =         0._R8P              ; d(5,5,1,S) =     90593._R8P /  40320._R8P
      ! stencil 2
      !             (i+2)*(i+2)              ;                (i+1)*(i+2)             ;                     i*(i+2)
      d(0,0,2,S) =    90593._R8P /  40320._R8P; d(1,0,2,S) =     -1240._R8P /     63._R8P; d(2,0,2,S) =    255397._R8P /   7560._R8P
      !             (i-1)*(i+2)              ;                (i-2)*(i+2)             ;                 (i-3)*(i+2)
      d(3,0,2,S) =  -288521._R8P /  10080._R8P; d(4,0,2,S) =    243127._R8P /  20160._R8P; d(5,0,2,S) =    -12281._R8P /   6048._R8P

      !                  /                   ;                (i+1)*(i+1)             ;                     i*(i+1)
      d(0,1,2,S) =        0._R8P              ; d(1,1,2,S) =   1884439._R8P /  40320._R8P; d(2,1,2,S) =  -5106971._R8P /  30240._R8P
      !             (i-1)*(i+1)              ;                (i-2)*(i+1)             ;                (i-3)*(i+1)
      d(3,1,2,S) =   248681._R8P /   1680._R8P; d(4,1,2,S) =   -643999._R8P /  10080._R8P; d(5,1,2,S) =    662503._R8P /  60480._R8P

      !                  /                   ;                     /                  ;                    i*i
      d(0,2,2,S) =        0._R8P              ; d(1,2,2,S) =         0._R8P              ; d(2,2,2,S) =   4877743._R8P /  30240._R8P
      !             (i-1)*i                  ;                (i-2)*i                 ;                (i-3)*i
      d(3,2,2,S) =  -559651._R8P /   1890._R8P; d(4,2,2,S) =   1991239._R8P /  15120._R8P; d(5,2,2,S) =   -139633._R8P /   6048._R8P

      !                  /                   ;                     /                  ;                     /
      d(0,3,2,S) =        0._R8P              ; d(1,3,2,S) =         0._R8P              ; d(2,3,2,S) =         0._R8P
      !             (i-1)*(i-1)              ;                (i-2)*(i-1)             ;                (i-3)*(i-1)
      d(3,3,2,S) =   159219._R8P /   1120._R8P; d(4,3,2,S) =  -1323367._R8P /  10080._R8P; d(5,3,2,S) =    178999._R8P /   7560._R8P

      !                  /                   ;                     /                  ;                     /
      d(0,4,2,S) =        0._R8P              ; d(1,4,2,S) =         0._R8P              ; d(2,4,2,S) =         0._R8P
      !                  /                   ;                (i-2)*(i-2)             ;                (i-3)*(i-2)
      d(3,4,2,S) =        0._R8P              ; d(4,4,2,S) =    141661._R8P /   4480._R8P; d(5,4,2,S) =   -178747._R8P /  15120._R8P

      !                  /                   ;                     /                  ;                     /
      d(0,5,2,S) =        0._R8P              ; d(1,5,2,S) =         0._R8P              ; d(2,5,2,S) =         0._R8P
      !                  /                   ;                     /                  ;                (i-3)*(i-3)
      d(3,5,2,S) =        0._R8P              ; d(4,5,2,S) =         0._R8P              ; d(5,5,2,S) =    139633._R8P / 120960._R8P
      ! stencil 3
      !             (i+3)*(i+3)              ;                (i+2)*(i+3)             ;                 (i+1)*(i+3)
      d(0,0,3,S) =   139633._R8P / 120960._R8P; d(1,0,3,S) =   -178747._R8P /  15120._R8P; d(2,0,3,S) =    178999._R8P /   7560._R8P
      !                 i*(i+3)              ;                (i-1)*(i+3)             ;                 (i-2)*(i+3)
      d(3,0,3,S) =  -139633._R8P /   6048._R8P; d(4,0,3,S) =    662503._R8P /  60480._R8P; d(5,0,3,S) =    -12281._R8P /   6048._R8P

      !                  /                   ;                (i+2)*(i+2)             ;                 (i+1)*(i+2)
      d(0,1,3,S) =        0._R8P              ; d(1,1,3,S) =    141661._R8P /   4480._R8P; d(2,1,3,S) =  -1323367._R8P /  10080._R8P
      !                 i*(i+2)              ;                (i-1)*(i+2)             ;                 (i-2)*(i+2)
      d(3,1,3,S) =  1991239._R8P /  15120._R8P; d(4,1,3,S) =   -643999._R8P /  10080._R8P; d(5,1,3,S) =    243127._R8P /  20160._R8P

      !                  /                   ;                     /                  ;                    i*(i+1)
      d(0,2,3,S) =        0._R8P              ; d(1,2,3,S) =         0._R8P              ; d(2,2,3,S) =    159219._R8P /   1120._R8P
      !                 i*(i+1)              ;                (i-1)*(i+1)             ;                (i-2)*(i+1)
      d(3,2,3,S) =  -559651._R8P /   1890._R8P; d(4,2,3,S) =    248681._R8P /   1680._R8P; d(5,2,3,S) =   -288521._R8P /  10080._R8P

      !                  /                   ;                     /                  ;                     /
      d(0,3,3,S) =        0._R8P              ; d(1,3,3,S) =         0._R8P              ; d(2,3,3,S) =         0._R8P
      !                 i*i                  ;                (i-1)*i                 ;                (i-2)*i
      d(3,3,3,S) =  4877743._R8P /  30240._R8P; d(4,3,3,S) =  -5106971._R8P /  30240._R8P; d(5,3,3,S) =    255397._R8P /   7560._R8P

      !                  /                   ;                     /                  ;                     /
      d(0,4,3,S) =        0._R8P              ; d(1,4,3,S) =         0._R8P              ; d(2,4,3,S) =         0._R8P
      !                  /                   ;                (i-1)*(i-1)             ;                (i-2)*(i-1)
      d(3,4,3,S) =        0._R8P              ; d(4,4,3,S) =   1884439._R8P /  40320._R8P; d(5,4,3,S) =     -1240._R8P /     63._R8P

      !                  /                   ;                     /                  ;                     /
      d(0,5,3,S) =        0._R8P              ; d(1,5,3,S) =         0._R8P              ; d(2,5,3,S) =         0._R8P
      !                  /                   ;                     /                  ;                (i-2)*(i-2)
      d(3,5,3,S) =        0._R8P              ; d(4,5,3,S) =         0._R8P              ; d(5,5,3,S) =     90593._R8P /  40320._R8P
      ! stencil 4
      !             (i+4)*(i+4)              ;                (i+3)*(i+4)             ;                 (i+2)*(i+4)
      d(0,0,4,S) =    90593._R8P /  40320._R8P; d(1,0,4,S) =   -188483._R8P /   7560._R8P; d(2,0,4,S) =    139471._R8P /   2520._R8P
      !             (i+1)*(i+4)              ;                    i*(i+4)             ;                 (i-1)*(i+4)
      d(3,0,4,S) =   -68601._R8P /   1120._R8P; d(4,0,4,S) =   2033509._R8P /  60480._R8P; d(5,0,4,S) =    -73379._R8P /  10080._R8P

      !                  /                   ;                (i+3)*(i+3)             ;                 (i+2)*(i+3)
      d(0,1,4,S) =        0._R8P              ; d(1,1,4,S) =   8449957._R8P / 120960._R8P; d(2,1,4,S) =  -9478331._R8P /  30240._R8P
      !             (i+1)*(i+3)              ;                    i*(i+3)             ;                 (i-1)*(i+3)
      d(3,1,4,S) =  5300629._R8P /  15120._R8P; d(4,1,4,S) =  -5877617._R8P /  30240._R8P; d(5,1,4,S) =   2567287._R8P /  60480._R8P

      !                  /                   ;                     /                  ;                 (i+2)*(i+2)
      d(0,2,4,S) =        0._R8P              ; d(1,2,4,S) =         0._R8P              ; d(2,2,4,S) =   1197047._R8P /   3360._R8P
      !             (i+1)*(i+2)              ;                    i*(i+2)             ;                 (i-1)*(i+2)
      d(3,2,4,S) =  -169859._R8P /    210._R8P; d(4,2,4,S) =   6881719._R8P /  15120._R8P; d(5,2,4,S) =  -1015303._R8P /  10080._R8P

      !                  /                   ;                     /                  ;                      /
      d(0,3,4,S) =        0._R8P              ; d(1,3,4,S) =         0._R8P              ; d(2,3,4,S) =         0._R8P
      !             (i+1)*(i+1)              ;                    i*(i+1)             ;                 (i-1)*(i+1)
      d(3,3,4,S) =  4721851._R8P /  10080._R8P; d(4,3,4,S) = -16306061._R8P /  30240._R8P; d(5,3,4,S) =     61427._R8P /    504._R8P

      !                  /                   ;                     /                  ;                      /
      d(0,4,4,S) =        0._R8P              ; d(1,4,4,S) =         0._R8P              ; d(2,4,4,S) =         0._R8P
      !                  /                   ;                    i*i                 ;                 (i-1)*i
      d(3,4,4,S) =        0._R8P              ; d(4,4,4,S) =  19365967._R8P / 120960._R8P; d(5,4,4,S) =  -1139749._R8P /  15120._R8P

      !                  /                   ;                     /                  ;                      /
      d(0,5,4,S) =        0._R8P              ; d(1,5,4,S) =         0._R8P              ; d(2,5,4,S) =         0._R8P
      !                  /                   ;                     /                  ;                 (i-1)*(i-1)
      d(3,5,4,S) =        0._R8P              ; d(4,5,4,S) =         0._R8P              ; d(5,5,4,S) =    384187._R8P /  40320._R8P
      ! stencil 5
      !             (i+5)*(i+5)              ;                (i+4)*(i+5)             ;                 (i+3)*(i+5)
      d(0,0,5,S) =   384187._R8P /  40320._R8P; d(1,0,5,S) =   -539591._R8P /   5040._R8P; d(2,0,5,S) =   1840141._R8P /   7560._R8P
      !             (i+2)*(i+5)              ;                (i+1)*(i+5)             ;                     i*(i+5)
      d(3,0,5,S) =  -314063._R8P /   1120._R8P; d(4,0,5,S) =    661145._R8P /   4032._R8P; d(5,0,5,S) =   -235637._R8P /   6048._R8P

      !                  /                   ;                (i+4)*(i+3)             ;                 (i+3)*(i+3)
      d(0,1,5,S) =        0._R8P              ; d(1,1,5,S) =  12160229._R8P /  40320._R8P; d(2,1,5,S) = -41615261._R8P /  30240._R8P
      !             (i+2)*(i+3)              ;                (i+1)*(i+3)             ;                     i*(i+3)
      d(3,1,5,S) =  2674951._R8P /   1680._R8P; d(4,1,5,S) =  -1048211._R8P /   1120._R8P; d(5,1,5,S) =   2706017._R8P /  12096._R8P

      !                  /                   ;                     /                  ;                 (i+3)*(i+2)
      d(0,2,5,S) =        0._R8P              ; d(1,2,5,S) =         0._R8P              ; d(2,2,5,S) =  47689393._R8P /  30240._R8P
      !             (i+2)*(i+2)              ;                (i+1)*(i+2)             ;                     i*(i+2)
      d(3,2,5,S) = -6937561._R8P /   1890._R8P; d(4,2,5,S) =  32862709._R8P /  15120._R8P; d(5,2,5,S) = -15848531._R8P /  30240._R8P

      !                  /                   ;                     /                  ;                      /
      d(0,3,5,S) =        0._R8P              ; d(1,3,5,S) =         0._R8P              ; d(2,3,5,S) =         0._R8P
      !             (i+2)*(i+1)              ;                (i+1)*(i+1)             ;                     i*(i+1)
      d(3,3,5,S) = 21703781._R8P /  10080._R8P; d(4,3,5,S) = -25980937._R8P /  10080._R8P; d(5,3,5,S) =   4762921._R8P /   7560._R8P

      !                  /                   ;                     /                  ;                      /
      d(0,4,5,S) =        0._R8P              ; d(1,4,5,S) =         0._R8P              ; d(2,4,5,S) =         0._R8P
      !                  /                   ;                (i+1)*i                 ;                     i*i
      d(3,4,5,S) =        0._R8P              ; d(4,4,5,S) =  31617079._R8P /  40320._R8P; d(5,4,5,S) =  -2966279._R8P /   7560._R8P

      !                  /                   ;                     /                  ;                      /
      d(0,5,5,S) =        0._R8P              ; d(1,5,5,S) =         0._R8P              ; d(2,5,5,S) =         0._R8P
      !                  /                   ;                     /                  ;                     i*(i-1)
      d(3,5,5,S) =        0._R8P              ; d(4,5,5,S) =         0._R8P              ; d(5,5,5,S) =   6150211._R8P / 120960._R8P
   endassociate
   endsubroutine initialize_S6

   pure subroutine compute_polynomials(self, S, V, VP)
   !< Compute WENO polynomials.
   class(weno_object), intent(in)  :: self             !< WENO object.
   integer(I4P),       intent(in)  :: S                !< Number of stencils used.
   real(R8P),          intent(in)  :: V (1:2,1-S:-1+S) !< Variable to be reconstructed.
   real(R8P),          intent(out) :: VP(1:2,0:S-1   ) !< Polynomial reconstructions.

   call weno_compute_polynomials(S=S, weno_p=self%p, V=V, VP=VP)
   endsubroutine compute_polynomials

   pure subroutine compute_weights(self, S, V, w)
   !< Compute WENO weights.
   class(weno_object), intent(in)  :: self                !< WENO object.
   integer(I4P),       intent(in)  :: S                   !< Number of stencils used.
   real(R8P),          intent(in)  :: V    (1:2,1-S:-1+S) !< Variable to be reconstructed.
   real(R8P),          intent(out) :: w    (1:2,0:S-1)    !< Weights of the stencils.

   call weno_compute_weights(S=S, weno_a=self%a, weno_d=self%d, weno_zeps=self%zeps, V=V, w=w)
   endsubroutine compute_weights

   ! non TBP
   pure subroutine weno_compute_convolution(S, VP, w, VR)
   !< Compute WENO convulution, non TBP.
   integer(I4P), intent(in)  :: S             !< Number of stencils used.
   real(R8P),    intent(in)  :: VP(1:2,0:S-1) !< Polynomial reconstructions.
   real(R8P),    intent(in)  :: w (1:2,0:S-1) !< Weights of the stencils.
   real(R8P),    intent(out) :: VR(1:2      ) !< Left and right (1,2) interface value of reconstructed V.
   integer(I4P)              :: k,f           !< Counter.

   VR = 0._R8P
   do k=0, S-1
      do f=1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
         VR(f) = VR(f) + w(f,k)*VP(f,k)
      enddo
   enddo
   endsubroutine weno_compute_convolution

   pure subroutine weno_compute_polynomials(S, weno_p, V, VP)
   !< Compute WENO polynomials, non TBP.
   integer(I4P), intent(in)  :: S                   !< Number of stencils used.
   real(R8P),    intent(in)  :: weno_p(1:,0:,0:,1:) !< Polinomials coefficients.
   real(R8P),    intent(in)  :: V (1:2,1-S:-1+S)    !< Variable to be reconstructed.
   real(R8P),    intent(out) :: VP(1:2,0:S-1   )    !< Polynomial reconstructions.
   integer(I4P)              :: s1,s2,f             !< Counter.

   ! computing the polynomials
   VP = 0._R8P
   do s1=0, S-1 ! stencil counter
      do s2=0, S-1 ! cell counter counter
         do f=1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
            VP(f,s1) = VP(f,s1) + weno_p(f,s2,s1,S)*V(f,-s2+s1)
         enddo
      enddo
   enddo
   endsubroutine weno_compute_polynomials

   pure subroutine weno_compute_weights(S, weno_a, weno_d, weno_zeps, V, w)
   !< Compute WENO weights, non TBP.
   integer(I4P), intent(in)  :: S                   !< Number of stencils used.
   real(R8P),    intent(in)  :: weno_a(1:,0:,1:)    !< Optimal weights.
   real(R8P),    intent(in)  :: weno_d(0:,0:,0:,1:) !< Smoothness indicators coefficients.
   real(R8P),    intent(in)  :: weno_zeps           !< Parameter for avoiding division by zero in computing IS.
   real(R8P),    intent(in)  :: V    (1:2,1-S:-1+S) !< Variable to be reconstructed.
   real(R8P),    intent(out) :: w    (1:2,0:S-1)    !< Weights of the stencils.
   real(R8P)                 :: IS   (1:2,0:S-1)    !< Smoothness indicators of the stencils.
   real(R8P)                 :: a    (1:2,0:S-1)    !< Alpha coifficients for the weights.
   real(R8P)                 :: a_tot(1:2)          !< Summ of the alpha coefficients.
   integer(I4P)              :: s1,s2,s3,f          !< Counter.

   ! computing smoothness indicators
   do s1=0,S-1 ! stencil counter
      do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
         IS(f,s1) = 0._R8P
         do s2=0,S-1
            do s3=0,S-1
               IS(f,s1) = IS(f,s1) + weno_d(s3,s2,s1,S)*V(f,s1-s3)*V(f,s1-s2)
            enddo
         enddo
      enddo
   enddo
   ! computing alfa coefficients
   a_tot = 0._R8P
   do s1=0,S-1
      do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
         a(f,s1) = weno_a(f,s1,S)*(1._R8P/(weno_zeps+IS(f,s1))**S) ; a_tot(f) = a_tot(f) + a(f,s1)
      enddo
   enddo
   ! computing weights
   do s1=0,S-1
      do f=1,2
         w(f,s1) = a(f,s1)/a_tot(f)
      enddo
   enddo
   endsubroutine weno_compute_weights

   pure subroutine weno_reconstruct_optimal(S, weno_a, weno_p, V, VR)
   !< Reconstruct by WENO method with optimal weigths (without smoothness indicators computations) of 2S-1 order.
   integer(I4P), intent(in)  :: S                   !< Number of stencils used.
   real(R8P),    intent(in)  :: weno_p(1:,0:,0:,1:) !< Polinomials coefficients.
   real(R8P),    intent(in)  :: weno_a(1:,0:,1:)    !< Optimal weights.
   real(R8P),    intent(in)  :: V (1:2,1-S:-1+S)    !< Variable to be reconstructed.
   real(R8P),    intent(out) :: VR(1:2         )    !< Left and right (1,2) interface value of reconstructed V.
   real(R8P)                 :: VP(1:2,0:S-1   )    !< Polynomial reconstructions.

   call weno_compute_polynomials(S=S, weno_p=weno_p, V=V(1:2,1-S:-1+S), VP=VP(1:2,0:S-1))
   call weno_compute_convolution(S=S, VP=VP(1:2,0:S-1), w=weno_a, VR=VR(1:2))
   endsubroutine weno_reconstruct_optimal

   pure subroutine weno_reconstruct_upwind(S, weno_a, weno_p, weno_d, weno_zeps, V, VR)
   !< Reconstruct by WENO upwind method of 2S-1 order, non TBP.
   integer(I4P), intent(in)  :: S                   !< Number of stencils used.
   real(R8P),    intent(in)  :: weno_a(1:,0:,1:)    !< Optimal weights.
   real(R8P),    intent(in)  :: weno_p(1:,0:,0:,1:) !< Polinomials coefficients.
   real(R8P),    intent(in)  :: weno_d(0:,0:,0:,1:) !< Smoothness indicators coefficients.
   real(R8P),    intent(in)  :: weno_zeps           !< Parameter for avoiding division by zero in computing IS.
   real(R8P),    intent(in)  :: V (1:2,1-S:-1+S)    !< Variables to be reconstructed.
   real(R8P),    intent(out) :: VR(1:2         )    !< Left and right (1,2) interface value of reconstructed V.
   real(R8P)                 :: VP(1:2,0:S-1   )    !< Polynomial reconstructions.
   real(R8P)                 :: w (1:2,0:S-1   )    !< Weights of the stencils.

   call weno_compute_polynomials(S=S, weno_p=weno_p, V=V(1:2,1-S:-1+S), VP=VP(1:2,0:S-1))
   call weno_compute_weights(S=S, weno_a=weno_a, weno_d=weno_d, weno_zeps=weno_zeps, V=V(1:2,1-S:-1+S), w=w(1:2,0:S-1))
   call weno_compute_convolution(S=S, VP=VP(1:2,0:S-1), w=w(1:2,0:S-1), VR=VR(1:2))
   endsubroutine weno_reconstruct_upwind
endmodule adam_weno_object
