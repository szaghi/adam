!< ADAM, WENO class definition, CPU backend.
module adam_weno_cpu_object
!< ADAM, WENO class definition, CPU backend.

use adam_mpih_object, only : mpih_object
use finer
use penf

implicit none
private
public :: weno_cpu_object

character(len=7), parameter :: INI_SECTION_NAME='schemes' !< INI (config) file section name containing RK configs.

type :: weno_cpu_object
   !< WENO class definition, CPU backend.
   type(mpih_object)         :: mpih            !< MPI handler.
   character(:), allocatable :: scheme          !< RK scheme, ['rk-11', 'rk-33', 'rk-43', 'rk-54']
   integer(I4P)              :: weno_s=1_I4P    !< Stencil number.
   real(R8P), allocatable    :: weno_c(:,:)     !< Central difference coefficients    [1:2,1:2*S].
   real(R8P), allocatable    :: weno_a(:,:)     !< Optimal weights                    [1:2,0:S-1].
   real(R8P), allocatable    :: weno_p(:,:,:)   !< Polinomials coefficients           [1:2,0:S-1,0:S-1].
   real(R8P), allocatable    :: weno_d(:,:,:)   !< Smoothness indicators coefficients [0:S-1,0:S-1,0:S-1].
   real(R8P)                 :: weno_eps=0._R8P !< Parameter for avoiding divided by zero when computing smoothness indicators.
   integer(I4P)              :: weno_exp=0_I4P  !< Exponent for growing the diffusive part of weights.
   contains
      ! public methods
      procedure, pass(self) :: initialize     !< Initialize RK.
      procedure, pass(self) :: load_from_file !< Load config from file.
      procedure, pass(self) :: reconstructed  !< Return WENO reconstruction of 2S-1 order.
      ! private methods
      procedure, pass(self), private :: convolution !< Return WENO convolution of the weighted polynomial recontructions.
      procedure, pass(self), private :: polynomials !< Return WENO polynomials
      procedure, pass(self), private :: weights     !< Return WENO weights of the polynomial reconstructions.
endtype weno_cpu_object

contains
   ! public methods
   subroutine initialize(self, file_parameters)
   !< Initialize the equation.
   class(weno_cpu_object), intent(inout) :: self            !< WENO.
   type(file_ini),         intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'weno_cpu_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   if (allocated(self%weno_c)) deallocate(self%weno_c) ; allocate(self%weno_c(1:2,1:2*self%weno_s))
   if (allocated(self%weno_a)) deallocate(self%weno_a) ; allocate(self%weno_a(1:2,0:self%weno_s-1))
   if (allocated(self%weno_p)) deallocate(self%weno_p) ; allocate(self%weno_p(1:2,0:self%weno_s-1,0:self%weno_s-1))
   if (allocated(self%weno_d)) deallocate(self%weno_d) ; allocate(self%weno_d(0:self%weno_s-1,0:self%weno_s-1,0:self%weno_s-1))
   associate(c=>self%weno_c, a=>self%weno_a, p=>self%weno_p, d=>self%weno_d)
   if (self%weno_s/=1) then
      self%weno_exp = self%weno_s
      ! if (self%weno_s>4) self%weno_exp = self%weno_s - 1
      self%weno_eps = 0.00000000001_R8P
      select case(self%weno_s)
      case(2) ! 3rd order WENO reconstruction
         ! central difference coefficients
         ! 1 => left interface (i-1/2)
         c(1,1) = -1._R8P/12._R8P ! cell -2
         c(1,2) =  7._R8P/12._R8P ! cell -1
         c(1,3) =  7._R8P/12._R8P ! cell  0
         c(1,4) = -1._R8P/12._R8P ! cell  1
         ! 2 => right interface (i+1/2)
         c(2,1) = -1._R8P/12._R8P ! cell -1
         c(2,2) =  7._R8P/12._R8P ! cell  0
         c(2,3) =  7._R8P/12._R8P ! cell  1
         c(2,4) = -1._R8P/12._R8P ! cell  2

         ! optimal weights
         ! 1 => left interface (i-1/2)
         a(1,0) = 2._R8P/3._R8P ! stencil 0
         a(1,1) = 1._R8P/3._R8P ! stencil 1
         ! 2 => right interface (i+1/2)
         a(2,0) = 1._R8P/3._R8P ! stencil 0
         a(2,1) = 2._R8P/3._R8P ! stencil 1

         ! polinomials coefficients
         ! 1 => left interface (i-1/2)
         !  cell  0          ;    cell  1
         p(1,0,0) =  0.5_R8P ; p(1,1,0) =  0.5_R8P ! stencil 0
         p(1,0,1) = -0.5_R8P ; p(1,1,1) =  1.5_R8P ! stencil 1
         ! 2 => right interface (i+1/2)
         !  cell  0          ;    cell  1
         p(2,0,0) =  1.5_R8P ; p(2,1,0) = -0.5_R8P ! stencil 0
         p(2,0,1) =  0.5_R8P ; p(2,1,1) =  0.5_R8P ! stencil 1

         ! smoothness indicators coefficients
         ! stencil 0
         !      i*i        ;       (i-1)*i
         d(0,0,0) = 1._R8P ; d(1,0,0) =-2._R8P
         !      /          ;       (i-1)*(i-1)
         d(0,1,0) = 0._R8P ; d(1,1,0) = 1._R8P
         ! stencil 1
         !     (i+1)*(i+1) ;       (i+1)*i
         d(0,0,1) = 1._R8P ; d(1,0,1) =-2._R8P
         !      /          ;        i*i
         d(0,1,1) = 0._R8P ; d(1,1,1) = 1._R8P
      case(3) ! 5th order WENO reconstruction
         ! central difference coefficients
         ! 1 => left interface (i-1/2)
         c(1,1) =  1._R8P/60._R8P ! cell -3
         c(1,2) = -7.5_R8P        ! cell -2
         c(1,3) = 37._R8P/60._R8P ! cell -1
         c(1,4) = 37._R8P/60._R8P ! cell  0
         c(1,5) = -7.5_R8P        ! cell  1
         c(1,6) =  1._R8P/60._R8P ! cell  2
         ! 2 => right interface (i+1/2)
         c(1,1) =  1._R8P/60._R8P ! cell -2
         c(1,2) = -7.5_R8P        ! cell -1
         c(1,3) = 37._R8P/60._R8P ! cell  0
         c(1,4) = 37._R8P/60._R8P ! cell  1
         c(1,5) = -7.5_R8P        ! cell  2
         c(1,6) =  1._R8P/60._R8P ! cell  3

         ! optimal weights
         ! 1 => left interface (i-1/2)
         a(1,0) = 0.3_R8P ! stencil 0
         a(1,1) = 0.6_R8P ! stencil 1
         a(1,2) = 0.1_R8P ! stencil 2
         ! 2 => right interface (i+1/2)
         a(2,0) = 0.1_R8P ! stencil 0
         a(2,1) = 0.6_R8P ! stencil 1
         a(2,2) = 0.3_R8P ! stencil 2

         ! polinomials coefficients
         ! 1 => left interface (i-1/2)
         !  cell  0                ;    cell  1                ;    cell  2
         p(1,0,0) =  1._R8P/3._R8P ; p(1,1,0) =  5._R8P/6._R8P ; p(1,2,0) = -1._R8P/6._R8P ! stencil 0
         p(1,0,1) = -1._R8P/6._R8P ; p(1,1,1) =  5._R8P/6._R8P ; p(1,2,1) =  1._R8P/3._R8P ! stencil 1
         p(1,0,2) =  1._R8P/3._R8P ; p(1,1,2) = -7._R8P/6._R8P ; p(1,2,2) = 11._R8P/6._R8P ! stencil 2
         ! 2 => right interface (i+1/2)
         !  cell  0                ;    cell  1                ;    cell  2
         p(2,0,0) = 11._R8P/6._R8P ; p(2,1,0) = -7._R8P/6._R8P ; p(2,2,0) =  1._R8P/3._R8P ! stencil 0
         p(2,0,1) =  1._R8P/3._R8P ; p(2,1,1) =  5._R8P/6._R8P ; p(2,2,1) = -1._R8P/6._R8P ! stencil 1
         p(2,0,2) = -1._R8P/6._R8P ; p(2,1,2) =  5._R8P/6._R8P ; p(2,2,2) =  1._R8P/3._R8P ! stencil 2

         ! smoothness indicators coefficients
         ! stencil 0
         !      i*i                 ;       (i-1)*i              ;       (i-2)*i
         d(0,0,0) =  10._R8P/3._R8P ; d(1,0,0) = -31._R8P/3._R8P ; d(2,0,0) =  11._R8P/3._R8P
         !      /                   ;       (i-1)*(i-1)          ;       (i-2)*(i-1)
         d(0,1,0) =   0._R8P        ; d(1,1,0) =  25._R8P/3._R8P ; d(2,1,0) = -19._R8P/3._R8P
         !      /                   ;        /                   ;       (i-2)*(i-2)
         d(0,2,0) =   0._R8P        ; d(1,2,0) =   0._R8P        ; d(2,2,0) =   4._R8P/3._R8P
         ! stencil 1
         !     (i+1)*(i+1)          ;        i*(i+1)             ;       (i-1)*(i+1)
         d(0,0,1) =   4._R8P/3._R8P ; d(1,0,1) = -13._R8P/3._R8P ; d(2,0,1) =   5._R8P/3._R8P
         !      /                   ;        i*i                 ;       (i-1)*i
         d(0,1,1) =   0._R8P        ; d(1,1,1) =  13._R8P/3._R8P ; d(2,1,1) = -13._R8P/3._R8P
         !      /                   ;        /                   ;       (i-1)*(i-1)
         d(0,2,1) =   0._R8P        ; d(1,2,1) =   0._R8P        ; d(2,2,1) =   4._R8P/3._R8P
         ! stencil 2
         !     (i+2)*(i+2)          ;       (i+1)*(i+2)          ;        i*(i+2)
         d(0,0,2) =   4._R8P/3._R8P ; d(1,0,2) = -19._R8P/3._R8P ; d(2,0,2) =  11._R8P/3._R8P
         !      /                   ;       (i+1)*(i+1)          ;        i*(i+1)
         d(0,1,2) =   0._R8P        ; d(1,1,2) =  25._R8P/3._R8P ; d(2,1,2) = -31._R8P/3._R8P
         !      /                   ;        /                   ;        i*i
         d(0,2,2) =   0._R8P        ; d(1,2,2) =   0._R8P        ; d(2,2,2) =  10._R8P/3._R8P
      case(4) ! 7th order WENO reconstruction
         ! optimal weights
         ! 1 => left interface (i-1/2)
         a(1,0) =  4._R8P/35._R8P ! stencil 0
         a(1,1) = 18._R8P/35._R8P ! stencil 1
         a(1,2) = 12._R8P/35._R8P ! stencil 2
         a(1,3) =  1._R8P/35._R8P ! stencil 3
         ! 2 => right interface (i+1/2)
         a(2,0) =  1._R8P/35._R8P ! stencil 0
         a(2,1) = 12._R8P/35._R8P ! stencil 1
         a(2,2) = 18._R8P/35._R8P ! stencil 2
         a(2,3) =  4._R8P/35._R8P ! stencil 3

         ! polinomials coefficients
         ! 1 => left interface (i-1/2)
         !  cell  0               ;    cell  1               ;    cell  2               ;    cell  3
         p(1,0,0)=  1._R8P/4._R8P ; p(1,1,0)= 13._R8P/12._R8P; p(1,2,0)= -5._R8P/12._R8P; p(1,3,0)=  1._R8P/12._R8P ! stencil 0
         p(1,0,1)= -1._R8P/12._R8P; p(1,1,1)=  7._R8P/12._R8P; p(1,2,1)=  7._R8P/12._R8P; p(1,3,1)= -1._R8P/12._R8P ! stencil 1
         p(1,0,2)=  1._R8P/12._R8P; p(1,1,2)= -5._R8P/12._R8P; p(1,2,2)= 13._R8P/12._R8P; p(1,3,2)=  1._R8P/4._R8P  ! stencil 2
         p(1,0,3)= -1._R8P/4._R8P ; p(1,1,3)= 13._R8P/12._R8P; p(1,2,3)=-23._R8P/12._R8P; p(1,3,3)= 25._R8P/12._R8P ! stencil 3
         ! 2 => right interface (i+1/2)
         !  cell  0               ;    cell  1               ;   cell  2                ;    cell  3
         p(2,0,0)= 25._R8P/12._R8P; p(2,1,0)=-23._R8P/12._R8P; p(2,2,0)= 13._R8P/12._R8P; p(2,3,0)= -1._R8P/4._R8P  ! stencil 0
         p(2,0,1)=  1._R8P/4._R8P ; p(2,1,1)= 13._R8P/12._R8P; p(2,2,1)= -5._R8P/12._R8P; p(2,3,1)=  1._R8P/12._R8P ! stencil 1
         p(2,0,2)= -1._R8P/12._R8P; p(2,1,2)=  7._R8P/12._R8P; p(2,2,2)=  7._R8P/12._R8P; p(2,3,2)= -1._R8P/12._R8P ! stencil 2
         p(2,0,3)=  1._R8P/12._R8P; p(2,1,3)= -5._R8P/12._R8P; p(2,2,3)= 13._R8P/12._R8P; p(2,3,3)=  1._R8P/4._R8P  ! stencil 3

         ! smoothness indicators coefficients
         ! stencil 0
         !              i*i              ;            (i-1)*i              ;            (i-2)*i
         d(0,0,0) =  2107._R8P / 240._R8P; d(1,0,0) = -1567._R8P /  40._R8P; d(2,0,0) =  3521._R8P / 120._R8P
         !          (i-3)*i
         d(3,0,0) =  -309._R8P /  40._R8P
         !               /               ;            (i-1)*(i-1)          ;            (i-2)*(i-1)
         d(0,1,0) =     0._R8P           ; d(1,1,0) = 11003._R8P / 240._R8P; d(2,1,0) = -8623._R8P / 120._R8P
         !          (i-3)*(i-1)
         d(3,1,0) =  2321._R8P / 120._R8P
         !               /               ;                 /               ;            (i-2)*(i-2)
         d(0,2,0) =     0._R8P           ; d(1,2,0) =     0._R8P           ; d(2,2,0) =  7043._R8P / 240._R8P
         !          (i-3)*(i-2)
         d(3,2,0) =  -647._R8P /  40._R8P
         !               /               ;                 /               ;                 /
         d(0,3,0) =     0._R8P           ; d(1,3,0) =     0._R8P           ; d(2,3,0) =     0._R8P
         !          (i-3)*(i-3)
         d(3,3,0) =   547._R8P / 240._R8P
         ! stencil 1
         !          (i+1)*(i+1)          ;                i*(i+1)          ;            (i-1)*(i+1)
         d(0,0,1) =   547._R8P / 240._R8P; d(1,0,1) = -1261._R8P / 120._R8P; d(2,0,1) =   961._R8P / 120._R8P
         !          (i-2)*(i+1)
         d(3,0,1) =  -247._R8P / 120._R8P
         !               /               ;                i*i              ;            (i-1)*i
         d(0,1,1) =     0._R8P           ; d(1,1,1) =  3443._R8P / 240._R8P; d(2,1,1) = -2983._R8P / 120._R8P
         !          (i-2)*i
         d(3,1,1) =   267._R8P /  40._R8P
         !               /               ;                 /               ;            (i-1)*(i-1)
         d(0,2,1) =     0._R8P           ; d(1,2,1) =     0._R8P           ; d(2,2,1) =  2843._R8P / 240._R8P
         !          (i-2)*(i-1)
         d(3,2,1) =  -821._R8P / 120._R8P
         !               /               ;                 /               ;                 /
         d(0,3,1) =     0._R8P           ; d(1,3,1) =     0._R8P           ; d(2,3,1) =     0._R8P
         !          (i-2)*(i-2)
         d(3,3,1) =    89._R8P /  80._R8P
         ! stencil 2
         !          (i+2)*(i+2)          ;            (i+1)*(i+2)          ;                i*(i+2)
         d(0,0,2) =    89._R8P /  80._R8P; d(1,0,2) =  -821._R8P / 120._R8P; d(2,0,2) =   267._R8P /  40._R8P
         !          (i-1)*(i+2)
         d(3,0,2) =  -247._R8P / 120._R8P
         !               /               ;            (i+1)*(i+1)          ;                i*(i+1)
         d(0,1,2) =     0._R8P           ; d(1,1,2) =  2843._R8P / 240._R8P; d(2,1,2) = -2983._R8P / 120._R8P
         !          (i-1)*(i+1)
         d(3,1,2) =   961._R8P / 120._R8P
         !               /               ;                 /               ;                i*i
         d(0,2,2) =     0._R8P           ; d(1,2,2) =     0._R8P           ; d(2,2,2) =  3443._R8P / 240._R8P
         !          (i-1)*i
         d(3,2,2) = -1261._R8P / 120._R8P
         !               /               ;                 /               ;                 /
         d(0,3,2) =     0._R8P           ; d(1,3,2) =     0._R8P           ; d(2,3,2) =     0._R8P
         !          (i-1)*(i-1)
         d(3,3,2) =   547._R8P / 240._R8P
         ! stencil 3
         !          (i+3)*(i+3)          ;            (i+2)*(i+3)          ;            (i+1)*(i+3)
         d(0,0,3) =   547._R8P / 240._R8P; d(1,0,3) =  -647._R8P /  40._R8P; d(2,0,3) =  2321._R8P / 120._R8P
         !              i*(i+3)
         d(3,0,3) =  -309._R8P /  40._R8P
         !               /               ;            (i+2)*(i+2)          ;      (i+1)*(i+2)
         d(0,1,3) =     0._R8P           ; d(1,1,3) =  7043._R8P / 240._R8P; d(2,1,3) = -8623._R8P / 120._R8P
         !              i*(i+2)
         d(3,1,3) =  3521._R8P / 120._R8P
         !               /               ;                 /               ;      (i+1)*(i+1)
         d(0,2,3) =     0._R8P           ; d(1,2,3) =     0._R8P           ; d(2,2,3) = 11003._R8P / 240._R8P
         !              i*(i+1)
         d(3,2,3) = -1567._R8P /  40._R8P
         !               /               ;                 /               ;           /
         d(0,3,3) =     0._R8P           ; d(1,3,3) =     0._R8P           ; d(2,3,3) =     0._R8P
         !              i*i
         d(3,3,3) =  2107._R8P / 240._R8P
      case(5) ! 9th order WENO reconstruction
         ! optimal weights
         ! 1 => left interface (i-1/2)
         a(1,0) =  5._R8P/126._R8P ! stencil 0
         a(1,1) = 20._R8P/63._R8P  ! stencil 1
         a(1,2) = 10._R8P/21._R8P  ! stencil 2
         a(1,3) = 10._R8P/63._R8P  ! stencil 3
         a(1,4) =  1._R8P/126._R8P ! stencil 4
         ! 2 => right interface (i+1/2)
         a(2,0) =  1._R8P/126._R8P ! stencil 0
         a(2,1) = 10._R8P/63._R8P  ! stencil 1
         a(2,2) = 10._R8P/21._R8P  ! stencil 2
         a(2,3) = 20._R8P/63._R8P  ! stencil 3
         a(2,4) =  5._R8P/126._R8P ! stencil 4

         ! polinomials coefficients
         ! 1 => left interface (i-1/2)
         !  cell  0                ;    cell  1                ;    cell  2                ;    cell  3
         p(1,0,0)=   1._R8P/5._R8P ; p(1,1,0)=  77._R8P/60._R8P; p(1,2,0)= -43._R8P/60._R8P; p(1,3,0)=  17._R8P/60._R8P  ! stencil 0
         p(1,0,1)=  -1._R8P/20._R8P; p(1,1,1)=   9._R8P/20._R8P; p(1,2,1)=  47._R8P/60._R8P; p(1,3,1)= -13._R8P/60._R8P  ! stencil 1
         p(1,0,2)=   1._R8P/30._R8P; p(1,1,2)= -13._R8P/60._R8P; p(1,2,2)=  47._R8P/60._R8P; p(1,3,2)=   9._R8P/20._R8P  ! stencil 2
         p(1,0,3)=  -1._R8P/20._R8P; p(1,1,3)=  17._R8P/60._R8P; p(1,2,3)= -43._R8P/60._R8P; p(1,3,3)=  77._R8P/60._R8P  ! stencil 3
         p(1,0,4)=   1._R8P/5._R8P ; p(1,1,4)= -21._R8P/20._R8P; p(1,2,4)= 137._R8P/60._R8P; p(1,3,4)=-163._R8P/60._R8P  ! stencil 4
         !  cell  4
         p(1,4,0)=  -1._R8P/20._R8P  ! stencil 0
         p(1,4,1)=   1._R8P/30._R8P  ! stencil 1
         p(1,4,2)=  -1._R8P/20._R8P  ! stencil 2
         p(1,4,3)=   1._R8P/5._R8P   ! stencil 3
         p(1,4,4)= 137._R8P/60._R8P  ! stencil 4
         ! 2 => right interface (i+1/2)
         !  cell  0               ;    cell  1               ;   cell  2                ;    cell  3
         p(2,0,0)= 137._R8P/60._R8P; p(2,1,0)=-163._R8P/60._R8P; p(2,2,0)= 137._R8P/60._R8P; p(2,3,0)= -21._R8P/20._R8P  ! stencil 0
         p(2,0,1)=   1._R8P/5._R8P ; p(2,1,1)=  77._R8P/60._R8P; p(2,2,1)= -43._R8P/60._R8P; p(2,3,1)=  17._R8P/60._R8P  ! stencil 1
         p(2,0,2)=  -1._R8P/20._R8P; p(2,1,2)=   9._R8P/20._R8P; p(2,2,2)=  47._R8P/60._R8P; p(2,3,2)= -13._R8P/60._R8P  ! stencil 2
         p(2,0,3)=   1._R8P/30._R8P; p(2,1,3)= -13._R8P/60._R8P; p(2,2,3)=  47._R8P/60._R8P; p(2,3,3)=   9._R8P/20._R8P  ! stencil 3
         p(2,0,4)=  -1._R8P/20._R8P; p(2,1,4)=  17._R8P/60._R8P; p(2,2,4)= -43._R8P/60._R8P; p(2,3,4)=  77._R8P/60._R8P  ! stencil 4
         !  cell  4
         p(2,4,0)=   1._R8P/5._R8P  ! stencil 0
         p(2,4,1)=  -1._R8P/20._R8P ! stencil 1
         p(2,4,2)=   1._R8P/30._R8P ! stencil 2
         p(2,4,3)=  -1._R8P/20._R8P ! stencil 3
         p(2,4,4)=   1._R8P/5._R8P  ! stencil 4

         ! smoothness indicators coefficients
         ! stencil 0
         !              i*i                 ;             (i-1)*i                ;             (i-2)*i
         d(0,0,0) =   53959._R8P / 2520._R8P; d(1,0,0) = -649501._R8P / 5040._R8P; d(2,0,0) =  252941._R8P / 1680._R8P
         !          (i-3)*i                 ;             (i-4)*i
         d(3,0,0) = -411487._R8P / 5040._R8P; d(4,0,0) =   86329._R8P / 5040._R8P
         !               /                  ;             (i-1)*(i-1)            ;             (i-2)*(i-1)
         d(0,1,0) =       0._R8P            ; d(1,1,0) = 1020563._R8P / 5040._R8P; d(2,1,0) =  -68391._R8P /  140._R8P
         !          (i-3)*(i-1)             ;             (i-4)*(i-1)
         d(3,1,0) =  679229._R8P / 2520._R8P; d(4,1,0) = -288007._R8P / 5040._R8P
         !               /                  ;                  /                 ;             (i-2)*(i-2)
         d(0,2,0) =       0._R8P            ; d(1,2,0) =       0._R8P            ; d(2,2,0) =  507131._R8P / 1680._R8P
         !          (i-3)*(i-2)             ;             (i-4)*(i-2)
         d(3,2,0) = -142033._R8P /  420._R8P; d(4,2,0) =  121621._R8P / 1680._R8P
         !               /                  ;                  /                 ;                  /
         d(0,3,0) =       0._R8P            ; d(1,3,0) =       0._R8P            ; d(2,3,0) =       0._R8P
         !          (i-3)*(i-3)             ;             (i-4)*(i-3)
         d(3,3,0) =  482963._R8P / 5040._R8P; d(4,3,0) = -208501._R8P / 5040._R8P
         !               /                  ;                  /                 ;                  /
         d(0,4,0) =       0._R8P            ; d(1,4,0) =       0._R8P            ; d(2,4,0) =       0._R8P
         !               /                  ;             (i-4)*(i-4)
         d(3,4,0) =       0._R8P            ; d(4,4,0) =   11329._R8P / 2520._R8P
         ! stencil 1
         !          (i+1)*(i+1)             ;                 i*(i+1)            ;             (i-1)*(i+1)
         d(0,0,1) =   11329._R8P / 2520._R8P; d(1,0,1) = -140251._R8P / 5040._R8P; d(2,0,1) =   55051._R8P / 1680._R8P
         !          (i-2)*(i+1)             ;             (i-3)*(i+1)
         d(3,0,1) =  -88297._R8P / 5040._R8P; d(4,0,1) =   18079._R8P / 5040._R8P
         !               /                  ;                 i*i                ;             (i-1)*i
         d(0,1,1) =       0._R8P            ; d(1,1,1) =  242723._R8P / 5040._R8P; d(2,1,1) =  -25499._R8P /  210._R8P
         !          (i-2)*i                 ;             (i-3)*i
         d(3,1,1) =  168509._R8P / 2520._R8P; d(4,1,1) =  -70237._R8P / 5040._R8P
         !               /                  ;                  /                 ;             (i-1)*(i-1)
         d(0,2,1) =       0._R8P            ; d(1,2,1) =       0._R8P            ; d(2,2,1) =  135431._R8P / 1680._R8P
         !          (i-2)*(i-1)             ;             (i-3)*(i-1)
         d(3,2,1) =   -3229._R8P /   35._R8P; d(4,2,1) =   33071._R8P / 1680._R8P
         !               /                  ;                  /                 ;                  /
         d(0,3,1) =       0._R8P            ; d(1,3,1) =       0._R8P            ; d(2,3,1) =       0._R8P
         !          (i-2)*(i-2)             ;             (i-3)*(i-2)
         d(3,3,1) =  138563._R8P / 5040._R8P; d(4,3,1) =  -60871._R8P / 5040._R8P
         !               /                  ;                  /                 ;                  /
         d(0,4,1) =       0._R8P            ; d(1,4,1) =       0._R8P            ; d(2,4,1) =       0._R8P
         !               /                  ;             (i-3)*(i-3)
         d(3,4,1) =       0._R8P            ; d(4,4,1) =    1727._R8P / 1260._R8P
         ! stencil 2
         !          (i+2)*(i+2)             ;             (i+1)*(i+2)            ;                 i*(i+2)
         d(0,0,2) =    1727._R8P / 1260._R8P; d(1,0,2) =  -51001._R8P / 5040._R8P; d(2,0,2) =    7547._R8P /  560._R8P
         !          (i-1)*(i+2)             ;             (i-2)*(i+2)
         d(3,0,2) =  -38947._R8P / 5040._R8P; d(4,0,2) =    8209._R8P / 5040._R8P
         !               /                  ;             (i+1)*(i+1)            ;                 i*(i+1)
         d(0,1,2) =       0._R8P            ; d(1,1,2) =  104963._R8P / 5040._R8P; d(2,1,2) =  -24923._R8P /  420._R8P
         !          (i-1)*(i+1)             ;             (i-2)*(i+1)
         d(3,1,2) =   89549._R8P / 2520._R8P; d(4,1,2) =  -38947._R8P / 5040._R8P
         !               /                  ;                  /                 ;                 i*i
         d(0,2,2) =       0._R8P            ; d(1,2,2) =       0._R8P            ; d(2,2,2) =   77051._R8P / 1680._R8P
         !          (i-1)*i                 ;             (i-2)*i
         d(3,2,2) =  -24923._R8P /  420._R8P; d(4,2,2) =    7547._R8P /  560._R8P
         !               /                  ;                  /                 ;                  /
         d(0,3,2) =       0._R8P            ; d(1,3,2) =       0._R8P            ; d(2,3,2) =       0._R8P
         !          (i-1)*(i-1)             ;             (i-2)*(i-1)
         d(3,3,2) =  104963._R8P / 5040._R8P; d(4,3,2) =  -51001._R8P / 5040._R8P
         !               /                  ;                  /                 ;                  /
         d(0,4,2) =       0._R8P            ; d(1,4,2) =       0._R8P            ; d(2,4,2) =       0._R8P
         !               /                  ;             (i-2)*(i-2)
         d(3,4,2) =       0._R8P            ; d(4,4,2) =    1727._R8P / 1260._R8P
         ! stencil 3
         !          (i+3)*(i+3)             ;             (i+2)*(i+3)            ;             (i+1)*(i+3)
         d(0,0,3) =    1727._R8P / 1260._R8P; d(1,0,3) =  -60871._R8P / 5040._R8P; d(2,0,3) =   33071._R8P / 1680._R8P
         !              i*(i+3)             ;             (i-1)*(i+3)
         d(3,0,3) =  -70237._R8P / 5040._R8P; d(4,0,3) =   18079._R8P / 5040._R8P
         !               /                  ;             (i+2)*(i+2)            ;             (i+1)*(i+2)
         d(0,1,3) =       0._R8P            ; d(1,1,3) =  138563._R8P / 5040._R8P; d(2,1,3) =   -3229._R8P /   35._R8P
         !              i*(i+2)             ;             (i-1)*(i+2)
         d(3,1,3) =  168509._R8P / 2520._R8P; d(4,1,3) =  -88297._R8P / 5040._R8P
         !               /                  ;                  /                 ;             (i+1)*(i+1)
         d(0,2,3) =       0._R8P            ; d(1,2,3) =       0._R8P            ; d(2,2,3) =  135431._R8P / 1680._R8P
         !              i*(i+1)             ;             (i-1)*(i+1)
         d(3,2,3) =  -25499._R8P /  210._R8P; d(4,2,3) =   55051._R8P / 1680._R8P
         !               /                  ;                  /                 ;                  /
         d(0,3,3) =       0._R8P            ; d(1,3,3) =       0._R8P            ; d(2,3,3) =       0._R8P
         !              i*i                 ;             (i-1)*i
         d(3,3,3) =  242723._R8P / 5040._R8P; d(4,3,3) = -140251._R8P / 5040._R8P
         !               /                  ;                  /                 ;                  /
         d(0,4,3) =       0._R8P            ; d(1,4,3) =       0._R8P            ; d(2,4,3) =       0._R8P
         !               /                  ;             (i-1)*(i-1)
         d(3,4,3) =       0._R8P            ; d(4,4,3) =   11329._R8P / 2520._R8P
         ! stencil 4
         !          (i+4)*(i+4)             ;             (i+3)*(i+4)            ;             (i+2)*(i+4)
         d(0,0,4) =   11329._R8P / 2520._R8P; d(1,0,4) = -208501._R8P / 5040._R8P; d(2,0,4) =  121621._R8P / 1680._R8P
         !          (i+1)*(i+4)             ;                 i*(i+4)
         d(3,0,4) = -288007._R8P / 5040._R8P; d(4,0,4) =   86329._R8P / 5040._R8P
         !               /                  ;             (i+3)*(i+3)            ;             (i+2)*(i+3)
         d(0,1,4) =       0._R8P            ; d(1,1,4) =  482963._R8P / 5040._R8P; d(2,1,4) = -142033._R8P /  420._R8P
         !          (i+1)*(i+3)             ;                 i*(i+3)
         d(3,1,4) =  679229._R8P / 2520._R8P; d(4,1,4) = -411487._R8P / 5040._R8P
         !               /                  ;                  /                 ;             (i+1)*(i+2)
         d(0,2,4) =       0._R8P            ; d(1,2,4) =       0._R8P            ; d(2,2,4) =  507131._R8P / 1680._R8P
         !          (i+1)*(i+2)             ;                 i*(i+2)
         d(3,2,4) =  -68391._R8P /  140._R8P; d(4,2,4) =  252941._R8P / 1680._R8P
         !               /                  ;                  /                 ;                  /
         d(0,3,4) =       0._R8P            ; d(1,3,4) =       0._R8P            ; d(2,3,4) =       0._R8P
         !          (i+1)*(i+1)             ;                 i*(i+1)
         d(3,3,4) = 1020563._R8P / 5040._R8P; d(4,3,4) = -649501._R8P / 5040._R8P
         !               /                  ;                  /                 ;                  /
         d(0,4,4) =       0._R8P            ; d(1,4,4) =       0._R8P            ; d(2,4,4) =       0._R8P
         !               /                  ;                 i*i
         d(3,4,4) =       0._R8P            ; d(4,4,4) =   53959._R8P / 2520._R8P
      endselect
   endif
   endassociate
   print '(A)', self%mpih%myrankstr//'weno_cpu_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(weno_cpu_object), intent(inout)        :: self            !< WENO.
   type(file_ini),         intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                      :: go_on_fail_     !< Go on if load fails.
   character(99)                                :: char_buff       !< Character buffer.
   integer(I4P)                                 :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='weno', val=char_buff, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(weno)')
   self%scheme = trim(adjustl(char_buff))
   select case(self%scheme)
   case('weno-1')
      self%weno_s = 1
   case('weno-3')
      self%weno_s = 2
   case('weno-5')
      self%weno_s = 3
   case('weno-7')
      self%weno_s = 4
   case('weno-9')
      self%weno_s = 5
   case default
      call self%mpih%error_stop(msg=': unknown weno scheme "'//self%scheme//'"')
   endselect
   endsubroutine load_from_file

   pure function reconstructed(self, s, v) result(vr)
   !< Return WENO reconstruction of 2S-1 order.
   class(weno_cpu_object), intent(in) :: self             !< The equation.
   integer(I4P),           intent(in) :: s                !< Number of stencils used.
   real(R8P),              intent(in) :: v (1:2,1-s:-1+s) !< Variable to be reconstructed.
   real(R8P)                          :: vr(1:2         ) !< Left and right (1,2) interface value of reconstructed V.
   real(R8P)                          :: vp(1:2,0:s-1   ) !< Polynomial reconstructions.
   real(R8P)                          :: w (1:2,0:s-1   ) !< Weights of the stencils.

   vp = self%polynomials(s=s, v=v       ) ! compute the polynomials
   w  = self%weights(    s=s, v=v       ) ! compute the weights associated to the polynomials
   vr = self%convolution(s=s, vp=vp, w=w) ! compute the convultion of reconstructing plynomials
   endfunction reconstructed

   ! private methods
   pure function convolution(self, s, vp, w) result(vr)
   !< Return WENO convolution of the weighted polynomial recontructions.
   class(weno_cpu_object), intent(in) :: self          !< The equation.
   integer(I4P),                     intent(in) :: s             !< Number of stencils used.
   real(R8P),                        intent(in) :: vp(1:2,0:s-1) !< Polynomial reconstructions.
   real(R8P),                        intent(in) :: w (1:2,0:s-1) !< Weights of the stencils.
   real(R8P)                                    :: vr(1:2      ) !< Left and right (1,2) interface value of reconstructed V.
   integer(I4P)                                 :: k, f          !< Counters.

   vr = 0._R8P
   do k=0,S-1
      do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
         vr(f) = vr(f) + w(f,k) * vp(f,k)
      enddo
   enddo
   endfunction convolution

   pure function polynomials(self, s, v) result(VP)
   !< Return WENO polynomials
   class(weno_cpu_object), intent(in) :: self             !< The equation.
   integer(I4P),           intent(in) :: s                !< Number of stencils used.
   real(R8P),              intent(in) :: v (1:2,1-s:-1+s) !< Variable to be reconstructed.
   real(R8P)                          :: vp(1:2,0:s-1   ) !< Polynomial reconstructions.
   integer(I4P)                       :: s1, s2, f        !< Counters.

   vp = 0._R8P
   do s1=0,s-1 ! stencil counter
      do s2=0,s-1 ! cell counter counter
         do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
            vp(f,s1) = vp(f,s1) + self%weno_p(f,s2,s1) * v(f,-s2+s1)
         enddo
      enddo
   enddo
   endfunction polynomials

   pure function weights(self, s, v) result(w)
   !< Return WENO weights of the polynomial reconstructions.
   class(weno_cpu_object), intent(in) :: self                !< The equation.
   integer(I4P),           intent(in) :: s                   !< Number of stencils used.
   real(R8P),              intent(in) :: v    (1:2,1-s:-1+s) !< Variable to be reconstructed.
   real(R8P)                          :: W    (1:2,0:s-1)    !< Weights of the stencils.
   real(R8P)                          :: IS   (1:2,0:s-1)    !< Smoothness indicators of the stencils.
   real(R8P)                          :: a    (1:2,0:s-1)    !< Alpha coefficients for the weights.
   real(R8P)                          :: a_tot(1:2)          !< Summ of the alpha coefficients.
   integer(I4P)                       :: s1, s2, s3, f       !< Counters.

   ! compute smoothness indicators
   do s1=0,S-1 ! stencil counter
      do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
         IS(f,s1) = 0._R8P
         do s2=0,S-1
            do s3=0,S-1
               IS(f,s1) = IS(f,s1) + self%weno_d(s3,s2,s1) * v(f,s1-s3) * v(f,s1-s2)
            enddo
         enddo
      enddo
   enddo
   ! compute alfa coefficients
   a_tot = 0._R8P
   do s1=0,S-1
      do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
         a(f,s1) = self%weno_a(f,s1) / ((self%weno_eps + IS(f,s1))**self%weno_exp)
         a_tot(f) = a_tot(f) + a(f,s1)
      enddo
   enddo
   ! compute the weights
   do s1=0,S-1
      do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
         w(f,s1) = a(f,s1) / a_tot(f)
      enddo
   enddo
   endfunction weights
endmodule adam_weno_cpu_object
