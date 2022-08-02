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
   integer(I4P)              :: weno_odd=1_I4P  !< Constant for distinguishing between odd and even number of stencils (mod(S,2)).
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
   associate(s=>self%weno_s, &
             weno_exp=>self%weno_exp, &
             weno_odd=>self%weno_odd, &
             weno_eps=>self%weno_eps, &
             weno_c=>self%weno_c,     &
             weno_a=>self%weno_a,     &
             weno_p=>self%weno_p,     &
             weno_d=>self%weno_d)
   if (s/=1) then
      ! initialize weno_exp
      weno_exp = s
      if (s>4) weno_exp = s - 1
      ! compute weno_odd
      weno_odd = mod(self%weno_s,2)
      weno_eps = 0.00000000001_R8P
      ! allocate variables
      if (allocated(self%weno_c)) deallocate(self%weno_c) ; allocate(self%weno_c(1:2,1:2*s))
      if (allocated(self%weno_a)) deallocate(self%weno_a) ; allocate(self%weno_a(1:2,0:s-1))
      if (allocated(self%weno_p)) deallocate(self%weno_p) ; allocate(self%weno_p(1:2,0:s-1,0:s-1))
      if (allocated(self%weno_d)) deallocate(self%weno_d) ; allocate(self%weno_d(0:s-1,0:s-1,0:s-1))
      ! inizialize the coefficients
      select case(s)
      case(2) ! 3rd order WENO reconstruction
        ! central difference coefficients
        ! 1 => left interface (i-1/2)
        weno_c(1,1) = -1._R8P/12._R8P ! cell -2
        weno_c(1,2) =  7._R8P/12._R8P ! cell -1
        weno_c(1,3) =  7._R8P/12._R8P ! cell  0
        weno_c(1,4) = -1._R8P/12._R8P ! cell  1
        ! 2 => right interface (i+1/2)
        weno_c(2,1) = -1._R8P/12._R8P ! cell -1
        weno_c(2,2) =  7._R8P/12._R8P ! cell  0
        weno_c(2,3) =  7._R8P/12._R8P ! cell  1
        weno_c(2,4) = -1._R8P/12._R8P ! cell  2

        ! optimal weights
        ! 1 => left interface (i-1/2)
        weno_a(1,0) = 2._R8P/3._R8P ! stencil 0
        weno_a(1,1) = 1._R8P/3._R8P ! stencil 1
        ! 2 => right interface (i+1/2)
        weno_a(2,0) = 1._R8P/3._R8P ! stencil 0
        weno_a(2,1) = 2._R8P/3._R8P ! stencil 1

        ! polinomials coefficients
        ! 1 => left interface (i-1/2)
        !  cell  0               ;    cell  1
        weno_p(1,0,0) =  0.5_R8P ; weno_p(1,1,0) =  0.5_R8P ! stencil 0
        weno_p(1,0,1) = -0.5_R8P ; weno_p(1,1,1) =  1.5_R8P ! stencil 1
        ! 2 => right interface (i+1/2)
        !  cell  0               ;    cell  1
        weno_p(2,0,0) =  1.5_R8P ; weno_p(2,1,0) = -0.5_R8P ! stencil 0
        weno_p(2,0,1) =  0.5_R8P ; weno_p(2,1,1) =  0.5_R8P ! stencil 1

        ! smoothness indicators coefficients
        ! stencil 0
        !      i*i             ;       (i-1)*i
        weno_d(0,0,0) = 1._R8P ; weno_d(1,0,0) =-2._R8P
        !      /               ;       (i-1)*(i-1)
        weno_d(0,1,0) = 0._R8P ; weno_d(1,1,0) = 1._R8P
        ! stencil 1
        !     (i+1)*(i+1)      ;       (i+1)*i
        weno_d(0,0,1) = 1._R8P ; weno_d(1,0,1) =-2._R8P
        !      /               ;        i*i
        weno_d(0,1,1) = 0._R8P ; weno_d(1,1,1) = 1._R8P
      case(3) ! 5th order WENO reconstruction
        ! central difference coefficients
        ! 1 => left interface (i-1/2)
        weno_c(1,1) =  1._R8P/60._R8P ! cell -3
        weno_c(1,2) = -7.5_R8P        ! cell -2
        weno_c(1,3) = 37._R8P/60._R8P ! cell -1
        weno_c(1,4) = 37._R8P/60._R8P ! cell  0
        weno_c(1,5) = -7.5_R8P        ! cell  1
        weno_c(1,6) =  1._R8P/60._R8P ! cell  2
        ! 2 => right interface (i+1/2)
        weno_c(1,1) =  1._R8P/60._R8P ! cell -2
        weno_c(1,2) = -7.5_R8P        ! cell -1
        weno_c(1,3) = 37._R8P/60._R8P ! cell  0
        weno_c(1,4) = 37._R8P/60._R8P ! cell  1
        weno_c(1,5) = -7.5_R8P        ! cell  2
        weno_c(1,6) =  1._R8P/60._R8P ! cell  3

        ! optimal weights
        ! 1 => left interface (i-1/2)
        weno_a(1,0) = 0.3_R8P ! stencil 0
        weno_a(1,1) = 0.6_R8P ! stencil 1
        weno_a(1,2) = 0.1_R8P ! stencil 2
        ! 2 => right interface (i+1/2)
        weno_a(2,0) = 0.1_R8P ! stencil 0
        weno_a(2,1) = 0.6_R8P ! stencil 1
        weno_a(2,2) = 0.3_R8P ! stencil 2

        ! polinomials coefficients
        ! 1 => left interface (i-1/2)
        !  cell  0                     ;    cell  1                     ;    cell  2
        weno_p(1,0,0) =  1._R8P/3._R8P ; weno_p(1,1,0) =  5._R8P/6._R8P ; weno_p(1,2,0) = -1._R8P/6._R8P ! stencil 0
        weno_p(1,0,1) = -1._R8P/6._R8P ; weno_p(1,1,1) =  5._R8P/6._R8P ; weno_p(1,2,1) =  1._R8P/3._R8P ! stencil 1
        weno_p(1,0,2) =  1._R8P/3._R8P ; weno_p(1,1,2) = -7._R8P/6._R8P ; weno_p(1,2,2) = 11._R8P/6._R8P ! stencil 2
        ! 2 => right interface (i+1/2)
        !  cell  0                     ;    cell  1                     ;    cell  2
        weno_p(2,0,0) = 11._R8P/6._R8P ; weno_p(2,1,0) = -7._R8P/6._R8P ; weno_p(2,2,0) =  1._R8P/3._R8P ! stencil 0
        weno_p(2,0,1) =  1._R8P/3._R8P ; weno_p(2,1,1) =  5._R8P/6._R8P ; weno_p(2,2,1) = -1._R8P/6._R8P ! stencil 1
        weno_p(2,0,2) = -1._R8P/6._R8P ; weno_p(2,1,2) =  5._R8P/6._R8P ; weno_p(2,2,2) =  1._R8P/3._R8P ! stencil 2

        ! smoothness indicators coefficients
        ! stencil 0
        !      i*i                      ;       (i-1)*i                   ;       (i-2)*i
        weno_d(0,0,0) =  10._R8P/3._R8P ; weno_d(1,0,0) = -31._R8P/3._R8P ; weno_d(2,0,0) =  11._R8P/3._R8P
        !      /                        ;       (i-1)*(i-1)               ;       (i-2)*(i-1)
        weno_d(0,1,0) =   0._R8P        ; weno_d(1,1,0) =  25._R8P/3._R8P ; weno_d(2,1,0) = -19._R8P/3._R8P
        !      /                        ;        /                        ;       (i-2)*(i-2)
        weno_d(0,2,0) =   0._R8P        ; weno_d(1,2,0) =   0._R8P        ; weno_d(2,2,0) =   4._R8P/3._R8P
        ! stencil 1
        !     (i+1)*(i+1)               ;        i*(i+1)                  ;       (i-1)*(i+1)
        weno_d(0,0,1) =   4._R8P/3._R8P ; weno_d(1,0,1) = -13._R8P/3._R8P ; weno_d(2,0,1) =   5._R8P/3._R8P
        !      /                        ;        i*i                      ;       (i-1)*i
        weno_d(0,1,1) =   0._R8P        ; weno_d(1,1,1) =  13._R8P/3._R8P ; weno_d(2,1,1) = -13._R8P/3._R8P
        !      /                        ;        /                        ;       (i-1)*(i-1)
        weno_d(0,2,1) =   0._R8P        ; weno_d(1,2,1) =   0._R8P        ; weno_d(2,2,1) =   4._R8P/3._R8P
        ! stencil 2
        !     (i+2)*(i+2)               ;       (i+1)*(i+2)               ;        i*(i+2)
        weno_d(0,0,2) =   4._R8P/3._R8P ; weno_d(1,0,2) = -19._R8P/3._R8P ; weno_d(2,0,2) =  11._R8P/3._R8P
        !      /                        ;       (i+1)*(i+1)               ;        i*(i+1)
        weno_d(0,1,2) =   0._R8P        ; weno_d(1,1,2) =  25._R8P/3._R8P ; weno_d(2,1,2) = -31._R8P/3._R8P
        !      /                        ;        /                        ;        i*i
        weno_d(0,2,2) =   0._R8P        ; weno_d(1,2,2) =   0._R8P        ; weno_d(2,2,2) =  10._R8P/3._R8P
      endselect
   endif
   print '(A)', self%mpih%myrankstr//'weno_cpu_object%initialize finish'
   endassociate
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
   if (.not.go_on_fail_.and.error>0) error stop self%mpih%myrankstr//'error: failed to load ['//INI_SECTION_NAME//'].(weno)'
   self%scheme = trim(adjustl(char_buff))
   select case(trim(adjustl(self%scheme)))
   case('weno-1')
      self%weno_s = 1
   case('weno-3')
      self%weno_s = 2
   case('weno-5')
      self%weno_s = 3
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

   vr = 0._R_P
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

   vp = 0._R_P
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
         IS(f,s1) = 0._R_P
         do s2=0,S-1
            do s3=0,S-1
               IS(f,s1) = IS(f,s1) + self%weno_d(s3,s2,s1) * v(f,s1-s3) * v(f,s1-s2)
            enddo
         enddo
      enddo
   enddo
   ! compute alfa coefficients
   a_tot = 0._R_P
   do s1=0,S-1
      do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
         a(f,s1) = self%weno_a(f,s1) * (1._R_P / (self%weno_eps + IS(f,s1))**s)
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
