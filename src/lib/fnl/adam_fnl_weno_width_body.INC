   ! Body of the FNL WENO upwind reconstruction for ONE stencil width, S a named constant: included by the routines
   ! `weno_reconstruct_upwind_s<N>_dev` and `weno_reconstruct_upwind_wratio_s<N>_dev` of `adam_fnl_weno_width_kernels`, after
   ! their declaration of `S` (the second flavour defines FNL_WENO_WRATIO: the same reconstruction plus the weights ratio).
   !
   ! It is `weno_reconstruct_upwind_dev` (`..._wratio_dev`) of `adam_fnl_weno_kernels` with its three helpers (polynomials,
   ! weights, convolution) merged, every operation in the same order, so that in the routine itself the stencil loops, the
   ! stencil index of the tables and the private arrays (`VP`, `w`, `IS`, `a`: sized by S, not by S_MAX) have constant bounds
   ! and the routine has no call left, whether or not the compiler inlines it into the face kernel of the same width (the
   ! AMD OpenMP-offload compiler does not: the generic routine stays a chain of calls with the stencil in scratch memory).
   ! The macro FNL_WENO_WRATIO only selects between two flavours of one routine (conditional compilation of a variant: Fortran
   ! has no template to share a body across a dummy list that differs by one argument).
   real(R8P),    intent(in)  :: weno_a(1:,0:,1:)    !< Optimal weights.
   real(R8P),    intent(in)  :: weno_p(1:,0:,0:,1:) !< Polinomials coefficients.
   real(R8P),    intent(in)  :: weno_d(0:,0:,0:,1:) !< Smoothness indicators coefficients.
   real(R8P),    intent(in)  :: weno_zeps           !< Parameter for avoiding division by zero in computing IS.
   real(R8P),    intent(in)  :: weno_rmu            !< Inverse of the descaler (1: Jiang-Shu weights).
   real(R8P),    intent(in)  :: V (1:2,1-S:-1+S)    !< Variables to be reconstructed.
   real(R8P),    intent(out) :: VR(1:2)             !< Left and right (1,2) interface value of reconstructed V.
#ifdef FNL_WENO_WRATIO
   real(R8P),    intent(out) :: WRATIO(1:2)         !< Left and right (1,2) interface smoothness ratio.
#endif
   real(R8P)                 :: VP  (1:2,0:S-1)     !< Polynomial reconstructions.
   real(R8P)                 :: w   (1:2,0:S-1)     !< Weights of the stencils.
   real(R8P)                 :: IS  (1:2,0:S-1)     !< Smoothness indicators of the stencils.
   real(R8P)                 :: a   (1:2,0:S-1)     !< Alpha coifficients for the weights.
   real(R8P)                 :: a_tot(1:2)          !< Summ of the alpha coefficients.
   integer(I4P)              :: wexp                !< Exponent of the smoothness indicators.
   integer(I4P)              :: s1,s2,s3,f,k        !< Counters.
   !$acc routine seq
   !$omp declare target

   ! polynomials
   VP = 0._R8P
   do s1=0, S-1 ! stencil counter
      do s2=0, S-1 ! cell counter counter
         do f=1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
            VP(f,s1) = VP(f,s1) + weno_p(f,s2,s1,S)*V(f,-s2+s1)
         enddo
      enddo
   enddo
   ! weights
   wexp = weno_weights_exponent(S=S)
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
   a_tot = 0._R8P
   do s1=0,S-1
      do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
         a(f,s1) = weno_a(f,s1,S)*(1._R8P/(weno_zeps+IS(f,s1)*weno_rmu*weno_rmu)**wexp) ; a_tot(f) = a_tot(f) + a(f,s1)
      enddo
   enddo
   do s1=0,S-1
      do f=1,2
         w(f,s1) = a(f,s1)/a_tot(f)
      enddo
   enddo
   ! convolution
   VR = 0._R8P
   do k=0, S-1
      do f=1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
         VR(f) = VR(f) + w(f,k)*VP(f,k)
      enddo
   enddo
#ifdef FNL_WENO_WRATIO
   ! smoothness ratio of each interface
   do f=1, 2
      WRATIO(f) = w(f,0) / weno_a(f,0,S)
      do k=1, S-1
         WRATIO(f) = min(WRATIO(f), w(f,k) / weno_a(f,k,S))
      enddo
   enddo
#endif
