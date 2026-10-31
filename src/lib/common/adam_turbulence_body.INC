   ! Body of the CPU eddy viscosity kernel of one model: ADAM_SGS_FORMULA names the pointwise routine (Fortran has no template to
   ! share a body across routines that differ by a callee, and a procedure dummy would be a call through a pointer per cell).
   integer(I4P), intent(in)    :: ir, iu                            !< Auxiliary indexes of the density and of u.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)    :: c                                 !< Model constant.
   real(R8P),    intent(in)    :: dxyz(1:,1:)                       !< Blocks space steps [3, nb].
   logical,      intent(in)    :: is_null(3)                        !< Null directions.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: mu_t(1-ngc:,1-ngc:,1-ngc:,1:)     !< Eddy viscosity.
   integer(I4P)                :: s(3)                              !< Steps along x, y, z: 1 active, 0 null.
   real(R8P)                   :: rdx(3)                            !< Inverse steps, 0 along null.
   real(R8P)                   :: dl(3)                             !< Steps, 0 along null.
   real(R8P)                   :: delta                             !< Filter width.
   real(R8P)                   :: g(3,3)                            !< Velocity gradient.
   real(R8P)                   :: nu                                !< Kinematic eddy viscosity.
   integer(I4P)                :: b, i, j, k, a                     !< Counters.

   s = merge(0_I4P, 1_I4P, is_null)
   !$omp parallel do collapse(4) default(firstprivate) shared(dxyz, q_aux, mu_t)
   do b=1, blocks_number
      do k=1-2*s(3), nk+2*s(3)
         do j=1-2*s(2), nj+2*s(2)
            do i=1-2*s(1), ni+2*s(1)
               dl = real(s, R8P) * dxyz(:,b)
               rdx = 0._R8P ; where (s == 1_I4P) rdx = 1._R8P / dxyz(:,b)
               delta = product(dxyz(:,b), mask=s == 1_I4P)**(1._R8P / real(max(1_I4P, sum(s)), R8P))
               do a=1, 3
                  g(a,1) = 0.5_R8P * rdx(1) * (q_aux(iu+a-1,i+s(1),j,k,b) - q_aux(iu+a-1,i-s(1),j,k,b))
                  g(a,2) = 0.5_R8P * rdx(2) * (q_aux(iu+a-1,i,j+s(2),k,b) - q_aux(iu+a-1,i,j-s(2),k,b))
                  g(a,3) = 0.5_R8P * rdx(3) * (q_aux(iu+a-1,i,j,k+s(3),b) - q_aux(iu+a-1,i,j,k-s(3),b))
               enddo
               call ADAM_SGS_FORMULA(c=c, dl=dl, delta=delta, g=g, nu=nu)
               mu_t(i,j,k,b) = q_aux(ir,i,j,k,b) * nu
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
