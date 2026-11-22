   ! Body of the PRISM CPU kernel `compute_fluxes_convective_weno_s<N>` for ONE WENO stencil width, S a named constant
   ! (issue #139, S5, decision D3): the face loops along `dir` with, per face, the Lax-Friedrichs split of the convective
   ! fluxes projected on the characteristic fields (the former `decompose_fluxes_convective`), the WENO reconstruction of
   ! each field and the back projection (the former `compute_fluxes_convective_ri_weno`), every operation in their order.
   ! The kernel is SETH's WENO-JS of width S (`weno_width`, used by the including routine), in SETH's stencil convention:
   ! `vp(m)` the + flux of cell m (left bias), `vm(m)` the - flux of cell m (right bias), cells by offset from the face.
   integer(I4P), intent(in)    :: dir                                !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)    :: blocks_number                      !< Number of blocks.
   integer(I4P), intent(in)    :: ni                                 !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                 !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                 !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                !< Ghost cells number.
   integer(I4P), intent(in)    :: nv_c                               !< Number of conservative varibales.
   real(R8P),    intent(in)    :: weno_eps                           !< Regulariser of the WENO smoothness indicators.
   real(R8P),    intent(in)    :: weno_sigma                         !< WENO descaler switch (0 Jiang-Shu, 1 scale-invariant).
   real(R8P),    intent(in)    :: evmax                              !< Maximum waves speed estimation.
   real(R8P),    intent(in)    :: erw(1:,1:,1:)                      !< Right eigenvectors for WENO reconstruction.
   real(R8P),    intent(in)    :: elw(1:,1:,1:)                      !< Left  eigenvectors for WENO reconstruction.
   real(R8P),    intent(in)    :: chi                                !< Speed parameter for divergence cleaning.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Field variables.
   real(R8P),    intent(inout) :: fluxes(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Fluxes.
   integer(I4P)                :: si(3), si_i, si_j, si_k            !< Directional (1=x,2=y,3=z) increment.
   real(R8P)                   :: sir(3)                             !< Directional (1=x,2=y,3=z) increment, real.
   real(R8P)                   :: vp(1-S:S,1:NV_MAX)                 !< Fluxes + of the characteristic fields (left bias).
   real(R8P)                   :: vm(1-S:S,1:NV_MAX)                 !< Fluxes - of the characteristic fields (right bias).
   real(R8P)                   :: rmu                                !< Inverse descaler of one field.
   real(R8P)                   :: fpmr(1:2,1:NV_MAX)                 !< Fluxes +- reconstructed.
   real(R8P)                   :: fmp(2)                             !< Fluxes -+ decomposition in each cell stencils.
   real(R8P)                   :: gc, wc                             !< Increments for fluxes decomposition.
   real(R8P)                   :: f(NV_MAX)                          !< Conservative fluxes.
   integer(I4P)                :: b, i, j, k                         !< Counter.
   integer(I4P)                :: v, vv, st, is, js, ks              !< Counter.

   select case(dir)
   case(1)
      si = [1,0,0]
   case(2)
      si = [0,1,0]
   case(3)
      si = [0,0,1]
   endselect
   sir = real(si,R8P)
   si_i = 1-si(1)
   si_j = 1-si(2)
   si_k = 1-si(3)
   vp(S,:) = 0._R8P ; vm(1-S,:) = 0._R8P ! not read by the kernel (SETH's stencil convention)

   do b=1, blocks_number
   do k=si_k, nk
   do j=si_j, nj
   do i=si_i, ni
      ! decomposition of the convective fluxes on the stencil
      do st=1-S, S
         is = i + (st) * si(1) ; js = j + (st) * si(2) ; ks = k + (st) * si(3)
         call compute_fluxes_maxwell(sir=sir,q=q(:,is,js,ks,b),f=f,chi=chi)
         do v=1, nv_c
            wc = 0._R8P
            gc = 0._R8P
            do vv=1, nv_c
               wc = wc + elw(vv,v,dir) * q(vv,is,js,ks,b)
               gc = gc + elw(vv,v,dir) * f(vv)
            enddo
            fmp(2) = 0.5_R8P * (gc + evmax * wc)
            fmp(1) = gc - fmp(2)
            if (st<S)   vp(st,v) = fmp(2)
            if (st>1-S) vm(st,v) = fmp(1)
         enddo
      enddo
      ! reconstruction
      do v=1, nv_c
         rmu = seth_weno_descaling(S=S, sigma=weno_sigma, vp=vp(:,v), vm=vm(:,v))
         call weno_width(vp=vp(:,v), vm=vm(:,v), eps=weno_eps, rmu=rmu, rp=fpmr(2,v), rm=fpmr(1,v))
      enddo
      ! back projection in conservative variables space
      do v=1, nv_c
         fluxes(v,i,j,k,b) = 0._R8P
         do vv=1,nv_c
            fluxes(v,i,j,k,b) = fluxes(v,i,j,k,b) + erw(vv,v,dir) * (fpmr(1,vv) + fpmr(2,vv))
         enddo
      enddo
   enddo
   enddo
   enddo
   enddo
