
   attributes(device) subroutine compute_eigenvectors_mhd_device(si,sir,b,i,j,k,ngc,nv,g,q_aux_gpu,el,er)
   ! Compute eigenvectors centered in inteface i,j,k/ip,jp,kp.
   integer(I4P), intent(in)         :: si(3)                                 !< Stencil increment.
   real(R8P)   , intent(in)         :: sir(3)                                !< Stencil increment, real cast.
   integer(I4P), intent(in), value  :: b, i, j, k                            !< Counter.
   integer(I4P), intent(in), value  :: ngc                                   !< Ghost cells number.
   integer(I4P), intent(in), value  :: nv                                    !< Number of conservative varibales.
   real(R8P),    intent(in), value  :: g                                     !< Specific heats ratio.
   real(R8P),    intent(in), device :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout)      :: el(1:9,1:9), er(1:9,1:9)              !< Left and right eigenvectors.
   real(R8P)                        :: uu, vv, ww, h, qq, c, ci, b1, b2      !< Roe average states.
   real(R8P)                        :: uvw, uvw_r1, uvw_r2                   !< Velocity rotation accordingly dir.

   call compute_roe_average_device(q_aux_gpu=q_aux_gpu, g=g, ngc=ngc, b=b, i=i, j=j, k=k, ip=i+si(1), jp=j+si(2), kp=k+si(3), &
                                   uu=uu, vv=vv, ww=ww, h=h, qq=qq, c=c, ci=ci, b1=b1, b2=b2)

   uvw    =  uu*sir(1)+vv*sir(2)+ww*sir(3)
   uvw_r1 =  uu*sir(3)+vv*sir(1)+ww*sir(2)
   uvw_r2 = -uu*sir(2)+vv*sir(3)+ww*sir(1)

   er(1,1)=1._R8P ; er(1,2)=uu-c*sir(1)   ; er(1,3)=vv-c*sir(2) ; er(1,4)=ww-c*sir(3)   ; er(1,5)=h-uvw*c
   er(2,1)=1._R8P ; er(2,2)=uu            ; er(2,3)=vv          ; er(2,4)=ww            ; er(2,5)=qq
   er(3,1)=1._R8P ; er(3,2)=uu+c*sir(1)   ; er(3,3)=vv+c*sir(2) ; er(3,4)=ww+c*sir(3)   ; er(3,5)=h+uvw*c
   er(4,1)=0._R8P ; er(4,2)=sir(2)+sir(3) ; er(4,3)=sir(1)      ; er(4,4)=0._R8P        ; er(4,5)=uvw_r1
   er(5,1)=0._R8P ; er(5,2)=0._R8P        ; er(5,3)=sir(3)      ; er(5,4)=sir(1)+sir(2) ; er(5,5)=uvw_r2

   el(1,1)= 0.5_R8P*(b1+uvw*ci)      ;el(1,2)= 1._R8P-b1;el(1,3)= 0.5_R8P*(b1-uvw*ci)      ;el(1,4)=-uvw_r1;el(1,5)=-uvw_r2
   el(2,1)=-0.5_R8P*(b2*uu+ci*sir(1));el(2,2)= b2*uu    ;el(2,3)=-0.5_R8P*(b2*uu-ci*sir(1));el(2,4)= sir(3);el(2,5)=-sir(2)
   el(3,1)=-0.5_R8P*(b2*vv+ci*sir(2));el(3,2)= b2*vv    ;el(3,3)=-0.5_R8P*(b2*vv-ci*sir(2));el(3,4)= sir(1);el(3,5)= sir(3)
   el(4,1)=-0.5_R8P*(b2*ww+ci*sir(3));el(4,2)= b2*ww    ;el(4,3)=-0.5_R8P*(b2*ww-ci*sir(3));el(4,4)= sir(2);el(4,5)= sir(1)
   el(5,1)= 0.5_R8P*b2               ;el(5,2)=-b2       ;el(5,3)= 0.5_R8P*b2               ;el(5,4)= 0._R8P;el(5,5)= 0._R8P
    ! x
   betax = 1._R8P
   betay = By / sqrt(By**2+Bz**2)
   betaz = Bz / sqrt(By**2+Bz**2)
   ! y
   betax = Bx / sqrt(Bx**2+Bz**2)
   betay = 1._R8P
   betaz = Bz / sqrt(Bx**2+Bz**2)
   ! z
   betax = Bx / sqrt(Bx**2+By**2)
   betay = By / sqrt(Bx**2+By**2)
   betaz = 1._R8P

   betax = sir(1)*1._R8P               + sir(2)*Bx/sqrt(Bx**2+Bz**2) + sir(3)*Bx/sqrt(Bx**2+By**2)
   betay = sir(1)*By/sqrt(By**2+Bz**2) + sir(2)*1._R8P               + sir(3)*By/sqrt(Bx**2+By**2)
   betaz = sir(1)*Bz/sqrt(By**2+Bz**2) + sir(2)*Bz/sqrt(Bx**2+Bz**2) + sir(3)*1._R8P

   Bmod2 = Bx**2+By**2+Bz**2
   gp_Bmod2 = g*p + Bmod2
   cf = sqrt(0.5_R8P*(gp_Bmod2 + sqrt(gp_Bmod2**2-4._R8P*g*p*(sir(1)*Bx**2+sir(1)*By**2+sir(3)*Bz**2)))/rho)
   cs = sqrt(0.5_R8P*(gp_Bmod2 - sqrt(gp_Bmod2**2-4._R8P*g*p*(sir(1)*Bx**2+sir(1)*By**2+sir(3)*Bz**2)))/rho)
   ca = abs(sir(1)*Bx**2+sir(1)*By**2+sir(3)*Bz**2)/sqrt(rho)
   ch = max(abs(uu)+cfx,abs(vv)+cfy,abs(ww)+cfz) ! cfx, cfy, cfz medie calcolate prima su tutte direzioni

   ! compute right and left eigenvectors matrices (at Roe state)
   er(1,1)=0._R8P ; er(1,2)=alphaF                ; er(1,3)=0._R8P             ; er(1,4)=alphaS                ; er(1,5)=1._R8P        ; er(1,6)=alphaS                ; er(1,7)=0._R8P             ; er(1,8)=alphaF                ; er(1,9)=0._R8P
   er(2,1)=0._R8P ; er(2,2)=alphaF*lambda2        ; er(2,3)=0._R8P             ; er(2,4)=alphaS*lamda4         ; er(2,5)=uu            ; er(2,6)=alphaS*lambda6        ; er(2,7)=0._R8P             ; er(2,8)=alphaF*lambda8        ; er(2,9)=0._R8P
   er(3,1)=0._R8P ; er(3,2)=alphaF*vv + Jf0*betaY ; er(3,3)=-betaZ*S           ; er(3,4)=alphaS*vv - Js0*betaY ; er(3,5)=vv            ; er(3,6)=alphaS*vv + Js0*betaY ; er(3,7)=-betaZ*S           ; er(3,8)=alphaF*vv - Jf0*betaY ; er(3,9)=0._R8P
   er(4,1)=0._R8P ; er(4,2)=alphaF*ww + Jf0*betaZ ; er(4,3)= betaY*S           ; er(4,4)=alphaS*ww - Js0*betaZ ; er(4,5)=ww            ; er(4,6)=alphaS*ww + Js0*betaZ ; er(4,7)= betaY*S           ; er(4,8)=alphaF*ww - Jf0*betaZ ; er(4,9)=0._R8P
   er(5,1)=1._R8P ; er(5,2)=0._R8P                ; er(5,3)=0._R8P             ; er(5,4)=0._R8P                ; er(5,5)=0._R8P        ; er(5,6)=0._R8P                ; er(5,7)=0._R8P             ; er(5,8)=0._R8P                ; er(5,9)=1._R8P
   er(6,1)=0._R8P ; er(6,2)=Jf1*betaY             ; er(6,3)=-betaZ*rho**(-0.5) ; er(6,4)=-Js1*BetaY            ; er(6,5)=0._R8P        ; er(6,6)=-Js1*betaY            ; er(6,7)=betaZ*rho**(-0.5)  ; er(6,8)=Jf1*betaY             ; er(6,9)=0._R8P
   er(7,1)=0._R8P ; er(7,2)=Jf1*betaZ             ; er(7,3)= betaY*rho**(-0.5) ; er(7,4)=-Js1*BetaZ            ; er(7,5)=0._R8P        ; er(7,6)=-Js1*betaZ            ; er(7,7)=-betaY*rho**(-0.5) ; er(7,8)=Jf1*betaZ             ; er(7,9)=0._R8P
   er(8,1)=0._R8P ; er(8,2)=Hf-GAMMAf             ; er(8,3)=-GAMMAa            ; er(8,4)=Hs-GAMMAs             ; er(8,5)=0.5_R8P*ni**2 ; er(8,6)=Hs+GAMMAs             ; er(8,7)=-GAMMAa            ; er(8,8)=Hf+GAMMAf             ; er(8,9)=0._R8P
   er(9,1)=-ch    ; er(9,2)=0._R8P                ; er(9,3)=0._R8P             ; er(9,4)=0._R8P                ; er(9,5)=0._R8P        ; er(9,6)=0._R8P                ; er(9,7)=0._R8P             ; er(9,8)=0._R8P                ; er(9,9)=ch
   !da scrivere la matrice autovettoriSX + paramteri
   el(1,1)=0._R8P                                      ; el(1,2)=0._R8P                          ; el(1,3)= 0._R8P                          ; el(1,4)=0._R8P                          ; el(1,5)=0._R8P ; el(1,6)= 0._R8P                              ; el(1,7)= 0._R8P                               ; el(1,8)= 0._R8P                               ; el(1,9)=-1._R8P/(2*ch)
   el(2,1)=0.5_R8P/(a**2)*(gamma1*alphaF*ni**2+GAMMAf) ; el(2,2)=0.5_R8P/(a**2)*(Ifuu-alphaF*cf) ; el(2,3)= 0.5_R8P/(a**2)*(Ifvv+Jf0*betaY) ; el(2,4)=0.5_R8P/(a**2)*(Ifww+Jf0*betaZ) ; el(2,5)=0._R8P ; el(2,6)= 0.5_R8P/(a**2)*(IfBY+Jf1*betaY)     ; el(2,7)= 0.5_R8P/(a**2)2*(IfBZ+Jf1*betaZ)     ; el(2,8)= 0.5_R8P/(a**2)*alphaF*(GAMMA-1._R8P) ; el(2,9)= 0._R8P
   el(3,1)=0.5_R8P*GAMMAa                              ; el(3,2)=0._R8P                          ; el(3,3)=-0.5_R8P*betaZ*S                 ; el(3,4)=0.5_R8P*betaY*S                 ; el(3,5)=0._R8P ; el(3,6)=-0.5_R8P*(rho**0.5)*BetaZ            ; el(3,7)= 0.5_R8P*(rho**0.5)*BetaY             ; el(3,8)= 0._R8P                               ; el(3,9)= 0._R8P
   el(4,1)=0.5_R8P/(a**2)*(gamma1*alphaS*ni**2+GAMMAs) ; el(4,2)=0.5_R8P/(a**2)*(Isuu-alphaS*cs) ; el(4,3)= 0.5_R8P/(a**2)*(Isvv-Js0*betaY) ; el(4,4)=0.5_R8P/(a**2)*(Isww-Js0*betaZ) ; el(4,5)=0._R8P ; el(4,6)= 0.5_R8P/(a**2)*(IsBY-Js1*rho*betaY) ; el(4,7)= 0.5_R8P/(a**2)2*(IsBZ-Jf1*rho*betaZ) ; el(4,8)= 0.5_R8P/(a**2)*alphaS*(GAMMA-1._R8P) ; el(4,9)= 0._R8P
   el(5,1)=1._R8P-0.5_R8P*GAMMAm*bi**2                 ; el(5,2)=GAMMAm*uu                       ; el(5,3)= GAMMAm*vv                       ; el(5,4)=GAMMAm*ww                       ; el(5,5)=0._R8P ; el(5,6)= GAMMA*By                            ; el(5,7)= GAMMAm*Bz                            ; el(5,8)=-GAMMAm                               ; el(5,9)= 0._R8P
   el(5,1)=0.5_R8P/(a**2)*(gamma1*alphaS*ni**2-GAMMAs) ; el(6,2)=0.5_R8P/(a**2)*(Isuu+alphaS*cs) ; el(6,3)= 0.5_R8P/(a**2)*(Isvv+Js0*betaY) ; el(6,4)=0.5_R8P/(a**2)*(Isww+Js0*betaZ) ; el(6,5)=0._R8P ; el(6,6)= 0.5_R8P/(a**2)*(IsBY-Js1*rho*betaY) ; el(6,7)= 0.5_R8P/(a**2)2*(IsBZ-Jf1*rho*betaZ) ; el(6,8)= 0.5_R8P/(a**2)*alphaS*(GAMMA-1._R8P) ; el(6,9)= 0._R8P
   el(7,1)=0.5_R8P*GAMMAa                              ; el(7,2)=0._R8P                          ; el(7,3)=-0.5_R8P*betaZ*S                 ; el(7,4)=0.5_R8P*betaY*S                 ; el(7,5)=0._R8P ; el(7,6)= 0.5_R8P*(rho**0.5)*BetaZ            ; el(7,7)=-0.5_R8P*(rho**0.5)*BetaY             ; el(7,8)= 0._R8P                               ; el(7,9)= 0._R8P
   el(8,1)=0.5_R8P/(a**2)*(gamma1*alphaF*ni**2-GAMMAf) ; el(8,2)=0.5_R8P/(a**2)*(Ifuu+alphaF*cf) ; el(8,3)= 0.5_R8P/(a**2)*(Ifvv-Jf0*betaY) ; el(8,4)=0.5_R8P/(a**2)*(Ifww-Jf0*betaZ) ; el(8,5)=0._R8P ; el(8,6)= 0.5_R8P/(a**2)*(IfBY+Jf1*rho*betaY) ; el(8,7)= 0.5_R8P/(a**2)*(IfBZ+Jf1*rho*betaZ)  ; el(8,8)= 0.5_R8P/(a**2)*alphaF*(GAMMA-1._R8P) ; el(8,9)= 0._R8P
   el(9,1)=0._R8P                                      ; el(9,2)=0._R8P                          ; el(9,3)= 0._R8P                          ; el(9,4)=0.5_R8P                         ; el(9,5)=0._R8P ; el(9,6)= 0._R8P                              ; el(9,7)= 0._R8P                               ; el(9,8)= 0._R8P                               ; el(9,9)= 1._R8P/(2*ch)
   endsubroutine compute_eigenvectors_mhd_device
     ! MHD
   ! conservative variables
   ! rho
   ! rho*u
   ! rho*v
   ! rho*w
   ! Bx
   ! By
   ! Bz
   ! rho*E
   ! psi
   ! auxiliary variable
   ! rho
   ! u
   ! v
   ! w
   attributes(device) subroutine compute_average_mhd_device(ngc, b, i, j, k, ip, jp, kp, g, q_aux_gpu, &
                                                            uu, vv, ww, h, qq, c, ci, b1, b2)
   !< Compute Roe averaged quantities.
   integer(I4P), intent(in)         :: ngc                                       !< Number of ghost cells.
   integer(I4P), intent(in)         :: b, i, j, k, ip, jp, kp                    !< Left/right cells indexes.
   real(R8P),    intent(in)         :: g                                         !< Specific heats ratio.
   real(R8P),    intent(in), device :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:) !< Auxiliary variables.
   real(R8P),    intent(out)        :: uu, vv, ww, h, qq, c, ci, b1, b2          !< Roe state average variables.
   real(R8P)                        :: cc                                        !< Local varbiables.

   r  =  0.5_R8P*(q_aux_gpu(b,i,j,k,1) + q_aux_gpu(b,ip,jp,kp,1))
   uu =  0.5_R8P*(q_aux_gpu(b,i,j,k,2) + q_aux_gpu(b,ip,jp,kp,2))
   vv =  0.5_R8P*(q_aux_gpu(b,i,j,k,3) + q_aux_gpu(b,ip,jp,kp,3))
   ww =  0.5_R8P*(q_aux_gpu(b,i,j,k,4) + q_aux_gpu(b,ip,jp,kp,4))
   h  =  0.5_R8P*(q_aux_gpu(b,i,j,k,8) + q_aux_gpu(b,ip,jp,kp,8))
   bx =  0.5_R8P*(q_aux_gpu(b,i,j,k,4) + q_aux_gpu(b,ip,jp,kp,4))
   qq =  0.5_R8P*(uu*uu+vv*vv+ww*ww)
   cc =  (g-1._R8P) * (h - qq)
   c  =  sqrt(cc)
   ci =  1._R8P/c
   b2 = (g-1)/cc  ! alias 1/(cp*theta)
   b1 = b2 * qq   ! alias q/(cp*theta)

   endsubroutine compute_average_mhd_device
