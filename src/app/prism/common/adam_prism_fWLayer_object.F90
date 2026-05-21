!< ADAM, PRISM (Plasma Research usIng Simulation Methods) fWLayer class definition, common backend.
module adam_prism_fWLayer_object

    !in input prendi il numero di celle che compone lo strato, chiamalo C e 3 flag per definire su quali lati ho lo strato
    !assegni la funzione f ad ognuna delle celle: essa dovrà avere tre elementi in modo tale da poter
    !tenere conto delle celle sulla diagonale

    !sicuramente devi dargli in pasto field o comunque un array di celle che contenga le coordinate delle celle stesse

    !ragiona su come scrivere (in cpu però, non qui) la funzione di aggiornamento dei campi. La puoi fare ricorsiva e senza if se
    !vale a livello matematico, altrimenti dovrai aggiungere un flag(3) per individuare se la cella appartiene a uno o più lati dello strato e a quali (in ogni elemento + o -1)

! ADAM singleton objects
use :: adam_mpih_global,  only : mpih
use :: adam_grid_global,  only : grid
use :: adam_field_object, only : field_object
! PRISM modules
use :: adam_prism_parameters
use :: adam_prism_physics_object, only : prism_physics_object
! third party modules
use :: finer, only : file_ini
use :: penf,  only : I4P, R8P, str
use :: stringifor

implicit none
private
public :: prism_fWLayer_object
public :: apply_fWL_correction_fun

character(len=7), parameter :: INI_SECTION_NAME='fWLayer' !< INI file section name containing flWLayer datas.

type :: prism_fWLayer_object
   !< PRISM fWLayer class definition.
   logical                   :: layer(6) = .false.                    !< Layer flags for each side (-x, +x, -y, +y, -z, +z).
   integer(I4P)              :: C        = 0_I4P                      !< Layer cell width.
   integer(I4P)              :: ni_fWL(2,6), nj_fWL(2,6), nk_fWL(2,6) !< Dimensions of FWL domain.
   real(R8P)                 :: s2(6)                                 !< Side coefficient.
   integer(I4P)              :: n(6)                                  !< FWL f function index.
   integer(I4P)              :: alfa_D(6), beta_D(6)                  !< Corrected var index of D (Barbas' notation).
   integer(I4P)              :: alfa_B(6), beta_B(6)                  !< Corrected var index of B (Barbas' notation).
   real(R8P), allocatable    :: f(:,:,:,:,:)                          !< fWLayer function values.
   type(string), allocatable :: f_name(:)                             !< fWLayer function values names.
contains
  ! public methods
  procedure, pass(self) :: description    !< Return pretty-printed object description.
  procedure, pass(self) :: initialize     !< Initialize physics.
  procedure, pass(self) :: load_from_file !< Load config from file.
endtype prism_fWLayer_object

contains
   ! public methods
   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_fWLayer_object), intent(in) :: self             !< Physics.
   character(len=:), allocatable           :: desc             !< Description.
   character(len=1), parameter             :: NL=new_line('a') !< New line character.

   if (self%C == 0_I4P) then
      desc = mpih%myrankstr//'   No fWLayer implemented'
   else
      desc = mpih%myrankstr//'   fWLayer datas:'
      desc = desc//NL//mpih%myrankstr//'      Layer cell width: '//trim(str(self%C))
      desc = desc//NL//mpih%myrankstr//'      Layer on -x side: '//trim(str(self%layer(1)))
      desc = desc//NL//mpih%myrankstr//'      Layer on +x side: '//trim(str(self%layer(2)))
      desc = desc//NL//mpih%myrankstr//'      Layer on -y side: '//trim(str(self%layer(3)))
      desc = desc//NL//mpih%myrankstr//'      Layer on +y side: '//trim(str(self%layer(4)))
      desc = desc//NL//mpih%myrankstr//'      Layer on -z side: '//trim(str(self%layer(5)))
      desc = desc//NL//mpih%myrankstr//'      Layer on +z side: '//trim(str(self%layer(6)))
    endif
   endfunction description

   subroutine initialize(self, field, file_parameters, physics)
   !< Initialize the fWLayer.
   class(prism_fWLayer_object), intent(inout) :: self            !< fWLayer.
   type(field_object), intent(in) :: field !< Field (sibling realm component, threaded in).
   type(file_ini),              intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   type(prism_physics_object),  intent(in)    :: physics         !< Physics.
   integer(I4P)                               :: i,j,k,b         !< Counters.
   real(R8P)                                  :: fi              !< Cell function
   real(R8P)                                  :: distance        !< Distance between cell and physical boundary
   real(R8P)                                  :: C_r             !< Number (real) of ghost cells in the layer
   real(R8P)                                  :: i_r, j_r, k_r   !< (Real) counters
   real(R8P)                                  :: ds              !< Cells distance in x, y or z.

   print '(A)', mpih%myrankstr//'prism_fWLayer_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   if (self%C == 0_I4P) return

   !Inizializzo funzione f nelle celle dello strato
   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, blocks_number=>field%blocks_number,          &
            ngc=>grid%ngc, nb=>field%nb, dx=>field%dxyz(1,:), dy=>field%dxyz(2,:), dz=>field%dxyz(3,:), &
            C=>self%C, ni_fWL=>self%ni_fWL, nj_fWL=>self%nj_fWL, nk_fWL=>self%nk_fWL, n=>self%n, s2=>self%s2, &
            alfa_D=>self%alfa_D, alfa_B=>self%alfa_B, beta_D=>self%beta_D, beta_B=>self%beta_B)

   allocate(self%f(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   allocate(self%f_name(1:3))
   self%f = 1.0_R8P
   self%f_name(1) = 'fWLayer_x'
   self%f_name(2) = 'fWLayer_y'
   self%f_name(3) = 'fWLayer_z'
   C_r = real(C, R8P)

   if (C >0) then
      ni_fWL(1,1)=1_I4P; ni_fWL(1,2)=ni-C+1_I4P; ni_fWL(1,3)=1_I4P; ni_fWL(1,4)=1_I4P     ; ni_fWL(1,5)=1_I4P; ni_fWL(1,6)=1_I4P
      ni_fWL(2,1)=C    ; ni_fWL(2,2)=ni        ; ni_fWL(2,3)=ni   ; ni_fWL(2,4)=ni        ; ni_fWL(2,5)=ni   ; ni_fWL(2,6)=ni
      nj_fWL(1,1)=1_I4P; nj_fWL(1,2)=1_I4P     ; nj_fWL(1,3)=1_I4P; nj_fWL(1,4)=nj-C+1_I4P; nj_fWL(1,5)=1_I4P; nj_fWL(1,6)=1_I4P
      nj_fWL(2,1)=nj   ; nj_fWL(2,2)=nj        ; nj_fWL(2,3)=C    ; nj_fWL(2,4)=nj        ; nj_fWL(2,5)=nj   ; nj_fWL(2,6)=nj
      nk_fWL(1,1)=1_I4P; nk_fWL(1,2)=1_I4P ; nk_fWL(1,3)=1_I4P; nk_fWL(1,4)=1_I4P     ; nk_fWL(1,5)=1_I4P; nk_fWL(1,6)=nk-C+1_I4P
      nk_fWL(2,1)=nk   ; nk_fWL(2,2)=nk        ; nk_fWL(2,3)=nk   ; nk_fWL(2,4)=nk        ; nk_fWL(2,5)=C    ; nk_fWL(2,6)=nk

      n(1)      =1_I4P  ; n(2)      =1_I4P   ; n(3)      =2_I4P  ; n(4)      =2_I4P   ; n(5)      =3_I4P  ; n(6)      =3_I4P
      s2(1)     =1.0_R8P; s2(2)     =-1.0_R8P; s2(3)     =1.0_R8P; s2(4)     =-1.0_R8P; s2(5)     =1.0_R8P; s2(6)     =-1.0_R8P
      alfa_D(1) =2_I4P  ; alfa_D(2) =2_I4P   ; alfa_D(3) =3_I4P  ; alfa_D(4) =3_I4P   ; alfa_D(5) =1_I4P  ; alfa_D(6) =1_I4P
      beta_D(1) =3_I4P  ; beta_D(2) =3_I4P   ; beta_D(3) =1_I4P  ; beta_D(4) =1_I4P   ; beta_D(5) =2_I4P  ; beta_D(6) =2_I4P
      alfa_B(1) =5_I4P  ; alfa_B(2) =5_I4P   ; alfa_B(3) =6_I4P  ; alfa_B(4) =6_I4P   ; alfa_B(5) =4_I4P  ; alfa_B(6) =4_I4P
      beta_B(1) =6_I4P  ; beta_B(2) =6_I4P   ; beta_B(3) =4_I4P  ; beta_B(4) =4_I4P   ; beta_B(5) =5_I4P  ; beta_B(6) =5_I4P

      if (C < 40_I4P) then
         fi = 1/150._R8P*(-7.0_R8P*C_r**2 + 255._R8P*C_r + 250._R8P)
      else
         fi = 25.0_R8P
      endif
      do b=1, blocks_number
         !x- side
         if (self%layer(1)) then
            ds = dx(b)
            do k=1,nk
               do j=1, nj
                  do i=1, C
                     i_r = real(i, R8P)
                     distance = i_r*ds-ds/2
                     self%f(1,i,j,k,b) = 1._R8P/fi*LOG10((distance)/(C_r*ds)*(10._R8P**fi-1._R8P)+1._R8P)
                  enddo
               enddo
            enddo
         endif
         !x+ side
         if (self%layer(2)) then
            ds = dx(b)
            do k=1,nk
               do j=1, nj
                  do i=ni-C+1_I4P, ni
                     i_r = real(i, R8P)
                     distance = (ni-i_r)*ds+ds/2
                     self%f(1,i,j,k,b) = 1._R8P/fi*LOG10((distance)/(C_r*ds)*(10._R8P**fi-1._R8P)+1._R8P)
                  enddo
               enddo
            enddo
         endif
         !y- side
         if (self%layer(3)) then
            ds = dy(b)
            do k=1,nk
               do i=1, ni
                  do j=1, C
                     j_r = real(j, R8P)
                     distance = j_r*ds-ds/2
                     self%f(2,i,j,k,b) = 1._R8P/fi*LOG10((distance)/(C_r*ds)*(10._R8P**fi-1._R8P)+1._R8P)
                  enddo
               enddo
            enddo
         endif
         !y+ side
         if (self%layer(4)) then
            ds = dy(b)
            do k=1,nk
               do i=1, ni
                  do j=nj-C+1_I4P, nj
                     j_r = real(j, R8P)
                     distance = (nj-j_r)*ds+ds/2
                     self%f(2,i,j,k,b) = 1._R8P/fi*LOG10((distance)/(C_r*ds)*(10._R8P**fi-1._R8P)+1._R8P)
                  enddo
               enddo
            enddo
         endif
            !z- side
         if(self%layer(5)) then
            ds = dz(b)
            do i=1,ni
               do j=1, nj
                  do k=1, C
                     k_r = real(k, R8P)
                     distance = k_r*ds-ds/2
                     self%f(3,i,j,k,b) = 1._R8P/fi*LOG10((distance)/(C_r*ds)*(10._R8P**fi-1._R8P)+1._R8P)
                  enddo
               enddo
            enddo
         endif
         !x+ side
         if (self%layer(6)) then
            ds = dz(b)
            do i=1,ni
               do j=1, nj
                  do k=nk-C+1_I4P, nk
                     k_r = real(k, R8P)
                     distance = (nk-k_r)*ds+ds/2
                     self%f(3,i,j,k,b) = 1._R8P/fi*LOG10((distance)/(C_r*ds)*(10._R8P**fi-1._R8P)+1._R8P)
                  enddo
               enddo
            enddo
         endif
      enddo
   endif
   endassociate
   print '(A)', mpih%myrankstr//'prism_fWLayer_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(prism_fWLayer_object), intent(inout)        :: self            !< Physics.
   type(file_ini),              intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                     intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                           :: go_on_fail_     !< Go on if load fails.
   integer(I4P)                                      :: error           !< Error status.
   character(99)                                     :: buff            !< Character buffer.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='C', val=self%C, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(C)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='x_minus_layer', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(x_minus_layer)')
   select case(trim(adjustl(buff)))
   case('NO', 'no', 'No', 'nO')
      self%layer(1) = .false.
   case('YES', 'yes', 'Yes', 'yES')
      self%layer(1) = .true.
   case default
      self%layer(1) = .false.
   endselect
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='x_plus_layer', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(x_plus_layer)')
   select case(trim(adjustl(buff)))
   case('NO', 'no', 'No', 'nO')
      self%layer(2) = .false.
   case('YES', 'yes', 'Yes', 'yES')
      self%layer(2) = .true.
   case default
      self%layer(2) = .false.
   endselect
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='y_minus_layer', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(y_minus_layer)')
   select case(trim(adjustl(buff)))
   case('NO', 'no', 'No', 'nO')
      self%layer(3) = .false.
   case('YES', 'yes', 'Yes', 'yES')
      self%layer(3) = .true.
   case default
      self%layer(3) = .false.
   endselect
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='y_plus_layer', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(y_plus_layer)')
   select case(trim(adjustl(buff)))
   case('NO', 'no', 'No', 'nO')
      self%layer(4) = .false.
   case('YES', 'yes', 'Yes', 'yES')
      self%layer(4) = .true.
   case default
      self%layer(4) = .false.
   endselect
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='z_minus_layer', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(z_minus_layer)')
   select case(trim(adjustl(buff)))
   case('NO', 'no', 'No', 'nO')
      self%layer(5) = .false.
   case('YES', 'yes', 'Yes', 'yES')
      self%layer(5) = .true.
   case default
      self%layer(5) = .false.
   endselect
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='z_plus_layer', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(z_plus_layer)')
   select case(trim(adjustl(buff)))
   case('NO', 'no', 'No', 'nO')
      self%layer(6) = .false.
   case('YES', 'yes', 'Yes', 'yES')
      self%layer(6) = .true.
   case default
      self%layer(6) = .false.
   endselect
   endsubroutine load_from_file

   subroutine apply_fWL_correction_fun(blocks_number,ngc,ni1,ni2,nj1,nj2,nk1,nk2,n,s2,alfa_D,beta_D,alfa_B,beta_B,f,q)
   !< Applay FWL correction, direction agnostic.
   integer(I4P), intent(in)    :: blocks_number                     !< Blocks number.
   integer(I4P), intent(in)    :: ngc                               !< Number of ghost cells.
   integer(I4P), intent(in)    :: ni1,ni2,nj1,nj2,nk1,nk2           !< Dimensions of FWL domain.
   integer(I4P), intent(in)    :: n                                 !< f component.
   real(R8P),    intent(in)    :: s2                                !< Side coefficient.
   integer(I4P), intent(in)    :: alfa_D, beta_D                    !< Corrected var index of D (Barbas' notation).
   integer(I4P), intent(in)    :: alfa_B, beta_B                    !< Corrected var index of D (Barbas' notation).
   real(R8P),    intent(in)    :: f(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< fWLayer function values.
   real(R8P),    intent(inout) :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Field variables.
   real(R8P)                   :: fm1, fp1                          !< fWLayer function values in -+ cell.
   integer(I4P)                :: b,i,j,k                           !< Counter.

   do b=1,blocks_number
      do k=nk1, nk2
         do j=nj1, nj2
            do i=ni1, ni2
               fm1 = f(n,i,j,k,b) - 1._R8P
               fp1 = f(n,i,j,k,b) + 1._R8P
               q(alfa_D,i,j,k,b) = MU0_SQ_I2  * ( s2*fm1*q(beta_B,i,j,k,b)*EPS0_SQ +    fp1*q(alfa_D,i,j,k,b)*MU0_SQ)
               q(beta_D,i,j,k,b) = MU0_SQ_I2  * (-s2*fm1*q(alfa_B,i,j,k,b)*EPS0_SQ +    fp1*q(beta_D,i,j,k,b)*MU0_SQ)
               q(alfa_B,i,j,k,b) = EPS0_SQ_I2 * (    fp1*q(alfa_B,i,j,k,b)*EPS0_SQ - s2*fm1*q(beta_D,i,j,k,b)*MU0_SQ)
               q(beta_B,i,j,k,b) = EPS0_SQ_I2 * (    fp1*q(beta_B,i,j,k,b)*EPS0_SQ + s2*fm1*q(alfa_D,i,j,k,b)*MU0_SQ)

   !                  q(alfa_D,i,j,k,b) = 1/(2*MU0**0.5_R8P)*(s2*(f(n,i,j,k,b)-1._R8P)*q(beta_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + (f(n,i,j,k,b)+1._R8P)*q(alfa_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(beta_D,i,j,k,b) = 1/(2*MU0**0.5_R8P)*(-s2*(f(n,i,j,k,b)-1._R8P)*q(alfa_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + (f(n,i,j,k,b)+1._R8P)*q(beta_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(alfa_B,i,j,k,b) = 1/(2*EPS0**0.5_R8P)*((f(n,i,j,k,b)+1._R8P)*q(alfa_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P - s2*(f(n,i,j,k,b)-1._R8P)*q(beta_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(beta_B,i,j,k,b) = 1/(2*EPS0**0.5_R8P)*((f(n,i,j,k,b)+1._R8P)*q(beta_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + s2*(f(n,i,j,k,b)-1._R8P)*q(alfa_D,i,j,k,b)*MU0**0.5_R8P)
            enddo
         enddo
      enddo
   enddo
   endsubroutine apply_fWL_correction_fun
endmodule adam_prism_fWLayer_object

   !subroutine apply_fWL_correction(self, q)
   !!< Apply correction if a fWL is present
   !class(prism_cpu_object), intent(inout) :: self                    !< The equation.
   !real(R8P),               intent(inout) :: q(1:,         &
   !                                            1-self%ngc:,&
   !                                            1-self%ngc:,&
   !                                            1-self%ngc:,&
   !                                            1:)                   !< Conservative variables.
   !real(R8P)                              :: s2                      !< Side coefficient
   !integer(I4P)                           :: i,j,k,b,n               !< Counters
   !integer(I4P)                           :: alfa_D, beta_D, gamma_D !< Indici alfa beta gamma come in Barbas.
   !integer(I4P)                           :: alfa_B, beta_B, gamma_B !< Indici alfa beta gamma come in Barbas.
   !associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, blocks_number=>self%blocks_number, &
   !   f=>self%fWLayer%f, layer=>self%fWLayer%layer, C=>self%fWLayer%C)
   !if (C>0) then
   !   !x- side
   !   do b=1,blocks_number
   !      if (layer(1)) then
   !         n = 1_I4P
   !         s2 = 1.0_R8P
   !         alfa_D = 2_I4P
   !         beta_D = 3_I4P
   !         gamma_D = 1_I4P
   !         alfa_B = 5_I4P
   !         beta_B = 6_I4P
   !         gamma_B = 4_I4P
   !         do k=1,nk
   !            do j=1, nj
   !               do i=1, C
   !                  q(alfa_D,i,j,k,b) = 1/(2*MU0**0.5_R8P)*(s2*(f(n,i,j,k,b)-1._R8P)*q(beta_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + (f(n,i,j,k,b)+1._R8P)*q(alfa_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(beta_D,i,j,k,b) = 1/(2*MU0**0.5_R8P)*(-s2*(f(n,i,j,k,b)-1._R8P)*q(alfa_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + (f(n,i,j,k,b)+1._R8P)*q(beta_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(alfa_B,i,j,k,b) = 1/(2*EPS0**0.5_R8P)*((f(n,i,j,k,b)+1._R8P)*q(alfa_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P - s2*(f(n,i,j,k,b)-1._R8P)*q(beta_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(beta_B,i,j,k,b) = 1/(2*EPS0**0.5_R8P)*((f(n,i,j,k,b)+1._R8P)*q(beta_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + s2*(f(n,i,j,k,b)-1._R8P)*q(alfa_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(gamma_D,i,j,k,b) = q(gamma_D,i,j,k,b)
!
   !                  q(gamma_B,i,j,k,b) = q(gamma_B,i,j,k,b)
   !               enddo
   !            enddo
   !         enddo
   !      endif
   !      !x+ side
   !      if(layer(2)) then
   !         n = 1_I4P
   !         s2 = -1.0_R8P
   !         alfa_D = 2_I4P
   !         beta_D = 3_I4P
   !         gamma_D = 1_I4P
   !         alfa_B = 5_I4P
   !         beta_B = 6_I4P
   !         gamma_B = 4_I4P
   !         do k=1,nk
   !            do j=1, nj
   !               do i=ni-C, ni
   !                  q(alfa_D,i,j,k,b) = 1/(2*MU0**0.5_R8P)*(s2*(f(n,i,j,k,b)-1._R8P)*q(beta_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + (f(n,i,j,k,b)+1._R8P)*q(alfa_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(beta_D,i,j,k,b) = 1/(2*MU0**0.5_R8P)*(-s2*(f(n,i,j,k,b)-1._R8P)*q(alfa_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + (f(n,i,j,k,b)+1._R8P)*q(beta_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(alfa_B,i,j,k,b) = 1/(2*EPS0**0.5_R8P)*((f(n,i,j,k,b)+1._R8P)*q(alfa_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P - s2*(f(n,i,j,k,b)-1._R8P)*q(beta_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(beta_B,i,j,k,b) = 1/(2*EPS0**0.5_R8P)*((f(n,i,j,k,b)+1._R8P)*q(beta_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + s2*(f(n,i,j,k,b)-1._R8P)*q(alfa_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(gamma_D,i,j,k,b) = q(gamma_D,i,j,k,b)
!
   !                  q(gamma_B,i,j,k,b) = q(gamma_B,i,j,k,b)
   !               enddo
   !            enddo
   !         enddo
   !      endif
   !      !y- side
   !      if (layer(3)) then
   !         n = 2_I4P
   !         s2 = 1.0_R8P
   !         alfa_D = 3_I4P
   !         beta_D = 1_I4P
   !         gamma_D = 2_I4P
   !         alfa_B = 6_I4P
   !         beta_B = 4_I4P
   !         gamma_B = 5_I4P
   !         do k=1,nk
   !            do i=1, ni
   !               do j=1, C
   !                  q(alfa_D,i,j,k,b) = 1/(2*MU0**0.5_R8P)*(s2*(f(n,i,j,k,b)-1._R8P)*q(beta_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + (f(n,i,j,k,b)+1._R8P)*q(alfa_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(beta_D,i,j,k,b) = 1/(2*MU0**0.5_R8P)*(-s2*(f(n,i,j,k,b)-1._R8P)*q(alfa_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + (f(n,i,j,k,b)+1._R8P)*q(beta_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(alfa_B,i,j,k,b) = 1/(2*EPS0**0.5_R8P)*((f(n,i,j,k,b)+1._R8P)*q(alfa_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P - s2*(f(n,i,j,k,b)-1._R8P)*q(beta_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(beta_B,i,j,k,b) = 1/(2*EPS0**0.5_R8P)*((f(n,i,j,k,b)+1._R8P)*q(beta_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + s2*(f(n,i,j,k,b)-1._R8P)*q(alfa_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(gamma_D,i,j,k,b) = q(gamma_D,i,j,k,b)
!
   !                  q(gamma_B,i,j,k,b) = q(gamma_B,i,j,k,b)
   !               enddo
   !            enddo
   !         enddo
   !      endif
   !      !y+ side
   !      if (layer(4)) then
   !         n = 2_I4P
   !         s2 = -1.0_R8P
   !         alfa_D = 3_I4P
   !         beta_D = 1_I4P
   !         gamma_D = 2_I4P
   !         alfa_B = 6_I4P
   !         beta_B = 4_I4P
   !         gamma_B = 5_I4P
   !         do k=1,nk
   !            do i=1, ni
   !               do j=nj-C, nj
   !                  q(alfa_D,i,j,k,b) = 1/(2*MU0**0.5_R8P)*(s2*(f(n,i,j,k,b)-1._R8P)*q(beta_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + (f(n,i,j,k,b)+1._R8P)*q(alfa_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(beta_D,i,j,k,b) = 1/(2*MU0**0.5_R8P)*(-s2*(f(n,i,j,k,b)-1._R8P)*q(alfa_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + (f(n,i,j,k,b)+1._R8P)*q(beta_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(alfa_B,i,j,k,b) = 1/(2*EPS0**0.5_R8P)*((f(n,i,j,k,b)+1._R8P)*q(alfa_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P - s2*(f(n,i,j,k,b)-1._R8P)*q(beta_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(beta_B,i,j,k,b) = 1/(2*EPS0**0.5_R8P)*((f(n,i,j,k,b)+1._R8P)*q(beta_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + s2*(f(n,i,j,k,b)-1._R8P)*q(alfa_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(gamma_D,i,j,k,b) = q(gamma_D,i,j,k,b)
!
   !                  q(gamma_B,i,j,k,b) = q(gamma_B,i,j,k,b)
   !               enddo
   !            enddo
   !         enddo
   !      endif
   !      !z- side
   !      if (layer(5)) then
   !         n = 3_I4P
   !         s2 = 1.0_R8P
   !         alfa_D = 1_I4P
   !         beta_D = 2_I4P
   !         gamma_D = 3_I4P
   !         alfa_B = 4_I4P
   !         beta_B = 5_I4P
   !         gamma_B = 6_I4P
   !         do i=1,ni
   !            do j=1, nj
   !               do k=1, C
   !                  q(alfa_D,i,j,k,b) = 1/(2*MU0**0.5_R8P)*(s2*(f(n,i,j,k,b)-1._R8P)*q(beta_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + (f(n,i,j,k,b)+1._R8P)*q(alfa_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(beta_D,i,j,k,b) = 1/(2*MU0**0.5_R8P)*(-s2*(f(n,i,j,k,b)-1._R8P)*q(alfa_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + (f(n,i,j,k,b)+1._R8P)*q(beta_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(alfa_B,i,j,k,b) = 1/(2*EPS0**0.5_R8P)*((f(n,i,j,k,b)+1._R8P)*q(alfa_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P - s2*(f(n,i,j,k,b)-1._R8P)*q(beta_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(beta_B,i,j,k,b) = 1/(2*EPS0**0.5_R8P)*((f(n,i,j,k,b)+1._R8P)*q(beta_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + s2*(f(n,i,j,k,b)-1._R8P)*q(alfa_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(gamma_D,i,j,k,b) = q(gamma_D,i,j,k,b)
!
   !                  q(gamma_B,i,j,k,b) = q(gamma_B,i,j,k,b)
   !               enddo
   !            enddo
   !         enddo
   !      endif
   !      !z+ side
   !      if (layer(6)) then
   !         n = 3_I4P
   !         s2 = -1.0_R8P
   !         alfa_D = 1_I4P
   !         beta_D = 2_I4P
   !         gamma_D = 3_I4P
   !         alfa_B = 4_I4P
   !         beta_B = 5_I4P
   !         gamma_B = 6_I4P
   !         do i=1,ni
   !            do j=1, nj
   !               do k=nk-C, nk
   !                  q(alfa_D,i,j,k,b) = 1/(2*MU0**0.5_R8P)*(s2*(f(n,i,j,k,b)-1._R8P)*q(beta_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + (f(n,i,j,k,b)+1._R8P)*q(alfa_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(beta_D,i,j,k,b) = 1/(2*MU0**0.5_R8P)*(-s2*(f(n,i,j,k,b)-1._R8P)*q(alfa_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + (f(n,i,j,k,b)+1._R8P)*q(beta_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(alfa_B,i,j,k,b) = 1/(2*EPS0**0.5_R8P)*((f(n,i,j,k,b)+1._R8P)*q(alfa_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P - s2*(f(n,i,j,k,b)-1._R8P)*q(beta_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(beta_B,i,j,k,b) = 1/(2*EPS0**0.5_R8P)*((f(n,i,j,k,b)+1._R8P)*q(beta_B,i,j,k,b) &
   !                  *EPS0**0.5_R8P + s2*(f(n,i,j,k,b)-1._R8P)*q(alfa_D,i,j,k,b)*MU0**0.5_R8P)
!
   !                  q(gamma_D,i,j,k,b) = q(gamma_D,i,j,k,b)
!
   !                  q(gamma_B,i,j,k,b) = q(gamma_B,i,j,k,b)
   !               enddo
   !            enddo
   !         enddo
   !      endif
   !   enddo
   !endif
   !endassociate
   !endsubroutine apply_fWL_correction
