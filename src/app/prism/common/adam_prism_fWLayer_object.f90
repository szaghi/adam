!< ADAM, PRISM (Plasma Research usIng Simulation Methods) fWLayer class definition, common backend.
module adam_prism_fWLayer_object

    !in input prendi il numero di celle che compone lo strato, chiamalo C e 3 flag per definire su quali lati ho lo strato
    !assegni la funzione f ad ognuna delle celle: essa dovrà avere tre elementi in modo tale da poter 
    !tenere conto delle celle sulla diagonale

    !sicuramente devi dargli in pasto field o comunque un array di celle che contenga le coordinate delle celle stesse

    !ragiona su come scrivere (in cpu però, non qui) la funzione di aggiornamento dei campi. La puoi fare ricorsiva e senza if se 
    !vale a livello matematico, altrimenti dovrai aggiungere un flag(3) per individuare se la cella appartiene a uno o più lati dello strato e a quali (in ogni elemento + o -1)

! ADAM modules
use :: adam_mpih_object, only : mpih_object
use :: adam_field_object, only : field_object
! PRISM modules
use :: adam_prism_parameters
use :: adam_prism_physics_object, only : prism_physics_object
! third party modules
use :: finer
use :: penf

implicit none
private
public :: prism_fWLayer_object

character(len=7), parameter :: INI_SECTION_NAME='fWLayer' !< INI file section name containing flWLayer datas.

type :: prism_fWLayer_object
   !< PRISM fWLayer class definition.
   type(mpih_object)           :: mpih                  !< MPI handler.
   integer(I4P)                :: C        = 0_I4P      !< Layer cell width.
   logical                     :: layer(6) = .false.    !< Layer flags for each side (-x, +x, -y, +y, -z, +z).
   real(R8P), allocatable      :: f(:,:,:,:,:)          !< fWLayer function values.

contains
  ! public methods
  procedure, pass(self) :: description    !< Return pretty-printed object description.
  procedure, pass(self) :: initialize     !< Initialize physics.
  procedure, pass(self) :: load_from_file !< Load config from file. 
endtype prism_fWLayer_object  

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_fWLayer_object), intent(in) :: self             !< Physics.
   character(len=:), allocatable           :: desc             !< Description.
   character(len=1), parameter             :: NL=new_line('a') !< New line character.

   if (self%C == 0_I4P) then
      desc = self%mpih%myrankstr//'No fWLayer implemented'
   else
      desc = self%mpih%myrankstr//'fWLayer datas:'//NL
      desc = desc//self%mpih%myrankstr//'  Layer cell width: '//trim(str(self%C))//NL
      desc = desc//self%mpih%myrankstr//'  Layer on -x side: '//trim(str(self%layer(1)))//NL
      desc = desc//self%mpih%myrankstr//'  Layer on +x side: '//trim(str(self%layer(2)))//NL
      desc = desc//self%mpih%myrankstr//'  Layer on -y side: '//trim(str(self%layer(3)))//NL
      desc = desc//self%mpih%myrankstr//'  Layer on +y side: '//trim(str(self%layer(4)))//NL
      desc = desc//self%mpih%myrankstr//'  Layer on -z side: '//trim(str(self%layer(5)))//NL
      desc = desc//self%mpih%myrankstr//'  Layer on +z side: '//trim(str(self%layer(6)))
    endif
   endfunction description

   subroutine initialize(self, file_parameters, physics, field)
   !< Initialize the fWLayer.
   class(prism_fWLayer_object), intent(inout) :: self            !< fWLayer.
   type(file_ini),              intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   type(prism_physics_object),  intent(in)    :: physics         !< Physics.
   type(field_object),          intent(in)    :: field           !< Field.
   integer(I4P)                               :: i,j,k,b         !< Counters.
   real(R8P)                                  :: fi              !< Cell function
   real(R8P)                                  :: distance        !< Distance between cell and physical boundary
   real(R8P)                                  :: C_r             !< Number (real) of ghost cells in the layer
   real(R8P)                                  :: i_r, j_r, k_r   !< (Real) counters
   real(R8P)                                  :: ds              !< Cells distance in x, y or z.


   print '(A)', self%mpih%myrankstr//'prism_fWLayer_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   
   !Inizializzo funzione f nelle celle dello strato
   associate(ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, blocks_number=>field%blocks_number, &
            ngc=>field%grid%ngc, nb=>field%nb, dx=>field%dxyz(1,:), dy=>field%dxyz(2,:), dz=>field%dxyz(3,:), &
            C=>self%C)

   allocate(self%f(1:3,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:nb))
   self%f = 1.0_R8P
   C_r = real(C, R8P)
   
   if (C >0) then
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
                  do i=ni-C, ni
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
                  do j=nj-C, nj
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
                  do k=nk-C, nk
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
   print '(A)', self%mpih%myrankstr//'prism_fWLayer_object%initialize finish'
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
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(C)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='x_minus_layer', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(x_minus_layer)')
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
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(x_plus_layer)')
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
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(y_minus_layer)')
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
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(y_plus_layer)')
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
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(z_minus_layer)')
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
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(z_plus_layer)')
   select case(trim(adjustl(buff)))
   case('NO', 'no', 'No', 'nO')
      self%layer(6) = .false.
   case('YES', 'yes', 'Yes', 'yES')
      self%layer(6) = .true.
   case default
      self%layer(6) = .false.
   endselect
   endsubroutine load_from_file

endmodule adam_prism_fWLayer_object