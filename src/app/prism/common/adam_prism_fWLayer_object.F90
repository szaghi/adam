!< ADAM, PRISM (Plasma Research usIng Simulation Methods) fWLayer class definition, common backend.
module adam_prism_fWLayer_object

    ! The fWLayer is configured through a physical width and the active
    ! boundary faces; the corresponding cell count is derived block by block
    ! from field%dxyz so the same setup remains meaningful across non-uniform
    ! grids and AMR levels.

! ADAM singleton objects
use :: adam_mpih_global,       only : mpih
use :: adam_grid_object,       only : grid_object
use :: adam_field_object,      only : field_object
use :: adam_tree_object,       only : tree_object, NODE_BOUNDARY_CONDITION
use :: adam_tree_node_object,  only : tree_node_object
! PRISM modules
use :: adam_prism_parameters
use :: adam_prism_physics_object, only : prism_physics_object
! third party modules
use :: finer, only : file_ini
use :: penf,  only : I4P, R8P, str

implicit none
private
public :: prism_fWLayer_object
public :: apply_fWL_correction_fun
public :: compute_fwl_factor

character(len=7), parameter :: INI_SECTION_NAME='fWLayer' !< INI file section name containing flWLayer datas.

type :: prism_fWLayer_object
   !< PRISM fWLayer class definition.
   logical                   :: layer(6) = .false.                       !< Layer flags for each side (-x, +x, -y, +y, -z, +z).
   real(R8P)                 :: width    = 0._R8P                        !< Requested physical layer width.
   integer(I4P), allocatable :: C(:,:)                                   !< Derived layer width in cells for each block/face [nb,6].
   integer(I4P), allocatable :: ni_fWL(:,:,:), nj_fWL(:,:,:), nk_fWL(:,:,:) !< FWL bounds for each block/face [2,nb,6].
   real(R8P)                 :: s2(6)                                    !< Side coefficient.
   integer(I4P)              :: n(6)                                     !< FWL f function index.
   integer(I4P)              :: alfa_D(6), beta_D(6)                     !< Corrected var index of D (Barbas' notation).
   integer(I4P)              :: alfa_B(6), beta_B(6)                     !< Corrected var index of B (Barbas' notation).
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

   if (self%width <= 0._R8P .or. .not. any(self%layer)) then
      desc = mpih%myrankstr//'   No fWLayer implemented'
   else
      desc = mpih%myrankstr//'   fWLayer datas:'
      desc = desc//NL//mpih%myrankstr//'      Layer physical width: '//trim(str(self%width))
      desc = desc//NL//mpih%myrankstr//'      Layer on -x side: '//trim(str(self%layer(1)))
      desc = desc//NL//mpih%myrankstr//'      Layer on +x side: '//trim(str(self%layer(2)))
      desc = desc//NL//mpih%myrankstr//'      Layer on -y side: '//trim(str(self%layer(3)))
      desc = desc//NL//mpih%myrankstr//'      Layer on +y side: '//trim(str(self%layer(4)))
      desc = desc//NL//mpih%myrankstr//'      Layer on -z side: '//trim(str(self%layer(5)))
      desc = desc//NL//mpih%myrankstr//'      Layer on +z side: '//trim(str(self%layer(6)))
    endif
   endfunction description

   subroutine initialize(self, field, grid, tree, file_parameters, physics)
   !< Initialize the fWLayer.
   class(prism_fWLayer_object), intent(inout) :: self            !< fWLayer.
   type(field_object), intent(in)             :: field           !< Field (sibling realm component, threaded in).
   type(grid_object),           intent(in)    :: grid            !< Grid (sibling realm component, threaded in).
   type(tree_object),           intent(in)    :: tree            !< Tree (sibling realm component, threaded in).
   type(file_ini),              intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   type(prism_physics_object),  intent(in)    :: physics         !< Physics.
   type(tree_node_object),      pointer       :: node_ptr        !< Pointer to current node.
   integer(I4P)                               :: b,fec           !< Counters.
   integer(I4P)                               :: neighbor_type   !< Neighbors type.
   integer(I4P)                               :: C_face          !< Number of cells covering the requested width on this face/block.
   real(R8P)                                  :: ds              !< Cells distance in x, y or z.
   integer(I4P)                               :: alloc_error     !< Allocation status.
   character(len=256)                         :: alloc_message   !< Allocation error message.

   print '(A)', mpih%myrankstr//'prism_fWLayer_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   if (allocated(self%C     )) deallocate(self%C)
   if (allocated(self%ni_fWL)) deallocate(self%ni_fWL)
   if (allocated(self%nj_fWL)) deallocate(self%nj_fWL)
   if (allocated(self%nk_fWL)) deallocate(self%nk_fWL)
   print '(A)', self%description()
   
   if (self%width <= 0._R8P .or. .not. any(self%layer)) return

   allocate(self%C(1:field%nb,1:6), stat=alloc_error, errmsg=alloc_message)
   if (alloc_error /= 0) call mpih%error_stop(msg='prism_fWLayer_object%initialize: failed to allocate C: '//trim(alloc_message))
   allocate(self%ni_fWL(1:2,1:field%nb,1:6), stat=alloc_error, errmsg=alloc_message)
   if (alloc_error /= 0) &
      call mpih%error_stop(msg='prism_fWLayer_object%initialize: failed to allocate ni_fWL: '//trim(alloc_message))
   allocate(self%nj_fWL(1:2,1:field%nb,1:6), stat=alloc_error, errmsg=alloc_message)
   if (alloc_error /= 0) &
      call mpih%error_stop(msg='prism_fWLayer_object%initialize: failed to allocate nj_fWL: '//trim(alloc_message))
   allocate(self%nk_fWL(1:2,1:field%nb,1:6), stat=alloc_error, errmsg=alloc_message)
   if (alloc_error /= 0) &
      call mpih%error_stop(msg='prism_fWLayer_object%initialize: failed to allocate nk_fWL: '//trim(alloc_message))
   self%C      = 0_I4P
   self%ni_fWL = 0_I4P
   self%nj_fWL = 0_I4P
   self%nk_fWL = 0_I4P

   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, blocks_number=>field%blocks_number,          &
            dx=>field%dxyz(1,:), dy=>field%dxyz(2,:), dz=>field%dxyz(3,:),                       &
            C=>self%C, ni_fWL=>self%ni_fWL, nj_fWL=>self%nj_fWL, nk_fWL=>self%nk_fWL, &
            n=>self%n, s2=>self%s2, alfa_D=>self%alfa_D, alfa_B=>self%alfa_B,         &
            beta_D=>self%beta_D, beta_B=>self%beta_B)
   n(1) =1_I4P; n(2)=1_I4P; n(3)=2_I4P
   n(4) =2_I4P; n(5)=3_I4P; n(6)=3_I4P

   s2(1)= 1.0_R8P; s2(2)=-1.0_R8P; s2(3)= 1.0_R8P
   s2(4)=-1.0_R8P; s2(5)= 1.0_R8P; s2(6)=-1.0_R8P

   alfa_D(1)=2_I4P; alfa_D(2)=2_I4P; alfa_D(3)=3_I4P
   alfa_D(4)=3_I4P; alfa_D(5)=1_I4P; alfa_D(6)=1_I4P

   beta_D(1)=3_I4P; beta_D(2)=3_I4P; beta_D(3)=1_I4P
   beta_D(4)=1_I4P; beta_D(5)=2_I4P; beta_D(6)=2_I4P

   alfa_B(1)=5_I4P; alfa_B(2)=5_I4P; alfa_B(3)=6_I4P
   alfa_B(4)=6_I4P; alfa_B(5)=4_I4P; alfa_B(6)=4_I4P

   beta_B(1)=6_I4P; beta_B(2)=6_I4P; beta_B(3)=4_I4P
   beta_B(4)=4_I4P; beta_B(5)=5_I4P; beta_B(6)=5_I4P

   do b=1, blocks_number
      node_ptr => tree%node(code=field%code(b))
      do fec=1, 6
         if (.not. self%layer(fec)) cycle
         neighbor_type = node_ptr%neighbor(fec)%ntype
         if (neighbor_type /= NODE_BOUNDARY_CONDITION) cycle

         select case(fec)
         case(1, 2)
            ds = dx(b)
            C_face = min(ni, ceiling(self%width / ds, kind=I4P))
            ni_fWL(1,b,fec) = merge(1_I4P, ni-C_face+1_I4P, fec == 1)
            ni_fWL(2,b,fec) = merge(C_face, ni, fec == 1)
            nj_fWL(1,b,fec) = 1_I4P
            nj_fWL(2,b,fec) = nj
            nk_fWL(1,b,fec) = 1_I4P
            nk_fWL(2,b,fec) = nk
         case(3, 4)
            ds = dy(b)
            C_face = min(nj, ceiling(self%width / ds, kind=I4P))
            ni_fWL(1,b,fec) = 1_I4P
            ni_fWL(2,b,fec) = ni
            nj_fWL(1,b,fec) = merge(1_I4P, nj-C_face+1_I4P, fec == 3)
            nj_fWL(2,b,fec) = merge(C_face, nj, fec == 3)
            nk_fWL(1,b,fec) = 1_I4P
            nk_fWL(2,b,fec) = nk
         case(5, 6)
            ds = dz(b)
            C_face = min(nk, ceiling(self%width / ds, kind=I4P))
            ni_fWL(1,b,fec) = 1_I4P
            ni_fWL(2,b,fec) = ni
            nj_fWL(1,b,fec) = 1_I4P
            nj_fWL(2,b,fec) = nj
            nk_fWL(1,b,fec) = merge(1_I4P, nk-C_face+1_I4P, fec == 5)
            nk_fWL(2,b,fec) = merge(C_face, nk, fec == 5)
         endselect

         C(b,fec) = C_face
      enddo
   enddo
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

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='width', val=self%width, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(width)')
   if (self%width < 0._R8P) call mpih%error_stop(msg=': invalid ['//INI_SECTION_NAME//'].(width), must be >= 0')
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

   pure real(R8P) function compute_fwl_factor(offset, cells_number, ds) result(f_value)
   !< Return the local fWLayer damping factor for a cell at integer offset from the boundary.
   integer(I4P), intent(in) :: offset       !< Cell offset from the boundary, 0 at the boundary cell.
   integer(I4P), intent(in) :: cells_number !< Layer thickness in cells.
   real(R8P),    intent(in) :: ds           !< Cell size along the layer-normal direction.
   real(R8P)                :: C_r          !< Layer thickness in reals.
   real(R8P)                :: fi           !< Profile-shape coefficient.
   real(R8P)                :: distance     !< Physical distance from the boundary.
   !$acc routine seq

   if (cells_number <= 0_I4P) then
      f_value = 1._R8P
      return
   endif

   C_r = real(cells_number, R8P)
   if (cells_number < 40_I4P) then
      fi = 1._R8P / 150._R8P * (-7._R8P * C_r**2 + 255._R8P * C_r + 250._R8P)
   else
      fi = 25._R8P
   endif
   distance = real(offset, R8P) * ds
   f_value = 1._R8P / fi * log10(distance / (C_r * ds) * (10._R8P**fi - 1._R8P) + 1._R8P)
   endfunction compute_fwl_factor

   subroutine apply_fWL_correction_fun(blocks_number,ngc,C,ni_fWL,nj_fWL,nk_fWL,n,s2,alfa_D,beta_D,alfa_B,beta_B,dxyz,q)
   !< Applay FWL correction, direction agnostic.
   integer(I4P), intent(in)    :: blocks_number                     !< Blocks number.
   integer(I4P), intent(in)    :: ngc                               !< Number of ghost cells.
   integer(I4P), intent(in)    :: C(1:)                             !< Layer width in cells for each block on the current face.
   integer(I4P), intent(in)    :: ni_fWL(1:,1:), nj_fWL(1:,1:), nk_fWL(1:,1:) !< FWL bounds for each block on the current face.
   integer(I4P), intent(in)    :: n                                 !< f component.
   real(R8P),    intent(in)    :: s2                                !< Side coefficient.
   integer(I4P), intent(in)    :: alfa_D, beta_D                    !< Corrected var index of D (Barbas' notation).
   integer(I4P), intent(in)    :: alfa_B, beta_B                    !< Corrected var index of D (Barbas' notation).
   real(R8P),    intent(in)    :: dxyz(1:,1:)                       !< Block mesh spacing [3,nb].
   real(R8P),    intent(inout) :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Field variables.
   real(R8P)                   :: f_value                           !< Local fWLayer factor.
   real(R8P)                   :: fm1, fp1                          !< fWLayer function values in -+ cell.
   real(R8P)                   :: D_alfa, D_beta                    !< components of tangential fields before correction
   real(R8P)                   :: B_alfa, B_beta                    !< components of tangential fields before correction
   integer(I4P)                :: b,i,j,k                           !< Counter.
   integer(I4P)                :: offset                            !< Cell offset from the active boundary.

   do b=1,blocks_number
      if (C(b) <= 0_I4P) cycle
      do k=nk_fWL(1,b), nk_fWL(2,b)
         do j=nj_fWL(1,b), nj_fWL(2,b)
            do i=ni_fWL(1,b), ni_fWL(2,b)
               select case(n)
               case(1)
                  if (s2 > 0._R8P) then
                     offset = i - ni_fWL(1,b)
                  else
                     offset = ni_fWL(2,b) - i
                  endif
               case(2)
                  if (s2 > 0._R8P) then
                     offset = j - nj_fWL(1,b)
                  else
                     offset = nj_fWL(2,b) - j
                  endif
               case default
                  if (s2 > 0._R8P) then
                     offset = k - nk_fWL(1,b)
                  else
                     offset = nk_fWL(2,b) - k
                  endif
               endselect
               f_value = compute_fwl_factor(offset=offset, cells_number=C(b), ds=dxyz(n,b))
               fm1 = f_value - 1._R8P
               fp1 = f_value + 1._R8P
               D_alfa = q(alfa_D,i,j,k,b)
               D_beta = q(beta_D,i,j,k,b)
               B_alfa = q(alfa_B,i,j,k,b)
               B_beta = q(beta_B,i,j,k,b)
               q(alfa_D,i,j,k,b) = MU0_SQ_I2  * ( s2*fm1*B_beta*EPS0_SQ +    fp1*D_alfa*MU0_SQ)
               q(beta_D,i,j,k,b) = MU0_SQ_I2  * (-s2*fm1*B_alfa*EPS0_SQ +    fp1*D_beta*MU0_SQ)
               q(alfa_B,i,j,k,b) = EPS0_SQ_I2 * (    fp1*B_alfa*EPS0_SQ - s2*fm1*D_beta*MU0_SQ)
               q(beta_B,i,j,k,b) = EPS0_SQ_I2 * (    fp1*B_beta*EPS0_SQ + s2*fm1*D_alfa*MU0_SQ)
            enddo
         enddo
      enddo
   enddo
   endsubroutine apply_fWL_correction_fun
endmodule adam_prism_fWLayer_object
