!< ADAM, PRISM perfectly matched layer class definition, common backend.
module adam_prism_pml_object

! ADAM singleton objects
use :: adam_mpih_global,      only : mpih
use :: adam_grid_object,      only : grid_object
use :: adam_field_object,     only : field_object
use :: adam_tree_object,      only : tree_object
use :: adam_prism_absorbing_layer_geometry, only : compute_absorbing_face_range
! third party modules
use :: finer, only : file_ini
use :: penf,  only : I4P, R8P, str

implicit none
private

public :: prism_pml_object

character(len=3), parameter :: INI_SECTION_NAME  = 'PML'
character(len=8), parameter :: PML_TYPE_NONE     = 'NONE'
character(len=8), parameter :: PML_TYPE_CLASSIC  = 'CLASSIC'
character(len=8), parameter :: PML_TYPE_CFS      = 'CFS'
character(len=8), parameter :: PML_TYPE_BERMUDEZ = 'BERMUDEZ'
integer(I4P),     parameter :: PML_VARS_PER_FACE = 4_I4P
integer(I4P),     parameter :: PML_FACE_X_M      = 1_I4P
integer(I4P),     parameter :: PML_FACE_X_P      = 2_I4P
integer(I4P),     parameter :: PML_FACE_Y_M      = 3_I4P
integer(I4P),     parameter :: PML_FACE_Y_P      = 4_I4P
integer(I4P),     parameter :: PML_FACE_Z_M      = 5_I4P
integer(I4P),     parameter :: PML_FACE_Z_P      = 6_I4P
character(len=2), parameter :: FACE_LABEL(6)     = ['-x', '+x', '-y', '+y', '-z', '+z']
character(len=8), parameter :: PML_X_VAR_NAME(4) = ['psi_Ey_x', 'psi_Ez_x', 'psi_Hy_x', 'psi_Hz_x']
character(len=8), parameter :: PML_Y_VAR_NAME(4) = ['psi_Ex_y', 'psi_Ez_y', 'psi_Hx_y', 'psi_Hz_y']
character(len=8), parameter :: PML_Z_VAR_NAME(4) = ['psi_Ex_z', 'psi_Ey_z', 'psi_Hx_z', 'psi_Hy_z']

type :: prism_pml_object
   !< PRISM PML class definition.
   !< Each face-local state array stores four split variables in the order
   !< "electric first, magnetic second", e.g. q_pml_x_* =
   !< [psi_Ey_x, psi_Ez_x, psi_Hy_x, psi_Hz_x]. Storage is reduced to
   !< boundary blocks only and compacted face by face.
   character(len=99)         :: pml_type = PML_TYPE_NONE
   logical                   :: enabled  = .false.
   logical                   :: layer(6) = .false.
   real(R8P)                 :: width     = 0._R8P
   real(R8P)                 :: gamma_max = 0._R8P
   real(R8P)                 :: gamma_exponent = 0._R8P
   real(R8P)                 :: alpha_max = 0._R8P
   real(R8P)                 :: k_max     = 1._R8P
   real(R8P)                 :: beta      = 0._R8P
   real(R8P)                 :: profile_span(6) = 0._R8P !< Face-wise maximum center distance inside the PML.
   integer(I4P), allocatable :: ni_pml(:,:,:)  !< Local i-range of active PML support [2,nb,6].
   integer(I4P), allocatable :: nj_pml(:,:,:)  !< Local j-range of active PML support [2,nb,6].
   integer(I4P), allocatable :: nk_pml(:,:,:)  !< Local k-range of active PML support [2,nb,6].
   integer(I4P), allocatable :: block_lid(:,:) !< Global-block to compact face-local slot [nb,6]; 0 if inactive.
   integer(I4P)              :: max_cells(6)     = 0_I4P !< Max local PML thickness on each face.
   integer(I4P)              :: active_blocks(6) = 0_I4P !< Boundary blocks carrying PML on each face.
   integer(I4P), allocatable :: blocks_x_m(:) !< Global block ids carrying x-minus PML.
   integer(I4P), allocatable :: blocks_x_p(:) !< Global block ids carrying x-plus PML.
   integer(I4P), allocatable :: blocks_y_m(:) !< Global block ids carrying y-minus PML.
   integer(I4P), allocatable :: blocks_y_p(:) !< Global block ids carrying y-plus PML.
   integer(I4P), allocatable :: blocks_z_m(:) !< Global block ids carrying z-minus PML.
   integer(I4P), allocatable :: blocks_z_p(:) !< Global block ids carrying z-plus PML.
   real(R8P), allocatable    :: q_pml_x_m(:,:,:,:,:) !< [4,Cmax_xm,nj,nk,nblocks_xm]
   real(R8P), allocatable    :: q_pml_x_p(:,:,:,:,:) !< [4,Cmax_xp,nj,nk,nblocks_xp]
   real(R8P), allocatable    :: q_pml_y_m(:,:,:,:,:) !< [4,ni,Cmax_ym,nk,nblocks_ym]
   real(R8P), allocatable    :: q_pml_y_p(:,:,:,:,:) !< [4,ni,Cmax_yp,nk,nblocks_yp]
   real(R8P), allocatable    :: q_pml_z_m(:,:,:,:,:) !< [4,ni,nj,Cmax_zm,nblocks_zm]
   real(R8P), allocatable    :: q_pml_z_p(:,:,:,:,:) !< [4,ni,nj,Cmax_zp,nblocks_zp]
contains
   procedure, pass(self) :: description
   procedure, pass(self) :: initialize
   procedure, pass(self) :: load_from_file
endtype prism_pml_object

contains
   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_pml_object), intent(in) :: self
   character(len=:), allocatable       :: desc
   character(len=1), parameter         :: NL = new_line('a')
   integer(I4P)                        :: face

   if (.not. self%enabled) then
      desc = mpih%myrankstr//'   No PML implemented'
      return
   endif

   desc = mpih%myrankstr//'   PML data:'
   desc = desc//NL//mpih%myrankstr//'      PML type: '//trim(self%pml_type)
   desc = desc//NL//mpih%myrankstr//'      Layer physical width: '//trim(str(self%width))
   select case (trim(self%pml_type))
   case (PML_TYPE_CLASSIC)
      desc = desc//NL//mpih%myrankstr//'      gamma_max: '//trim(str(self%gamma_max))
      desc = desc//NL//mpih%myrankstr//'      gamma exponent: '//trim(str(self%gamma_exponent))
   case (PML_TYPE_CFS)
      desc = desc//NL//mpih%myrankstr//'      gamma_max: '//trim(str(self%gamma_max))
      desc = desc//NL//mpih%myrankstr//'      alfa_max:  '//trim(str(self%alpha_max))
      desc = desc//NL//mpih%myrankstr//'      k_max:     '//trim(str(self%k_max))
   case (PML_TYPE_BERMUDEZ)
      desc = desc//NL//mpih%myrankstr//'      beta:      '//trim(str(self%beta))
      desc = desc//NL//mpih%myrankstr//'      gamma exponent: '//trim(str(self%gamma_exponent))
   endselect

   do face=1, 6
      desc = desc//NL//mpih%myrankstr//'      Layer on '//FACE_LABEL(face)//' side: '//trim(str(self%layer(face)))
      if (allocated(self%block_lid)) then
         desc = desc//NL//mpih%myrankstr//'         active blocks: '//trim(str(self%active_blocks(face)))// &
                      ', max cells: '//trim(str(self%max_cells(face)))
      endif
   enddo
   endfunction description

   subroutine initialize(self, field, grid, tree, file_parameters)
   !< Initialize the PML support and allocate reduced face-local states.
   class(prism_pml_object), intent(inout) :: self
   type(field_object),      intent(in)    :: field
   type(grid_object),       intent(in)    :: grid
   type(tree_object),       intent(in)    :: tree
   type(file_ini),          intent(in)    :: file_parameters
   integer(I4P)                           :: alloc_error
   integer(I4P)                           :: b
   integer(I4P)                           :: face
   integer(I4P)                           :: face_cells
   integer(I4P)                           :: face_counter(6)
   integer(I4P)                           :: range_i(2), range_j(2), range_k(2)
   character(len=256)                     :: alloc_message
   real(R8P)                              :: face_last_center_distance
   real(R8P)                              :: face_profile_extent

   print '(A)', mpih%myrankstr//'prism_pml_object%initialize start'
   call reset_pml_object(self=self)
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()

   if (.not. self%enabled) then
      print '(A)', mpih%myrankstr//'prism_pml_object%initialize finish'
      return
   endif

   allocate(self%ni_pml(1:2,1:field%nb,1:6), stat=alloc_error, errmsg=alloc_message)
   if (alloc_error /= 0_I4P) &
      call mpih%error_stop(msg='prism_pml_object%initialize: failed to allocate ni_pml: '//trim(alloc_message))
   allocate(self%nj_pml(1:2,1:field%nb,1:6), stat=alloc_error, errmsg=alloc_message)
   if (alloc_error /= 0_I4P) &
      call mpih%error_stop(msg='prism_pml_object%initialize: failed to allocate nj_pml: '//trim(alloc_message))
   allocate(self%nk_pml(1:2,1:field%nb,1:6), stat=alloc_error, errmsg=alloc_message)
   if (alloc_error /= 0_I4P) &
      call mpih%error_stop(msg='prism_pml_object%initialize: failed to allocate nk_pml: '//trim(alloc_message))
   allocate(self%block_lid(1:field%nb,1:6), stat=alloc_error, errmsg=alloc_message)
   if (alloc_error /= 0_I4P) &
      call mpih%error_stop(msg='prism_pml_object%initialize: failed to allocate block_lid: '//trim(alloc_message))
   self%ni_pml       = 0_I4P
   self%nj_pml       = 0_I4P
   self%nk_pml       = 0_I4P
   self%block_lid    = 0_I4P
   self%max_cells    = 0_I4P
   self%active_blocks = 0_I4P
   self%profile_span = 0._R8P

   do b=1, field%blocks_number
      do face=1, 6
         if (.not. self%layer(face)) cycle
         call compute_absorbing_face_range(face=face, width=self%width, ni=grid%ni, nj=grid%nj, nk=grid%nk,          &
                                           domain_emin=grid%domain_emin, domain_emax=grid%domain_emax,               &
                                           block_emin=field%emin(:,b), block_emax=field%emax(:,b),                   &
                                           dxyz=field%dxyz(:,b), ni_range=range_i, nj_range=range_j, nk_range=range_k, &
                                           cells=face_cells, last_center_distance=face_last_center_distance,         &
                                           profile_extent=face_profile_extent)
         if (face_cells <= 0_I4P) cycle

         self%ni_pml(:,b,face) = range_i
         self%nj_pml(:,b,face) = range_j
         self%nk_pml(:,b,face) = range_k

         self%active_blocks(face) = self%active_blocks(face) + 1_I4P
         self%max_cells(face) = max(self%max_cells(face), face_cells)
         self%profile_span(face) = max(self%profile_span(face), face_last_center_distance)
      enddo
   enddo

   call allocate_face_metadata(self=self)

   face_counter = 0_I4P
   do b=1, field%blocks_number
      do face=1, 6
         if (.not. face_is_active(self=self, block_id=b, face=face)) cycle
         face_counter(face) = face_counter(face) + 1_I4P
         self%block_lid(b,face) = face_counter(face)
         select case (face)
         case (PML_FACE_X_M)
            self%blocks_x_m(face_counter(face)) = b
         case (PML_FACE_X_P)
            self%blocks_x_p(face_counter(face)) = b
         case (PML_FACE_Y_M)
            self%blocks_y_m(face_counter(face)) = b
         case (PML_FACE_Y_P)
            self%blocks_y_p(face_counter(face)) = b
         case (PML_FACE_Z_M)
            self%blocks_z_m(face_counter(face)) = b
         case (PML_FACE_Z_P)
            self%blocks_z_p(face_counter(face)) = b
         endselect
      enddo
   enddo

   call allocate_face_storage(self=self, grid=grid)
   print '(A)', mpih%myrankstr//'prism_pml_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load PML configuration from file.
   class(prism_pml_object), intent(inout)        :: self
   type(file_ini),          intent(in)           :: file_parameters
   logical,                 intent(in), optional :: go_on_fail
   logical                                      :: go_on_fail_
   integer(I4P)                                 :: error
   character(99)                                :: buff

   call reset_pml_configuration(self=self)
   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='PML_type', val=buff, error=error)
   if (error > 0_I4P) return

   select case (trim(adjustl(buff)))
   case ('NO', 'no', 'No', 'nO', 'NONE', 'none', 'None')
      self%pml_type = PML_TYPE_NONE
   case ('PML', 'pml', 'Pml', 'CLASSIC', 'classic', 'Classic', 'STANDARD', 'standard', 'Standard')
      self%pml_type = PML_TYPE_CLASSIC
   case ('CFS', 'cfs', 'Cfs', 'CFS_PML', 'cfs_pml', 'CPML', 'cpml')
      self%pml_type = PML_TYPE_CFS
   case ('BERMUDEZ', 'bermudez', 'Bermudez')
      self%pml_type = PML_TYPE_BERMUDEZ
   case default
      call mpih%print_message(msg='warning: PML type "'//trim(adjustl(buff))//'" unknown, PML disabled')
      self%pml_type = PML_TYPE_NONE
   endselect

   if (trim(self%pml_type) == PML_TYPE_NONE) return

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='width', val=self%width, error=error)
   if (.not. go_on_fail_ .and. error > 0_I4P) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(width)')
   if (self%width <= 0._R8P) call mpih%error_stop(msg=': invalid ['//INI_SECTION_NAME//'].(width), must be > 0')

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='x_minus_layer', val=self%layer(PML_FACE_X_M), error=error)
   if (.not. go_on_fail_ .and. error > 0_I4P) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(x_minus_layer)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='x_plus_layer', val=self%layer(PML_FACE_X_P), error=error)
   if (.not. go_on_fail_ .and. error > 0_I4P) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(x_plus_layer)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='y_minus_layer', val=self%layer(PML_FACE_Y_M), error=error)
   if (.not. go_on_fail_ .and. error > 0_I4P) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(y_minus_layer)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='y_plus_layer', val=self%layer(PML_FACE_Y_P), error=error)
   if (.not. go_on_fail_ .and. error > 0_I4P) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(y_plus_layer)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='z_minus_layer', val=self%layer(PML_FACE_Z_M), error=error)
   if (.not. go_on_fail_ .and. error > 0_I4P) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(z_minus_layer)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='z_plus_layer', val=self%layer(PML_FACE_Z_P), error=error)
   if (.not. go_on_fail_ .and. error > 0_I4P) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(z_plus_layer)')

   if (.not. any(self%layer)) then
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] active but all face flags are false; '// &
                               'enable at least one among x/y/z +/- layer')
   endif

   select case (trim(self%pml_type))
   case (PML_TYPE_CLASSIC)
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='gamma_max', val=self%gamma_max, error=error)
      if (.not. go_on_fail_ .and. error > 0_I4P) &
         call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(gamma_max)')
      if (self%gamma_max < 0._R8P) call mpih%error_stop(msg=': invalid ['//INI_SECTION_NAME//'].(gamma_max), must be >= 0')
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='gamma_exponent', val=self%gamma_exponent, error=error)
      if (.not. go_on_fail_ .and. error > 0_I4P) &
         call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(gamma_exponent)')
      if (self%gamma_exponent <= 0._R8P) &
         call mpih%error_stop(msg=': invalid ['//INI_SECTION_NAME//'].(gamma_exponent), must be > 0')
   case (PML_TYPE_CFS)
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='gamma_max', val=self%gamma_max, error=error)
      if (.not. go_on_fail_ .and. error > 0_I4P) &
         call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(gamma_max)')
      if (self%gamma_max < 0._R8P) call mpih%error_stop(msg=': invalid ['//INI_SECTION_NAME//'].(gamma_max), must be >= 0')

      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='alfa_max', val=self%alpha_max, error=error)
      if (.not. go_on_fail_ .and. error > 0_I4P) &
         call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(alfa_max)')
      if (self%alpha_max < 0._R8P) call mpih%error_stop(msg=': invalid ['//INI_SECTION_NAME//'].(alfa_max), must be >= 0')

      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='k_max', val=self%k_max, error=error)
      if (.not. go_on_fail_ .and. error > 0_I4P) &
         call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(k_max)')
      if (self%k_max <= 0._R8P) call mpih%error_stop(msg=': invalid ['//INI_SECTION_NAME//'].(k_max), must be > 0')
   case (PML_TYPE_BERMUDEZ)
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='beta', val=self%beta, error=error)
      if (.not. go_on_fail_ .and. error > 0_I4P) &
         call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(beta)')
      if (self%beta < 0._R8P) call mpih%error_stop(msg=': invalid ['//INI_SECTION_NAME//'].(beta), must be >= 0')
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='gamma_exponent', val=self%gamma_exponent, error=error)
      if (.not. go_on_fail_ .and. error > 0_I4P) &
         call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(gamma_exponent)')
      if (self%gamma_exponent <= 0._R8P) &
         call mpih%error_stop(msg=': invalid ['//INI_SECTION_NAME//'].(gamma_exponent), must be > 0')
   endselect

   self%enabled = .true.
   endsubroutine load_from_file

   subroutine allocate_face_metadata(self)
   !< Allocate per-face block lists only where the local rank owns boundary blocks.
   class(prism_pml_object), intent(inout) :: self
   integer(I4P)                           :: alloc_error
   character(len=256)                     :: alloc_message

   if (self%active_blocks(PML_FACE_X_M) > 0_I4P) then
      allocate(self%blocks_x_m(1:self%active_blocks(PML_FACE_X_M)), stat=alloc_error, errmsg=alloc_message)
      if (alloc_error /= 0_I4P) &
         call mpih%error_stop(msg='prism_pml_object%allocate_face_metadata: failed to allocate blocks_x_m: '// &
                                  trim(alloc_message))
      self%blocks_x_m = 0_I4P
   endif
   if (self%active_blocks(PML_FACE_X_P) > 0_I4P) then
      allocate(self%blocks_x_p(1:self%active_blocks(PML_FACE_X_P)), stat=alloc_error, errmsg=alloc_message)
      if (alloc_error /= 0_I4P) &
         call mpih%error_stop(msg='prism_pml_object%allocate_face_metadata: failed to allocate blocks_x_p: '// &
                                  trim(alloc_message))
      self%blocks_x_p = 0_I4P
   endif
   if (self%active_blocks(PML_FACE_Y_M) > 0_I4P) then
      allocate(self%blocks_y_m(1:self%active_blocks(PML_FACE_Y_M)), stat=alloc_error, errmsg=alloc_message)
      if (alloc_error /= 0_I4P) &
         call mpih%error_stop(msg='prism_pml_object%allocate_face_metadata: failed to allocate blocks_y_m: '// &
                                  trim(alloc_message))
      self%blocks_y_m = 0_I4P
   endif
   if (self%active_blocks(PML_FACE_Y_P) > 0_I4P) then
      allocate(self%blocks_y_p(1:self%active_blocks(PML_FACE_Y_P)), stat=alloc_error, errmsg=alloc_message)
      if (alloc_error /= 0_I4P) &
         call mpih%error_stop(msg='prism_pml_object%allocate_face_metadata: failed to allocate blocks_y_p: '// &
                                  trim(alloc_message))
      self%blocks_y_p = 0_I4P
   endif
   if (self%active_blocks(PML_FACE_Z_M) > 0_I4P) then
      allocate(self%blocks_z_m(1:self%active_blocks(PML_FACE_Z_M)), stat=alloc_error, errmsg=alloc_message)
      if (alloc_error /= 0_I4P) &
         call mpih%error_stop(msg='prism_pml_object%allocate_face_metadata: failed to allocate blocks_z_m: '// &
                                  trim(alloc_message))
      self%blocks_z_m = 0_I4P
   endif
   if (self%active_blocks(PML_FACE_Z_P) > 0_I4P) then
      allocate(self%blocks_z_p(1:self%active_blocks(PML_FACE_Z_P)), stat=alloc_error, errmsg=alloc_message)
      if (alloc_error /= 0_I4P) &
         call mpih%error_stop(msg='prism_pml_object%allocate_face_metadata: failed to allocate blocks_z_p: '// &
                                  trim(alloc_message))
      self%blocks_z_p = 0_I4P
   endif
   endsubroutine allocate_face_metadata

   subroutine allocate_face_storage(self, grid)
   !< Allocate reduced face-local PML state arrays.
   class(prism_pml_object), intent(inout) :: self
   type(grid_object),       intent(in)    :: grid
   integer(I4P)                           :: alloc_error
   character(len=256)                     :: alloc_message

   if (self%active_blocks(PML_FACE_X_M) > 0_I4P) then
      allocate(self%q_pml_x_m(1:PML_VARS_PER_FACE, 1:self%max_cells(PML_FACE_X_M), 1:grid%nj, 1:grid%nk, &
                              1:self%active_blocks(PML_FACE_X_M)), stat=alloc_error, errmsg=alloc_message)
      if (alloc_error /= 0_I4P) &
         call mpih%error_stop(msg='prism_pml_object%allocate_face_storage: failed to allocate q_pml_x_m: '// &
                                  trim(alloc_message))
      self%q_pml_x_m = 0._R8P
   endif
   if (self%active_blocks(PML_FACE_X_P) > 0_I4P) then
      allocate(self%q_pml_x_p(1:PML_VARS_PER_FACE, 1:self%max_cells(PML_FACE_X_P), 1:grid%nj, 1:grid%nk, &
                              1:self%active_blocks(PML_FACE_X_P)), stat=alloc_error, errmsg=alloc_message)
      if (alloc_error /= 0_I4P) &
         call mpih%error_stop(msg='prism_pml_object%allocate_face_storage: failed to allocate q_pml_x_p: '// &
                                  trim(alloc_message))
      self%q_pml_x_p = 0._R8P
   endif
   if (self%active_blocks(PML_FACE_Y_M) > 0_I4P) then
      allocate(self%q_pml_y_m(1:PML_VARS_PER_FACE, 1:grid%ni, 1:self%max_cells(PML_FACE_Y_M), 1:grid%nk, &
                              1:self%active_blocks(PML_FACE_Y_M)), stat=alloc_error, errmsg=alloc_message)
      if (alloc_error /= 0_I4P) &
         call mpih%error_stop(msg='prism_pml_object%allocate_face_storage: failed to allocate q_pml_y_m: '// &
                                  trim(alloc_message))
      self%q_pml_y_m = 0._R8P
   endif
   if (self%active_blocks(PML_FACE_Y_P) > 0_I4P) then
      allocate(self%q_pml_y_p(1:PML_VARS_PER_FACE, 1:grid%ni, 1:self%max_cells(PML_FACE_Y_P), 1:grid%nk, &
                              1:self%active_blocks(PML_FACE_Y_P)), stat=alloc_error, errmsg=alloc_message)
      if (alloc_error /= 0_I4P) &
         call mpih%error_stop(msg='prism_pml_object%allocate_face_storage: failed to allocate q_pml_y_p: '// &
                                  trim(alloc_message))
      self%q_pml_y_p = 0._R8P
   endif
   if (self%active_blocks(PML_FACE_Z_M) > 0_I4P) then
      allocate(self%q_pml_z_m(1:PML_VARS_PER_FACE, 1:grid%ni, 1:grid%nj, 1:self%max_cells(PML_FACE_Z_M), &
                              1:self%active_blocks(PML_FACE_Z_M)), stat=alloc_error, errmsg=alloc_message)
      if (alloc_error /= 0_I4P) &
         call mpih%error_stop(msg='prism_pml_object%allocate_face_storage: failed to allocate q_pml_z_m: '// &
                                  trim(alloc_message))
      self%q_pml_z_m = 0._R8P
   endif
   if (self%active_blocks(PML_FACE_Z_P) > 0_I4P) then
      allocate(self%q_pml_z_p(1:PML_VARS_PER_FACE, 1:grid%ni, 1:grid%nj, 1:self%max_cells(PML_FACE_Z_P), &
                              1:self%active_blocks(PML_FACE_Z_P)), stat=alloc_error, errmsg=alloc_message)
      if (alloc_error /= 0_I4P) &
         call mpih%error_stop(msg='prism_pml_object%allocate_face_storage: failed to allocate q_pml_z_p: '// &
                                  trim(alloc_message))
      self%q_pml_z_p = 0._R8P
   endif
   endsubroutine allocate_face_storage

   subroutine reset_pml_configuration(self)
   !< Reset the input-driven configuration, leaving storage handling to reset_pml_object.
   class(prism_pml_object), intent(inout) :: self

   self%pml_type  = PML_TYPE_NONE
   self%enabled   = .false.
   self%layer     = .false.
   self%width     = 0._R8P
   self%gamma_max = 0._R8P
   self%gamma_exponent = 0._R8P
   self%alpha_max = 0._R8P
   self%k_max     = 1._R8P
   self%beta      = 0._R8P
   self%profile_span = 0._R8P
   endsubroutine reset_pml_configuration

   subroutine reset_pml_object(self)
   !< Reset configuration and release all derived storage.
   class(prism_pml_object), intent(inout) :: self

   call reset_pml_configuration(self=self)
   self%max_cells     = 0_I4P
   self%active_blocks = 0_I4P
   if (allocated(self%ni_pml))     deallocate(self%ni_pml)
   if (allocated(self%nj_pml))     deallocate(self%nj_pml)
   if (allocated(self%nk_pml))     deallocate(self%nk_pml)
   if (allocated(self%block_lid))  deallocate(self%block_lid)
   if (allocated(self%blocks_x_m)) deallocate(self%blocks_x_m)
   if (allocated(self%blocks_x_p)) deallocate(self%blocks_x_p)
   if (allocated(self%blocks_y_m)) deallocate(self%blocks_y_m)
   if (allocated(self%blocks_y_p)) deallocate(self%blocks_y_p)
   if (allocated(self%blocks_z_m)) deallocate(self%blocks_z_m)
   if (allocated(self%blocks_z_p)) deallocate(self%blocks_z_p)
   if (allocated(self%q_pml_x_m))  deallocate(self%q_pml_x_m)
   if (allocated(self%q_pml_x_p))  deallocate(self%q_pml_x_p)
   if (allocated(self%q_pml_y_m))  deallocate(self%q_pml_y_m)
   if (allocated(self%q_pml_y_p))  deallocate(self%q_pml_y_p)
   if (allocated(self%q_pml_z_m))  deallocate(self%q_pml_z_m)
   if (allocated(self%q_pml_z_p))  deallocate(self%q_pml_z_p)
   endsubroutine reset_pml_object

   pure logical function face_is_active(self, block_id, face) result(is_active)
   !< Return true if the given block carries PML support on the selected face.
   class(prism_pml_object), intent(in) :: self
   integer(I4P),            intent(in) :: block_id
   integer(I4P),            intent(in) :: face

   is_active = self%ni_pml(1,block_id,face) > 0_I4P .or. self%nj_pml(1,block_id,face) > 0_I4P .or. &
               self%nk_pml(1,block_id,face) > 0_I4P
   endfunction face_is_active

endmodule adam_prism_pml_object
