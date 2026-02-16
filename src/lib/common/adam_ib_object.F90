!< ADAM, IB class definition, CPU backend.
module adam_ib_object
!< ADAM, IB class definition, CPU backend.

use adam_field_object, only : field_object
use adam_grid_object, only : grid_object
use adam_mpih_object, only : mpih_object
use finer
use penf

implicit none
private
public :: ib_object
public :: BCS_VISCOUS
public :: BCS_EULER
public :: IB_ANALYTICAL_SPHERE

character(len=5), parameter :: INI_SECTION_NAME="solid" !< INI (config) file section name containing IB configs.

character(len=20), parameter :: IB_ANALYTICAL_SPHERE   = 'analytical_sphere   ' !< Analytical sphere solid.
character(len=20), parameter :: IB_ANALYTICAL_CIRCLE   = 'analytical_circle   ' !< Analytical circle solid.
character(len=20), parameter :: IB_ANALYTICAL_RECTANGLE= 'analytical_rectangle' !< Analytical rectangle solid.
character(len=20), parameter :: IB_DEFINITIONS(4)=[IB_ANALYTICAL_SPHERE,    &
                                                   IB_ANALYTICAL_CIRCLE,    &
                                                   IB_ANALYTICAL_RECTANGLE, &
                                                   'file.off            ']      !< Available solid definitions.

integer(I4P), parameter :: BCS_VISCOUS = 1_I4P !< Visous wall.
integer(I4P), parameter :: BCS_EULER   = 2_I4P !< Inviscid wall.

type :: analytical_sphere_object
   !< Analytical sphere (or circle) solid class.
   real(R8P)    :: center(3) = [0._R8P,0._R8P,0._R8P] !< Sphere center.
   real(R8P)    :: radius    =  0._R8P                !< Sphere radius.
   character(1) :: axis      =  ' '                   !< Axis (x,y,z) normal in case of circle solid.
endtype analytical_sphere_object

type :: analytical_rectangle_object
   !< Analytical rectangle solid class.
   real(R8P)    :: center(3) = [0._R8P,0._R8P,0._R8P] !< Sphere center.
   real(R8P)    :: edge(2)   =  0._R8P                !< Major/minor edge length.
   character(1) :: axis      =  ' '                   !< Axis (x,y,z) normal.
endtype analytical_rectangle_object

type :: ib_object
   !< IB class definition, CPU backend.
   type(mpih_object)                              :: mpih            !< MPI handler.
   integer(I4P)                                   :: solids_number=0 !< Number of solids (only 1 supported now).
   character(99),                     allocatable :: s_name(:)       !< Solid name.
   integer(I4P),                      allocatable :: bc_type(:)      !< Boundary condition type.
   real(R8P),                         allocatable :: q(:, :)         !< Variables array for boundary conditions.
   character(99),                     allocatable :: definition(:)   !< (Type of) Solid definition.
   type(analytical_sphere_object),    allocatable :: sphere(:)       !< Analytical sphere/circle solid.
   type(analytical_rectangle_object), allocatable :: rectangle(:)    !< Analytical rectangle solid.
   integer(I4P)                                   :: n_eikonal=2     !< Number of eikonal integration steps.
   ! Pointers to ADAM data for easy handling.
   type(field_object), pointer :: field=>null() !< The field.
   type(grid_object),  pointer :: grid =>null() !< The grid.
   ! Large arrays.
   real(R8P), allocatable ::  phi(:,:,:,:,:) !< IB distance function.
   contains
      ! public methods
      procedure, pass(self) :: compute_phi            !< Compute distance function.
      procedure, pass(self) :: compute_phi_all_solids !< Compute last phi index, all solids summary.
      procedure, pass(self) :: description            !< Return pretty-printed object description.
      procedure, pass(self) :: evolve_eikonal         !< Evolve eikonal equation.
      procedure, pass(self) :: initialize             !< Initialize IB.
      procedure, pass(self) :: invert_eikonal         !< Invert eikonal equation over q inside IB.
      procedure, pass(self) :: load_from_file         !< Load config from file.
      procedure, pass(self) :: move_phi               !< Move phi and the actual ptree representation.
      procedure, pass(self) :: sphere_to_array        !< Convert analytical sphere class data to array data.
      ! private methods
      procedure, pass(self), private :: compute_phi_analytical_sphere    !< Compute distance for analytical sphere solids.
      procedure, pass(self), private :: compute_phi_analytical_circle    !< Compute distance for analytical circle solids.
      procedure, pass(self), private :: compute_phi_analytical_rectangle !< Compute distance for analytical rectangle solids.
endtype ib_object
contains
   ! public methods
   subroutine compute_phi(self, verbose)
   !< Compute phi, distance from IB solid.
   class(ib_object), intent(inout)        :: self     !< IB.
   logical,          intent(in), optional :: verbose  !< Flag to trigger verbose prints.
   logical                                :: verbose_ !< Flag to trigger verbose prints, local variable.
   integer(I4P)                           :: ib       !< Counter.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (self%solids_number > 0) then
      if (verbose_) call self%mpih%print_message('ib_object%compute_phi start')
      do ib=1, self%solids_number
         select case(trim(adjustl(self%definition(ib))))
         case(trim(IB_ANALYTICAL_SPHERE))
            call self%compute_phi_analytical_sphere(solid=ib, sphere=self%sphere(ib))
         case(trim(IB_ANALYTICAL_CIRCLE))
            call self%compute_phi_analytical_circle(solid=ib, sphere=self%sphere(ib))
         case(trim(IB_ANALYTICAL_RECTANGLE))
            call self%compute_phi_analytical_rectangle(solid=ib, rectangle=self%rectangle(ib))
         endselect
      enddo
      call self%compute_phi_all_solids(verbose=verbose)
      if (verbose_) call self%mpih%print_message('ib_object%compute_phi finish')
   endif
   endsubroutine compute_phi

   subroutine compute_phi_all_solids(self, verbose)
   !< Compute phi, distance from IB solid.
   class(ib_object), intent(inout)        :: self          !< IB.
   logical,          intent(in), optional :: verbose       !< Flag to trigger verbose prints.
   logical                                :: verbose_      !< Flag to trigger verbose prints, local variable.
   integer(I4P)                           :: b, i, j, k    !< Counter.
   integer(I4P)                           :: all_solids    !< Last phi index, all solids summary.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (self%solids_number > 0) then
      if (verbose_) call self%mpih%print_message('ib_object%compute_phi_all_solids start')
      all_solids = ubound(self%phi, dim=1)
      do b=1, self%field%blocks_number
         do k=1, self%grid%nk
            do j=1, self%grid%nj
               do i=1, self%grid%ni
                  self%phi(all_solids,i,j,k,b) = maxval(self%phi(1:all_solids-1,i,j,k,b))
               enddo
            enddo
         enddo
      enddo
      if (verbose_) call self%mpih%print_message('ib_object%compute_phi_all_solids finish')
   endif
   endsubroutine compute_phi_all_solids

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(ib_object), intent(in)  :: self             !< IB.
   character(len=:), allocatable :: desc             !< Description.
   character(len=1), parameter   :: NL=new_line('a') !< New line character.
   integer(I4P)                  :: s                !< Counter.

   desc =       self%mpih%myrankstr//'IB main data'//NL
   desc = desc//self%mpih%myrankstr//'  solids number: '//trim(str(self%solids_number))
   desc = desc//self%mpih%myrankstr//'  n_eikonal:     '//trim(str(self%n_eikonal))
   do s=1, self%solids_number
      desc = desc//NL//self%mpih%myrankstr//'  Solid '//trim(str(s,.true.))
      desc = desc//NL//self%mpih%myrankstr//'    BC type:    '//trim(str(self%bc_type(s)))
      desc = desc//NL//self%mpih%myrankstr//'    definition: '//trim(self%definition(s))
      select case(trim(adjustl(self%definition(s))))
      case(trim(IB_ANALYTICAL_SPHERE))
      desc = desc//NL//self%mpih%myrankstr//'    sphere_center_x: '//trim(str(self%sphere(s)%center(1)))
      desc = desc//NL//self%mpih%myrankstr//'    sphere_center_y: '//trim(str(self%sphere(s)%center(2)))
      desc = desc//NL//self%mpih%myrankstr//'    sphere_center_z: '//trim(str(self%sphere(s)%center(3)))
      desc = desc//NL//self%mpih%myrankstr//'    sphere_radius  : '//trim(str(self%sphere(s)%radius   ))
      case(trim(IB_ANALYTICAL_CIRCLE))
      case(trim(IB_ANALYTICAL_RECTANGLE))
      endselect
   enddo
   endfunction description

   subroutine initialize(self, grid, field, file_parameters)
   !< Initialize the equation.
   class(ib_object),   intent(inout)      :: self            !< IB.
   type(grid_object),  intent(in), target :: grid            !< The grid.
   type(field_object), intent(in), target :: field           !< The field.
   type(file_ini),     intent(inout)      :: file_parameters !< INI file handler.

   call self%mpih%initialize
   call self%mpih%print_message('ib_object%initialize start')

   ! associate ADAM main data
   self%field => field
   self%grid  => grid

   call self%load_from_file(file_parameters=file_parameters)

   ! allocate large arrays
   associate(ngc=>self%grid%ngc, ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, nb=>self%field%nb, &
             solids_number=>self%solids_number)
   if (solids_number > 0) then
      ! the phi array is allocated with solids_number + 1 elements: in the last element there is the all-solids-summary, namely it
      ! contains the information if the cell is outside all IB solids or inside of at least one solid
      allocate(self%phi(1:solids_number+1, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
      self%phi = -1._R8P
   endif
   endassociate

   print '(A)', self%description()
   call self%mpih%print_message('ib_object%initialize finish')
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(ib_object), intent(inout)        :: self            !< IB.
   type(file_ini),   intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,          intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                :: go_on_fail_     !< Go on if load fails.
   character(:), allocatable              :: sname           !< Section name.
   integer(I4P)                           :: i               !< Counter.
   integer(I4P)                           :: err             !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
   call file_parameters%get(section_name=INI_SECTION_NAME//'s', option_name='number', val=self%solids_number, error=err)
   if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'s].(number)')
   call file_parameters%get(section_name=INI_SECTION_NAME//'s', option_name='n_eikonal', val=self%n_eikonal, error=err)
   if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'s].(n_eikonal)')

   if (self%solids_number>=1) then
      allocate(self%s_name(self%solids_number))
      allocate(self%bc_type(self%solids_number))
      allocate(self%q(6,self%solids_number))
      allocate(self%definition(self%solids_number))
      allocate(self%sphere(self%solids_number))
      allocate(self%rectangle(self%solids_number))
      do i=1, self%solids_number
         sname = INI_SECTION_NAME//'_'//trim(str(i,.true.))
         call file_parameters%get(section_name=sname, option_name='name', val=self%s_name(i), error=err)
         if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(name)')
         call file_parameters%get(section_name=sname, option_name='bc_type', val=self%bc_type(i), error=err)
         if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(bc_type)')
         call file_parameters%get(section_name=sname, option_name='definition', val=self%definition(i), error=err)
         if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(definition)')
         select case(trim(adjustl(self%definition(i))))
         case(trim(IB_ANALYTICAL_SPHERE))
            call file_parameters%get(section_name=sname, option_name='sphere_center_x', val=self%sphere(i)%center(1), error=err)
            if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(sphere_center_x)')
            call file_parameters%get(section_name=sname, option_name='sphere_center_y', val=self%sphere(i)%center(2), error=err)
            if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(sphere_center_y)')
            call file_parameters%get(section_name=sname, option_name='sphere_center_z', val=self%sphere(i)%center(3), error=err)
            if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(sphere_center_z)')
            call file_parameters%get(section_name=sname, option_name='sphere_radius', val=self%sphere(i)%radius, error=err)
            if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(sphere_radius)')
         case(trim(IB_ANALYTICAL_CIRCLE))
            call file_parameters%get(section_name=sname, option_name='circle_center_x', val=self%sphere(i)%center(1), error=err)
            if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(circle_center_x)')
            call file_parameters%get(section_name=sname, option_name='circle_center_y', val=self%sphere(i)%center(2), error=err)
            if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(circle_center_y)')
            call file_parameters%get(section_name=sname, option_name='circle_center_z', val=self%sphere(i)%center(3), error=err)
            if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(circle_center_z)')
            call file_parameters%get(section_name=sname, option_name='circle_radius', val=self%sphere(i)%radius, error=err)
            if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(circle_radius)')
            call file_parameters%get(section_name=sname, option_name='circle_axis', val=self%sphere(i)%axis, error=err)
            if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(circle_axis)')
         case(trim(IB_ANALYTICAL_RECTANGLE))
            call file_parameters%get(section_name=sname,option_name='rectangle_center_x',val=self%rectangle(i)%center(1),error=err)
            if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(rectangle_center_x)')
            call file_parameters%get(section_name=sname,option_name='rectangle_center_y',val=self%rectangle(i)%center(2),error=err)
            if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(rectangle_center_y)')
            call file_parameters%get(section_name=sname,option_name='rectangle_center_z',val=self%rectangle(i)%center(3),error=err)
            if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(rectangle_center_z)')
            call file_parameters%get(section_name=sname, option_name='rectangle_edge_1', val=self%rectangle(i)%edge(1), error=err)
            if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(rectangle_edge_1)')
            call file_parameters%get(section_name=sname, option_name='rectangle_edge_2', val=self%rectangle(i)%edge(2), error=err)
            if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(rectangle_edge_2)')
            call file_parameters%get(section_name=sname, option_name='rectangle_axis', val=self%rectangle(i)%axis, error=err)
            if (.not.go_on_fail_.and.err>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(rectangle_axis)')
         endselect
      enddo
   endif
   endsubroutine load_from_file

   subroutine move_phi(self, velocity, s)
   !< Move phi.
   class(ib_object), intent(inout) :: self        !< IB.
   real(R8P),        intent(in)    :: velocity(3) !< Velocity of the movement.
   integer(I4P),     intent(in)    :: s           !< Solid index.

   ! to be implemented
   endsubroutine move_phi

   function sphere_to_array(self, ib) result(array)
   !< Convert analytical sphere class data to array data.
   class(ib_object), intent(in) :: self       !< IB.
   integer(I4P),     intent(in) :: ib         !< Index of IB solid.
   real(R8P)                    :: array(1:4) !< Data array.

   array = [self%sphere(ib)%center(1), self%sphere(ib)%center(2), self%sphere(ib)%center(3), self%sphere(ib)%radius]
   endfunction sphere_to_array

   subroutine evolve_eikonal(self, q)
   !< Evolve eikonal equation.
   class(ib_object), intent(in)    :: self                             !< IB.
   real(R8P),        intent(inout) ::  q(1:,               &
                                         1-self%grid%ngc:, &
                                         1-self%grid%ngc:, &
                                         1-self%grid%ngc:, &
                                         1:)                           !< Conservative variables.
   real(R8P)                       :: dq(1:self%field%nv)              !< Conservative variables differences.
   real(R8P)                       :: n_phi_x, n_phi_y, n_phi_z, n_phi !< Eikonal directions.
   integer(I4P)                    :: i, j, k, b, s                    !< Counter.

   associate(blocks_number=>self%field%blocks_number, ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, ngc=>self%grid%ngc, &
             nv=>self%field%nv, solids_number=>self%solids_number)
   !$omp parallel do collapse(4) default(firstprivate) shared(q,self)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1,ni
               solids_loop : do s=1, solids_number
                  if (self%phi(s,i,j,k,b) > 0._R8P) then
                     ! compute dq
                     n_phi_x = (self%phi(s,i+1,j,k,b) - self%phi(s,i-1,j,k,b))
                     n_phi_y = (self%phi(s,i,j+1,k,b) - self%phi(s,i,j-1,k,b))
                     n_phi_z = (self%phi(s,i,j,k+1,b) - self%phi(s,i,j,k-1,b))
                     n_phi = abs(n_phi_x) + abs(n_phi_y) + abs(n_phi_z) + 10e-12
                     n_phi = 0.9_R8P / n_phi
                     n_phi_x = n_phi_x * n_phi
                     n_phi_y = n_phi_y * n_phi
                     n_phi_z = n_phi_z * n_phi
                     dq(:) = 0._R8P
                     if (n_phi_x > 0._R8P) then
                        dq(:) = dq(:) + abs(n_phi_x) * (q(:,i,j,k,b) - q(:,i-1,j,k,b))
                     else
                        dq(:) = dq(:) + abs(n_phi_x) * (q(:,i,j,k,b) - q(:,i+1,j,k,b))
                     endif
                     if (n_phi_y > 0._R8P) then
                        dq(:) = dq(:) + abs(n_phi_y) * (q(:,i,j,k,b) - q(:,i,j-1,k,b))
                     else
                        dq(:) = dq(:) + abs(n_phi_y) * (q(:,i,j,k,b) - q(:,i,j+1,k,b))
                     endif
                     if (n_phi_z > 0._R8P) then
                        dq(:) = dq(:) + abs(n_phi_z) * (q(:,i,j,k,b) - q(:,i,j,k-1,b))
                     else
                        dq(:) = dq(:) + abs(n_phi_z) * (q(:,i,j,k,b) - q(:,i,j,k+1,b))
                     endif
                     ! evolve q
                     q(:,i,j,k,b) = q(:,i,j,k,b) - dq(:)
                     exit solids_loop
                  endif
               enddo solids_loop
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endassociate
   endsubroutine evolve_eikonal

   subroutine invert_eikonal(self, q)
   !< Invert eikonal equation over q inside IB.
   class(ib_object), intent(in)    :: self                  !< IB.
   real(R8P),        intent(inout) ::  q(1:,               &
                                         1-self%grid%ngc:, &
                                         1-self%grid%ngc:, &
                                         1-self%grid%ngc:, &
                                         1:)                !< Conservative variables.
   integer(I4P)                :: i, j, k, b, s             !< Counter.
   real(R8P)                   :: n_phi_x, n_phi_y, n_phi_z !< Distance function normals.
   real(R8P)                   :: n_phi_mod, un_mod         !< Distance abs normal and normal velocity.

   associate(blocks_number=>self%field%blocks_number, ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, ngc=>self%grid%ngc, &
             solids_number=>self%solids_number)
   !$omp parallel do collapse(4) default(firstprivate) shared(q,self)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1,ni
               solids_loop : do s=1, solids_number
                  if     (self%bc_type(s) == BCS_VISCOUS) then
                     if (self%phi(s,i,j,k,b) > 0) then
                        q(2,i,j,k,b) = - q(2,i,j,k,b)
                        q(3,i,j,k,b) = - q(3,i,j,k,b)
                        q(4,i,j,k,b) = - q(4,i,j,k,b)
                     endif
                  elseif (self%bc_type(s) == BCS_EULER  ) then
                     if (self%phi(s,i,j,k,b) > 0) then
                        n_phi_x = self%phi(s,i+1,j,k,b) - self%phi(s,i-1,j,k,b)
                        n_phi_y = self%phi(s,i,j+1,k,b) - self%phi(s,i,j-1,k,b)
                        n_phi_z = self%phi(s,i,j,k+1,b) - self%phi(s,i,j,k-1,b)
                        n_phi_mod = sqrt(n_phi_x**2 + n_phi_y**2 + n_phi_z**2)
                        n_phi_x = n_phi_x/n_phi_mod
                        n_phi_y = n_phi_y/n_phi_mod
                        n_phi_z = n_phi_z/n_phi_mod
                        un_mod = q(2,i,j,k,b)*n_phi_x + q(3,i,j,k,b)*n_phi_y + q(4,i,j,k,b)*n_phi_z

                        q(2,i,j,k,b) = q(2,i,j,k,b) - 2*un_mod*n_phi_x
                        q(3,i,j,k,b) = q(3,i,j,k,b) - 2*un_mod*n_phi_y
                        q(4,i,j,k,b) = q(4,i,j,k,b) - 2*un_mod*n_phi_z
                     endif
                  endif
               enddo solids_loop
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endassociate
   endsubroutine invert_eikonal

   ! private methods
   subroutine compute_phi_analytical_sphere(self, solid, sphere)
   !< Compute distance for analytical sphere solid.
   class(ib_object),               intent(inout) :: self       !< IB.
   integer(I4P),                   intent(in)    :: solid      !< Solid index.
   type(analytical_sphere_object), intent(in)    :: sphere     !< Analytical sphere solid.
   integer(I4P)                                  :: b, i, j, k !< Counter.

   associate(blocks_number=>self%field%blocks_number, ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, ngc=>self%grid%ngc)
   !$omp parallel do collapse(4) default(firstprivate) shared(self)
   do b=1, blocks_number
      do i=1-ngc, ni+ngc
         do j=1-ngc, nj+ngc
            do k=1-ngc, nk+ngc
               self%phi(solid,i,j,k,b) = - (sqrt((self%field%x_cell(i,b)-sphere%center(1))**2 + &
                                                 (self%field%y_cell(j,b)-sphere%center(2))**2 + &
                                                 (self%field%z_cell(k,b)-sphere%center(3))**2) - sphere%radius)
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endassociate
   endsubroutine compute_phi_analytical_sphere

   subroutine compute_phi_analytical_circle(self, solid, sphere)
   !< Compute distance function for analytical circle (2D) solid.
   class(ib_object),               intent(inout) :: self       !< IB.
   integer(I4P),                   intent(in)    :: solid      !< Solid index.
   type(analytical_sphere_object), intent(in)    :: sphere     !< Analytical circle solid.
   integer(I4P)                                  :: b, i, j, k !< Counter.

   associate(blocks_number=>self%field%blocks_number, ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, ngc=>self%grid%ngc)
   select case(sphere%axis)
   case('x')
      !$omp parallel do collapse(4) default(firstprivate) shared(self)
      do b=1, blocks_number
         do i=1-ngc, ni+ngc
            do j=1-ngc, nj+ngc
               do k=1-ngc, nk+ngc
                  self%phi(solid,i,j,k,b) = - (sqrt((self%field%y_cell(j,b)-sphere%center(2))**2 + &
                                                    (self%field%z_cell(k,b)-sphere%center(3))**2) - sphere%radius)
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   case('y')
      !$omp parallel do collapse(4) default(firstprivate) shared(self)
      do b=1, blocks_number
         do i=1-ngc, ni+ngc
            do j=1-ngc, nj+ngc
               do k=1-ngc, nk+ngc
                  self%phi(solid,i,j,k,b) = - (sqrt((self%field%x_cell(i,b)-sphere%center(1))**2 + &
                                                    (self%field%z_cell(k,b)-sphere%center(3))**2) - sphere%radius)
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   case('z')
      !$omp parallel do collapse(4) default(firstprivate) shared(self)
      do b=1, blocks_number
         do i=1-ngc, ni+ngc
            do j=1-ngc, nj+ngc
               do k=1-ngc, nk+ngc
                  self%phi(solid,i,j,k,b) = - (sqrt((self%field%x_cell(i,b)-sphere%center(1))**2 + &
                                                    (self%field%y_cell(j,b)-sphere%center(2))**2) - sphere%radius)
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   endselect
   endassociate
   endsubroutine compute_phi_analytical_circle

   subroutine compute_phi_analytical_rectangle(self, solid, rectangle)
   !< Compute distance function for analytical rectangle (2D) solid.
   class(ib_object),                  intent(inout) :: self         !< IB.
   integer(I4P),                      intent(in)    :: solid        !< Solid index.
   type(analytical_rectangle_object), intent(in)    :: rectangle    !< Analytical rectangle solid.
   integer(I4P)                                     :: b, i, j, k   !< Counter.
   real(R8P)                                        :: edges_min(2) !< Rectangle min edges position.
   real(R8P)                                        :: edges_max(2) !< Rectangle max edges position.

   associate(blocks_number=>self%field%blocks_number, ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, ngc=>self%grid%ngc, &
             x_cell=>self%field%x_cell, y_cell=>self%field%y_cell, z_cell=>self%field%z_cell, phi=>self%phi)
   select case(rectangle%axis)
   case('x')
      edges_min = rectangle%center(2:3) - rectangle%edge * 0.5_R8P
      edges_max = rectangle%center(2:3) + rectangle%edge * 0.5_R8P
      do b=1, blocks_number
         do i=1-ngc, ni+ngc
            do j=1-ngc, nj+ngc
               do k=1-ngc, nk+ngc
                  phi(solid,i,j,k,b) = edges_distance(point=[y_cell(j,b),z_cell(k,b)])
               enddo
            enddo
         enddo
      enddo
   case('y')
      edges_min = [rectangle%center(1),rectangle%center(3)] - rectangle%edge * 0.5_R8P
      edges_max = [rectangle%center(1),rectangle%center(3)] + rectangle%edge * 0.5_R8P
      do b=1, blocks_number
         do i=1-ngc, ni+ngc
            do j=1-ngc, nj+ngc
               do k=1-ngc, nk+ngc
                  phi(solid,i,j,k,b) = edges_distance(point=[x_cell(i,b),z_cell(k,b)])
               enddo
            enddo
         enddo
      enddo
   case('z')
      edges_min = rectangle%center(1:2) - rectangle%edge * 0.5_R8P
      edges_max = rectangle%center(1:2) + rectangle%edge * 0.5_R8P
      do b=1, blocks_number
         do i=1-ngc, ni+ngc
            do j=1-ngc, nj+ngc
               do k=1-ngc, nk+ngc
                  phi(solid,i,j,k,b) = edges_distance(point=[x_cell(i,b),y_cell(j,b)])
               enddo
            enddo
         enddo
      enddo
   endselect
   endassociate
   contains
      function edges_distance(point) result(distance)
      !< Return point-rectangle distance with sign.
      real(R8P), intent(in) :: point(2)           !< Point coordinate.
      real(R8P)             :: distance           !< Point-rectangle distance.
      real(R8P)             :: distances(2), d(2) !< Point-edges distances.

      ! direction 1
      distances(1) = edges_min(1) - point(1)
      distances(2) = point(1)     - edges_max(1)
      if     (distances(1)<0._R8P.and.distances(2)>0._R8P) then
         d(1) = distances(2)
      elseif (distances(1)<=0._R8P.and.distances(2)<=0._R8P) then
         d(1) = minval(distances, dim=1)
      elseif (distances(1)>0._R8P.and.distances(2)<0._R8P) then
         d(1) = distances(1)
      endif
      ! direction 2
      distances(1) = edges_min(2) - point(2)
      distances(2) = point(2)     - edges_max(2)
      if     (distances(1)<0._R8P.and.distances(2)>0._R8P) then
         d(2) = distances(2)
      elseif (distances(1)<=0._R8P.and.distances(2)<=0._R8P) then
         d(2) = minval(distances, dim=1)
      elseif (distances(1)>0._R8P.and.distances(2)<0._R8P) then
         d(2) = distances(1)
      endif
      ! correct sign in case point falls at corners quarter
      if (d(1)>=0._R8P.and.d(2)>=0._R8P) d(1) = - d(1)
      distance = sign(sqrt(d(1)*d(1) + d(2)*d(2)), d(1)*d(2))
      endfunction edges_distance
   endsubroutine compute_phi_analytical_rectangle
endmodule adam_ib_object
