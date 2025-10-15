!< ADAM, PRISM coil source definition, CPU backend.
module adam_prism_coil_object
    !< ADAM, PRISM coil source definition, CPU backend.

! ADAM modules
use adam_field_object, only : field_object
use adam_mpih_object, only : mpih_object
! PRISM modules
use adam_prism_physics_object, only : prism_physics_object
use adam_prism_parameters
! third party modules
use finer
use penf

implicit none
private
public :: INI_SECTION_NAME
public :: COIL_TYPE_RECTANGULAR
public :: COIL_TYPE_CIRCULAR
public :: CURRENT_TYPE_AC
public :: CURRENT_TYPE_DC
public :: INF_COIL_REPRESENTATION
public :: FIN_COIL_REPRESENTATION
public :: CONST_CURRENT_DISTRIBUTION
public :: GAUSS_CURRENT_DISTRIBUTION
public :: prism_coil_object

character(len=11), parameter :: INI_SECTION_NAME="coils_input"        !< INI (config) file section name containing coils configs.
character(len=11), parameter :: COIL_TYPE_RECTANGULAR="rectangular"   !< Rectangular shape coil.
character(len=8),  parameter :: COIL_TYPE_CIRCULAR="circular"         !< Circular shape coil.
character(len=10), parameter :: CURRENT_TYPE_DC="DC_current"          !< DC current.
character(len=10), parameter :: CURRENT_TYPE_AC="AC_current"          !< AC current
character(len=8),  parameter :: INF_COIL_REPRESENTATION="Infinite"    !< Coil composed by infinite wires.
character(len=6),  parameter :: FIN_COIL_REPRESENTATION="Finite"      !< Coil composed by finite wires.
character(len=8),  parameter :: CONST_CURRENT_DISTRIBUTION="Constant"      !< Coil composed by finite wires.
character(len=8),  parameter :: GAUSS_CURRENT_DISTRIBUTION="Gaussian"      !< Coil composed by finite wires.

type :: prism_coil_object
   !< ADAM, PRISM coil source definition, CPU backend.
   type(mpih_object)              :: mpih                                  !< MPI handler.
   character(len=99), allocatable :: coil_type(:)                          !< Coil type.
   character(len=99), allocatable :: current_type(:)                       !< Current type.
   character(len=99), allocatable :: coil_rep(:)                           !< Coil representation.
   character(len=99), allocatable :: current_distribution(:)               !< Coil representation.
   real(R8P), allocatable         :: A(:)                                  !< Current amplitude (A)
   real(R8P), allocatable         :: f(:)                                  !< Current frequency, if AC (Hz)
   real(R8P), allocatable         :: phi(:,:,:,:,:)                        !< Distance function from the coils
   real(R8P), allocatable         :: phase(:)                              !< Current initial phase, if AC
   real(R8P), allocatable         :: d(:)                                  !< Coil wire diameter
   real(R8P), allocatable         :: x_center(:), y_center(:), z_center(:) !< Coil center
   real(R8P), allocatable         :: lx(:), ly(:)                          !< Rectangle's sizes (if rectangular coil)
   real(R8P), allocatable         :: r_c(:)                                !< Rectangle's radius of curvature (if rectangular coil)
   real(R8P), allocatable         :: r_coil(:)                             !< Circle's radius (if circular coil)
   real(R8P), allocatable         :: normal(:,:)                           !< Versore normale alla spira, che identifica anche verso
   real(R8P), pointer             :: J_vec(:,:,:,:,:)                      !< Matrice contenente versori corrente spire (se assente
   real(R8P)                      :: td                                    !< Delay di accensione della spira
   integer(I4P), pointer          :: coil_flag(:,:,:,:)                    !< Matrice contenente informazioni su quale spira pass pe
   integer(I4P)                   :: circular_coils_number=0_I4P           !< Number of circular coils
   integer(I4P)                   :: rectangular_coils_number=0_I4P        !< Number of rectangular coils
   integer(I4P)                   :: total_coils_number=0_I4P              !< Number of coils
   contains
      ! public methods
      procedure, pass(self) :: compute_distance_naive            !< Compute distance between point and wire (naive way).
      procedure, pass(self) :: curved_wire                       !< Set coils current on PRISM fields.
      procedure, pass(self) :: description                       !< Return pretty-printed object description.
      procedure, pass(self) :: initialize                        !< Initialize IC.
      procedure, pass(self) :: load_from_file                    !< Load config from file.
      procedure, pass(self) :: set_coils                         !< Set coil_object on PRISM fields.
      procedure, pass(self) :: set_circular_coil                 !< Set circular coils on PRISM fields.
      procedure, pass(self) :: set_rectangular_coil              !< Set rectangular coils on PRISM fields.
      procedure, pass(self) :: set_rectangular_coil_quad_section !< Set rectangular coils on PRISM fields with quadratic section
      procedure, pass(self) :: set_rectangular_coil_junction     !< Set rectangular coils on PRISM fields with junction
      procedure, pass(self) :: set_rectangular_coil_infinite     !< Set rectangular coils on PRISM fields with infinite wires
      procedure, pass(self) :: straight_wire                     !< Straight wire function.

endtype prism_coil_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_coil_object), intent(in) :: self             !< IC.
   character(len=:), allocatable        :: desc             !< Description.
   character(len=1), parameter          :: NL=new_line('a') !< New line character.
   integer(I4P)                         :: r                !< Counter.

   desc =       self%mpih%myrankstr//'Coils main data'//NL
   if (self%total_coils_number > 0_I4P) then
      do r=1, self%total_coils_number
         desc = desc//self%mpih%myrankstr//'  Coil('//trim(str(r,.true.))//')'
         select case(self%coil_type(r))
         case(COIL_TYPE_CIRCULAR)
         !desc = desc//NL//self%mpih%myrankstr//'    Coil type: '//trim(str(self%coil_type(r)))
         !desc = desc//NL//self%mpih%myrankstr//'    Current type:'//trim(str(self%current_type(r)))
         desc = desc//NL//self%mpih%myrankstr//'    Radius: '//trim(str(self%r_coil(r)))
         desc = desc//NL//self%mpih%myrankstr//'    Section diameter: '//trim(str(self%d(r)))
         desc = desc//NL//self%mpih%myrankstr//'    Normal: '//trim(str(self%normal(:,r)))
         desc = desc//NL//self%mpih%myrankstr//'    X_center: '//trim(str(self%x_center(r)))
         desc = desc//NL//self%mpih%myrankstr//'    Y_center: '//trim(str(self%y_center(r)))
         desc = desc//NL//self%mpih%myrankstr//'    Z_center: '//trim(str(self%z_center(r)))
         desc = desc//NL//self%mpih%myrankstr//'    Amplitude: '//trim(str(self%A(r)))
         desc = desc//NL//self%mpih%myrankstr//'    requency: '//trim(str(self%f(r)))
         desc = desc//NL//self%mpih%myrankstr//'    Phase: '//trim(str(self%phase(r)))
         case(COIL_TYPE_RECTANGULAR)
         !desc = desc//NL//self%mpih%myrankstr//'    Coil type: '//trim(str(self%coil_type(r)))
         !desc = desc//NL//self%mpih%myrankstr//'    Current type:'//trim(str(self%current_type(r)))
         desc = desc//NL//self%mpih%myrankstr//'    L1: '//trim(str(self%lx(r)))
         desc = desc//NL//self%mpih%myrankstr//'    L2: '//trim(str(self%ly(r)))
         desc = desc//NL//self%mpih%myrankstr//'    Section diameter: '//trim(str(self%d(r)))
         desc = desc//NL//self%mpih%myrankstr//'    Normal: '//trim(str(self%normal(:,r)))
         desc = desc//NL//self%mpih%myrankstr//'    X_center: '//trim(str(self%x_center(r)))
         desc = desc//NL//self%mpih%myrankstr//'    Y_center: '//trim(str(self%y_center(r)))
         desc = desc//NL//self%mpih%myrankstr//'    Z_center: '//trim(str(self%z_center(r)))
         desc = desc//NL//self%mpih%myrankstr//'    Amplitude: '//trim(str(self%A(r)))
         desc = desc//NL//self%mpih%myrankstr//'    requency: '//trim(str(self%f(r)))
         desc = desc//NL//self%mpih%myrankstr//'    Phase: '//trim(str(self%phase(r)))
         endselect
      enddo
   else
      desc = desc//self%mpih%myrankstr//'  No coils defined.'
   endif
   endfunction description

   subroutine initialize(self, file_parameters, field) !Cfr ic%initialize, ma commentata parte descrizione perchè da implementare
   !< Initialize the equation.
   class(prism_coil_object), intent(inout) :: self            !< Coils.
   type(file_ini),         intent(in)      :: file_parameters !< Simulation parameters ini file handler.
   type(field_object),     intent(in)      :: field           !< The field.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'prism_coil_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters, field=field)
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'prism_coil_object%initialize finish'
   endsubroutine initialize

   subroutine compute_distance_naive(self,field)

      class(prism_coil_object), intent(inout) :: self
      type(field_object), intent(in)          :: field
      real(R8P)                               :: x_cell(1-field%grid%ngc:field%grid%ni+field%grid%ngc), &
                                                 y_cell(1-field%grid%ngc:field%grid%nj+field%grid%ngc), &
                                                 z_cell(1-field%grid%ngc:field%grid%nk+field%grid%ngc)  !< Vettori posizione centro celle del blocco b
      real(R8P)                               :: distance
      integer(I4P)                            :: i,j,k,ii,jj,kk,b,r,n_sweep

      associate(ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, blocks_number=>field%blocks_number, &
                ngc=>field%grid%ngc, nb=>field%nb, dx=>field%dxyz(1,:), dy=>field%dxyz(2,:), dz=>field%dxyz(3,:))
      self%phi = 1.0E4_R8P
      n_sweep = 1_I4P !numero di celle da considerare per il calcolo della distanza minima
      !n_sweep = max(ni/2,nj/2,nk/2)
      do r = 1, self%total_coils_number
         do b=1, blocks_number
            call field%grid%cell_xyz(coordinates = field%coordinates(:,b), x_cell = x_cell, y_cell = y_cell, z_cell = z_cell)
            do k=1, nk
               do j=1, nj
                  do i=1, ni
                     if (self%coil_flag(i,j,k,b) == 0_I4P) then
                        do kk= k-n_sweep, k+n_sweep
                           if (kk > nk + ngc .or. kk < 1 - ngc) cycle
                           do jj= j-n_sweep, j+n_sweep
                              if (jj > nj + ngc .or. jj < 1 - ngc) cycle
                              do ii= i-n_sweep, i+n_sweep
                                 if (ii > ni + ngc .or. ii < 1 - ngc) cycle
                                 if (self%coil_flag(ii,jj,kk,b) /= 0_I4P) then
                                    distance = sqrt((x_cell(i)-x_cell(ii))**2 + (y_cell(j)-y_cell(jj))**2 &
                                     + (z_cell(k)-z_cell(kk))**2) - dx(b)/2
                                    self%phi(r,i,j,k,b) = min(self%phi(r,i,j,k,b), distance)
                                 endif
                              enddo
                           enddo
                        enddo
                     else
                        self%phi(r,i,j,k,b) = -dx(b)/2
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
      self%phi = -self%phi
      endassociate

   endsubroutine compute_distance_naive


   subroutine load_from_file(self, file_parameters, field, go_on_fail)
   !< Load config from file.
   class(prism_coil_object), intent(inout)      :: self            !< coils.
   type(file_ini),         intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   type(field_object),     intent(in)           :: field           !< The field.
   logical,                intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                      :: go_on_fail_     !< Go on if load fails.
   character(:), allocatable                    :: sname           !< Section name.
   integer(I4P)                                 :: i               !< Counter.
   integer(I4P)                                 :: error           !< Error status.
   character(99)                                :: buff_char       !< Option character buffer.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='circular_coils_number', &
                            val=self%circular_coils_number, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load &
                                                               ['//INI_SECTION_NAME//'].(circular_coils_number)')

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='rectangular_coils_number', &
                            val=self%rectangular_coils_number, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load &
                                                               ['//INI_SECTION_NAME//'].(rectangular_coils_number)')

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='time_delay', val=self%td, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(time_delay)')

   self%total_coils_number = self%circular_coils_number + self%rectangular_coils_number

      ! Alloczione variabili dell'oggetto spira
      allocate(self%r_coil              (0:self%total_coils_number))
      allocate(self%ly                  (0:self%total_coils_number))
      allocate(self%lx                  (0:self%total_coils_number))
      allocate(self%d                   (0:self%total_coils_number))
      allocate(self%r_c                 (0:self%total_coils_number))
      allocate(self%normal            (3,0:self%total_coils_number))
      allocate(self%x_center            (0:self%total_coils_number))
      allocate(self%y_center            (0:self%total_coils_number))
      allocate(self%z_center            (0:self%total_coils_number))
      allocate(self%coil_type           (0:self%total_coils_number))
      allocate(self%current_type        (0:self%total_coils_number))
      allocate(self%coil_rep            (0:self%total_coils_number))
      allocate(self%current_distribution(0:self%total_coils_number))
      allocate(self%A                   (0:self%total_coils_number))
      allocate(self%f                   (0:self%total_coils_number))
      allocate(self%phase               (0:self%total_coils_number))
      self%r_coil = 0.0_R8P
      self%ly = 0.0_R8P
      self%lx = 0.0_R8P
      self%d = 0.0_R8P
      self%r_c = 0.0_R8P
      self%normal = 0.0_R8P
      self%x_center = 0.0_R8P
      self%y_center = 0.0_R8P
      self%z_center = 0.0_R8P
      self%coil_type = ' '
      self%current_type = ' '
      self%coil_rep = ' '
      self%current_distribution = ' '
      self%A = 0.0_R8P
      self%f = 0.0_R8P
      self%phase = 0.0_R8P

      !Allocazione matrice identificazione spire nelle celle e matrice versori corrente spire nelle celle
      associate(ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, blocks_number=>field%blocks_number, &
                ngc=>field%grid%ngc, nb=>field%nb)

      allocate(self%coil_flag(1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
      self%coil_flag = 0_I4P

      allocate(self%J_vec(4, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
      self%J_vec = 0._R8P
      
      allocate(self%phi(self%total_coils_number,1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
      self%phi = 0._R8P

      endassociate

   if (self%total_coils_number>=1_I4P) then
      
      do i=1, self%total_coils_number
         sname = INI_SECTION_NAME//'_coil_'//trim(str(i,.true.))

         call file_parameters%get(section_name=sname, option_name='coil_type', val=buff_char, error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(coil_type)')
         self%coil_type(i) = trim(buff_char)
         self%coil_type(i) = trim(self%coil_type(i))

         call file_parameters%get(section_name=sname, option_name='current_type', val=buff_char, error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(current_type)')
         self%current_type(i) = trim(buff_char)
         self%current_type(i) = trim(self%current_type(i))

         call file_parameters%get(section_name=sname, option_name='current_distribution', val=buff_char, error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(current_distribution)')
         self%current_distribution(i) = trim(buff_char)
         self%current_distribution(i) = trim(self%current_distribution(i))

         select case(self%coil_type(i))
         case(COIL_TYPE_CIRCULAR)
            call file_parameters%get(section_name=sname, option_name='r_coil', val=self%r_coil(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(r_coil)')

            call file_parameters%get(section_name=sname, option_name='d', val=self%d(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(d)')

            call file_parameters%get(section_name=sname, option_name='x_center', val=self%x_center(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(x_center)')

            call file_parameters%get(section_name=sname, option_name='y_center', val=self%y_center(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(y_center)')

            call file_parameters%get(section_name=sname, option_name='z_center', val=self%z_center(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(z_center)')

            call file_parameters%get(section_name=sname, option_name='nx', val=self%normal(1,i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(nx)')

            call file_parameters%get(section_name=sname, option_name='ny', val=self%normal(2,i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(ny)')

            call file_parameters%get(section_name=sname, option_name='nz', val=self%normal(3,i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(nz)')

            self%lx(i) = 0.0_R8P

            self%ly(i) = 0.0_R8P

         case(COIL_TYPE_RECTANGULAR)
            call file_parameters%get(section_name=sname, option_name='coil_representation', val=buff_char, error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop (msg=': failed to load ['//sname//'].(coil_rep)')
            self%coil_rep(i) = trim(buff_char)
            self%coil_rep(i) = trim(self%coil_rep(i))

            call file_parameters%get(section_name=sname, option_name='lx', val=self%lx(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(lx)')

            call file_parameters%get(section_name=sname, option_name='ly', val=self%ly(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(ly)')

            call file_parameters%get(section_name=sname, option_name='d', val=self%d(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(d)')

            call file_parameters%get(section_name=sname, option_name='r_c', val=self%r_c(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(r_c)')

            call file_parameters%get(section_name=sname, option_name='x_center', val=self%x_center(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(x_center)')

            call file_parameters%get(section_name=sname, option_name='y_center', val=self%y_center(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(y_center)')

            call file_parameters%get(section_name=sname, option_name='z_center', val=self%z_center(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(z_center)')

            call file_parameters%get(section_name=sname, option_name='nx', val=self%normal(1,i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(nx)')

            call file_parameters%get(section_name=sname, option_name='ny', val=self%normal(2,i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(ny)')

            call file_parameters%get(section_name=sname, option_name='nz', val=self%normal(3,i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(nz)')

            self%r_coil(i) = 0.0_R8P

         endselect

         select case(self%current_type(i))
         case(CURRENT_TYPE_DC)

            call file_parameters%get(section_name=sname, option_name='Amplitude', val=self%A(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Amplitude)')

            self%f(i) = 0.0_R8P

            self%phase(i) = 0.0_R8P

         case(CURRENT_TYPE_AC)

            call file_parameters%get(section_name=sname, option_name='Amplitude', val=self%A(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Amplitude)')

            call file_parameters%get(section_name=sname, option_name='Frequency', val=self%f(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Frequency)')

            call file_parameters%get(section_name=sname, option_name='Phase', val=self%phase(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Phase)')

         endselect
         enddo
   endif
   endsubroutine load_from_file

   !Schema: -subroutine set coils: ciclo do con lettura di numero di coil
   !                                    -lettura coil type e conseguente chiamata alla subroutine di riferimento

   subroutine set_coils(self, physics, field)
   !< Set initial conditions on PRISM fields.
   class(prism_coil_object),     intent(inout) :: self    !< Coils
   type(field_object),           intent(inout) :: field   !< Field object.
   type(prism_physics_object),   intent(in)    :: physics !< Fluids physics.
   integer(I4P)                                :: i       !< Counter.

   if (self%total_coils_number >= 1_I4P) then
      do i=1, self%total_coils_number

         select case(self%coil_type(i))
            
         case(COIL_TYPE_CIRCULAR) !Caso spire circolari
            call self%set_circular_coil (physics = physics, field = field, n = i)

         case(COIL_TYPE_RECTANGULAR) !Caso spire rettangolari
            selectcase(self%coil_rep(i))
            case(INF_COIL_REPRESENTATION) !Caso spire rettangolari infinite
               call self%set_rectangular_coil_infinite(physics = physics, field = field, n = i)
            case(FIN_COIL_REPRESENTATION) !Caso spire rettangolari finite
               call self%set_rectangular_coil(physics = physics, field = field, n = i)
            endselect

         endselect

      enddo
      call self%compute_distance_naive(field=field)
   endif
   endsubroutine set_coils

   subroutine set_circular_coil(self, physics, field, n) !agli input aggiungo n del contatore per sapere a quale
                                                         !spira faccio riferimento
      !< Set coils on PRISM fields. La subroutine restituirà il vettore q contenuto in fields
      !< completo anche dei valori normalizzati delle correnti che passano per le celle (elementi 7,8,9)
      !< da calcolare poi tramite la funzione che assegna il valore della corrente compute_coils_current
      class(prism_coil_object),     intent(inout) :: self                                                                !< Coils
      type(field_object),           intent(inout) :: field                                                               !< Field object.
      type(prism_physics_object),   intent(in)    :: physics                                                             !< Fluids physiscs.
      integer(I4P),                 intent(in)    :: n                                                                   !< Coil number.
      !real(R8P),                    allocatable   :: flag(:,:,:,:)                                                      !< Flag per identificare se la spira passa per la cella
      real(R8P),                    allocatable   :: Gaussian(:,:,:,:)                                                   !< Matrice gaussiana per distribuzione corrente 
      real(R8P)                                   :: dmax, dist                                                                !< Vincolo distanza massima dalla spira.
      real(R8P)                                   :: c_c(3)                                                              !< Vettore posizione centro spira
      real(R8P)                                   :: cell_coord(3)                                                       !< Vettore posizione centro cella
      real(R8P)                                   :: x_cell(1-field%grid%ngc:field%grid%ni+field%grid%ngc), &
                                                     y_cell(1-field%grid%ngc:field%grid%nj+field%grid%ngc), &
                                                     z_cell(1-field%grid%ngc:field%grid%nk+field%grid%ngc)               !< Vettori posizione centro celle del blocco b
      integer(I4P)                                :: b,i,j,k                                                             !< Counter.
      !associo per dati su posizioni delle celle e contatori
      associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, &
                x_c => self%x_center(n), y_c => self%y_center(n), z_c => self%z_center(n), &
                dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:), r_coil => self%r_coil(n), &
                normal => self%normal(:,n), d => self%d(n), nb=>field%nb, &
                current_distribution => self%current_distribution(n))

      c_c = [ x_c, y_c, z_c ] !Vettore posizione centro spira

      !allocate(flag(ni,nj,nk,blocks_number))
      do b=1, blocks_number
         ! chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
         call field%grid%cell_xyz(coordinates = field%coordinates(:,b), x_cell = x_cell, y_cell = y_cell, z_cell = z_cell)
         !calcolo distanza massima dall'asse del filo della spira: somma di raggio del filo e metà della dimensione
         !massima della cella associata ai vettori dx dy e dz contenuti in field
         dmax = d/2
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  cell_coord = [x_cell(i), y_cell(j), z_cell(k)]
                  !if (sq_norm(cell_coord-c_c) <= (r_coil+dmax)**2 .and. (r_coil-dmax)**2 <= sq_norm(cell_coord-c_c) .and. &
                      !PI/2-acos(abs((dotproduct(a=(cell_coord-c_c),b=normal)))/sqrt(sq_norm(cell_coord-c_c))) <= asin(d/(2*r_coil)) & 
                      !.and. self%coil_flag(i,j,k,b) == 0_I4P ) then
                     !q(7:9,i,j,k,b) = crossproduct(a=normal,b=(cell_coord-c_c))
                  if ((dotproduct(a=(cell_coord-c_c),b=normal))**2 + (sqrt(sq_norm(cell_coord-c_c) - &
                     (dotproduct(a=(cell_coord-c_c),b=normal))**2) - r_coil)**2 <= (d/2)**2 .and. &
                     self%coil_flag(i,j,k,b) == 0_I4P) then

                     self%J_vec(1:3,i,j,k,b) = crossproduct(a=normal,b=(cell_coord-c_c))

                     !normalizzo per ottenere, alla fine il versore della corrente nella cella
                     !q(7:9,i,j,k,b) = q(7:9,i,j,k,b)/sqrt(sq_norm(q(7:9,i,j,k,b)))
                     self%J_vec(1:3,i,j,k,b) = self%J_vec(1:3,i,j,k,b)/sqrt(sq_norm(self%J_vec(1:3,i,j,k,b)))

                     !metto flag su quale spira passa per la cella
                     self%coil_flag(i,j,k,b) = n

                     selectcase (current_distribution)
                     case (GAUSS_CURRENT_DISTRIBUTION)
                        dist = sqrt((dotproduct(a=(cell_coord-c_c),b=normal))**2 + (sqrt(sq_norm(cell_coord-c_c) - &
                               (dotproduct(a=(cell_coord-c_c),b=normal))**2) - r_coil)**2)
                        self%J_vec(4,i,j,k,b) = gaussian_2D_ind(sigma = d/6, r = dist)
                     case (CONST_CURRENT_DISTRIBUTION)
                        self%J_vec(4,i,j,k,b) = 4.0/(PI*d**2)
                     endselect

                  endif
               enddo
            enddo
         enddo
      enddo
      endassociate
   endsubroutine set_circular_coil

   subroutine set_rectangular_coil(self, physics, field, n) !modificata per avere input equivalente a Filippo, ossia spira a
      !quadrata con dimensione pari a quella della cella. OSS se la spira passa su un'interfaccia non la trova!
      !DA CORREGGERE TUTTA LA SUBROUTINE
   class(prism_coil_object),     intent(inout) :: self                                                            !< Coils
   type(field_object),           intent(inout) :: field                                                           !< Field object.
   type(prism_physics_object),   intent(in)    :: physics                                                         !< Fluids physiscs.
   integer(I4P),                 intent(in)    :: n                                                               !< Coil number.
   integer(I4P),                 allocatable   :: flag(:,:,:,:)                                                   !< Flag per identificare se la spira passa per la cella
   real(R8P)                                   :: dmax                                                            !< Vincolo distanza massima dalla spira.
   real(R8P)                                   :: c_c(3)                                                          !< Vettore posizione centro spira
   real(R8P)                                   :: cell_coord(3)                                                   !< Vettore posizione centro cella
   real(R8P)                                   :: x_cell(1-field%grid%ngc:field%grid%ni+field%grid%ngc), &
                                                  y_cell(1-field%grid%ngc:field%grid%nj+field%grid%ngc), &
                                                  z_cell(1-field%grid%ngc:field%grid%nk+field%grid%ngc)           !< Vettori posizione centro celle del blocco b
   real(R8P)                                   :: vx(3), vy(3), vz(3)                                             !< Versori assi cartesiani
   real(R8P)                                   :: V1(3), V2(3), V3(3), V4(3), V(4,3)                              !< Vertici rettangolo e relativa matrice
   real(R8P)                                   :: v_l1(3), v_l2(3), vec(4,3)                                      !< Versori lati rettangolo (vale regola mano dx) e relativa matrice
   real(R8P)                                   :: n1(3), d1, n2(3), d2                                            !< Parametri piani diagonali perpendicolari a spira, per evitare sovrapposizioni
   real(R8P)                                   :: kappa(3), K_rot(3,3), Id(3,3), theta                            !< Vettore, angolo e matrice di appoggio per formula Rodrigues + matrice identità
   real(R8P)                                   :: Kquad(3,3), R(3,3)                                              !< Matrice K^2 e matrice rotazione tra vz e normale alla spira
   real(R8P)                                   :: dist, prj_v(3)                                                  !< Distanza punto retta e proiezione del punto sulla retta                                                                          !< Variabile utilizzata per definire direzione corrente
   integer(I4P)                                :: b,i,j,k,w                                                       !< Counter.
   !associo per dati su posizioni delle celle e contatori
   associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, &
            x_c => self%x_center(n), y_c => self%y_center(n), z_c => self%z_center(n),                                         &
            dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:), lx => self%lx(n),                             &
            ly => self%ly(n), normal => self%normal(:,n), d => self%d(n), nb =>field%nb)
   c_c = [ x_c, y_c, z_c ] !Vettore posizione centro spira
   !vertici del rettangolo, lo costruisco come se avesse normale asse z e fosse centrato nell'origine;
   !vertici in senso antiorario, partendo da in basso a sinistra; si ipotizza
   !rotazione rispetto al centro, quindi, ad esempio, se la vuoi con normale lungo x
   !con lato lungo lungo y AB sarà il lato corto (pensa a rotazione 3D). Con
   !normale lungo y il problema non si pone
   vx = [1._R8P,0._R8P,0._R8P]
   vy = [0._R8P,1._R8P,0._R8P]
   vz = [0._R8P,0._R8P,1._R8P]
   V1 = [-lx/2, -ly/2, 0._R8P]
   V2 = [+lx/2, -ly/2, 0._R8P]
   V3 = [+lx/2, +ly/2, 0._R8P]
   V4 = [-lx/2, +ly/2, 0._R8P]

   !            V4 ________________ V3     A
   !              |                |       A
   !              |                |       |
   !              |                |       |
   !              |                |       y
   !              |________________|
   !            V1                  V2     x --->>>

   !calcolo rotazione tra i due vettori normali %OSS una rotazione di 180° dà problemi, qua ci metti un go to.
   !isnan per il caso rotazione nulla rispetto
   !a n // z
   kappa = crossproduct(a=vz,b=normal)

   if (any(normal /= vz )) then
       kappa = kappa/sqrt(sq_norm(kappa))
   endif

   theta = acos(dotproduct(a=vz,b=normal))

   K_rot(1,1) = 0._R8P
   K_rot(1,2) = -kappa(3)
   K_rot(1,3) = kappa(2)
   K_rot(2,1) = kappa(3)
   K_rot(2,2) = 0._R8P
   K_rot(2,3) = -kappa(1)
   K_rot(3,1) = -kappa(2)
   K_rot(3,2) = kappa(1)
   K_rot(3,3) = 0._R8P
   Kquad = matmul(K_rot,K_rot) !matrice K^2
   Id(:,1) = vx
   Id(:,2) = vy
   Id(:,3) = vz

   R = Id+sin(theta)*K_rot+(1-cos(theta))*Kquad !costruisco matrice di rotazione, ruoto i vettori posizione e poi traslo
                                               ! a seconda della posizione del centro richiesta

   V1 = matmul(R,V1)+c_c
   V2 = matmul(R,V2)+c_c
   V3 = matmul(R,V3)+c_c
   V4 = matmul(R,V4)+c_c
   V(1,:) = V1
   V(2,:) = V2
   V(3,:) = V3
   V(4,:) = V4
   !genero matrice [V1; V2; V3; V4], quindi le colonne sono le coordinate x y z e le righe i vari vertici nell'ordine descritto prima

   v_l1 = matmul(R,vx)
   v_l1 = v_l1/sqrt(sq_norm(v_l1))
   v_l2 = matmul(R,vy);
   v_l2 = v_l2/sqrt(sq_norm(v_l2))

   !matrice dei versori dei lati, generata ruotando i versori del rettangolo tramite la matrice di rotazione precedentemente calcolata
   vec(1,:) = v_l1
   vec(2,:) = v_l2
   vec(3,:) = -v_l1
   vec(4,:) = -v_l2

   !calcolo piani perpendicolari alla spira su cui giacciono le due diagonali
   !del rettangolo. Per convenzione, le normali puntano "verso i vertici" D e C
   !piano 1, diagonale V1-V3 e normale verso V4
   n1 = crossproduct(a=(V1-c_c),b=normal)/sqrt(sq_norm(V1-c_c)) !normale al piano
   d1 = -dotproduct(a=n1,b=c_c) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale

   !print *, n1, d1, 'normale e d del piano 1'

   !piano 2, diagonale V2-V4 e normale verso V3
   n2 = crossproduct(a=(V4-c_c),b=normal)/sqrt(sq_norm(V4-c_c)) !normale al piano
   d2 = -dotproduct(a=n2,b=c_c) !parametro d dell'equazione ax + by + cz + d2 = 0, con a b c coseni direttori della normale

   !print *, n2, d2, 'normale e d del piano 2'

   allocate(flag(1:ni,1:nj,1:nk,1:nb))
   flag(:,:,:,:) = 0_I4P !inizializzo matrice flag a zero, per indicare che nessun lato passa per le celle

   do w = 1, 4 !per ogni lato del rettangolo
      do b=1, blocks_number
         ! chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
         call field%grid%cell_xyz(coordinates = field%coordinates(:,b), x_cell = x_cell, y_cell = y_cell, z_cell = z_cell)
         !print *, x_cell, 'xcell'
         !print *, y_cell, 'ycell'
         !print *, z_cell, 'zcell'
         !calcolo distanza massima dall'asse del filo della spira: somma di raggio del filo e metà della dimensione massima della cella
         !associata ai vettori dx dy e dz contenuti in field
         !dmax = d/2 + maxval([dx(b),dy(b),dz(b)])/2

         !modificata per avere termine sorgente come Filippo, se infittiamo o aumentiamo sezione spira torna la precedente
         dmax = d/2

         do k=1, nk
            do j=1, nj
               do i=1, ni
                  cell_coord = [x_cell(i), y_cell(j), z_cell(k)]
                  dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:)))) !Distanza punto retta |(P-A) x v| / |v| con A punto sulla retta e v versore della retta
                  prj_v = V(w,:)+dotproduct(a=(cell_coord-V(w,:)),b=vec(w,:))*vec(w,:) !Formula proiezione di un punto su una retta con A punto della retta A+[(P-A)*v]*v
                  !primo if: Il centro cella deve avere distanza dalla retta passante per il lato inferiore alla distanza massima e deve essere all'interno
                  !della proiezione dei lati, altrimenti prendo i punti su tutta la retta. Metto inoltre if sul flag, altrimenti ho sovrapposizioni dal secondo loop

                  !if ( flag(i,j,k,b) == 0_I4P .and. dist <= dmax .and. prj_v(1) <= maxval(V(:,1)) + dmax .and. &
                  !     prj_v(2) <= maxval(V(:,2)) + dmax .and. prj_v(3) <= maxval(V(:,3)) + dmax .and. &
                  !     minval(V(:,1))-dmax <= prj_v(1) .and. minval(V(:,2))-dmax <= prj_v(2) .and. &
                  !     minval(V(:,3))-dmax <= prj_v(3) ) then

                  !if ( flag(i,j,k,b) /= 1_I4P .and. flag(i,j,k,b) /= 3_I4P  .and. &
                  if (    dist <= dmax .and. prj_v(1) <= maxval(V(:,1)) + dmax .and. &
                     prj_v(2) <= maxval(V(:,2)) + dmax .and. prj_v(3) <= maxval(V(:,3)) + dmax .and. &
                     minval(V(:,1))-dmax <= prj_v(1) .and. minval(V(:,2))-dmax <= prj_v(2) .and. &
                     minval(V(:,3))-dmax <= prj_v(3) ) then

                     flag(i,j,k,b) = w

                     !if (flag(i,j,k,b) == w .and. self%coil_flag(i,j,k,b) == 0_I4P) then

                        !q(7:9,i,j,k,b) = vec(flag(i,j,k,b),:)
                     self%J_vec(1:3,i,j,k,b) = self%J_vec(1:3,i,j,k,b) + vec(w,:)
                     self%coil_flag(i,j,k,b) = n
                     self%J_vec(4,i,j,k,b) = 1.0/(d**2) !Distribuzione uniforme, poi normalizzo in compute_coils_current
                        !print*, cell_coord
                        !print*, w
                     !endif

                  endif
                   !secondo if: per ogni lato verifico di essere dal "lato giusto" dei piani definiti dalle diagonali, al fine di
                   !non avere sovrapposizioni in prossimità dei vertici
                  !if (w == 1) then
                  !   if ((dotproduct(a=n1,b=cell_coord)+d1 >= 0 .or. dotproduct(a=n2,b=cell_coord)+d2 >= 0) .and. &
                  !        flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
                  !      flag(i,j,k,b) = 0_I4P
                  !      print *, w
                  !   endif
                  !elseif (w == 2) then
                  !   if ((dotproduct(a=n1,b=cell_coord)+d1 >= 0 .or. dotproduct(a=n2,b=cell_coord)+d2 <= 0) .and. &
                  !        flag(i,j,k,b) == w) then!aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
                  !      flag(i,j,k,b) = 0_I4P
                  !      print *, w
                  !   endif
                  !elseif (w == 3) then
                  !   if ((dotproduct(a=n1,b=cell_coord)+d1 <= 0 .or. dotproduct(a=n2,b=cell_coord)+d2 <= 0) .and. &
                  !        flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
                  !      flag(i,j,k,b) = 0_I4P
                  !      print *, w
                  !   endif
                  !elseif (w == 4) then
                  !   if ((dotproduct(a=n1,b=cell_coord)+d1 <= 0 .or. dotproduct(a=n2,b=cell_coord)+d2 >= 0) .and. &
                  !        flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
                  !      flag(i,j,k,b) = 0_I4P
                  !      print *, w
                  !   endif
                  !endif
               enddo
            enddo
         enddo
      enddo
   enddo
   !Ho un flag pari a 1 2 3 4 nelle celle per cui passa uno dei dati della spira. La direzione della corrente è
   !Ceorente con quella dei versori dei lati precedentemente descritti
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               !if (flag(i,j,k,b) /= 0 .and. self%coil_flag(i,j,k,b) == 0_I4P) then

                  !self%J_vec(1:3,i,j,k,b) = self%J_vec(1:3,i,j,k,b) + vec(flag(i,j,k,b),:)

                  if (abs(self%J_vec(1,i,j,k,b)) < 1.0e-10_R8P) then
                     self%J_vec(1,i,j,k,b) = 0._R8P
                  endif
                  if (abs(self%J_vec(2,i,j,k,b)) < 1.0e-10_R8P) then
                     self%J_vec(2,i,j,k,b) = 0._R8P
                  endif
                  if (abs(self%J_vec(3,i,j,k,b)) < 1.0e-10_R8P) then
                     self%J_vec(3,i,j,k,b) = 0._R8P
                  endif

                  !metto flag su quale spira passa per la cella
                  !self%coil_flag(i,j,k,b) = n

                  !Azzero gli spigoli
               if (sq_norm(self%J_vec(1:3,i,j,k,b)) > 1.1_R8P) then
                  self%J_vec(1:4,i,j,k,b) = 0._R8P
                  self%coil_flag(i,j,k,b) = 0_I4P     
                  print *, 'spigolo azzerato'             
               endif
                  
               !endif
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine set_rectangular_coil

   subroutine set_rectangular_coil_quad_section(self, physics, field, n)

   class(prism_coil_object),     intent(inout) :: self                                                            !< Coils
   type(field_object),           intent(inout) :: field                                                           !< Field object.
   type(prism_physics_object),   intent(in)    :: physics                                                         !< Fluids physiscs.
   integer(I4P),                 intent(in)    :: n                                                               !< Coil number.
   integer(I4P),                 allocatable   :: flag(:,:,:,:)                                                   !< Flag per identificare se la spira passa per la cella
   real(R8P),                    allocatable   :: Gaussian(:,:,:,:)                                               !< Matrice gaussiana per distribuzione corrente 
   real(R8P)                                   :: dmax                                                            !< Vincolo distanza massima dalla spira.
   real(R8P)                                   :: c_c(3)                                                          !< Vettore posizione centro spira
   real(R8P)                                   :: cell_coord(3)                                                   !< Vettore posizione centro cella
   real(R8P)                                   :: x_cell(1-field%grid%ngc:field%grid%ni+field%grid%ngc), &
                                                  y_cell(1-field%grid%ngc:field%grid%nj+field%grid%ngc), &
                                                  z_cell(1-field%grid%ngc:field%grid%nk+field%grid%ngc)           !< Vettori posizione centro celle del blocco b
   real(R8P)                                   :: vx(3), vy(3), vz(3)                                             !< Versori assi cartesiani
   real(R8P)                                   :: V1(3), V2(3), V3(3), V4(3), V(4,3)                              !< Vertici rettangolo e relativa matrice
   real(R8P)                                   :: v_l1(3), v_l2(3), vec(4,3)                                      !< Versori lati rettangolo (vale regola mano dx) e relativa matrice
   real(R8P)                                   :: n1(3), d1, n2(3), d2                                            !< Parametri piani diagonali perpendicolari a spira, per evitare sovrapposizioni
   real(R8P)                                   :: kappa(3), K_rot(3,3), Id(3,3), theta                            !< Vettore, angolo e matrice di appoggio per formula Rodrigues + matrice identità
   real(R8P)                                   :: Kquad(3,3), R(3,3)                                              !< Matrice K^2 e matrice rotazione tra vz e normale alla spira
   real(R8P)                                   :: dist, prj_v(3)                                                  !< Distanza punto retta e proiezione del punto sulla retta                                                                          !< Variabile utilizzata per definire direzione corrente
   integer(I4P)                                :: b,i,j,k,w                                                       !< Counter.
   !associo per dati su posizioni delle celle e contatori
   associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, &
             x_c => self%x_center(n), y_c => self%y_center(n), z_c => self%z_center(n),                                        &
             dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:), lx => self%lx(n),                            &
             ly => self%ly(n), normal => self%normal(:,n), d => self%d(n), nb =>field%nb,                                      &
             current_distribution => self%current_distribution(n))
   c_c = [ x_c, y_c, z_c ] !Vettore posizione centro spira
   !vertici del rettangolo, lo costruisco come se avesse normale asse z e fosse centrato nell'origine;
   !vertici in senso antiorario, partendo da in basso a sinistra; si ipotizza
   !rotazione rispetto al centro, quindi, ad esempio, se la vuoi con normale lungo x
   !con lato lungo lungo y AB sarà il lato corto (pensa a rotazione 3D). Con
   !normale lungo y il problema non si pone
   vx = [1._R8P,0._R8P,0._R8P]
   vy = [0._R8P,1._R8P,0._R8P]
   vz = [0._R8P,0._R8P,1._R8P]
   V1 = [-lx/2, -ly/2, 0._R8P]
   V2 = [+lx/2, -ly/2, 0._R8P]
   V3 = [+lx/2, +ly/2, 0._R8P]
   V4 = [-lx/2, +ly/2, 0._R8P]

   !            V4 ________________ V3     A
   !              |                |       A
   !              |                |       |
   !              |                |       |
   !              |                |       y
   !              |________________|
   !            V1                  V2     x --->>>

   !calcolo rotazione tra i due vettori normali %OSS una rotazione di 180° dà problemi, qua ci metti un go to.
   !isnan per il caso rotazione nulla rispetto
   !a n // z
   kappa = crossproduct(a=vz,b=normal)

   if (any(normal /= vz )) then
       kappa = kappa/sqrt(sq_norm(kappa))
   endif

   theta = acos(dotproduct(a=vz,b=normal))

   K_rot(1,1) = 0._R8P
   K_rot(1,2) = -kappa(3)
   K_rot(1,3) = kappa(2)
   K_rot(2,1) = kappa(3)
   K_rot(2,2) = 0._R8P
   K_rot(2,3) = -kappa(1)
   K_rot(3,1) = -kappa(2)
   K_rot(3,2) = kappa(1)
   K_rot(3,3) = 0._R8P
   Kquad = matmul(K_rot,K_rot) !matrice K^2
   Id(:,1) = vx
   Id(:,2) = vy
   Id(:,3) = vz

   R = Id+sin(theta)*K_rot+(1-cos(theta))*Kquad !costruisco matrice di rotazione, ruoto i vettori posizione e poi traslo
                                               ! a seconda della posizione del centro richiesta

   V1 = matmul(R,V1)+c_c
   V2 = matmul(R,V2)+c_c
   V3 = matmul(R,V3)+c_c
   V4 = matmul(R,V4)+c_c
   V(1,:) = V1
   V(2,:) = V2
   V(3,:) = V3
   V(4,:) = V4 !genero matrice [V1; V2; V3; V4], quindi le colonne sono le coordinate x y z e le righe i vari vertici nell'ordine descritto prima

   v_l1 = matmul(R,vx)
   v_l1 = v_l1/sqrt(sq_norm(v_l1))
   v_l2 = matmul(R,vy);
   v_l2 = v_l2/sqrt(sq_norm(v_l2))

   !matrice dei versori dei lati, generata ruotando i versori del rettangolo tramite la matrice di rotazione precedentemente calcolata
   vec(1,:) = v_l1
   vec(2,:) = v_l2
   vec(3,:) = -v_l1
   vec(4,:) = -v_l2

   !calcolo piani perpendicolari alla spira su cui giacciono le due diagonali
   !del rettangolo. Per convenzione, le normali puntano "verso i vertici" D e C
   !piano 1, diagonale V1-V3 e normale verso V4
   n1 = crossproduct(a=(V1-c_c),b=normal)/sqrt(sq_norm(V1-c_c)) !normale al piano
   d1 = -dotproduct(a=n1,b=c_c) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale

   !print *, n1, d1, 'normale e d del piano 1'

   !piano 2, diagonale V2-V4 e normale verso V3
   n2 = crossproduct(a=(V4-c_c),b=normal)/sqrt(sq_norm(V4-c_c)) !normale al piano
   d2 = -dotproduct(a=n2,b=c_c) !parametro d dell'equazione ax + by + cz + d2 = 0, con a b c coseni direttori della normale

   !print *, n2, d2, 'normale e d del piano 2'

   allocate(flag(1:ni,1:nj,1:nk,1:nb))
   flag(:,:,:,:) = 0_I4P !inizializzo matrice flag a zero, per indicare che nessun lato passa per le celle
   allocate(Gaussian(1:ni,1:nj,1:nk,1:nb))
   Gaussian(:,:,:,:) = 0_I4P !inizializzo matrice gaussiana a zero

   !modificata per avere termine sorgente come Filippo, se infittiamo o aumentiamo sezione spira torna la precedente
   do w = 1, 4 !per ogni lato del rettangolo
      do b=1, blocks_number
            ! If che mi serve per fare spire quadrate avendo dimensioni celle comparabili a sezione spira
            if (abs(d-dx(b))>1.0e-10_R8P) then
               dmax = (d+dx(b))/2
            else
               dmax = d/2
            endif
         ! chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
         call field%grid%cell_xyz(coordinates = field%coordinates(:,b), x_cell = x_cell, y_cell = y_cell, z_cell = z_cell)
         !print *, x_cell, 'xcell'
         !print *, y_cell, 'ycell'
         !print *, z_cell, 'zcell'
         !calcolo distanza massima dall'asse del filo della spira: somma di raggio del filo e metà della dimensione massima della cella
         !associata ai vettori dx dy e dz contenuti in field
         !dmax = d/2 + maxval([dx(b),dy(b),dz(b)])/2

         do k=1, nk
            do j=1, nj
               do i=1, ni
                  cell_coord = [x_cell(i), y_cell(j), z_cell(k)]
                  dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:)))) !Distanza punto retta |(P-A) x v| / |v| con A punto sulla retta e v versore della retta
                  prj_v = V(w,:)+dotproduct(a=(cell_coord-V(w,:)),b=vec(w,:))*vec(w,:); !Formula proiezione di un punto su una retta con A punto della retta A+[(P-A)*v]*v
                  !primo if: Il centro cella deve avere distanza dalla retta passante per il lato inferiore alla distanza massima e deve essere all'interno
                  !della proiezione dei lati, altrimenti prendo i punti su tutta la retta. Metto inoltre if sul flag, altrimenti ho sovrapposizioni dal secondo loop

                  if ( flag(i,j,k,b) == 0_I4P .and. dist <= dmax .and. prj_v(1) <= maxval(V(:,1)) + dmax .and. &
                       prj_v(2) <= maxval(V(:,2)) + dmax .and. prj_v(3) <= maxval(V(:,3)) + dmax .and. &
                       minval(V(:,1))-dmax <= prj_v(1) .and. minval(V(:,2))-dmax <= prj_v(2) .and. &
                       minval(V(:,3))-dmax <= prj_v(3) ) then
                     
                     flag(i,j,k,b) = w
                  endif
                   !secondo if: per ogni lato verifico di essere dal "lato giusto" dei piani definiti dalle diagonali, al fine di
                   !non avere sovrapposizioni in prossimità dei vertici
                  if (w == 1) then
                     if ((dotproduct(a=n1,b=cell_coord)+d1 >= 0 .or. dotproduct(a=n2,b=cell_coord)+d2 >= 0) .and. &
                          flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
                        flag(i,j,k,b) = 0_I4P
                        !print *, w
                     endif
                  elseif (w == 2) then
                     if ((dotproduct(a=n1,b=cell_coord)+d1 >= 0 .or. dotproduct(a=n2,b=cell_coord)+d2 <= 0) .and. &
                          flag(i,j,k,b) == w) then!aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
                        flag(i,j,k,b) = 0_I4P
                        !print *, w
                     endif
                  elseif (w == 3) then
                     if ((dotproduct(a=n1,b=cell_coord)+d1 <= 0 .or. dotproduct(a=n2,b=cell_coord)+d2 <= 0) .and. &
                          flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
                        flag(i,j,k,b) = 0_I4P
                        !print *, w
                     endif
                  elseif (w == 4) then
                     if ((dotproduct(a=n1,b=cell_coord)+d1 <= 0 .or. dotproduct(a=n2,b=cell_coord)+d2 >= 0) .and. &
                          flag(i,j,k,b) == w) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
                        flag(i,j,k,b) = 0_I4P
                        !print *, w
                     endif
                  endif
                  if (flag(i,j,k,b) == w) then
                     selectcase (current_distribution)
                     case (GAUSS_CURRENT_DISTRIBUTION)
                        Gaussian(i,j,k,b) = gaussian_2D_ind(sigma = d/6, r = dist)
                     case (CONST_CURRENT_DISTRIBUTION)
                        Gaussian(i,j,k,b) = 1/d**2
                     endselect
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   !Ho un flag pari a 1 2 3 4 nelle celle per cui passa uno dei dati della spira. La direzione della corrente è
   !Ceorente con quella dei versori dei lati precedentemente descritti
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               if (flag(i,j,k,b) /= 0 .and. self%coil_flag(i,j,k,b) == 0_I4P) then

                  self%J_vec(1:3,i,j,k,b) = vec(flag(i,j,k,b),:)
                  self%J_vec(4,i,j,k,b) = Gaussian(i,j,k,b)

                  if (abs(self%J_vec(1,i,j,k,b)) < 1.0e-10_R8P) then
                     self%J_vec(1,i,j,k,b) = 0._R8P
                  endif
                  if (abs(self%J_vec(2,i,j,k,b)) < 1.0e-10_R8P) then
                     self%J_vec(2,i,j,k,b) = 0._R8P
                  endif
                  if (abs(self%J_vec(3,i,j,k,b)) < 1.0e-10_R8P) then
                     self%J_vec(3,i,j,k,b) = 0._R8P
                  endif

                  self%coil_flag(i,j,k,b) = n

               endif
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine set_rectangular_coil_quad_section

   subroutine set_rectangular_coil_infinite(self, physics, field, n)
   !< Set rectangular coil infinite, i.e. with infinite wires
   class(prism_coil_object),     intent(inout) :: self                                                            !< Coils
   type(field_object),           intent(inout) :: field                                                           !< Field object.
   type(prism_physics_object),   intent(in)    :: physics                                                         !< Fluids physiscs.
   integer(I4P),                 intent(in)    :: n                                                               !< Coil number.
   integer(I4P),                 allocatable   :: flag(:,:,:,:)                                                   !< Flag per identificare se la spira passa per la cella
   real(R8P)                                   :: dmax                                                            !< Vincolo distanza massima dalla spira.
   real(R8P)                                   :: c_c(3)                                                          !< Vettore posizione centro spira
   real(R8P)                                   :: cell_coord(3)                                                   !< Vettore posizione centro cella
   real(R8P)                                   :: x_cell(1-field%grid%ngc:field%grid%ni+field%grid%ngc), &
                                                  y_cell(1-field%grid%ngc:field%grid%nj+field%grid%ngc), &
                                                  z_cell(1-field%grid%ngc:field%grid%nk+field%grid%ngc)           !< Vettori posizione centro celle del blocco b
   real(R8P)                                   :: vx(3), vy(3), vz(3)                                             !< Versori assi cartesiani
   real(R8P)                                   :: V1(3), V2(3), V3(3), V4(3), V(4,3)                              !< Vertici rettangolo e relativa matrice
   real(R8P)                                   :: v_l1(3), v_l2(3), vec(4,3)                                      !< Versori lati rettangolo (vale regola mano dx) e relativa matrice
   real(R8P)                                   :: n1(3), d1, n2(3), d2                                            !< Parametri piani diagonali perpendicolari a spira, per evitare sovrapposizioni
   real(R8P)                                   :: kappa(3), K_rot(3,3), Id(3,3), theta                            !< Vettore, angolo e matrice di appoggio per formula Rodrigues + matrice identità
   real(R8P)                                   :: Kquad(3,3), R(3,3)                                              !< Matrice K^2 e matrice rotazione tra vz e normale alla spira
   real(R8P)                                   :: dist, prj_v(3)                                                  !< Distanza punto retta e proiezione del punto sulla retta                                                                          !< Variabile utilizzata per definire direzione corrente
   integer(I4P)                                :: b,i,j,k,w                                                       !< Counter.
   !associo per dati su posizioni delle celle e contatori
   associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, &
             x_c => self%x_center(n), y_c => self%y_center(n), z_c => self%z_center(n),                                        &
             dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:), lx => self%lx(n),                            &
             ly => self%ly(n), normal => self%normal(:,n), d => self%d(n), nb =>field%nb)
   c_c = [ x_c, y_c, z_c ] !Vettore posizione centro spira
   !vertici del rettangolo, lo costruisco come se avesse normale asse z e fosse centrato nell'origine;
   !vertici in senso antiorario, partendo da in basso a sinistra; si ipotizza
   !rotazione rispetto al centro, quindi, ad esempio, se la vuoi con normale lungo x
   !con lato lungo lungo y AB sarà il lato corto (pensa a rotazione 3D). Con
   !normale lungo y il problema non si pone
   vx = [1._R8P,0._R8P,0._R8P]
   vy = [0._R8P,1._R8P,0._R8P]
   vz = [0._R8P,0._R8P,1._R8P]
   V1 = [-lx/2, -ly/2, 0._R8P]
   V2 = [+lx/2, -ly/2, 0._R8P]
   V3 = [+lx/2, +ly/2, 0._R8P]
   V4 = [-lx/2, +ly/2, 0._R8P]

   !            V4 ________________ V3     A
   !              |                |       A
   !              |                |       |
   !              |                |       |
   !              |                |       y
   !              |________________|
   !            V1                  V2     x --->>>

   !calcolo rotazione tra i due vettori normali %OSS una rotazione di 180° dà problemi, qua ci metti un go to.
   !isnan per il caso rotazione nulla rispetto
   !a n // z
   kappa = crossproduct(a=vz,b=normal)

   if (any(normal /= vz )) then
       kappa = kappa/sqrt(sq_norm(kappa))
   endif

   theta = acos(dotproduct(a=vz,b=normal))

   K_rot(1,1) = 0._R8P
   K_rot(1,2) = -kappa(3)
   K_rot(1,3) = kappa(2)
   K_rot(2,1) = kappa(3)
   K_rot(2,2) = 0._R8P
   K_rot(2,3) = -kappa(1)
   K_rot(3,1) = -kappa(2)
   K_rot(3,2) = kappa(1)
   K_rot(3,3) = 0._R8P
   Kquad = matmul(K_rot,K_rot) !matrice K^2
   Id(:,1) = vx
   Id(:,2) = vy
   Id(:,3) = vz

   R = Id+sin(theta)*K_rot+(1-cos(theta))*Kquad !costruisco matrice di rotazione, ruoto i vettori posizione e poi traslo
                                               ! a seconda della posizione del centro richiesta

   V1 = matmul(R,V1)+c_c
   V2 = matmul(R,V2)+c_c
   V3 = matmul(R,V3)+c_c
   V4 = matmul(R,V4)+c_c
   V(1,:) = V1
   V(2,:) = V2
   V(3,:) = V3
   V(4,:) = V4 !genero matrice [V1; V2; V3; V4], quindi le colonne sono le coordinate x y z e le righe i vari vertici nell'ordine descritto prima

   v_l1 = matmul(R,vx)
   v_l1 = v_l1/sqrt(sq_norm(v_l1))
   v_l2 = matmul(R,vy);
   v_l2 = v_l2/sqrt(sq_norm(v_l2))

   !matrice dei versori dei lati, generata ruotando i versori del rettangolo tramite la matrice di rotazione precedentemente calcolata
   vec(1,:) = v_l1
   vec(2,:) = v_l2
   vec(3,:) = -v_l1
   vec(4,:) = -v_l2

   !calcolo piani perpendicolari alla spira su cui giacciono le due diagonali
   !del rettangolo. Per convenzione, le normali puntano "verso i vertici" D e C
   !piano 1, diagonale V1-V3 e normale verso V4
   n1 = crossproduct(a=(V1-c_c),b=normal)/sqrt(sq_norm(V1-c_c)) !normale al piano
   d1 = -dotproduct(a=n1,b=c_c) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale

   !print *, n1, d1, 'normale e d del piano 1'

   !piano 2, diagonale V2-V4 e normale verso V3
   n2 = crossproduct(a=(V4-c_c),b=normal)/sqrt(sq_norm(V4-c_c)) !normale al piano
   d2 = -dotproduct(a=n2,b=c_c) !parametro d dell'equazione ax + by + cz + d2 = 0, con a b c coseni direttori della normale

   !print *, n2, d2, 'normale e d del piano 2'

   allocate(flag(1:ni,1:nj,1:nk,1:nb))
   flag(:,:,:,:) = 0_I4P !inizializzo matrice flag a zero, per indicare che nessun lato passa per le celle

   do w = 1,4 !per ogni lato del rettangolo
      do b=1, blocks_number
         ! chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
         call field%grid%cell_xyz(coordinates = field%coordinates(:,b), x_cell = x_cell, y_cell = y_cell, z_cell = z_cell)
         !print *, x_cell, 'xcell'
         !print *, y_cell, 'ycell'
         !print *, z_cell, 'zcell'
         !calcolo distanza massima dall'asse del filo della spira: somma di raggio del filo e metà della dimensione massima della cella
         !associata ai vettori dx dy e dz contenuti in field
         !dmax = d/2 + maxval([dx(b),dy(b),dz(b)])/2

         !modificata per avere termine sorgente come Filippo, se infittiamo o aumentiamo sezione spira torna la precedente
                  !modificata per avere termine sorgente come Filippo, se infittiamo o aumentiamo sezione spira torna la precedente
         dmax =  maxval([dx(b),dy(b),dz(b)])/2
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  cell_coord = [x_cell(i), y_cell(j), z_cell(k)]
                  dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:)))) !Distanza punto retta |(P-A) x v| / |v| con A punto sulla retta e v versore della retta
                  prj_v = V(w,:)+dotproduct(a=(cell_coord-V(w,:)),b=vec(w,:))*vec(w,:); !Formula proiezione di un punto su una retta con A punto della retta A+[(P-A)*v]*v
                  !primo if: Il centro cella deve avere distanza dalla retta passante per il lato inferiore alla distanza massima e deve essere all'interno
                  !della proiezione dei lati, altrimenti prendo i punti su tutta la retta. Metto inoltre if sul flag, altrimenti ho sovrapposizioni dal secondo loop

                  if (  dist <= dmax ) then

                     flag(i,j,k,b) = w
                     self%J_vec(1:3,i,j,k,b) = self%J_vec(1:3,i,j,k,b)+vec(flag(i,j,k,b),:)
                     self%coil_flag(i,j,k,b) = n

                     if (abs(self%J_vec(2,i,j,k,b)) < 1.0e-10_R8P) then
                        self%J_vec(1,i,j,k,b) = 0._R8P 
                     endif                 
                     if (abs(self%J_vec(2,i,j,k,b)) < 1.0e-10_R8P) then
                        self%J_vec(2,i,j,k,b) = 0._R8P
                     endif
                     if (abs(self%J_vec(3,i,j,k,b)) < 1.0e-10_R8P) then
                        self%J_vec(3,i,j,k,b) = 0._R8P
                     endif

                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine set_rectangular_coil_infinite

   subroutine set_rectangular_coil_junction(self, physics, field, n)
   class(prism_coil_object),     intent(inout) :: self                                           !< Coils
   type(field_object),           intent(inout) :: field                                          !< Field object.
   type(prism_physics_object),   intent(in)    :: physics                                        !< Fluids physics.
   real(R8P)                                   :: x_cell(1-field%grid%ngc:field%grid%ni+field%grid%ngc), &
                                                  y_cell(1-field%grid%ngc:field%grid%nj+field%grid%ngc), &
                                                  z_cell(1-field%grid%ngc:field%grid%nk+field%grid%ngc)           !< Vettori posizione centro celle del blocco b
   real(R8P)                                   :: c_c(3)                                         !< Vettore posizione centro spira
   real(R8P)                                   :: dmax                                           !< Vincolo distanza massima dalla spira.
   real(R8P)                                   :: kappa(3), K_rot(3,3), Id(3,3), theta           !< Vettore, angolo e matrice di appoggio per formula Rodrigues + matrice identità
   real(R8P)                                   :: Kquad(3,3), R(3,3)                             !< Matrice K^2 e matrice rotazione tra vz e normale alla spira
   real(R8P)                                   :: vx(3), vy(3), vz(3)                            !< Versori assi cartesiani
   real(R8P)                                   :: V1(3), V2(3), V3(3), V4(3), V(4,3)             !< Vertici rettangolo e relativa matrice
   real(R8P)                                   :: C1(3), C2(3), C3(3), C4(3), C(4,3)             !< Centri di curvatura e relativa matrice
   real(R8P)                                   :: prj_C1_l1(3), prj_C2_l1(3), C_l1(2,3) 
   real(R8P)                                   :: prj_C2_l2(3), prj_C3_l2(3), C_l2(2,3)                      
   real(R8P)                                   :: prj_C3_l3(3), prj_C4_l3(3), C_l3(2,3) 
   real(R8P)                                   :: prj_C4_l4(3), prj_C1_l4(3), C_l4(2,3)          !< Proiezioni dei centri di curvatura sui lati
   real(R8P)                                   :: v_l1(3), v_l2(3), vec(4,3)                     !< Versori lati rettangolo (vale regola mano dx) e relativa matrice
   integer(I4P),                 intent(in)    :: n                                              !< Coil number.
   integer(I4P)                                :: b                                               !< Counter.

   
   !associo per dati su posizioni delle celle e contatori
   associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk,      &
               ngc=>field%grid%ngc, x_c => self%x_center(n), y_c => self%y_center(n), z_c => self%z_center(n), &
               dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:), lx => self%lx(n),          &
               ly => self%ly(n), normal => self%normal(:,n), d => self%d(n), nb =>field%nb,                    &
               r_c => self%r_c(n))

   c_c = [ x_c, y_c, z_c ] !Vettore posizione centro spira
   !vertici del rettangolo, lo costruisco come se avesse normale asse z e fosse centrato nell'origine;
   !vertici in senso antiorario, partendo da in basso a sinistra; si ipotizza
   !rotazione rispetto al centro, quindi, ad esempio, se la vuoi con normale lungo x
   !con lato lungo lungo y AB sarà il lato corto (pensa a rotazione 3D). Con
   !normale lungo y il problema non si pone
   vx = [1._R8P,0._R8P,0._R8P]
   vy = [0._R8P,1._R8P,0._R8P]
   vz = [0._R8P,0._R8P,1._R8P]
   V1 = [-lx/2, -ly/2, 0._R8P]
   V2 = [+lx/2, -ly/2, 0._R8P]
   V3 = [+lx/2, +ly/2, 0._R8P]
   V4 = [-lx/2, +ly/2, 0._R8P]

   !Posizione dei centri di curvatura dei vari raccordi e raggio di curvatura (unico per tutti i vertici)
   !C1 = V1 + [r_c, r_c, 0]
   C1(1) = V1(1) + r_c
   C1(2) = V1(2) + r_c
   C1(3) = V1(3)
   !C2 = V2 + [-r_c, r_c, 0]
   C2(1) = V2(1) - r_c
   C2(2) = V2(2) + r_c
   C2(3) = V2(3)
   !C3 = V3 + [-r_c, -r_c, 0]
   C3(1) = V3(1) - r_c
   C3(2) = V3(2) - r_c
   C3(3) = V3(3)
   !C4 = V4 + [r_c, -r_c, 0]
   C4(1) = V4(1) + r_c
   C4(2) = V4(2) - r_c
   C4(3) = V4(3)

   !            V4 ________________ V3     A
   !              |                |       A
   !              |                |       |
   !              |                |       |
   !              |                |       y
   !              |________________|
   !            V1                  V2     x --->>>

   !calcolo rotazione tra i due vettori normali %OSS una rotazione di 180° dà problemi, qua ci metti un go to.
   !isnan per il caso rotazione nulla rispetto
   !a n // z
   kappa = crossproduct(a=vz,b=normal)

   if (any(normal /= vz )) then
       kappa = kappa/sqrt(sq_norm(kappa))
   endif

   theta = acos(dotproduct(a=vz,b=normal))

   K_rot(1,1) = 0._R8P
   K_rot(1,2) = -kappa(3)
   K_rot(1,3) = kappa(2)
   K_rot(2,1) = kappa(3)
   K_rot(2,2) = 0._R8P
   K_rot(2,3) = -kappa(1)
   K_rot(3,1) = -kappa(2)
   K_rot(3,2) = kappa(1)
   K_rot(3,3) = 0._R8P
   Kquad = matmul(K_rot,K_rot) !matrice K^2
   Id(:,1) = vx
   Id(:,2) = vy
   Id(:,3) = vz

   R = Id+sin(theta)*K_rot+(1-cos(theta))*Kquad !costruisco matrice di rotazione, ruoto i vettori posizione e poi traslo
                                               ! a seconda della posizione del centro richiesta

   V1 = matmul(R,V1)+c_c
   V2 = matmul(R,V2)+c_c
   V3 = matmul(R,V3)+c_c
   V4 = matmul(R,V4)+c_c
   V(1,:) = V1
   V(2,:) = V2
   V(3,:) = V3
   V(4,:) = V4 !genero matrice [V1; V2; V3; V4], quindi le colonne sono le coordinate x y z e le righe i vari vertici nell'ordine descritto prima

   C1 = matmul(R,C1)+c_c
   C2 = matmul(R,C2)+c_c
   C3 = matmul(R,C3)+c_c
   C4 = matmul(R,C4)+c_c
   C(1,:) = C1
   C(2,:) = C2
   C(3,:) = C3
   C(4,:) = C4 !genero matrice [C1; C2; C3; C4], quindi le colonne sono le coordinate x y z e le righe i vari vertici nell'ordine descritto prima

   v_l1 = matmul(R,vx)
   v_l1 = v_l1/sqrt(sq_norm(v_l1))
   v_l2 = matmul(R,vy);
   v_l2 = v_l2/sqrt(sq_norm(v_l2))

   !matrice dei versori dei lati, generata ruotando i versori del rettangolo tramite la matrice di rotazione precedentemente calcolata
   vec(1,:) = v_l1
   vec(2,:) = v_l2
   vec(3,:) = -v_l1
   vec(4,:) = -v_l2   

   !proiezione dei centri di curvatura sui lati del rettangolo
   prj_C1_l1 = V1+dotproduct(a=(C1-V1),b=vec(1,:))*vec(1,:)
   prj_C2_l1 = V1+dotproduct(a=(C2-V1),b=vec(1,:))*vec(1,:)
   prj_C2_l2 = V2+dotproduct(a=(C2-V2),b=vec(2,:))*vec(2,:)
   prj_C3_l2 = V2+dotproduct(a=(C3-V2),b=vec(2,:))*vec(2,:)
   prj_C3_l3 = V3+dotproduct(a=(C3-V3),b=vec(3,:))*vec(3,:)
   prj_C4_l3 = V3+dotproduct(a=(C4-V3),b=vec(3,:))*vec(3,:)
   prj_C4_l4 = V4+dotproduct(a=(C4-V4),b=vec(4,:))*vec(4,:)
   prj_C1_l4 = V4+dotproduct(a=(C1-V4),b=vec(4,:))*vec(4,:)

   C_l1(1,:) = prj_C1_l1
   C_l1(2,:) = prj_C2_l1
   C_l2(1,:) = prj_C2_l2
   C_l2(2,:) = prj_C3_l2
   C_l3(1,:) = prj_C3_l3
   C_l3(2,:) = prj_C4_l3
   C_l4(1,:) = prj_C4_l4
   C_l4(2,:) = prj_C1_l4


   dmax = 0.001

   do b=1, blocks_number
      !chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
      call field%grid%cell_xyz(coordinates = field%coordinates(:,b), x_cell = x_cell, y_cell = y_cell, z_cell = z_cell)
      call self%straight_wire(x_cell=x_cell, y_cell=y_cell, z_cell=z_cell, V=C_l1, dmax=dmax, vec=vec(1,:), &
                              ni=ni, nj=nj, nk=nk, blocks_number=blocks_number, ngc=ngc, b=b, numb=n)
      call self%straight_wire(x_cell=x_cell, y_cell=y_cell, z_cell=z_cell, V=C_l2, dmax=dmax, vec=vec(2,:), &
                              ni=ni, nj=nj, nk=nk, blocks_number=blocks_number, ngc=ngc, b=b, numb=n)
      call self%straight_wire(x_cell=x_cell, y_cell=y_cell, z_cell=z_cell, V=C_l3, dmax=dmax, vec=vec(3,:), &
                              ni=ni, nj=nj, nk=nk, blocks_number=blocks_number, ngc=ngc, b=b, numb=n)
      call self%straight_wire(x_cell=x_cell, y_cell=y_cell, z_cell=z_cell, V=C_l4, dmax=dmax, vec=vec(4,:), &
                              ni=ni, nj=nj, nk=nk, blocks_number=blocks_number, ngc=ngc, b=b, numb=n)

   enddo
   endassociate
   endsubroutine set_rectangular_coil_junction

   subroutine curved_wire(self,x_cell,y_cell,z_cell,C,r,dmax,alfa1,alfa2,t,normal,ni,nj,nk,blocks_number,ngc,b,numb)
   !< Compute the current density vector for a curved wire.
   class(prism_coil_object),     intent(inout) :: self                                    !< Coils
   real(R8P),    intent(in)   :: x_cell(1-ngc:ni+ngc), &
                                 y_cell(1-ngc:nj+ngc), &
                                 z_cell(1-ngc:nk+ngc)                                     !< Block b cell center coordinates.
   real(R8P),    intent(in)   :: C(3)                                                     !< Center of curvature of the wire.
   real(R8P),    intent(in)   :: r                                                        !< Radius of curvature of the wire.
   real(R8P),    intent(in)   :: dmax                                                     !< Maximum distance from the wire axis.
   real(R8P),    intent(in)   :: alfa1                                                    !< Angle of the first end of the wire.
   real(R8P),    intent(in)   :: alfa2                                                    !< Angle of the second end of the wire.
   real(R8P),    intent(in)   :: t(3)     !< Vettore tangente al piano per prendere metà corretta della semicirconferenza
   real(R8P),    intent(in)   :: normal(3)                                                !< Coil normal.
   real(R8P)                  :: cell_coord(3)                                            !< Cell center coordinates.
   real(R8P)                  :: n1(3), d1                                                !< Normal and d parameter of the  plane.
   integer(I4P), intent(in)   :: ni,nj,nk,blocks_number,ngc,b                             !< Grid dimensions & block number
   integer(I4P), intent(in)   :: numb                                                        !< Coil number.
   integer(I4P)               :: flag(1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:blocks_number)!< Flag for the current density vector.
   integer(I4P)               :: i,j,k                                                       !< Counters.

   flag(:,:,:,:) = 0_I4P
   n1 = crossproduct(normal,t)
   d1 = -dotproduct(n1,C) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale

   do k=1, nk
      do j=1, nj
         do i=1, ni
            cell_coord = [x_cell(i), y_cell(j), z_cell(k)]
            if (sq_norm(cell_coord-C) <= (r+dmax)**2 .and. (r-dmax)**2 <= sq_norm(cell_coord-C) .and. &
               PI/2 - acos(abs(dotproduct(a=(cell_coord-C),b=normal))) <= asin(dmax/r)) then

               flag(i,j,k,b) = numb
            endif
            if (acos((dotproduct(cell_coord-C,t)))/sq_norm(cell_coord-C)**0.5 >= alfa2 .or. &
               acos((dotproduct(cell_coord-C,t)))/sq_norm(cell_coord-C)**0.5 <= alfa1 .or. &
               dotproduct(n1,cell_coord-C)+d1 >= 0) then

               flag(i,j,k,b) = 0_I4P
            endif
            if (flag(i,j,k,b)/= 0_I4P) then
               self%coil_flag(i,j,k,b) = numb
               self%J_vec(1:3,i,j,k,b) = crossproduct(a=normal,b=(cell_coord-C))/sq_norm(cell_coord-C)**0.5
            endif
         enddo
      enddo
   enddo
   endsubroutine curved_wire

   subroutine straight_wire(self,x_cell,y_cell,z_cell,V,dmax,vec,ni,nj,nk,blocks_number,ngc,b,numb)
   !< Compute the current density vector for a straight wire.
   class(prism_coil_object),     intent(inout) :: self                                       !< Coils
   real(R8P),    intent(in)   :: x_cell(1-ngc:ni+ngc), &
                                 y_cell(1-ngc:nj+ngc), &
                                 z_cell(1-ngc:nk+ngc)                                        !< Block b cell center coordinates.
   real(R8P),    intent(in)   :: V(2,3)                                                      !< Extremal points of the wire.
   real(R8P),    intent(in)   :: dmax                                                        !< Maximum distance from the wire axis.
   real(R8P),    intent(in)   :: vec(3)                                                      !< Direction vector of the wire. 
   real(R8P)                  :: cell_coord(3), dist, prj_v(3)                               !< Cell center coordinates.
   integer(I4P), intent(in)   :: ni,nj,nk,blocks_number,ngc,b                                !< Grid dimensions & block number
   integer(I4P), intent(in)   :: numb                                                        !< Coil number.
   integer(I4P)               :: i,j,k                                                       !< Counters.
   do k=1, nk
      do j=1, nj
         do i=1, ni
            cell_coord = [x_cell(i), y_cell(j), z_cell(k)]
            dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(1,:)),b=vec))) !Distanza punto retta |(P-A) x v| / |v| con A punto sulla retta e v versore della retta
            prj_v = V(1,:)+dotproduct(a=(cell_coord-V(1,:)),b=vec)*vec !Formula proiezione di un punto su una retta con A punto della retta A+[(P-A)*v]*v

            if (dist <= dmax .and. prj_v(1) <= (maxval(V(:,1)) + dmax) .and. &
               prj_v(2) <= (maxval(V(:,2)) + dmax) .and. prj_v(3) <= (maxval(V(:,3)) + dmax) .and. &
               (minval(V(:,1))-dmax) <= prj_v(1) .and. (minval(V(:,2))-dmax) <= prj_v(2) .and. &
               (minval(V(:,3))-dmax) <= prj_v(3)) then

               self%coil_flag(i,j,k,b) = numb
               self%J_vec(1:3,i,j,k,b) = vec
            endif
         enddo
      enddo
   enddo
   endsubroutine straight_wire

   function dotproduct(a, b) result(dot)
   !< Compute the scalar (dot) product.
   real(R8P), intent(in) :: a(3) !< Left hand side.
   real(R8P), intent(in) :: b(3) !< Left hand side.
   real(R8P)             :: dot  !< Dot product.

   dot = (a(1) * b(1)) + (a(2) * b(2)) + (a(3) * b(3))
   endfunction dotproduct

   function crossproduct(a, b) result(cross)
   !< Compute the cross product.
   !<
   !< $$ \vec V=\left({y_1 z_2 - z_1 y_2}\right)\vec i +
   !<           \left({z_1 x_2 - x_1 z_2}\right)\vec j +
   !<           \left({x_1 y_2 - y_1 x_2}\right)\vec k $$
   !< where \( x_i \), \( y_i \) and \( z_i \) \( i=1,2 \) are the components of the vectors.
   !<
   !<```fortran
   !< type(vector_RPP) :: pt(0:2)
   !< pt(1) = 2 * ex_RPP
   !< pt(2) = ex_RPP
   !< pt(0) = pt(1).cross.pt(2)
   !< print "(3(F3.1,1X))", abs(pt(0)%x), abs(pt(0)%y), abs(pt(0)%z)
   !<```
   !=> 0.0 0.0 0.0 <<<
   real(R8P), intent(in) :: a(3)     !< Left hand side.
   real(R8P), intent(in) :: b(3)     !< Left hand side.
   real(R8P)             :: cross(3) !< Cross product.

   cross(1) = (a(2) * b(3)) - (a(3) * b(2))
   cross(2) = (a(3) * b(1)) - (a(1) * b(3))
   cross(3) = (a(1) * b(2)) - (a(2) * b(1))
   endfunction crossproduct

   function sq_norm(a) result(sq)
   !< Return the square of the norm of vector.
   real(R8P), intent(in)  :: a(3)     !< Input vector
   real(R8P)              :: sq       !< Square norm of input

   sq = (a(1) * a(1)) + (a(2) * a(2)) + (a(3) * a(3))
   endfunction sq_norm

   function gaussian_2D_ind(sigma, r) result(f)
   ! Compute the 2D Gaussian function for indipendent variables x,y linked with the distance from the cell centre to the
   ! coil wire in the plain perpendicular to the coil direction. Mu is 0 (obviously)
   real(R8P), intent(in)  :: sigma, r  !< Standard deviation and distance from the coil wire
   real(R8P)              :: f         !< Density of probability function value

   f = exp(-0.5 * (r / sigma)**2) / (2.0 * PI * sigma**2)
   endfunction gaussian_2D_ind

endmodule adam_prism_coil_object
