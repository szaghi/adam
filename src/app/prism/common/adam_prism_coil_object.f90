!< ADAM, PRISM coil source definition, CPU backend.
module adam_prism_coil_object
    !< ADAM, PRISM coil source definition, CPU backend.

use adam_field_object, only : field_object
use adam_mpih_object, only : mpih_object
use adam_prism_physics_object, only : prism_physics_object
use adam_prism_parameters
use finer
use penf
use stringifor !da  controllare

implicit none
private
public :: INI_SECTION_NAME
public :: COIL_TYPE_RECTANGULAR
public :: COIL_TYPE_CIRCULAR
public :: CURRENT_TYPE_AC
public :: CURRENT_TYPE_DC
public :: prism_coil_object

character(len=11), parameter :: INI_SECTION_NAME="coils_input"        !< INI (config) file section name containing coils configs.
character(len=11), parameter :: COIL_TYPE_RECTANGULAR="rectangular"   !< Rectangular shape coil.
character(len=8),  parameter :: COIL_TYPE_CIRCULAR="circular"         !< Circular shape coil.
character(len=15), parameter :: CURRENT_TYPE_DC="DC_current"          !< DC current.
character(len=15), parameter :: CURRENT_TYPE_AC="AC_current"          !< AC current

type :: prism_coil_object
   !< ADAM, PRISM coil source definition, CPU backend.
   type(mpih_object)               :: mpih                                    !< MPI handler.
   !integer(I4P)              :: amr_iterations=1_I4P !< Number of AMR iterations for coils.
   character(len=99), allocatable  :: coil_type(:)                            !< Coil type.
   character(len=99), allocatable  :: current_type(:)                         !< Current type.
   real(R8P), allocatable          :: A(:)                                    !< Current amplitude (A)
   real(R8P), allocatable          :: f(:)                                    !< Current frequency, if AC (Hz)
   real(R8P), allocatable          :: phase(:)                                !< Current initial phase, if AC
   real(R8P), allocatable          :: d(:)                                    !< Coil wire diameter
   real(R8P), allocatable          :: x_center(:), y_center(:), z_center(:)   !< Coil center
   real(R8P), allocatable          :: lx(:), ly(:)                            !< Rectangle's sizes (if rectangular coil)
   real(R8P), allocatable          :: r_coil(:)                               !< Circle's radius (if circular coil)
   real(R8P), allocatable          :: normal(:,:)                             !< Versore normale alla spira, che identifica anche verso della corrente con regola mano dx
   real(R8P), allocatable          :: J_vec(:,:,:,:,:)                        !< Matrice contenente versori corrente spire (se assente = 0)
   real(R8P)                       :: td                                      !< Delay di accensione della spira
   integer(I4P), allocatable       :: coil_flag(:,:,:,:)                      !< Matrice contenente informazioni su quale spira pass per una certa cella
   integer(I4P)                    :: circular_coils_number=0_I4P             !< Number of circular coils
   integer(I4P)                    :: rectangular_coils_number=0_I4P          !< Number of rectangular coils
   integer(I4P)                    :: total_coils_number=0_I4P                !< Number of coils

   !integer(I4P)              :: regions_number=1_I4P !< Number of IC regions.
   !real(R8P), allocatable    :: q(:,:)               !< Primitive variables (r,u,v,w,p), s fluid specie index at IC for each region.
   !real(R8P), allocatable    :: emin(:,:), emax(:,:) !< IC regions bounding box.
   contains
      ! public methods
      !procedure, pass(self) :: description                        !< Return pretty-printed object description.
      procedure, pass(self) :: initialize                          !< Initialize IC.
      procedure, pass(self) :: load_from_file                      !< Load config from file.
      procedure, pass(self) :: set_coils                           !< Set coil_object on PRISM fields.
      procedure, pass(self) :: set_circular_coil                   !< Set circular coils on PRISM fields.
      procedure, pass(self) :: set_rectangular_coil                !< Set rectangular coils on PRISM fields.
      procedure, pass(self) :: set_rectangular_coil_quad_section   !< Set rectangular coils on PRISM fields with quadratic section
endtype prism_coil_object

contains
   ! public methods

   subroutine initialize(self, file_parameters, field) !Cfr ic%initialize, ma commentata parte descrizione perchè da implementare
   !< Initialize the equation.
   class(prism_coil_object), intent(inout) :: self            !< Coils.
   type(file_ini),         intent(in)      :: file_parameters !< Simulation parameters ini file handler.
   type(field_object),     intent(in)      :: field           !< The field.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'prism_coil_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters, field=field)
   !print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'prism_coil_object%initialize finish'
   endsubroutine initialize

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

   if (self%total_coils_number>=1_I4P) then

      ! Alloczione variabili dell'oggetto spira
      allocate(self%r_coil(1:self%total_coils_number))
      allocate(self%ly(1:self%total_coils_number))
      allocate(self%lx(1:self%total_coils_number))
      allocate(self%d(1:self%total_coils_number))
      allocate(self%normal(3,1:self%total_coils_number))
      allocate(self%x_center(1:self%total_coils_number))
      allocate(self%y_center(1:self%total_coils_number))
      allocate(self%z_center(1:self%total_coils_number))
      allocate(self%coil_type(1:self%total_coils_number))
      allocate(self%current_type(1:self%total_coils_number))

      allocate(self%A(1:self%total_coils_number))
      allocate(self%f(1:self%total_coils_number))
      allocate(self%phase(1:self%total_coils_number))

      !Allocazione matrice identificazione spire nelle celle e matrice versori corrente spire nelle celle
      associate(ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, blocks_number=>field%blocks_number, &
                  ngc=>field%grid%ngc, nb=>field%nb)

      allocate(self%coil_flag(1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
      self%coil_flag = 0_I4P

      allocate(self%J_vec(3, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
      self%J_vec = 0._R8P

      endassociate

      do i=1, self%total_coils_number
         sname = INI_SECTION_NAME//'_coil_'//trim(str(i,.true.))

         call file_parameters%get(section_name=sname, option_name='coil_type', val=buff_char, error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(coil_type)')
         self%coil_type(i) = trim(buff_char)
         self%coil_type(i) =trim(self%coil_type(i))

         call file_parameters%get(section_name=sname, option_name='current_type', val=buff_char, error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(current_type)')
         self%current_type(i) = trim(buff_char)
         self%current_type(i) =trim(self%current_type(i))

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

            call file_parameters%get(section_name=sname, option_name='lx', val=self%lx(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(lx)')

            call file_parameters%get(section_name=sname, option_name='ly', val=self%ly(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(ly)')

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
   class(prism_coil_object),     intent(inout) :: self                 !< Coils
   type(field_object),           intent(inout) :: field                !< Field object.
   type(prism_physics_object),   intent(in)    :: physics              !< Fluids physics.
   integer(I4P)                                :: i                    !< Counter.

   if (self%total_coils_number >= 1_I4P) then
      do i=1, self%total_coils_number

         select case(self%coil_type(i))

         case(COIL_TYPE_CIRCULAR) !Caso spire circolari

            call self%set_circular_coil (physics = physics, field = field, n = i)

         case(COIL_TYPE_RECTANGULAR) !Caso spire rettangolari

            call self%set_rectangular_coil_quad_section(physics = physics, field = field, n = i)

         endselect

      enddo
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
      !real(R8P),                    allocatable   :: flag(:,:,:,:)                                                       !< Flag per identificare se la spira passa per la cella
      real(R8P)                                   :: dmax                                                                !< Vincolo distanza massima dalla spira.
      real(R8P)                                   :: c_c(3)                                                              !< Vettore posizione centro spira
      real(R8P)                                   :: cell_coord(3)                                                       !< Vettore posizione centro cella
      real(R8P)                                   :: x_cell(1-field%grid%ngc:field%grid%ni+field%grid%ngc), &
                                                     y_cell(1-field%grid%ngc:field%grid%nj+field%grid%ngc), &
                                                     z_cell(1-field%grid%ngc:field%grid%nk+field%grid%ngc)               !< Vettori posizione centro celle del blocco b
      integer(I4P)                                :: b,i,j,k                                                             !< Counter.
      !associo per dati su posizioni delle celle e contatori
      associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, &
                q=>field%q, x_c => self%x_center(n), y_c => self%y_center(n), z_c => self%z_center(n), &
                dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:), r_coil => self%r_coil(n), &
                normal => self%normal(:,n), d => self%d(n), nb=>field%nb)

      c_c = [ x_c, y_c, z_c ] !Vettore posizione centro spira

      !allocate(flag(ni,nj,nk,blocks_number))
      do b=1, blocks_number
         ! chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
         call field%grid%cell_xyz(coordinates = field%coordinates(:,b), x_cell = x_cell, y_cell = y_cell, z_cell = z_cell)
         !calcolo distanza massima dall'asse del filo della spira: somma di raggio del filo e metà della dimensione
         !massima della cella associata ai vettori dx dy e dz contenuti in field
         dmax = maxval([dx(b),dy(b),dz(b)])/2
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  !Scrivi qua vettore posizione cella b i j k :: cell_coord
                  cell_coord = [x_cell(i), y_cell(j), z_cell(k)]
                  if ( sq_norm(cell_coord-c_c) <= (r_coil+dmax)**2 .and. (r_coil-dmax)**2 <= sq_norm(cell_coord-c_c) .and. &
                        abs(dotproduct(a=(cell_coord-c_c),b=normal)) <= d/r_coil .and. self%coil_flag(i,j,k,b) /= 0_I4P ) then

                     !q(7:9,i,j,k,b) = crossproduct(a=normal,b=(cell_coord-c_c))
                     self%J_vec(:,i,j,k,b) = crossproduct(a=normal,b=(cell_coord-c_c))

                     !normalizzo per ottenere, alla fine il versore della corrente nella cella
                     !q(7:9,i,j,k,b) = q(7:9,i,j,k,b)/sqrt(sq_norm(q(7:9,i,j,k,b)))
                     self%J_vec(:,i,j,k,b) = self%J_vec(:,i,j,k,b)/sqrt(sq_norm(self%J_vec(:,i,j,k,b)))

                     !metto flag su quale spira passa per la cella
                     self%coil_flag(i,j,k,b) = n

                  endif
               enddo
            enddo
         enddo
      enddo
      endassociate
   endsubroutine set_circular_coil

   subroutine set_rectangular_coil(self, physics, field, n) !modificata per avere input equivalente a Filippo, ossia spira a
      !quadrata con dimensione pari a quella della cella. OSS se la spira passa su un'interfaccia non la trova!
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
            q=>field%q, x_c => self%x_center(n), y_c => self%y_center(n), z_c => self%z_center(n), &
            dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:), lx => self%lx(n),  &
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
         dmax =  maxval([dx(b),dy(b),dz(b)])/2
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  cell_coord = [x_cell(i), y_cell(j), z_cell(k)]
                  dist = sqrt(sq_norm(crossproduct(a=(cell_coord-V(w,:)),b=vec(w,:)))) !Distanza punto retta |(P-A) x v| / |v| con A punto sulla retta e v versore della retta
                  prj_v = V(w,:)+dotproduct(a=(cell_coord-V(w,:)),b=vec(w,:))*vec(w,:); !Formula proiezione di un punto su una retta con A punto della retta A+[(P-A)*v]*v
                  !primo if: Il centro cella deve avere distanza dalla retta passante per il lato inferiore alla distanza massima e deve essere all'interno
                  !della proiezione dei lati, altrimenti prendo i punti su tutta la retta. Metto inoltre if sul flag, altrimenti ho sovrapposizioni dal secondo loop
                  !if ( flag(i,j,k,b) == 0_I4P .and. dist <= dmax .and. prj_v(1) <= maxval(V(:,1)) + dmax .and. &
                  !     prj_v(2) <= maxval(V(:,2)) + dmax .and. prj_v(3) <= maxval(V(:,3)) + dmax .and. &
                  !     minval(V(:,1))-dmax <= prj_v(1) .and. minval(V(:,2))-dmax <= prj_v(2) .and. &
                  !     minval(V(:,3))-dmax <= prj_v(3) ) then
                  if ( flag(i,j,k,b) /= 1_I4P .and. flag(i,j,k,b) /= 3_I4P  .and. &
                       dist <= dmax .and. prj_v(1) <= maxval(V(:,1)) + dmax .and. &
                       prj_v(2) <= maxval(V(:,2)) + dmax .and. prj_v(3) <= maxval(V(:,3)) + dmax .and. &
                       minval(V(:,1))-dmax <= prj_v(1) .and. minval(V(:,2))-dmax <= prj_v(2) .and. &
                       minval(V(:,3))-dmax <= prj_v(3) ) then

                        flag(i,j,k,b) = w

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
               if (flag(i,j,k,b) /= 0 .and. self%coil_flag(i,j,k,b) == 0_I4P) then

                  !q(7:9,i,j,k,b) = vec(flag(i,j,k,b),:)
                  self%J_vec(:,i,j,k,b) = vec(flag(i,j,k,b),:)

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

                  self%coil_flag(i,j,k,b) = n
                  !print*, self%coil_flag(i,j,k,b)

               !elseif (flag(i,j,k,b) == 0) then

                  !q(7:9,i,j,k,b) = 0._R8P

               endif
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
             q=>field%q, x_c => self%x_center(n), y_c => self%y_center(n), z_c => self%z_center(n),                            &
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
         dmax =  d/2
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

                  !q(7:9,i,j,k,b) = vec(flag(i,j,k,b),:)
                  self%J_vec(:,i,j,k,b) = vec(flag(i,j,k,b),:)

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

                  self%coil_flag(i,j,k,b) = n

               !elseif (flag(i,j,k,b) == 0) then

                  !q(7:9,i,j,k,b) = 0._R8P
                  !print *, self%coil_flag(i,j,k,b), 'coil flag'
                  !print *, self%J_vec(:,i,j,k,b), 'J vec'

               endif
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine set_rectangular_coil_quad_section
   !Funzioni per prodotto scalare, prodotto vettoriale e norma^2. Per ora sono in set_coils, vediamo se
   !poi metterle in punti più congeniali per evitare riscritture

   function dotproduct(a, b) result(dot)

   !< Compute the scalar (dot) product.

   real(R8P), intent(in) :: a(3)     !< Left hand side.
   real(R8P), intent(in) :: b(3)     !< Left hand side.
   real(R8P)             :: dot       !< Dot product.

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
   !<
   !< The square norm if defined as \( N = x^2  + y^2  + z^2 \).
   !<
   !<```fortran
   !< type(vector_RPP) :: pt
   !< pt = ex_RPP + ey_RPP
   !< print "(F3.1)", pt%sq_norm()
   !<```
   !=> 2.0 <<<
   !<
   !<```fortran
   !< type(vector_RPP) :: pt
   !< pt = ex_RPP + ey_RPP
   !< print "(F3.1)", sq_norm_RPP(pt)
   !<```
   !=> 2.0 <<<
   real(R8P), intent(in)  :: a(3)     !< Input vector
   real(R8P)              :: sq       !< Square norm of input

   sq = (a(1) * a(1)) + (a(2) * a(2)) + (a(3) * a(3))
   endfunction sq_norm

endmodule adam_prism_coil_object
