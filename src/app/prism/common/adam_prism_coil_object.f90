!< ADAM, PRISM coil source definition, CPU backend.
module adam_prism_coil_object
    !< ADAM, PRISM coil source definition, CPU backend.

use adam_field_object, only : field_object
use adam_mpih_object, only : mpih_object
use finer
use penf
use stringifor !da  controllare

implicit none
private
public :: COIL_TYPE_RECTANGULAR
public :: COIL_TYPE_CIRCULAR
public :: CURRENT_TYPE_AC
public :: CURRENT_TYPE_DC
public :: prism_coil_object

character(len=11), parameter :: INI_SECTION_NAME="coils_input"        !< INI (config) file section name containing coils configs.
character(len=11), parameter :: COIL_TYPE_RECTANGULAR="rectangular"   !< Rectangular shape coil.
character(len=8),  parameter :: COIL_TYPE_CIRCULAR="circular"         !< Circular shape coil.
character(len=10), parameter :: CURRENT_TYPE_DC="DC_current"          !< DC current.
character(len=15), parameter :: CURRENT_TYPE_AC="AC_current"          !< AC current


type :: prism_coil_object
   !< ADAM, PRISM coil source definition, CPU backend.
   !type(mpih_object)         :: mpih                 !< MPI handler.
   !integer(I4P)              :: amr_iterations=1_I4P !< Number of AMR iterations imposing IC.
   type(string), allocatable  :: coil_type(:)                            !< Coil type.
   type(string), allocatable  :: current_type(:)                         !< Current type.
   real(R8P), allocatable     :: A(:)                                    !< Current amplitude (A)
   real(R8P), allocatable     :: f(:)                                    !< Current frequency, if AC (Hz)
   real(R8P), allocatable     :: phase(:)                                !< Current initial phase, if AC
   real(R8P), allocatable     :: d(:)                                    !< Coil wire diameter
   real(R8P), allocatable     :: x_center(:), y_center(:), z_center(:)   !< Coil center 
   real(R8P), allocatable     :: lx(:), ly(:)                            !< Rectangle's sizes (if rectangular coil)
   real(R8P), allocatable     :: r_coil(:)                               !< Circle's radius (if circular coil)
   real(R8P), allocatable     :: normal(3,:)                             !< Versore normale alla spira, che identifica anche verso della corrente con regola mano dx
   integer(I4P)               :: circular_coils_number=0_I4P             !< Number of circular coils
   integer(I4P)               :: rectangular_coils_number=0_I4P          !< Number of rectangular coils
   integer(I4P)               :: total_coils_number=0_I4P                !< Number of coils


   !integer(I4P)              :: regions_number=1_I4P !< Number of IC regions.
   !real(R8P), allocatable    :: q(:,:)               !< Primitive variables (r,u,v,w,p), s fluid specie index at IC for each region.
   !real(R8P), allocatable    :: emin(:,:), emax(:,:) !< IC regions bounding box.
   contains
      ! public methods
      !procedure, pass(self) :: description            !< Return pretty-printed object description.
      !procedure, pass(self) :: initialize             !< Initialize IC.
      procedure, pass(self) :: load_from_file         !< Load config from file.
      procedure, pass(self) :: set_coils              !< Set coils on PRISM fields.
endtype prism_coil_object


contains
   ! public methods

    subroutine load_from_file(self, file_parameters, go_on_fail)
    !< Load config from file.
    class(prism_coil_object), intent(inout)      :: self            !< coils.
    type(file_ini),         intent(in)           :: file_parameters !< Simulation parameters ini file handler.
    logical,                intent(in), optional :: go_on_fail      !< Go on if load fails.
    logical                                      :: go_on_fail_     !< Go on if load fails.
    character(:), allocatable                    :: sname           !< Section name.
    integer(I4P)                                 :: i               !< Counter.
    integer(I4P)                                 :: error           !< Error status.
    character(99)                                :: buff_char       !< Option character buffer.
 
    go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
     
    call file_parameters%get(section_name=INI_SECTION_NAME, option_name='circular_coils_number', val=self%circular_coils_number, error=error)
    if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(circular_coils_number)')

    call file_parameters%get(section_name=INI_SECTION_NAME, option_name='rectangular_coils_number', val=self%rectangular_coils_number, error=error)
    if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(rectangular_coils_number)')

    self%total_coils_number = self%circular_coils_number + self%rectangular_coils_number
 
    if (self%total_coils_number>=1) then


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

       do i=1, self%total_coils_number
        sname = INI_SECTION_NAME//'_coil_'//trim(str(i,.true.))

        call file_parameters%get(section_name=sname, option_name='coil_type', val=buff_char, error=error)
        if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(coil_type)') 
        self%coil_type(i) = trim(buff_char)

        call file_parameters%get(section_name=sname, option_name='current_type', val=buff_char, error=error)
        if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(current_type)') 
        self%current_type(i) = trim(buff_char)
        
        select case(self%coil_type(i))
        case(COIL_TYPE_CIRCULAR)

            call file_parameters%get(section_name=sname, option_name='r_coil', val=self%lmax(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(r_coil)')

            call file_parameters%get(section_name=sname, option_name='d', val=self%d(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(d)') 

            call file_parameters%get(section_name=sname, option_name='x_center', val=self%x_center(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(x_center)')

            call file_parameters%get(section_name=sname, option_name='y_center', val=self%y_center(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(y_center)') 

            call file_parameters%get(section_name=sname, option_name='z_center', val=self%z_center(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(z_center)')

            call file_parameters%get(section_name=sname, option_name='normal', val=self%normal(:,i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(normal)')

            self%lx(i) = 0.0_R8P

            self%ly(i) = 0.0_R8P

        case(COIL_TYPE_RECTANGULAR)

            call file_parameters%get(section_name=sname, option_name='lx', val=self%lx(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(lx)')
            
            call file_parameters%get(section_name=sname, option_name='ly', val=self%ly(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(lmin)') 

            call file_parameters%get(section_name=sname, option_name='d', val=self%d(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(d)') 

            call file_parameters%get(section_name=sname, option_name='x_center', val=self%x_center(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(x_center)')

            call file_parameters%get(section_name=sname, option_name='y_center', val=self%y_center(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(y_center)') 

            call file_parameters%get(section_name=sname, option_name='z_center', val=self%z_center(i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(z_center)') 

            call file_parameters%get(section_name=sname, option_name='normal', val=self%normal(:,i), error=error)
            if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(normal)')

            self%r_coil(i) = 0.0_R8P
            
        endselect

        select case(self%current_type(1,i))
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

    !Di sicuro va aggiunta la procedura di inzializzazione, quantomeno per leggere gli input
    !Probabilmente semplicemente sulla falsariga di quella delle IC


    !Schema: -subroutine set coils: ciclo do con lettura di numero di coil
    !                                    -lettura coil type e conseguente chiamata alla subroutine di riferimento


    subroutine set_coils(self, physics, field)  
    !< Set initial conditions on PRISM fields.
    class(prism_coil_object),     intent(in)    :: self                 !< IC.
    type(prism_physics_object),   intent(in)    :: physics              !< Fluids physiscs.
    type(field_object),           intent(inout) :: field                !< Field object.
    integer(I4P)                                :: i       !< Counter.      


    do i=1, self%total_coils_number

        select case(self%coil_type(i))

        case(COIL_TYPE_CIRCULAR) !Caso spire circolari

            call set_circular_coil (physics = physics, field = field, n = i)
        

        case(COIL_TYPE_RECTANGULAR) !Caso spire rettangolari

            call set_rectangular_coil(physics = physics, field = field, n = i)

        endselect

    enddo

    contains
        subroutine set_circular_coil(self, physics, field, n)  !agli input aggiungo n del contatore per sapere a quale 
                                                               !spira faccio riferimento

            !< Set coils on PRISM fields. La subroutine restituirà il vettore q contenuto in fields
            !< completo anche dei valori normalizzati delle correnti che passano per le celle (elementi 7,8,9) 
            !< da calcolare poi tramite la funziona che assegna il valore della corrente, da implementare forse
            !< in CPU (o comunque serve la info del tempo avendo anche AC) 

            class(prism_coil_object),     intent(in)    :: self                                                                !< Coils
            type(prism_physics_object),   intent(in)    :: physics                                                             !< Fluids physiscs.
            type(field_object),           intent(inout) :: field                                                               !< Field object.
            integer(I4P),                 intent(in)    :: n                                                                   !< Coil number.
            !real(R8P),                    allocatable   :: flag(:,:,:,:)                                                       !< Flag per identificare se la spira passa per la cella
            real(R8P)                                   :: dmax                                                                !< Vincolo distanza massima dalla spira.
            real(R8P)                                   :: c_c(3)                                                              !< Vettore posizione centro spira
            real(R8P)                                   :: cell_coord(3)                                                       !< Vettore posizione centro cella
            real(R8P)                                   :: x_cell(1-field%ngc:field%ni+field%ngc), y_cell(1-field%ngc:field%nj+field%ngc), &
                                                           z_cell(1-field%ngc:field%nk+field%ngc)                              !< Vettori posizione centro celle del blocco b
            real(R8P)                                   :: dx_b, dy_b, dz_b                                                    !< dimensione lati celle del blocco b-esimo
            integer(I4P)                                :: b,i,j,k                                                             !< Counter.


            !associo per dati su posizioni delle celle e contatori
            associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, &
                q=>field%q, x_c => self%x_center(n), y_c => self%y_center(n), z_c => self%z_center(n), &
                dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:), r_coil => self%r_coil(n), normal => self%normal(n), d => self%d(n))

            c_c = [ x_c, y_c, z_c ] !Vettore posizione centro spira
            
            !allocate(flag(ni,nj,nk,blocks_number))

            do b=1, blocks_number

                ! chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
                call grid%cell_xyz(coordinates = field%coordinates(:,b), x_cell = x_cell, y_cell = y_cell, z_cell = z_cell)

                !calcolo distanza massima dall'asse del filo della spira: somma di raggio del filo e metà della dimensione massima della cella
                !associata ai vettori dx dy e dz contenuti in field
                dmax = d/2 + max([dx(b),dy(b),dz(b)])

                do k=1, nk
                   do j=1, nj
                      do i=1, ni
                        !Scrivi qua vettore posizione cella b i j k :: cell_coord
                        cell_coord = [x_cell(i), y_cell(j), z_cell(k)]

                        if ( sq_norm(cell_coord-c_c) <= sq_norm(r_coil+dmax) .and. sq_norm(r_coil-dmax) <= sq_norm(cell_coord-c_c) .and. abs(dotproduct((cell_coord-c_c),normal)) <= d/r_coil ) then

                            !flag(i,j,k,b) = 1._R8P
                            q(7:9,i,j,k,b) = crossproduct(normal,(cell_coord-c_c))

                            !normalizzo per ottenere, alla fine il versore della corrente nella cella
                            q(7:9,i,j,k,b) = q(7:9,i,j,k,b)/sqrt(sq_norm(q(7:9,i,j,k,b)))

                        endif

                      enddo
                    enddo
                enddo
            enddo
            endassociate

        endsubroutine set_circular_coil


        subroutine set_rectangular_coil(self, physics, field, n)

        class(prism_coil_object),     intent(in)    :: self                                                                                !< Coils
        type(prism_physics_object),   intent(in)    :: physics                                                                             !< Fluids physiscs.
        type(field_object),           intent(inout) :: field                                                                               !< Field object.
        integer(I4P),                 intent(in)    :: n                                                                                   !< Coil number.
        integer(I4P),                 allocatable   :: flag(:,:,:,:)                                                                       !< Flag per identificare se la spira passa per la cella
        real(R8P)                                   :: dmax                                                                                !< Vincolo distanza massima dalla spira.
        real(R8P)                                   :: c_c(3)                                                                              !< Vettore posizione centro spira
        real(R8P)                                   :: cell_coord(3)                                                                       !< Vettore posizione centro cella
        real(R8P)                                   :: x_cell(1-field%ngc:field%ni+field%ngc), y_cell(1-field%ngc:field%nj+field%ngc), &
                                                       z_cell(1-field%ngc:field%nk+field%ngc)                                              !< Vettori posizione centro celle del blocco b
        real(R8P)                                   :: dx_b, dy_b, dz_b                                                                    !< dimensione lati celle del blocco b-esimo
        real(R8P)                                   :: vx(3), vy(3), vz(3)                                                                 !< Versori assi cartesiani
        real(R8P)                                   :: V1(3), V2(3), V3(3), V4(3), V(4,3)                                                  !< Vertici rettangolo e relativa matrice 
        real(R8P)                                   :: v_l1(3), v_l2(3), vec(4,3)                                                          !< Versori lati rettangolo (vale regola mano dx) e relativa matrice 
        real(R8P)                                   :: n1(3), d1, n2(3), d2                                                                !< Parametri piani diagonali perpendicolari a spira, per evitare sovrapposizioni
        real(R8P)                                   :: kappa(3), K(3,3), Id(3,3), theta                                                    !< Vettore, angolo e matrice di appoggio per formula Rodrigues + matrice identità
        real(R8P)                                   :: Kquad(3,3), R(3,3)                                                                  !< Matrice K^2 e matrice rotazione tra vz e normale alla spira
        real(R8P)                                   :: dist, prj_v(3)                                                                      !< Distanza punto retta e proiezione del punto sulla retta
        real(R8P)                                   :: flag_real                                                                           !< Variabile utilizzata per definire direzione corrente
        integer(I4P),                               :: b,i,j,k,w                                                                           !< Counter.

        !associo per dati su posizioni delle celle e contatori
        associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, &
            q=>field%q, x_c => self%x_center(n), y_c => self%y_center(n), z_c => self%z_center(n), &
            dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:), lx => self%lx(n),  &
            ly => self%ly(n), normal => self%normal(n), d => self%d(n))

        c_c = [ x_c, y_c, z_c ] !Vettore posizione centro spira

        !vertici del rettangolo, lo costruisco come se avesse normale asse z e fosse centrato nell'origine;
        !vertici in senso antiorario, partendo da in basso a sinistra; si ipotizza
        !rotazione rispetto al centro, quindi, ad esempio, se la vuoi con normale lungo x
        !con lato lungo lungo y AB sarà il lato corto (pensa a rotazione 3D). Con
        !normale lungo y il problema non si pone

        vx = [1,0,0]
        vy = [0,1,0]
        vz = [0,0,1]

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

        kappa = crossproduct(vz,normal) 

        if ( normal /= vz ) then

            kappa = kappa/sqrt(sq_norm(kappa))

        endif

        theta = acos(dotproduct(vz,normal)); 
        
        K(1,1) = 0._R8P
        K(1,2) = -kappa(3)
        K(1,3) = kappa(2)
        K(2,1) = kappa(3)
        K(2,2) = 0._R8P
        K(2,3) = -kappa(1) 
        K(3,1) = -kappa(2) 
        K(3,2) = kappa(1) 
        K(3,3) = 0._R8P
        Kquad = K**2
        I = [vx, vy, vz]
        
        R = I+sin(theta)*K+(1-cos(theta))*Kquad
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
        
        !matrice dei versori dei lati
        vec(1,:) = v_l1
        vec(2,:) = v_l2
        vec(3,:) = -v_l1
        vec(4,:) = -v_l2


        !calcolo piani perpendicolari alla spira su cui giacciono le due diagonali
        !del rettangolo. Per convenzione, le normali puntano "verso i vertici" D e C


        !CONTROLLA PRODOTTI SCALARE E VETTORIALE
        !piano 1, diagonale V1-V3 e normale verso V4
        n1 = crossproduct((V1-c_c),normal)/sqrt(sq_norm(V1-c_c)) !normale al piano 
        d1 = -dotproduct(n1,c_c) !parametro d dell'equazione ax + by + cz + d1 = 0, con a b c coseni direttori della normale
        
        !piano 2, diagonale V2-V4 e normale verso V3
        n2 = crossproduct((V4-c_c),normal)/sqrt(sq_norm(V4-c_c)) !normale al piano 
        d2 = -dotproduct(n2,c_c) !parametro d dell'equazione ax + by + cz + d2 = 0, con a b c coseni direttori della normale


        allocate(flag(ni,nj,nk,blocks_number))

        do w = 1, 4
            do b=1, blocks_number

                ! chiamo funzione che restituisce le coordinate delle varie celle che compongono il blocco
                call grid%cell_xyz(coordinates = field%coordinates(:,b), x_cell = x_cell, y_cell = y_cell, z_cell = z_cell)

                !calcolo distanza massima dall'asse del filo della spira: somma di raggio del filo e metà della dimensione massima della cella
                !associata ai vettori dx dy e dz contenuti in field
                dmax = d/2 + max([dx(b),dy(b),dz(b)])


                do k=1, nk
                   do j=1, nj
                        do i=1, ni

                        !Scrivi qua vettore posizione cella b i j k :: cell_coord
                        cell_coord = [x_cell(i), y_cell(j), z_cell(k)]

                        dist = sqrt(sq_norm(crossproduct((cell_coord-V(w,:)),vec(w,:)))) !Distanza punto rettta ||P-A|| x v / |v| con A punto sulla retta e v versore della retta

                        !correggi questo prodotto. 
                        prj_v = V(w,:)+dotproduct((cell_coord-V(w,:)),vec(w,:))*vec(w,:); !Formula proiezione di un punto su una retta con A punto della retta (da cambiare con ciclo for avendo messo V1 ora)
                                                                              !A+[(P-A)*v]*v. 

                        !Verifica la scrittura di 0 "intero"
                        if ( flag(b,i,j,k) == 0 .and. d <= dmax .and. prj_v(1) <= max(V(:,1)) + dmax .and. prj_v(2) <= max(V(:,2)) + dmax .and. prj_v(3) <= max(V(:,3)) & 
                             + dmax .and. min(V(:,1))-dmax <= prj_v(1) .and. min(V(:,2))-dmax <= prj_v(2) .and. min(V(:,3))-dmax <= prj_v(3) ) then
                            flag(b,i,j,k) = j
                        endif

                        !Scrivi bene questi prodotti scalari
                        if (j == 1) then
                            if ((dotproduct(n1,cell_coord)+d1 >= 0 .or. dotproduct(n2,cell_coord)+d2 >= 0) .and. flag(b,i,j,k) == j) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
                                flag(b,i,j,k) = 0_I4P
                            endif
                        elseif (j == 2) then
                            if ((dotproduct(n1,cell_coord)+d1 >= 0 .or. dotproduct(n2,cell_coord)+d2 <= 0) .and. flag(b,i,j,k) == j) then!aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
                                flag(b,i,j,k) = 0_I4P
                            endif
                        elseif (j == 3) then
                            if ((dotproduct(n1,cell_coord)+d1 <= 0 .or. dotproduct(n2,cell_coord)+d2 <= 0) .and. flag(b,i,j,k) == j) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
                                flag(b,i,j,k) = 0_I4P
                            endif                
                        elseif (j == 4) then 
                            if ((dotproduct(n1,cell_coord)+d1 <= 0 .or. dotproduct(n2,cell_coord)+d2 >= 0) .and. flag(b,i,j,k) == j) then !aggiungo secondo if per evitare sovrapposizioni tra celle per i vari lati
                                flag(b,i,j,k) = 0_I4P
                            endif
                        endif
                        enddo
                    enddo
                enddo
            enddo
        enddo

        !Ho un flag pari a 1 2 3 4 nelle celle per cui passa uno dei dati della spira. La direzione della corrente è 
        !Ceorente con quella dei versori dei lati precedentemente descritti
        do k=1, nk
            do j=1, nj
                 do i=1, ni
                    if (flag(b,i,j,k) ~= 0) then
                        q(7:9,i,j,k,b) = vec(flag(b,i,j,k),:)
                    elseif (flag(b,i,j,k) == 0) then
                        q(7:9,i,j,k,b) = 0._R8P
                    endif
                enddo
            enddo
        enddo
        
        endassociate
        endsubroutine set_rectangular_coil

        !Funzioni per prodotto scalare, prodotto vettoriale e norma^2. Per ora sono in set_coils, vediamo se poi metterle in punti più congeniali per evitare riscritture
    
        elemental function dotproduct(a, b) result(dot)
    
        !< Compute the scalar (dot) product.
    
        real(R8P), intent(in) :: a(3)     !< Left hand side.
        real(R8P), intent(in) :: b(3)     !< Left hand side.
        real(R8P), intent(out):: dot       !< Dot product.
     
        dot = (a(1) * b(1)) + (a(2) * b(2)) + (a(3) * b(3))
    
        endfunction dotproduct


        elemental function crossproduct(a, b) result(cross)
    
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
        real(R8P), intent(out):: cross(3) !< Cross product.
    
        cross(1) = (a(2) * b(3)) - (a(3) * b(2))
        cross(2) = (a(3) * b(1)) - (a(1) * b(3))
        cross(3) = (a(1) * b(2)) - (a(2) * b(1))
    
        endfunction crossproduct


        elemental function sq_norm(a) result(b)
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
        real(R8P), intent(out) :: sq       !< Square norm of input
     
        sq = (a(1) * a(1)) + (a(2) * a(2)) + (a(3) * a(3))
        endfunction sq_norm


    endsubroutine set_coils

endmodule adam_prism_coil_object


   







