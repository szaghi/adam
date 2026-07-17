!< ADAM, FLAIL Fortran Linear Algebra Interface Library class definition, CPU common backend.
module adam_flail_object
!< ADAM, FLAIL Fortran Linear Algebra Interface Library class definition, CPU common backend.

! ADAM singleton objects
use :: adam_mpih_global, only : mpih
use :: adam_field_object
use :: adam_parameters, only : FEC_1_6_ARRAY
! third party modules
use :: finer
use :: penf

implicit none
private
public :: flail_object
public :: compute_smoothing_interface
public :: compute_smoothing_gauss_seidel
public :: compute_smoothing_gauss_seidel_centered_dg
public :: compute_smoothing_gauss_seidel_2nd
public :: compute_smoothing_gauss_seidel_4th
public :: compute_smoothing_gauss_seidel_6th
public :: compute_smoothing_gauss_seidel_8th
public :: compute_smoothing_multigrid
public :: compute_smoothing_sor
public :: compute_smoothing_sor_omp

character(9),  parameter, public :: SMOOTHING_MULTIGRID   ='MULTIGRID'    !< Smoothing multigrid parameter.
character(12), parameter, public :: SMOOTHING_GAUSS_SEIDEL='GAUSS-SEIDEL' !< Smoothing Gauss-Seidel parameter.
character(3),  parameter, public :: SMOOTHING_SOR         ='SOR'          !< Smoothing SOR parameter.
character(7),  parameter, public :: SMOOTHING_SOR_OMP     ='SOR-OMP'      !< Smoothing SOR-OpenMP parameter.
integer(I4P),  parameter         :: ELL_BC_DIRICHLET      = 1_I4P         !< Elliptic Dirichlet BC.
integer(I4P),  parameter         :: ELL_BC_PERIODIC       = 2_I4P         !< Elliptic periodic BC.
integer(I4P),  parameter         :: ELL_BC_EXACT_OPEN     = 3_I4P         !< Elliptic exact/open BC.

character(len=14), parameter :: INI_SECTION_NAME="linear-algebra" !< INI (config) file section name containing FLAIL configs.

type :: flail_object
   !< FLAIL Fortran Linear Algebra Interface Library class definition, CPU common backend.
   integer(I4P)              :: iterations=1_I4P        !< Number of iterations, general.
   integer(I4P)              :: iterations_init=1_I4P   !< Number of iterations on fine grid for initialize guess.
   integer(I4P)              :: iterations_fine=1_I4P   !< Number of iterations on fine grid.
   integer(I4P)              :: iterations_coarse=1_I4P !< Number of iterations on coarse grid.
   real(R8P)                 :: tolerance=0.000001_R8P  !< Tolerance on maximum residuals.
   character(:), allocatable :: smoothing               !< Iterative smoothing method.
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize time handler.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype flail_object

interface
   subroutine compute_smoothing_interface(ni, nj, nk, ngc, nv, blocks_number, dxyz, f, q, dq, dq_max, &
                                          iterations_init, iterations_fine, iterations_coarse)
   !< Compute smoothing, abstract interface.
   import :: I4P, R8P
   integer(I4P), intent(in)              :: ni,nj,nk,ngc      !< Grid dimensions.
   integer(I4P), intent(in)              :: nv                !< Number of q variables.
   integer(I4P), intent(in)              :: blocks_number     !< Number of current blocks.
   real(R8P),    intent(in)              :: dxyz(1:,1:)       !< Space steps.
   real(R8P),    intent(in)              :: f(1:,    &
                                              1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:,&
                                              1:)             !< Forcing distribution.
   real(R8P),    intent(inout)           :: q(1:,    &
                                              1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:,&
                                              1:)             !< Field variables.
   real(R8P),    intent(inout), optional :: dq(1:,    &
                                               1-ngc:,&
                                               1-ngc:,&
                                               1-ngc:,&
                                               1:)            !< Residuals.
   real(R8P),    intent(inout), optional :: dq_max            !< Maximum residual.
   integer(I4P), intent(in),    optional :: iterations_init   !< Smoothing iterations to initialize guess.
   integer(I4P), intent(in),    optional :: iterations_fine   !< Smoothing iterations for fine grid.
   integer(I4P), intent(in),    optional :: iterations_coarse !< Smoothing iterations for coarse grid.
   endsubroutine compute_smoothing_interface
endinterface
contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(flail_object), intent(in) :: self             !< Time handler.
   character(len=:), allocatable   :: desc             !< Description.
   character(len=1), parameter     :: NL=new_line('a') !< New line character.

   desc =       mpih%myrankstr//'Linear Algebra methods data'//NL
   desc = desc//mpih%myrankstr//'  smoothing method:  '//self%smoothing                   //NL
   desc = desc//mpih%myrankstr//'  iterations:        '//trim(str(self%iterations       ))//NL
   desc = desc//mpih%myrankstr//'  iterations init:   '//trim(str(self%iterations_init  ))//NL
   desc = desc//mpih%myrankstr//'  iterations fine:   '//trim(str(self%iterations_fine  ))//NL
   desc = desc//mpih%myrankstr//'  tolerance:         '//trim(str(self%tolerance        ))//NL
   desc = desc//mpih%myrankstr//'  iterations coarse: '//trim(str(self%iterations_coarse))
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize time handler.
   class(flail_object), intent(inout) :: self            !< Time handler.
   type(file_ini),      intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   print '(A)', mpih%myrankstr//'flail_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', mpih%myrankstr//'flail_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(flail_object), intent(inout)        :: self            !< Time handler.
   type(file_ini),      intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,             intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                   :: go_on_fail_     !< Go on if load fails.
   character(99)                             :: buffer          !< Buffer.
   integer(I4P)                              :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='smoothing', val=buffer, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(smoothing)')
   select case(trim(adjustl(buffer)))
   case('MULTIGRID','multigrid','Multigrid')
      self%smoothing = SMOOTHING_MULTIGRID
   case('GAUSS-SEIDEL','gauss-seidel','Gauss-Seidel')
      self%smoothing = SMOOTHING_GAUSS_SEIDEL
   case('SOR','sor','Sor')
      self%smoothing = SMOOTHING_SOR
   case('SOR-OMP','sor-omp','Sor-omp')
      self%smoothing = SMOOTHING_SOR_OMP
   endselect
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='iterations', val=self%iterations, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(iterations)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='iterations_init', val=self%iterations_init, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(iterations_init)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='iterations_fine', val=self%iterations_fine, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(iterations_fine)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='iterations_coarse', val=self%iterations_coarse, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(iterations_coarse)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='tolerance', val=self%tolerance, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(tolerance)')
   endsubroutine load_from_file

   ! non TBP
   subroutine apply_bc_dirichlet(ni, nj, nk, ngc, blocks_number, q)
   integer(I4P), intent(in)    :: ni,nj,nk,ngc  !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number !< Number of current blocks.
   real(R8P),    intent(inout) :: q(1:,    &
                                    1-ngc:,&
                                    1-ngc:,&
                                    1-ngc:,&
                                    1:)         !< Field variables.
   integer(I4P)                :: b, gc         !< Counter.

   do b=1, blocks_number
      do gc = 1, ngc
      q(:,1-gc,:,:,b) = 0._R8P ; q(:,ni+gc,:   ,:   ,b) = 0._R8P
      q(:,:,1-gc,:,b) = 0._R8P ; q(:,:   ,nj+gc,:   ,b) = 0._R8P
      q(:,:,:,1-gc,b) = 0._R8P ; q(:,:   ,:   ,nk+gc,b) = 0._R8P
      enddo
   enddo
   !do b=1, blocks_number
   !   do gc = 1, ngc
   !   q(:,1-gc,:,:,b) = -q(:,gc,:,:,b) ; q(:,ni+gc,:   ,:   ,b) = -q(:,ni-gc+1,:   ,:   ,b)
   !   q(:,:,1-gc,:,b) = -q(:,:,gc,:,b) ; q(:,:   ,nj+gc,:   ,b) = -q(:,:   ,nj-gc+1,:   ,b)
   !   q(:,:,:,1-gc,b) = -q(:,:,:,gc,b) ; q(:,:   ,:   ,nk+gc,b) = -q(:,:   ,:   ,nk-gc+1,b)
   !   enddo
   !enddo
   endsubroutine apply_bc_dirichlet

   subroutine apply_bc_neumann(ni, nj, nk, ngc, blocks_number, q)
   integer(I4P), intent(in)    :: ni,nj,nk,ngc  !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number !< Number of current blocks.
   real(R8P),    intent(inout) :: q(1:,    &
                                    1-ngc:,&
                                    1-ngc:,&
                                    1-ngc:,&
                                    1:)         !< Field variables.
   integer(I4P)                :: b, gc         !< Counter.

   do b=1, blocks_number
      do gc = 1, ngc
      q(:,1-gc,:,:,b) = q(:,gc,:,:,b) ; q(:,ni+gc,:   ,:   ,b) = q(:,ni-gc+1,:   ,:   ,b)
      q(:,:,1-gc,:,b) = q(:,:,gc,:,b) ; q(:,:   ,nj+gc,:   ,b) = q(:,:   ,nj-gc+1,:   ,b)
      q(:,:,:,1-gc,b) = q(:,:,:,gc,b) ; q(:,:   ,:   ,nk+gc,b) = q(:,:   ,:   ,nk-gc+1,b)
      enddo
   enddo
   endsubroutine apply_bc_neumann

   subroutine apply_bc_analytic(ni, nj, nk, ngc, blocks_number, ivar, mu, eps, field, rho, current, q)
   integer(I4P),       intent(in)              :: ni,nj,nk,ngc                 !< Grid dimensions.
   integer(I4P),       intent(in)              :: blocks_number                !< Number of current blocks.
   integer(I4P),       intent(in)              :: ivar                         !< Variable (start) index in q.
   real(R8P),          intent(in), optional    :: mu, eps
   type(field_object), intent(in)              :: field                        !< Field (sibling realm component, threaded in).
   real(R8P),          intent(in), optional    :: rho(1:,     &
                                                      1-ngc:, &
                                                      1-ngc:, &
                                                      1-ngc:, &
                                                      1:)                      !< Source term for scalar potential.
   real(R8P),          intent(in), optional    :: current(1:, &
                                                      1-ngc:, &
                                                      1-ngc:, &
                                                      1-ngc:, &
                                                      1:)                      !< Source term for vector potential.
   real(R8P),          intent(inout)           :: q(1:,       &
                                                    1-ngc:,   &
                                                    1-ngc:,   &
                                                    1-ngc:,   &
                                                    1:)                        !< Field variables.
   integer(I4P)                                :: i,j,k,b,f                    !< Counter.
   integer(I4P)                                :: i_f,j_f,k_f,b_f              !< Counter.
   real(R8P)                                   :: gc_coord(3)                  !< Ghost cell coordinates
   real(R8P)                                   :: cell_coord(3)                !< Cell coordinates
   integer(I4P)                                :: i_dir_n, i_dir_a, i_dir_b    !< Normal and tangential directions indices for ghost cell reconstruction
   integer(I4P), parameter                     :: n_faces = 6_I4P              !< Number of ghost-cell face slabs.
   integer(I4P)                                :: i1_f(n_faces), i2_f(n_faces)
   integer(I4P)                                :: j1_f(n_faces), j2_f(n_faces)
   integer(I4P)                                :: k1_f(n_faces), k2_f(n_faces)

   associate(x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell, &
               dx=>field%dxyz(1,:), dy=>field%dxyz(2,:), dz=>field%dxyz(3,:))

   i_dir_n = 1_I4P
   i_dir_a = 2_I4P
   i_dir_b = 3_I4P

   ! Faccia -x
   i1_f(1) = 1_I4P - ngc
   i2_f(1) = 0_I4P
   j1_f(1) = 1_I4P
   j2_f(1) = nj
   k1_f(1) = 1_I4P
   k2_f(1) = nk

   ! Faccia +x
   i1_f(2) = ni    + 1_I4P
   i2_f(2) = ni    + ngc
   j1_f(2) = 1_I4P
   j2_f(2) = nj
   k1_f(2) = 1_I4P
   k2_f(2) = nk

   ! Faccia -y
   i1_f(3) = 1_I4P
   i2_f(3) = ni
   j1_f(3) = 1_I4P - ngc
   j2_f(3) = 0_I4P
   k1_f(3) = 1_I4P
   k2_f(3) = nk

   ! Faccia +y
   i1_f(4) = 1_I4P
   i2_f(4) = ni
   j1_f(4) = nj    + 1_I4P
   j2_f(4) = nj    + ngc
   k1_f(4) = 1_I4P
   k2_f(4) = nk

   ! Faccia -z
   i1_f(5) = 1_I4P
   i2_f(5) = ni
   j1_f(5) = 1_I4P
   j2_f(5) = nj
   k1_f(5) = 1_I4P - ngc
   k2_f(5) = 0_I4P

   ! Faccia +z
   i1_f(6) = 1_I4P
   i2_f(6) = ni
   j1_f(6) = 1_I4P
   j2_f(6) = nj
   k1_f(6) = nk    + 1_I4P
   k2_f(6) = nk    + ngc
  
   if(ivar==1_I4P) then
      do f = 1_I4P, n_faces
         do b_f = 1, blocks_number
            do k_f = k1_f(f), k2_f(f)
               do j_f = j1_f(f), j2_f(f)
                  do i_f = i1_f(f), i2_f(f)
                     gc_coord = [x_cell(i_f,b_f), y_cell(j_f,b_f), z_cell(k_f,b_f)]
                     q(1, i_f, j_f, k_f, b_f) = 0.0_R8P
                     do b=1, blocks_number
                     do k = 1, nk
                     do j = 1, nj
                     do i = 1, ni
                        cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
                        q(1, i_f, j_f, k_f, b_f) = q(1, i_f, j_f, k_f, b_f) + 1/(4* acos(-1.0)*eps)*(rho(1,i,j,k,b)* &
                                                   (dx(b)*dy(b)*dz(b)))/sqrt((gc_coord(1)-cell_coord(1))**2 +        & 
                                                   (gc_coord(2)-cell_coord(2))**2 + (gc_coord(3)-cell_coord(3))**2)
                     enddo
                     enddo
                     enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   elseif(ivar==4_I4P) then
      do f = 1_I4P, n_faces
         do b_f = 1, blocks_number
            do k_f = k1_f(f), k2_f(f)
               do j_f = j1_f(f), j2_f(f)
                  do i_f = i1_f(f), i2_f(f)
                     gc_coord = [x_cell(i_f,b_f), y_cell(j_f,b_f), z_cell(k_f,b_f)]
                     q(1, i_f, j_f, k_f, b_f) = 0.0_R8P
                     q(2, i_f, j_f, k_f, b_f) = 0.0_R8P
                     q(3, i_f, j_f, k_f, b_f) = 0.0_R8P
                     do b=1, blocks_number
                     do k = 1, nk
                     do j = 1, nj
                     do i = 1, ni
                        cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
                        q(1, i_f, j_f, k_f, b_f) = q(1, i_f, j_f, k_f, b_f) + mu/(4* acos(-1.0))*current(1,i,j,k,b) * &
                                                   (dy(b)*dz(b))*dx(b)/sqrt((gc_coord(1)-cell_coord(1))**2 +          & 
                                                   (gc_coord(2)-cell_coord(2))**2 + (gc_coord(3)-cell_coord(3))**2)
                        q(2, i_f, j_f, k_f, b_f) = q(2, i_f, j_f, k_f, b_f) + mu/(4* acos(-1.0))*current(2,i,j,k,b) * &
                                                   (dx(b)*dz(b))*dy(b)/sqrt((gc_coord(1)-cell_coord(1))**2 +          & 
                                                   (gc_coord(2)-cell_coord(2))**2 + (gc_coord(3)-cell_coord(3))**2)
                        q(3, i_f, j_f, k_f, b_f) = q(3, i_f, j_f, k_f, b_f) + mu/(4* acos(-1.0))*current(3,i,j,k,b) * &
                                                   (dx(b)*dy(b))*dz(b)/sqrt((gc_coord(1)-cell_coord(1))**2 +          & 
                                                   (gc_coord(2)-cell_coord(2))**2 + (gc_coord(3)-cell_coord(3))**2)
                     enddo
                     enddo
                     enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   endif
   endassociate
   endsubroutine apply_bc_analytic

   subroutine apply_bc_elliptic_from_faces(ni, nj, nk, ngc, blocks_number, ell_bc_type, local_map_bc_crown, ivar, mu, eps, field, &
                                           f, q, rebuild_exact_open)
   !< Apply elliptic BCs using the face/crown map built for the hyperbolic BC machinery.
   integer(I4P),       intent(in)    :: ni,nj,nk,ngc          !< Grid dimensions.
   integer(I4P),       intent(in)    :: blocks_number         !< Number of current blocks.
   integer(I4P),       intent(in)    :: ell_bc_type(6)        !< Elliptic BC types on the 6 faces.
   integer(I8P),       intent(in)    :: local_map_bc_crown(:,:,:) !< Face BC crown map.
   integer(I4P),       intent(in)    :: ivar                  !< Variable (start) index in q.
   real(R8P),          intent(in), optional :: mu, eps        !< Electromagnetic constants.
   type(field_object), intent(in)    :: field                 !< Field (sibling realm component, threaded in).
   real(R8P),          intent(in)    :: f(1:,     &
                                          1-ngc:, &
                                          1-ngc:, &
                                          1-ngc:, &
                                          1:)                 !< Forcing distribution.
   real(R8P),          intent(inout) :: q(1:,     &
                                          1-ngc:, &
                                          1-ngc:, &
                                          1-ngc:, &
                                          1:)                 !< Field variables.
   logical,            intent(in), optional :: rebuild_exact_open !< Rebuild exact/open BC values.
   integer(I4P)                      :: crown                  !< Crown counter.
   integer(I4P)                      :: c                      !< Map entry counter.
   integer(I4P)                      :: b, i, j, k            !< Grid indexes.
   integer(I4P)                      :: face                   !< Face index in 1:6 numbering.
   integer(I4P)                      :: ip, jp, kp            !< Periodic source indexes.
   logical                           :: rebuild_exact_open_    !< Rebuild exact/open BC values, local var.

   rebuild_exact_open_ = .true.
   if (present(rebuild_exact_open)) rebuild_exact_open_ = rebuild_exact_open

   do crown=1, ngc
      do c=1, size(local_map_bc_crown, dim=1)
         b = int(local_map_bc_crown(c, 1, crown), I4P)
         if (b <= 0_I4P) cycle
         i = int(local_map_bc_crown(c, 2, crown), I4P)
         j = int(local_map_bc_crown(c, 3, crown), I4P)
         k = int(local_map_bc_crown(c, 4, crown), I4P)
         face = FEC_1_6_ARRAY(int(local_map_bc_crown(c, 9, crown), I4P))
         select case(ell_bc_type(face))
         case(ELL_BC_DIRICHLET)
            q(:,i,j,k,b) = 0._R8P
         case(ELL_BC_PERIODIC)
            ip = i ; jp = j ; kp = k
            if (ip < 1_I4P ) ip = ip + ni
            if (ip > ni   ) ip = ip - ni
            if (jp < 1_I4P ) jp = jp + nj
            if (jp > nj   ) jp = jp - nj
            if (kp < 1_I4P ) kp = kp + nk
            if (kp > nk   ) kp = kp - nk
            q(:,i,j,k,b) = q(:,ip,jp,kp,b)
         case(ELL_BC_EXACT_OPEN)
            if (rebuild_exact_open_) then
               call apply_bc_elliptic_face_exact_open(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, face=face, &
                                                      ivar=ivar, mu=mu, eps=eps, field=field, f=f, i_gc=i, j_gc=j, k_gc=k, &
                                                      b_gc=b, q=q)
            endif
         case default
            call mpih%error_stop(msg='apply_bc_elliptic_from_faces: unsupported elliptic face BC '//trim(str(ell_bc_type(face))))
         endselect
      enddo
   enddo
   endsubroutine apply_bc_elliptic_from_faces

   subroutine apply_bc_elliptic_face_exact_open(ni, nj, nk, ngc, blocks_number, face, ivar, mu, eps, field, f, i_gc, j_gc, k_gc, &
                                                b_gc, q)
   !< Reconstruct one exact/open elliptic ghost cell from the volumetric source distribution.
   integer(I4P),       intent(in)    :: ni,nj,nk,ngc          !< Grid dimensions.
   integer(I4P),       intent(in)    :: blocks_number         !< Number of current blocks.
   integer(I4P),       intent(in)    :: face                  !< Face index in 1:6 numbering.
   integer(I4P),       intent(in)    :: ivar                  !< Variable (start) index in q.
   real(R8P),          intent(in), optional :: mu, eps        !< Electromagnetic constants.
   type(field_object), intent(in)    :: field                 !< Field (sibling realm component, threaded in).
   real(R8P),          intent(in)    :: f(1:,     &
                                          1-ngc:, &
                                          1-ngc:, &
                                          1-ngc:, &
                                          1:)                 !< Forcing distribution.
   integer(I4P),       intent(in)    :: i_gc, j_gc, k_gc, b_gc !< Ghost-cell indexes.
   real(R8P),          intent(inout) :: q(1:,     &
                                          1-ngc:, &
                                          1-ngc:, &
                                          1-ngc:, &
                                          1:)                 !< Field variables.
   integer(I4P)                      :: i,j,k,b              !< Counters.
   real(R8P)                         :: gc_coord(3)          !< Ghost-cell coordinates.
   real(R8P)                         :: cell_coord(3)        !< Interior-cell coordinates.
   real(R8P)                         :: distance             !< Distance between source and target cell centers.

   if (face < 1_I4P .or. face > 6_I4P) then
      call mpih%error_stop(msg='apply_bc_elliptic_face_exact_open: invalid face index '//trim(str(face)))
   endif

   associate(x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell, &
             dx=>field%dxyz(1,:), dy=>field%dxyz(2,:), dz=>field%dxyz(3,:))
   gc_coord = [x_cell(i_gc,b_gc), y_cell(j_gc,b_gc), z_cell(k_gc,b_gc)]
   q(:,i_gc,j_gc,k_gc,b_gc) = 0._R8P
   if (ivar == 1_I4P) then
      if (.not. present(eps)) call mpih%error_stop(msg='apply_bc_elliptic_face_exact_open: eps missing for scalar potential')
      do b=1, blocks_number
         do k=1, nk
         do j=1, nj
         do i=1, ni
            cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
            distance = sqrt((gc_coord(1)-cell_coord(1))**2 + (gc_coord(2)-cell_coord(2))**2 + (gc_coord(3)-cell_coord(3))**2)
            q(1,i_gc,j_gc,k_gc,b_gc) = q(1,i_gc,j_gc,k_gc,b_gc) - &
                                        (f(1,i,j,k,b) * dx(b) * dy(b) * dz(b)) / (4._R8P * acos(-1._R8P) * distance)
         enddo
         enddo
         enddo
      enddo
   elseif (ivar == 4_I4P) then
      if (.not. present(mu)) call mpih%error_stop(msg='apply_bc_elliptic_face_exact_open: mu missing for vector potential')
      do b=1, blocks_number
         do k=1, nk
         do j=1, nj
         do i=1, ni
            cell_coord = [x_cell(i,b), y_cell(j,b), z_cell(k,b)]
            distance = sqrt((gc_coord(1)-cell_coord(1))**2 + (gc_coord(2)-cell_coord(2))**2 + (gc_coord(3)-cell_coord(3))**2)
            q(1,i_gc,j_gc,k_gc,b_gc) = q(1,i_gc,j_gc,k_gc,b_gc) - f(1,i,j,k,b) * dx(b) * dy(b) * dz(b) / &
                                        (4._R8P * acos(-1._R8P) * distance)
            q(2,i_gc,j_gc,k_gc,b_gc) = q(2,i_gc,j_gc,k_gc,b_gc) - f(2,i,j,k,b) * dx(b) * dy(b) * dz(b) / &
                                        (4._R8P * acos(-1._R8P) * distance)
            q(3,i_gc,j_gc,k_gc,b_gc) = q(3,i_gc,j_gc,k_gc,b_gc) - f(3,i,j,k,b) * dx(b) * dy(b) * dz(b) / &
                                        (4._R8P * acos(-1._R8P) * distance)
         enddo
         enddo
         enddo
      enddo
   else
      call mpih%error_stop(msg='apply_bc_elliptic_face_exact_open: unsupported ivar '//trim(str(ivar)))
   endif
   endassociate
   endsubroutine apply_bc_elliptic_face_exact_open

   subroutine apply_bc_elliptic(ni, nj, nk, ngc, blocks_number, f, q, bc_type, ell_bc_type, local_map_bc_crown, ivar, mu, eps, &
                                field, rebuild_exact_open)
   !< Apply face-based elliptic BCs, or fall back to the legacy scalar BC selector when needed.
   integer(I4P),       intent(in)              :: ni,nj,nk,ngc      !< Grid dimensions.
   integer(I4P),       intent(in)              :: blocks_number     !< Number of current blocks.
   real(R8P),          intent(in)              :: f(1:,     &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1:)             !< Forcing distribution.
   real(R8P),          intent(inout)           :: q(1:,     &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1:)             !< Field variables.
   character(len=*),   intent(in),    optional :: bc_type           !< Legacy scalar BC selector.
   integer(I4P),       intent(in),    optional :: ell_bc_type(6)    !< Elliptic BC types on the 6 faces.
   integer(I8P),       intent(in),    optional :: local_map_bc_crown(:,:,:) !< BC crown map.
   integer(I4P),       intent(in),    optional :: ivar              !< Variable (start) index in q.
   real(R8P),          intent(in),    optional :: mu, eps           !< Electromagnetic constants.
   type(field_object), intent(in),    optional :: field             !< Field (sibling realm component, threaded in).
   logical,            intent(in),    optional :: rebuild_exact_open !< Rebuild exact/open BC values.
   logical                                     :: rebuild_exact_open_ !< Rebuild exact/open BC values, local var.

   rebuild_exact_open_ = .true.
   if (present(rebuild_exact_open)) rebuild_exact_open_ = rebuild_exact_open

   if (present(ell_bc_type)) then
      if (.not. present(local_map_bc_crown)) call mpih%error_stop(msg='apply_bc_elliptic: missing local_map_bc_crown')
      if (.not. present(ivar)) call mpih%error_stop(msg='apply_bc_elliptic: missing ivar for face-based elliptic BCs')
      if (.not. present(field)) call mpih%error_stop(msg='apply_bc_elliptic: missing field for face-based elliptic BCs')
      call apply_bc_elliptic_from_faces(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, ell_bc_type=ell_bc_type, &
                                        local_map_bc_crown=local_map_bc_crown, ivar=ivar, mu=mu, eps=eps, field=field, f=f, q=q, &
                                        rebuild_exact_open=rebuild_exact_open_)
   elseif (present(bc_type)) then
      if (bc_type == 'analytic') then
         if (rebuild_exact_open_) then
            if (.not. present(ivar)) call mpih%error_stop(msg='apply_bc_elliptic: missing ivar for analytic elliptic BCs')
            if (.not. present(field)) call mpih%error_stop(msg='apply_bc_elliptic: missing field for analytic elliptic BCs')
            if (ivar == 1_I4P) then
               call apply_bc_analytic(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, ivar=ivar, eps=eps, &
                                      field=field, rho=-f*eps, q=q)
            elseif (ivar == 4_I4P) then
               call apply_bc_analytic(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, ivar=ivar, mu=mu, &
                                      field=field, current=-f/mu, q=q)
            endif
         endif
      elseif (bc_type == 'dirichlet') then
         call apply_bc_dirichlet(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, q=q)
      elseif (bc_type == 'neumann') then
         call apply_bc_neumann(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, q=q)
      endif
   else
      call apply_bc_dirichlet(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, q=q)
   endif
   endsubroutine apply_bc_elliptic

   subroutine compute_prolongation(ngc, nic, njc, nkc, nv, blocks_number, coarse, fine)
   !< Compute trilinear prolongation.
   integer(I4P), intent(in)    :: ngc            !< Number of ghost cells.
   integer(I4P), intent(in)    :: nic,njc,nkc    !< Coarse grid dimensions.
   integer(I4P), intent(in)    :: nv             !< Number of q variables.
   integer(I4P), intent(in)    :: blocks_number  !< Number of current blocks.
   real(R8P),    intent(in)    :: coarse(1:,    &
                                         1-ngc:,&
                                         1-ngc:,&
                                         1-ngc:,&
                                         1:)     !< Coarse grid field variable.
   real(R8P),    intent(inout) :: fine(1:,    &
                                       1-ngc:,&
                                       1-ngc:,&
                                       1-ngc:,&
                                       1:)       !< Fine grid field variable.
   integer(I4P)                :: i,j,k,b,v      !< Counter.
   integer(I4P)                :: ii,jj,kk       !< Counter.

   do b=1, blocks_number
      fine(:,:,:,:,b) = 0._R8P
      do k = 1, nkc
      do j = 1, njc
      do i = 1, nic
         ii = 2*i - 1
         jj = 2*j - 1
         kk = 2*k - 1
         do v=1, nv
            fine(v,ii,  jj,  kk  ,b) =              coarse(v,i,  j,  k  ,b)
            fine(v,ii+1,jj,  kk  ,b) = 0.5_R8P   * (coarse(v,i,  j,  k  ,b) + coarse(v,i+1,j,  k  ,b))
            fine(v,ii,  jj+1,kk  ,b) = 0.5_R8P   * (coarse(v,i,  j,  k  ,b) + coarse(v,i,  j+1,k  ,b))
            fine(v,ii,  jj,  kk+1,b) = 0.5_R8P   * (coarse(v,i,  j,  k  ,b) + coarse(v,i,  j,  k+1,b))
            fine(v,ii+1,jj+1,kk  ,b) = 0.25_R8P  * (coarse(v,i,  j,  k  ,b) + coarse(v,i+1,j,  k  ,b) + &
                                                    coarse(v,i,  j+1,k  ,b) + coarse(v,i+1,j+1,k  ,b))
            fine(v,ii+1,jj,  kk+1,b) = 0.25_R8P  * (coarse(v,i,  j,  k  ,b) + coarse(v,i+1,j,  k  ,b) + &
                                                    coarse(v,i,  j,  k+1,b) + coarse(v,i+1,j,  k+1,b))
            fine(v,ii,  jj+1,kk+1,b) = 0.25_R8P  * (coarse(v,i,  j,  k  ,b) + coarse(v,i,  j+1,k  ,b) + &
                                                    coarse(v,i,  j,  k+1,b) + coarse(v,i,  j+1,k+1,b))
            fine(v,ii+1,jj+1,kk+1,b) = 0.125_R8P * (coarse(v,i,  j,  k  ,b) + coarse(v,i+1,j,  k  ,b) + &
                                                    coarse(v,i,  j+1,k  ,b) + coarse(v,i,  j,  k+1,b) + &
                                                    coarse(v,i+1,j+1,k  ,b) + coarse(v,i+1,j,  k+1,b) + &
                                                    coarse(v,i,  j+1,k+1,b) + coarse(v,i+1,j+1,k+1,b))
         enddo
      enddo
      enddo
      enddo
   enddo
   endsubroutine compute_prolongation

   subroutine compute_laplacian_residuals(ni, nj, nk, ngc, nv, blocks_number, dxyz, f, q, dq)
   !< Compute residuals of `f-laplacian(q)`.
   integer(I4P), intent(in)    :: ni,nj,nk,ngc  !< Grid dimensions.
   integer(I4P), intent(in)    :: nv            !< Number of q variables.
   integer(I4P), intent(in)    :: blocks_number !< Number of current blocks.
   real(R8P),    intent(in)    :: dxyz(1:,1:)   !< Space steps.
   real(R8P),    intent(in)    :: f(1:,    &
                                    1-ngc:,&
                                    1-ngc:,&
                                    1-ngc:,&
                                    1:)         !< Forcing distribution.
   real(R8P),    intent(in)    :: q(1:,    &
                                    1-ngc:,&
                                    1-ngc:,&
                                    1-ngc:,&
                                    1:)         !< Field variables.
   real(R8P),    intent(inout) :: dq(1:,    &
                                     1-ngc:,&
                                     1-ngc:,&
                                     1-ngc:,&
                                     1:)        !< Residuals.
   real(R8P)                   :: laplacian     !< Laplacian value.
   real(R8P)                   :: dx2,dy2,dz2   !< Square space steps.
   integer(I4P)                :: i,j,k,b,v     !< Counter.

   do b=1, blocks_number
      dx2 = dxyz(1,b)*dxyz(1,b)
      dy2 = dxyz(2,b)*dxyz(2,b)
      dz2 = dxyz(3,b)*dxyz(3,b)
      do k = 1, nk
      do j = 1, nj
      do i = 1, ni
         do v=1, nv
            laplacian = (q(v,i+1,j,  k  ,b) - 2._R8P * q(v,i,j,k,b) + q(v,i-1,j,  k  ,b)) / dx2 + &
                        (q(v,i,  j+1,k  ,b) - 2._R8P * q(v,i,j,k,b) + q(v,i,  j-1,k  ,b)) / dy2 + &
                        (q(v,i,  j,  k+1,b) - 2._R8P * q(v,i,j,k,b) + q(v,i,  j,  k-1,b)) / dz2
            dq(v,i,j,k,b) = f(v,i,j,k,b) - laplacian
         enddo
      enddo
      enddo
      enddo
   enddo
   call apply_bc_dirichlet(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, q=dq)
   endsubroutine compute_laplacian_residuals

   subroutine compute_restriction(ngc, nic, njc, nkc, nv, blocks_number, fine, coarse)
   !< Compute full weighting restriction.
   integer(I4P), intent(in)    :: ngc            !< Number of ghost cells.
   integer(I4P), intent(in)    :: nic,njc,nkc    !< Coarse grid dimensions.
   integer(I4P), intent(in)    :: nv             !< Number of q variables.
   integer(I4P), intent(in)    :: blocks_number  !< Number of current blocks.
   real(R8P),    intent(in)    :: fine(1:,    &
                                       1-ngc:,&
                                       1-ngc:,&
                                       1-ngc:,&
                                       1:)       !< Fine grid field variable.
   real(R8P),    intent(inout) :: coarse(1:,    &
                                         1-ngc:,&
                                         1-ngc:,&
                                         1-ngc:,&
                                         1:)     !< Coarse grid field variable.
   integer(I4P)                :: i,j,k,b,v      !< Counter.
   integer(I4P)                :: ii,jj,kk       !< Counter.

   do b=1, blocks_number
      do k = 1, nkc
      do j = 1, njc
      do i = 1, nic
         ii = 2*i - 1
         jj = 2*j - 1
         kk = 2*k - 1
         do v=1, nv
            coarse(v,i,j,k,b) = 0.125_R8P * (fine(v,ii,  jj,  kk  ,b) + fine(v,ii+1,jj,  kk  ,b) + &
                                             fine(v,ii,  jj+1,kk  ,b) + fine(v,ii,  jj,  kk+1,b) + &
                                             fine(v,ii+1,jj+1,kk  ,b) + fine(v,ii+1,jj,  kk+1,b) + &
                                             fine(v,ii,  jj+1,kk+1,b) + fine(v,ii+1,jj+1,kk+1,b))
         enddo
      enddo
      enddo
      enddo
   enddo
   call apply_bc_dirichlet(ni=nic, nj=njc, nk=nkc, ngc=ngc, blocks_number=blocks_number, q=coarse)
   endsubroutine compute_restriction

   subroutine compute_smoothing_gauss_seidel_centered_dg(ni, nj, nk, ngc, nv, blocks_number, order, dxyz, f, q, dq, dq_max, &
                                                         iterations_init, iterations_fine, iterations_coarse, bc_type, ivar,  &
                                                         mu, eps, field, ell_bc_type, local_map_bc_crown, progress_label,     &
                                                         progress_counter, progress_total, progress_last_percent)
   !< Compute smoothing by Gauss-Seidel for the centered `D(G)` elliptic operator.
   !<
   !< On the uniform Cartesian blocks used here, centered FD and centered FV share
   !< the same composed `D(G)` coefficients. The dispatch therefore depends only
   !< on the formal order, not on `fdv_scheme`.
   integer(I4P),       intent(in)              :: ni,nj,nk,ngc      !< Grid dimensions.
   integer(I4P),       intent(in)              :: nv                !< Number of q variables.
   integer(I4P),       intent(in)              :: blocks_number     !< Number of current blocks.
   integer(I4P),       intent(in)              :: order             !< Formal order of the centered `D(G)` operator.
   real(R8P),          intent(in)              :: dxyz(1:,1:)       !< Space steps.
   real(R8P),          intent(in)              :: f(1:,     &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1:)             !< Forcing distribution.
   real(R8P),          intent(inout)           :: q(1:,     &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1:)             !< Field variables.
   real(R8P),          intent(inout), optional :: dq(1:,     &
                                                     1-ngc:, &
                                                     1-ngc:, &
                                                     1-ngc:, &
                                                     1:)            !< Residuals.
   real(R8P),          intent(inout), optional :: dq_max            !< Maximum residual.
   integer(I4P),       intent(in),    optional :: iterations_init   !< Smoothing iterations to initialize guess.
   integer(I4P),       intent(in),    optional :: iterations_fine   !< Smoothing iterations for fine grid.
   integer(I4P),       intent(in),    optional :: iterations_coarse !< Smoothing iterations for coarse grid.
   character(len=*),   intent(in),    optional :: bc_type           !< Boundary condition type.
   integer(I4P),       intent(in),    optional :: ivar              !< Variable (start) index in q.
   real(R8P),          intent(in),    optional :: mu, eps           !< Electromagnetic constants.
   type(field_object), intent(in),    optional :: field             !< Field (sibling realm component, threaded in).
   integer(I4P),       intent(in),    optional :: ell_bc_type(6)    !< Elliptic BC type for each face.
   integer(I8P),       intent(in),    optional :: local_map_bc_crown(:,:,:) !< BC crown map.
   character(len=*),   intent(in),    optional :: progress_label    !< Progress message label.
   integer(I4P),       intent(inout), optional :: progress_counter  !< Completed smoothing sweeps.
   integer(I4P),       intent(in),    optional :: progress_total    !< Planned smoothing sweeps.
   integer(I4P),       intent(inout), optional :: progress_last_percent !< Last printed progress percentage.

   select case(order)
   case(2_I4P)
      call compute_smoothing_gauss_seidel_2nd(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                              dxyz=dxyz, f=f, q=q, dq=dq, dq_max=dq_max,                         &
                                              iterations_init=iterations_init, iterations_fine=iterations_fine,   &
                                              iterations_coarse=iterations_coarse, bc_type=bc_type, ivar=ivar,   &
                                              mu=mu, eps=eps, field=field, ell_bc_type=ell_bc_type,              &
                                              local_map_bc_crown=local_map_bc_crown, progress_label=progress_label, &
                                              progress_counter=progress_counter, progress_total=progress_total,     &
                                              progress_last_percent=progress_last_percent)
   case(4_I4P)
      call compute_smoothing_gauss_seidel_4th(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                              dxyz=dxyz, f=f, q=q, dq=dq, dq_max=dq_max,                         &
                                              iterations_init=iterations_init, iterations_fine=iterations_fine,   &
                                              iterations_coarse=iterations_coarse, bc_type=bc_type, ivar=ivar,   &
                                              mu=mu, eps=eps, field=field, ell_bc_type=ell_bc_type,              &
                                              local_map_bc_crown=local_map_bc_crown, progress_label=progress_label, &
                                              progress_counter=progress_counter, progress_total=progress_total,     &
                                              progress_last_percent=progress_last_percent)
   case(6_I4P)
      call compute_smoothing_gauss_seidel_6th(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                              dxyz=dxyz, f=f, q=q, dq=dq, dq_max=dq_max,                         &
                                              iterations_init=iterations_init, iterations_fine=iterations_fine,   &
                                              iterations_coarse=iterations_coarse, bc_type=bc_type, ivar=ivar,   &
                                              mu=mu, eps=eps, field=field, ell_bc_type=ell_bc_type,              &
                                              local_map_bc_crown=local_map_bc_crown, progress_label=progress_label, &
                                              progress_counter=progress_counter, progress_total=progress_total,     &
                                              progress_last_percent=progress_last_percent)
   case(8_I4P)
      call compute_smoothing_gauss_seidel_8th(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                              dxyz=dxyz, f=f, q=q, dq=dq, dq_max=dq_max,                         &
                                              iterations_init=iterations_init, iterations_fine=iterations_fine,   &
                                              iterations_coarse=iterations_coarse, bc_type=bc_type, ivar=ivar,   &
                                              mu=mu, eps=eps, field=field, ell_bc_type=ell_bc_type,              &
                                              local_map_bc_crown=local_map_bc_crown, progress_label=progress_label, &
                                              progress_counter=progress_counter, progress_total=progress_total,     &
                                              progress_last_percent=progress_last_percent)
   case default
      call mpih%error_stop(msg='compute_smoothing_gauss_seidel_centered_dg: unsupported FDV order '//trim(str(order)))
   endselect
   endsubroutine compute_smoothing_gauss_seidel_centered_dg

   subroutine compute_smoothing_gauss_seidel(ni, nj, nk, ngc, nv, blocks_number, dxyz, f, q, dq, dq_max, &
                                             iterations_init, iterations_fine, iterations_coarse, bc_type, mu, eps, ivar, field, &
                                             ell_bc_type, local_map_bc_crown)
   !< Compute smoothing by Gauss-Seidel method.
   integer(I4P),       intent(in)              :: ni,nj,nk,ngc      !< Grid dimensions.
   integer(I4P),       intent(in)              :: nv                !< Number of q variables.
   integer(I4P),       intent(in)              :: blocks_number     !< Number of current blocks.
   real(R8P),          intent(in)              :: dxyz(1:,1:)       !< Space steps.
   real(R8P),          intent(in)              :: f(1:,     &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1:)             !< Forcing distribution.
   real(R8P),          intent(inout)           :: q(1:,     &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1:)             !< Field variables.
   real(R8P),          intent(inout), optional :: dq(1:,     &
                                                     1-ngc:, &
                                                     1-ngc:, &
                                                     1-ngc:, &
                                                     1:)            !< Residuals.
   real(R8P),          intent(inout), optional :: dq_max            !< Maximum residual.
   integer(I4P),       intent(in),    optional :: iterations_init   !< Smoothing iterations to initialize guess.
   integer(I4P),       intent(in),    optional :: iterations_fine   !< Smoothing iterations for fine grid.
   integer(I4P),       intent(in),    optional :: iterations_coarse !< Smoothing iterations for coarse grid.
   character(len=*),   intent(in),    optional :: bc_type           !< Boundary condition type
   real(R8P),          intent(in),    optional :: mu, eps           !< Constant of electromagnetism, =1 in adimensional case
   integer(I4P),       intent(in),    optional :: ivar              !< Variable (start) index in q.
   type(field_object), intent(in),    optional :: field             !< Field (sibling realm component, threaded in).
   integer(I4P),       intent(in),    optional :: ell_bc_type(6)    !< Elliptic BC type for each face.
   integer(I8P),       intent(in),    optional :: local_map_bc_crown(:,:,:) !< BC crown map.
   integer(I4P)                                :: iterations_       !< Smoothing iterations, local var.
   real(R8P)                                   :: dq_max_           !< Maximum residual, local var.
   real(R8P)                                   :: dx2,dy2,dz2       !< Square space steps.
   real(R8P)                                   :: q_old             !< Previous q.
   real(R8P)                                   :: factor            !< Jacobi relaxation factor.
   integer(I4P)                                :: i,j,k,b,v,iter    !< Counter.

   iterations_ = 1 ; if (present(iterations_fine)) iterations_ = iterations_fine
   dq_max_ = 0._R8P
   call apply_bc_elliptic(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, f=f, q=q, bc_type=bc_type, &
                          ell_bc_type=ell_bc_type, local_map_bc_crown=local_map_bc_crown, ivar=ivar, mu=mu, eps=eps, &
                          field=field, rebuild_exact_open=.true.)
   do iter=1, iterations_
      call apply_bc_elliptic(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, f=f, q=q, bc_type=bc_type, &
                             ell_bc_type=ell_bc_type, local_map_bc_crown=local_map_bc_crown, ivar=ivar, mu=mu, eps=eps, &
                             field=field, rebuild_exact_open=.false.)
      do b=1, blocks_number
         dx2 = dxyz(1,b)*dxyz(1,b)
         dy2 = dxyz(2,b)*dxyz(2,b)
         dz2 = dxyz(3,b)*dxyz(3,b)
         factor = 1._R8P / (2._R8P * (1._R8P/dx2 + 1._R8P/dy2 + 1._R8P/dz2))
         do k = 1, nk
         do j = 1, nj
         do i = 1, ni
            do v=1, nv
               q_old = q(v,i,j,k,b)
               q(v,i,j,k,b) = factor * ((q(v,i+1,j,  k  ,b) + q(v,i-1,j,  k  ,b)) / dx2 + &
                                        (q(v,i,  j+1,k  ,b) + q(v,i,  j-1,k  ,b)) / dy2 + &
                                        (q(v,i,  j,  k+1,b) + q(v,i,  j,  k-1,b)) / dz2 - &
                                         f(v,i,  j,  k  ,b))
               dq_max_ = max(dq_max_, abs(q(v,i,j,k,b) - q_old))
            enddo
         enddo
         enddo
         enddo
      enddo
      !if (present(bc_type)) then
      !   if (bc_type == 'dirichlet') then
      !      call apply_bc_dirichlet(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, q=q)
      !   elseif (bc_type == 'neumann') then
      !      call apply_bc_neumann(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, q=q)
      !   endif
      !else
      !   call apply_bc_dirichlet(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, q=q)
      !endif
   enddo

   if (present(dq_max)) dq_max = dq_max_
   endsubroutine compute_smoothing_gauss_seidel

   subroutine compute_smoothing_gauss_seidel_2nd(ni, nj, nk, ngc, nv, blocks_number, dxyz, f, q, dq, dq_max, &
                                             iterations_init, iterations_fine, iterations_coarse, bc_type, ivar, mu, eps, field, &
                                             ell_bc_type, local_map_bc_crown, progress_label, progress_counter, progress_total, &
                                             progress_last_percent)
   !< Compute smoothing by Gauss-Seidel using L = D(G), with second-order
   !< centered first-derivative operators.
   integer(I4P),       intent(in)              :: ni,nj,nk,ngc               !< Grid dimensions.
   integer(I4P),       intent(in)              :: nv                         !< Number of q variables.
   integer(I4P),       intent(in)              :: blocks_number              !< Number of current blocks.
   real(R8P),          intent(in)              :: dxyz(1:,1:)                !< Space steps.
   real(R8P),          intent(in)              :: f(1:,     &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1:)                      !< Forcing distribution.
   real(R8P),          intent(inout)           :: q(1:,     &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1:)                      !< Field variables.
   real(R8P),          intent(inout), optional :: dq(1:,     &
                                                     1-ngc:, &
                                                     1-ngc:, &
                                                     1-ngc:, &
                                                     1:)                     !< Residuals.
   real(R8P),          intent(inout), optional :: dq_max                     !< Maximum residual.
   integer(I4P),       intent(in),    optional :: iterations_init            !< Smoothing iterations to initialize guess.
   integer(I4P),       intent(in),    optional :: iterations_fine            !< Smoothing iterations for fine grid.
   integer(I4P),       intent(in),    optional :: iterations_coarse          !< Smoothing iterations for coarse grid.
   character(len=*),   intent(in),    optional :: bc_type                    !< Boundary condition type
   real(R8P),          intent(in),    optional :: mu, eps                    !< Constant of electromagnetism, =1 in adimensional case
   integer(I4P),       intent(in),    optional :: ivar                       !< Variable (start) index in q.
   type(field_object), intent(in),    optional :: field                      !< Field (sibling realm component, threaded in).
   integer(I4P),       intent(in),    optional :: ell_bc_type(6)             !< Elliptic BC type for each face.
   integer(I8P),       intent(in),    optional :: local_map_bc_crown(:,:,:)  !< BC crown map.
   character(len=*),   intent(in),    optional :: progress_label             !< Progress message label.
   integer(I4P),       intent(inout), optional :: progress_counter           !< Completed smoothing sweeps.
   integer(I4P),       intent(in),    optional :: progress_total             !< Planned smoothing sweeps.
   integer(I4P),       intent(inout), optional :: progress_last_percent      !< Last printed progress percentage.
   real(R8P), parameter                        :: c2 =   1._R8P /    4._R8P  !< Constant for laplacian computation
   integer(I4P)                                :: iterations_                !< Smoothing iterations, local var.
   real(R8P)                                   :: dq_max_                    !< Maximum residual, local var.
   real(R8P)                                   :: dx2,dy2,dz2                !< Square space steps.
   real(R8P)                                   :: idx2,idy2,idz2             !< Inverse square space steps.
   real(R8P)                                   :: q_old                      !< Previous q.
   real(R8P)                                   :: factor                     !< Relaxation factor.
   integer(I4P)                                :: i,j,k,b,v,iter             !< Counter.

   if (ngc < 2_I4P) error stop 'compute_smoothing_gauss_seidel_2nd: ngc must be >= 2 for second-order D(G) stencil'

   iterations_ = 1_I4P; if (present(iterations_fine)) iterations_ = iterations_fine
   dq_max_ = 0._R8P
   call apply_bc_elliptic(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, f=f, q=q, bc_type=bc_type, &
                          ell_bc_type=ell_bc_type, local_map_bc_crown=local_map_bc_crown, ivar=ivar, mu=mu, eps=eps, &
                          field=field, rebuild_exact_open=.true.)
   do iter=1, iterations_
      call apply_bc_elliptic(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, f=f, q=q, bc_type=bc_type, &
                             ell_bc_type=ell_bc_type, local_map_bc_crown=local_map_bc_crown, ivar=ivar, mu=mu, eps=eps, &
                             field=field, rebuild_exact_open=.false.)
      do b=1, blocks_number
         dx2 = dxyz(1,b)*dxyz(1,b)
         dy2 = dxyz(2,b)*dxyz(2,b)
         dz2 = dxyz(3,b)*dxyz(3,b)
         idx2 = 1._R8P / dx2
         idy2 = 1._R8P / dy2
         idz2 = 1._R8P / dz2
         factor = 2._R8P / (idx2 + idy2 + idz2)
         do k=1, nk
         do j=1, nj
         do i=1, ni
            do v=1, nv
               q_old = q(v,i,j,k,b)
               q(v,i,j,k,b) = factor * (idx2 * (c2 * (q(v,i+2,j,  k  ,b) + q(v,i-2,j,  k  ,b))) + &
                                        idy2 * (c2 * (q(v,i,  j+2,k  ,b) + q(v,i,  j-2,k  ,b))) + &
                                        idz2 * (c2 * (q(v,i,  j,  k+2,b) + q(v,i,  j,  k-2,b))) - f(v,i,j,k,b))
               dq_max_ = max(dq_max_, abs(q(v,i,j,k,b) - q_old))
            enddo
         enddo
         enddo
         enddo
      enddo
      call apply_bc_elliptic(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, f=f, q=q, bc_type=bc_type, &
                             ell_bc_type=ell_bc_type, local_map_bc_crown=local_map_bc_crown, ivar=ivar, mu=mu, eps=eps, &
                             field=field, rebuild_exact_open=.false.)
      call update_smoothing_progress(progress_label=progress_label, progress_counter=progress_counter,            &
                                     progress_total=progress_total, progress_last_percent=progress_last_percent)
   enddo

   if (present(dq_max)) dq_max = dq_max_
   endsubroutine compute_smoothing_gauss_seidel_2nd

   subroutine compute_smoothing_gauss_seidel_4th(ni, nj, nk, ngc, nv, blocks_number, dxyz, f, q, dq, dq_max, &
                                             iterations_init, iterations_fine, iterations_coarse, bc_type, ivar, mu, eps, field, &
                                             ell_bc_type, local_map_bc_crown, progress_label, progress_counter, progress_total, &
                                             progress_last_percent)
   !< Compute smoothing by Gauss-Seidel using L = D(G), with fourth-order
   !< centered first-derivative operators.
   integer(I4P),       intent(in)              :: ni,nj,nk,ngc               !< Grid dimensions.
   integer(I4P),       intent(in)              :: nv                         !< Number of q variables.
   integer(I4P),       intent(in)              :: blocks_number              !< Number of current blocks.
   real(R8P),          intent(in)              :: dxyz(1:,1:)                !< Space steps.
   real(R8P),          intent(in)              :: f(1:,     &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1:)                      !< Forcing distribution.
   real(R8P),          intent(inout)           :: q(1:,     &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1:)                      !< Field variables.
   real(R8P),          intent(inout), optional :: dq(1:,     &
                                                     1-ngc:, &
                                                     1-ngc:, &
                                                     1-ngc:, &
                                                     1:)                     !< Residuals.
   real(R8P),          intent(inout), optional :: dq_max                     !< Maximum residual.
   integer(I4P),       intent(in),    optional :: iterations_init            !< Smoothing iterations to initialize guess.
   integer(I4P),       intent(in),    optional :: iterations_fine            !< Smoothing iterations for fine grid.
   integer(I4P),       intent(in),    optional :: iterations_coarse          !< Smoothing iterations for coarse grid.
   character(len=*),   intent(in),    optional :: bc_type                    !< Boundary condition type
   real(R8P),          intent(in),    optional :: mu, eps                    !< Constant of electromagnetism, =1 in adimensional case
   integer(I4P),       intent(in),    optional :: ivar                       !< Variable (start) index in q.
   type(field_object), intent(in),    optional :: field                      !< Field (sibling realm component, threaded in).
   integer(I4P),       intent(in),    optional :: ell_bc_type(6)             !< Elliptic BC type for each face.
   integer(I8P),       intent(in),    optional :: local_map_bc_crown(:,:,:)  !< BC crown map.
   character(len=*),   intent(in),    optional :: progress_label             !< Progress message label.
   integer(I4P),       intent(inout), optional :: progress_counter           !< Completed smoothing sweeps.
   integer(I4P),       intent(in),    optional :: progress_total             !< Planned smoothing sweeps.
   integer(I4P),       intent(inout), optional :: progress_last_percent      !< Last printed progress percentage.
   real(R8P), parameter                        :: c1 =   1._R8P /   9._R8P  !< Constant for laplacian computation
   real(R8P), parameter                        :: c2 =   4._R8P /   9._R8P  !< Constant for laplacian computation
   real(R8P), parameter                        :: c3 =  -1._R8P /   9._R8P  !< Constant for laplacian computation
   real(R8P), parameter                        :: c4 =   1._R8P / 144._R8P  !< Constant for laplacian computation
   integer(I4P)                                :: iterations_                !< Smoothing iterations, local var.
   real(R8P)                                   :: dq_max_                    !< Maximum residual, local var.
   real(R8P)                                   :: dx2,dy2,dz2                !< Square space steps.
   real(R8P)                                   :: idx2,idy2,idz2             !< Inverse square space steps.
   real(R8P)                                   :: q_old                      !< Previous q.
   real(R8P)                                   :: factor                     !< Relaxation factor.
   integer(I4P)                                :: i,j,k,b,v,iter             !< Counter.

   if (ngc < 4_I4P) error stop 'compute_smoothing_gauss_seidel_4th: ngc must be >= 4 for fourth-order D(G) stencil'

   iterations_ = 1_I4P; if (present(iterations_fine)) iterations_ = iterations_fine
   dq_max_ = 0._R8P
   call apply_bc_elliptic(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, f=f, q=q, bc_type=bc_type, &
                          ell_bc_type=ell_bc_type, local_map_bc_crown=local_map_bc_crown, ivar=ivar, mu=mu, eps=eps, &
                          field=field, rebuild_exact_open=.true.)
   do iter=1, iterations_
      call apply_bc_elliptic(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, f=f, q=q, bc_type=bc_type, &
                             ell_bc_type=ell_bc_type, local_map_bc_crown=local_map_bc_crown, ivar=ivar, mu=mu, eps=eps, &
                             field=field, rebuild_exact_open=.false.)
      do b=1, blocks_number
         dx2 = dxyz(1,b)*dxyz(1,b)
         dy2 = dxyz(2,b)*dxyz(2,b)
         dz2 = dxyz(3,b)*dxyz(3,b)
         idx2 = 1._R8P / dx2
         idy2 = 1._R8P / dy2
         idz2 = 1._R8P / dz2
         factor = 72._R8P / (65._R8P * (idx2 + idy2 + idz2))
         do k=1, nk
         do j=1, nj
         do i=1, ni
            do v=1, nv
               q_old = q(v,i,j,k,b)
               q(v,i,j,k,b) = factor * (idx2 * (c1 * (q(v,i+1,j,  k  ,b) + q(v,i-1,j,  k  ,b))  + &
                                                c2 * (q(v,i+2,j,  k  ,b) + q(v,i-2,j,  k  ,b))  + &
                                                c3 * (q(v,i+3,j,  k  ,b) + q(v,i-3,j,  k  ,b))  + &
                                                c4 * (q(v,i+4,j,  k  ,b) + q(v,i-4,j,  k  ,b))) + &
                                        idy2 * (c1 * (q(v,i,  j+1,k  ,b) + q(v,i,  j-1,k  ,b))  + &
                                                c2 * (q(v,i,  j+2,k  ,b) + q(v,i,  j-2,k  ,b))  + &
                                                c3 * (q(v,i,  j+3,k  ,b) + q(v,i,  j-3,k  ,b))  + &
                                                c4 * (q(v,i,  j+4,k  ,b) + q(v,i,  j-4,k  ,b))) + &
                                        idz2 * (c1 * (q(v,i,  j,  k+1,b) + q(v,i,  j,  k-1,b))  + &
                                                c2 * (q(v,i,  j,  k+2,b) + q(v,i,  j,  k-2,b))  + &
                                                c3 * (q(v,i,  j,  k+3,b) + q(v,i,  j,  k-3,b))  + &
                                                c4 * (q(v,i,  j,  k+4,b) + q(v,i,  j,  k-4,b))) - f(v,i,j,k,b))
               dq_max_ = max(dq_max_, abs(q(v,i,j,k,b) - q_old))
            enddo
         enddo
         enddo
         enddo
      enddo
      call apply_bc_elliptic(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, f=f, q=q, bc_type=bc_type, &
                             ell_bc_type=ell_bc_type, local_map_bc_crown=local_map_bc_crown, ivar=ivar, mu=mu, eps=eps, &
                             field=field, rebuild_exact_open=.false.)
      call update_smoothing_progress(progress_label=progress_label, progress_counter=progress_counter,            &
                                     progress_total=progress_total, progress_last_percent=progress_last_percent)
   enddo

   if (present(dq_max)) dq_max = dq_max_
   endsubroutine compute_smoothing_gauss_seidel_4th

   subroutine compute_smoothing_gauss_seidel_6th(ni, nj, nk, ngc, nv, blocks_number, dxyz, f, q, dq, dq_max, &
                                             iterations_init, iterations_fine, iterations_coarse, bc_type, ivar, mu, eps, field, &
                                             ell_bc_type, local_map_bc_crown, progress_label, progress_counter, progress_total, &
                                             progress_last_percent)
   !< Compute smoothing by Gauss-Seidel using L = D(G), with sixth-order
   !< centered first-derivative operators.
   integer(I4P),       intent(in)              :: ni,nj,nk,ngc               !< Grid dimensions.
   integer(I4P),       intent(in)              :: nv                         !< Number of q variables.
   integer(I4P),       intent(in)              :: blocks_number              !< Number of current blocks.
   real(R8P),          intent(in)              :: dxyz(1:,1:)                !< Space steps.
   real(R8P),          intent(in)              :: f(1:,     &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1:)                      !< Forcing distribution.
   real(R8P),          intent(inout)           :: q(1:,     &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1:)                      !< Field variables.
   real(R8P),          intent(inout), optional :: dq(1:,     &
                                                     1-ngc:, &
                                                     1-ngc:, &
                                                     1-ngc:, &
                                                     1:)                     !< Residuals.
   real(R8P),          intent(inout), optional :: dq_max                     !< Maximum residual.
   integer(I4P),       intent(in),    optional :: iterations_init            !< Smoothing iterations to initialize guess.
   integer(I4P),       intent(in),    optional :: iterations_fine            !< Smoothing iterations for fine grid.
   integer(I4P),       intent(in),    optional :: iterations_coarse          !< Smoothing iterations for coarse grid.
   character(len=*),   intent(in),    optional :: bc_type                    !< Boundary condition type
   real(R8P),          intent(in),    optional :: mu, eps                    !< Constant of electromagnetism, =1 in adimensional case
   integer(I4P),       intent(in),    optional :: ivar                       !< Variable (start) index in q.
   type(field_object), intent(in),    optional :: field                      !< Field (sibling realm component, threaded in).
   integer(I4P),       intent(in),    optional :: ell_bc_type(6)             !< Elliptic BC type for each face.
   integer(I8P),       intent(in),    optional :: local_map_bc_crown(:,:,:)  !< BC crown map.
   character(len=*),   intent(in),    optional :: progress_label             !< Progress message label.
   integer(I4P),       intent(inout), optional :: progress_counter           !< Completed smoothing sweeps.
   integer(I4P),       intent(in),    optional :: progress_total             !< Planned smoothing sweeps.
   integer(I4P),       intent(inout), optional :: progress_last_percent      !< Last printed progress percentage.
   real(R8P), parameter                        :: c1 =  23._R8P /  100._R8P  !< Constant for laplacian computation
   real(R8P), parameter                        :: c2 =  43._R8P /   80._R8P  !< Constant for laplacian computation
   real(R8P), parameter                        :: c3 =  -9._R8P /   40._R8P  !< Constant for laplacian computation
   real(R8P), parameter                        :: c4 =  19._R8P /  400._R8P  !< Constant for laplacian computation
   real(R8P), parameter                        :: c5 =  -1._R8P /  200._R8P  !< Constant for laplacian computation
   real(R8P), parameter                        :: c6 =   1._R8P / 3600._R8P  !< Constant for laplacian computation
   integer(I4P)                                :: iterations_                !< Smoothing iterations, local var.
   real(R8P)                                   :: dq_max_                    !< Maximum residual, local var.
   real(R8P)                                   :: dx2,dy2,dz2                !< Square space steps.
   real(R8P)                                   :: idx2,idy2,idz2             !< Inverse square space steps.
   real(R8P)                                   :: q_old                      !< Previous q.
   real(R8P)                                   :: factor                     !< Relaxation factor.
   integer(I4P)                                :: i,j,k,b,v,iter             !< Counter.

   if (ngc < 6_I4P) error stop 'compute_smoothing_gauss_seidel_6th: ngc must be >= 6 for sixth-order D(G) stencil'

   iterations_ = 1_I4P; if (present(iterations_fine)) iterations_ = iterations_fine
   dq_max_ = 0._R8P
   call apply_bc_elliptic(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, f=f, q=q, bc_type=bc_type, &
                          ell_bc_type=ell_bc_type, local_map_bc_crown=local_map_bc_crown, ivar=ivar, mu=mu, eps=eps, &
                          field=field, rebuild_exact_open=.true.)
   do iter=1, iterations_
      call apply_bc_elliptic(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, f=f, q=q, bc_type=bc_type, &
                             ell_bc_type=ell_bc_type, local_map_bc_crown=local_map_bc_crown, ivar=ivar, mu=mu, eps=eps, &
                             field=field, rebuild_exact_open=.false.)
      do b=1, blocks_number
         dx2 = dxyz(1,b)*dxyz(1,b)
         dy2 = dxyz(2,b)*dxyz(2,b)
         dz2 = dxyz(3,b)*dxyz(3,b)
         idx2 = 1._R8P / dx2
         idy2 = 1._R8P / dy2
         idz2 = 1._R8P / dz2
         factor = 1800._R8P / (2107._R8P * (idx2 + idy2 + idz2))
         do k=1, nk
         do j=1, nj
         do i=1, ni
            do v=1, nv
               q_old = q(v,i,j,k,b)
               q(v,i,j,k,b) = factor * (idx2 * (c1 * (q(v,i+1,j,  k  ,b) + q(v,i-1,j,  k  ,b))  + &
                                                c2 * (q(v,i+2,j,  k  ,b) + q(v,i-2,j,  k  ,b))  + &
                                                c3 * (q(v,i+3,j,  k  ,b) + q(v,i-3,j,  k  ,b))  + &
                                                c4 * (q(v,i+4,j,  k  ,b) + q(v,i-4,j,  k  ,b))  + &
                                                c5 * (q(v,i+5,j,  k  ,b) + q(v,i-5,j,  k  ,b))  + &
                                                c6 * (q(v,i+6,j,  k  ,b) + q(v,i-6,j,  k  ,b))) + &
                                        idy2 * (c1 * (q(v,i,  j+1,k  ,b) + q(v,i,  j-1,k  ,b))  + &
                                                c2 * (q(v,i,  j+2,k  ,b) + q(v,i,  j-2,k  ,b))  + &
                                                c3 * (q(v,i,  j+3,k  ,b) + q(v,i,  j-3,k  ,b))  + &
                                                c4 * (q(v,i,  j+4,k  ,b) + q(v,i,  j-4,k  ,b))  + &
                                                c5 * (q(v,i,  j+5,k  ,b) + q(v,i,  j-5,k  ,b))  + &
                                                c6 * (q(v,i,  j+6,k  ,b) + q(v,i,  j-6,k  ,b))) + &
                                        idz2 * (c1 * (q(v,i,  j,  k+1,b) + q(v,i,  j,  k-1,b))  + &
                                                c2 * (q(v,i,  j,  k+2,b) + q(v,i,  j,  k-2,b))  + &
                                                c3 * (q(v,i,  j,  k+3,b) + q(v,i,  j,  k-3,b))  + &
                                                c4 * (q(v,i,  j,  k+4,b) + q(v,i,  j,  k-4,b))  + &
                                                c5 * (q(v,i,  j,  k+5,b) + q(v,i,  j,  k-5,b))  + &
                                                c6 * (q(v,i,  j,  k+6,b) + q(v,i,  j,  k-6,b))) - f(v,i,j,k,b))      
               dq_max_ = max(dq_max_, abs(q(v,i,j,k,b) - q_old))
            enddo
         enddo
         enddo
         enddo
      enddo
      call apply_bc_elliptic(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, f=f, q=q, bc_type=bc_type, &
                             ell_bc_type=ell_bc_type, local_map_bc_crown=local_map_bc_crown, ivar=ivar, mu=mu, eps=eps, &
                             field=field, rebuild_exact_open=.false.)
      call update_smoothing_progress(progress_label=progress_label, progress_counter=progress_counter,            &
                                     progress_total=progress_total, progress_last_percent=progress_last_percent)
   enddo

   if (present(dq_max)) dq_max = dq_max_
   endsubroutine compute_smoothing_gauss_seidel_6th

   subroutine compute_smoothing_gauss_seidel_8th(ni, nj, nk, ngc, nv, blocks_number, dxyz, f, q, dq, dq_max, &
                                                 iterations_init, iterations_fine, iterations_coarse, bc_type, &
                                                 ivar, mu, eps, field, ell_bc_type, local_map_bc_crown,       &
                                                 progress_label, progress_counter, progress_total,             &
                                                 progress_last_percent)
   !< Compute smoothing by Gauss-Seidel using L = D(G), with eighth-order
   !< centered first-derivative operators.
   integer(I4P),       intent(in)              :: ni,nj,nk,ngc               !< Grid dimensions.
   integer(I4P),       intent(in)              :: nv                         !< Number of q variables.
   integer(I4P),       intent(in)              :: blocks_number              !< Number of current blocks.
   real(R8P),          intent(in)              :: dxyz(1:,1:)                !< Space steps.
   real(R8P),          intent(in)              :: f(1:,     &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1:)                      !< Forcing distribution.
   real(R8P),          intent(inout)           :: q(1:,     &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1-ngc:, &
                                                    1:)                      !< Field variables.
   real(R8P),          intent(inout), optional :: dq(1:,     &
                                                     1-ngc:, &
                                                     1-ngc:, &
                                                     1-ngc:, &
                                                     1:)                     !< Residuals.
   real(R8P),          intent(inout), optional :: dq_max                     !< Maximum residual.
   integer(I4P),       intent(in),    optional :: iterations_init            !< Smoothing iterations to initialize guess.
   integer(I4P),       intent(in),    optional :: iterations_fine            !< Smoothing iterations for fine grid.
   integer(I4P),       intent(in),    optional :: iterations_coarse          !< Smoothing iterations for coarse grid.
   character(len=*),   intent(in),    optional :: bc_type                    !< Boundary condition type.
   integer(I4P),       intent(in),    optional :: ivar                       !< Variable (start) index in q.
   real(R8P),          intent(in),    optional :: mu, eps                    !< Constant of electromagnetism, =1 in adimensional case.
   type(field_object), intent(in),    optional :: field                      !< Field (sibling realm component, threaded in).
   integer(I4P),       intent(in),    optional :: ell_bc_type(6)             !< Elliptic BC type for each face.
   integer(I8P),       intent(in),    optional :: local_map_bc_crown(:,:,:)  !< BC crown map.
   character(len=*),   intent(in),    optional :: progress_label             !< Progress message label.
   integer(I4P),       intent(inout), optional :: progress_counter           !< Completed smoothing sweeps.
   integer(I4P),       intent(in),    optional :: progress_total             !< Planned smoothing sweeps.
   integer(I4P),       intent(inout), optional :: progress_last_percent      !< Last printed progress percentage.
   real(R8P), parameter                        :: c1 =  411._R8P /   1225._R8P   !< Constant for laplacian computation.
   real(R8P), parameter                        :: c2 = 1213._R8P /   2100._R8P   !< Constant for laplacian computation.
   real(R8P), parameter                        :: c3 =  -11._R8P /     35._R8P   !< Constant for laplacian computation.
   real(R8P), parameter                        :: c4 =   53._R8P /    525._R8P   !< Constant for laplacian computation.
   real(R8P), parameter                        :: c5 =  -11._R8P /    525._R8P   !< Constant for laplacian computation.
   real(R8P), parameter                        :: c6 =  127._R8P /  44100._R8P   !< Constant for laplacian computation.
   real(R8P), parameter                        :: c7 =   -1._R8P /   3675._R8P   !< Constant for laplacian computation.
   real(R8P), parameter                        :: c8 =    1._R8P /  78400._R8P   !< Constant for laplacian computation.
   integer(I4P)                                :: iterations_                   !< Smoothing iterations, local var.
   real(R8P)                                   :: dq_max_                       !< Maximum residual, local var.
   real(R8P)                                   :: dx2,dy2,dz2                   !< Square space steps.
   real(R8P)                                   :: idx2,idy2,idz2                !< Inverse square space steps.
   real(R8P)                                   :: q_old                         !< Previous q.
   real(R8P)                                   :: factor                        !< Relaxation factor.
   integer(I4P)                                :: i,j,k,b,v,iter                !< Counter.

   if (ngc < 8_I4P) error stop 'compute_smoothing_gauss_seidel_8th: ngc must be >= 8 for eighth-order D(G) stencil'

   iterations_ = 1_I4P; if (present(iterations_fine)) iterations_ = iterations_fine
   dq_max_ = 0._R8P
   call apply_bc_elliptic(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, f=f, q=q, bc_type=bc_type, &
                          ell_bc_type=ell_bc_type, local_map_bc_crown=local_map_bc_crown, ivar=ivar, mu=mu, eps=eps, &
                          field=field, rebuild_exact_open=.true.)
   do iter=1, iterations_
      call apply_bc_elliptic(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, f=f, q=q, bc_type=bc_type, &
                             ell_bc_type=ell_bc_type, local_map_bc_crown=local_map_bc_crown, ivar=ivar, mu=mu, eps=eps, &
                             field=field, rebuild_exact_open=.false.)
      do b=1, blocks_number
         dx2 = dxyz(1,b)*dxyz(1,b)
         dy2 = dxyz(2,b)*dxyz(2,b)
         dz2 = dxyz(3,b)*dxyz(3,b)
         idx2 = 1._R8P / dx2
         idy2 = 1._R8P / dy2
         idz2 = 1._R8P / dz2
         factor = 352800._R8P / (480841._R8P * (idx2 + idy2 + idz2))
         do k=1, nk
         do j=1, nj
         do i=1, ni
            do v=1, nv
               q_old = q(v,i,j,k,b)
               q(v,i,j,k,b) = factor * (idx2 * (c1 * (q(v,i+1,j,  k  ,b) + q(v,i-1,j,  k  ,b))  + &
                                                c2 * (q(v,i+2,j,  k  ,b) + q(v,i-2,j,  k  ,b))  + &
                                                c3 * (q(v,i+3,j,  k  ,b) + q(v,i-3,j,  k  ,b))  + &
                                                c4 * (q(v,i+4,j,  k  ,b) + q(v,i-4,j,  k  ,b))  + &
                                                c5 * (q(v,i+5,j,  k  ,b) + q(v,i-5,j,  k  ,b))  + &
                                                c6 * (q(v,i+6,j,  k  ,b) + q(v,i-6,j,  k  ,b))  + &
                                                c7 * (q(v,i+7,j,  k  ,b) + q(v,i-7,j,  k  ,b))  + &
                                                c8 * (q(v,i+8,j,  k  ,b) + q(v,i-8,j,  k  ,b))) + &
                                        idy2 * (c1 * (q(v,i,  j+1,k  ,b) + q(v,i,  j-1,k  ,b))  + &
                                                c2 * (q(v,i,  j+2,k  ,b) + q(v,i,  j-2,k  ,b))  + &
                                                c3 * (q(v,i,  j+3,k  ,b) + q(v,i,  j-3,k  ,b))  + &
                                                c4 * (q(v,i,  j+4,k  ,b) + q(v,i,  j-4,k  ,b))  + &
                                                c5 * (q(v,i,  j+5,k  ,b) + q(v,i,  j-5,k  ,b))  + &
                                                c6 * (q(v,i,  j+6,k  ,b) + q(v,i,  j-6,k  ,b))  + &
                                                c7 * (q(v,i,  j+7,k  ,b) + q(v,i,  j-7,k  ,b))  + &
                                                c8 * (q(v,i,  j+8,k  ,b) + q(v,i,  j-8,k  ,b))) + &
                                        idz2 * (c1 * (q(v,i,  j,  k+1,b) + q(v,i,  j,  k-1,b))  + &
                                                c2 * (q(v,i,  j,  k+2,b) + q(v,i,  j,  k-2,b))  + &
                                                c3 * (q(v,i,  j,  k+3,b) + q(v,i,  j,  k-3,b))  + &
                                                c4 * (q(v,i,  j,  k+4,b) + q(v,i,  j,  k-4,b))  + &
                                                c5 * (q(v,i,  j,  k+5,b) + q(v,i,  j,  k-5,b))  + &
                                                c6 * (q(v,i,  j,  k+6,b) + q(v,i,  j,  k-6,b))  + &
                                                c7 * (q(v,i,  j,  k+7,b) + q(v,i,  j,  k-7,b))  + &
                                                c8 * (q(v,i,  j,  k+8,b) + q(v,i,  j,  k-8,b))) - f(v,i,j,k,b))
               dq_max_ = max(dq_max_, abs(q(v,i,j,k,b) - q_old))
            enddo
         enddo
         enddo
         enddo
      enddo
      call apply_bc_elliptic(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, f=f, q=q, bc_type=bc_type, &
                             ell_bc_type=ell_bc_type, local_map_bc_crown=local_map_bc_crown, ivar=ivar, mu=mu, eps=eps, &
                             field=field, rebuild_exact_open=.false.)
      call update_smoothing_progress(progress_label=progress_label, progress_counter=progress_counter,            &
                                     progress_total=progress_total, progress_last_percent=progress_last_percent)
   enddo

   if (present(dq_max)) dq_max = dq_max_
   endsubroutine compute_smoothing_gauss_seidel_8th

   subroutine compute_smoothing_multigrid(ni, nj, nk, ngc, nv, blocks_number, dxyz, f, q, dq, dq_max, &
                                          iterations_init, iterations_fine, iterations_coarse, bc_type, ivar, mu, eps, field, &
                                          ell_bc_type, local_map_bc_crown)
   integer(I4P), intent(in)              :: ni,nj,nk,ngc         !< Grid dimensions.
   integer(I4P), intent(in)              :: nv                   !< Number of q variables.
   integer(I4P), intent(in)              :: blocks_number        !< Number of current blocks.
   real(R8P),    intent(in)              :: dxyz(1:,1:)          !< Space steps.
   real(R8P),    intent(in)              :: f(1:,    &
                                              1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:,&
                                              1:)                !< Forcing distribution.
   real(R8P),    intent(inout)           :: q(1:,    &
                                              1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:,&
                                              1:)                !< Field variables.
   real(R8P),    intent(inout), optional :: dq(1:,    &
                                               1-ngc:,&
                                               1-ngc:,&
                                               1-ngc:,&
                                               1:)               !< Residuals.
   real(R8P),    intent(inout), optional :: dq_max               !< Maximum residual.
   integer(I4P), intent(in),    optional :: iterations_init      !< Smoothing iterations to initialize guess.
   integer(I4P), intent(in),    optional :: iterations_fine      !< Smoothing iterations for fine grid.
   integer(I4P), intent(in),    optional :: iterations_coarse    !< Smoothing iterations for coarse grid.
   character(len=*),   intent(in),    optional :: bc_type        !< Legacy scalar BC selector.
   integer(I4P),       intent(in),    optional :: ivar           !< Variable (start) index in q.
   real(R8P),          intent(in),    optional :: mu, eps        !< Electromagnetic constants.
   type(field_object), intent(in),    optional :: field          !< Field (sibling realm component, threaded in).
   integer(I4P),       intent(in),    optional :: ell_bc_type(6) !< Elliptic BC type for each face.
   integer(I8P),       intent(in),    optional :: local_map_bc_crown(:,:,:) !< BC crown map.
   real(R8P)                             :: f_c(1:nv,          &
                                                1-ngc:ni/2+ngc,&
                                                1-ngc:nj/2+ngc,&
                                                1-ngc:nk/2+ngc,&
                                                1:blocks_number) !< Forcing distribution, coarse grid.
   real(R8P)                             :: dq_c(1:nv,          &
                                                 1-ngc:ni/2+ngc,&
                                                 1-ngc:nj/2+ngc,&
                                                 1-ngc:nk/2+ngc,&
                                                 1:blocks_number)!< Residuals, coarse grid.
   real(R8P)                             :: q_c(1:nv,          &
                                                1-ngc:ni/2+ngc,&
                                                1-ngc:nj/2+ngc,&
                                                1-ngc:nk/2+ngc,&
                                                1:blocks_number) !< Field variables, coarse grid.
   integer(I4P)                          :: b                    !< Counter.

   ! V-cicle
   call compute_smoothing_gauss_seidel(ni=ni,nj=nj,nk=nk,ngc=ngc,nv=nv,blocks_number=blocks_number,dxyz=dxyz,f=f,q=q,&
                                       iterations_fine=iterations_init, bc_type=bc_type, ivar=ivar, mu=mu, eps=eps, &
                                       field=field, ell_bc_type=ell_bc_type, local_map_bc_crown=local_map_bc_crown)
   call compute_laplacian_residuals(ni=ni,nj=nj,nk=nk,ngc=ngc,nv=nv,blocks_number=blocks_number,dxyz=dxyz,f=f,q=q,dq=dq)
   call compute_restriction(ngc=ngc,nic=ni/2,njc=nj/2,nkc=nk/2,nv=nv,blocks_number=blocks_number,fine=dq,coarse=f_c)
   do b=1, blocks_number
      q_c(:,:,:,:,b) = 0._R8P
   enddo
   call compute_smoothing_gauss_seidel(ni=ni/2,nj=nj/2,nk=nk/2,ngc=ngc,nv=nv,blocks_number=blocks_number,dxyz=dxyz*2_R8P,&
                                       f=f_c,q=q_c,iterations_fine=iterations_coarse, bc_type=bc_type, ivar=ivar, mu=mu, &
                                       eps=eps, field=field, ell_bc_type=ell_bc_type, local_map_bc_crown=local_map_bc_crown)
   call compute_prolongation(ngc=ngc,nic=ni/2,njc=nj/2,nkc=nk/2,nv=nv,blocks_number=blocks_number,coarse=q_c,fine=dq)
   do b=1, blocks_number
      q(:,:,:,:,b) = q(:,:,:,:,b) + dq(:,:,:,:,b)
   enddo
   call compute_smoothing_gauss_seidel(ni=ni,nj=nj,nk=nk,ngc=ngc,nv=nv,blocks_number=blocks_number,dxyz=dxyz,f=f,q=q,dq_max=dq_max,&
                                       iterations_fine=iterations_fine, bc_type=bc_type, ivar=ivar, mu=mu, eps=eps, &
                                       field=field, ell_bc_type=ell_bc_type, local_map_bc_crown=local_map_bc_crown)
   endsubroutine compute_smoothing_multigrid

   subroutine compute_smoothing_sor(ni, nj, nk, ngc, nv, blocks_number, dxyz, f, q, dq, dq_max, &
                                    iterations_init, iterations_fine, iterations_coarse)
   !< Compute smoothing by SOR, Successive Overrelaxation.
   integer(I4P), intent(in)              :: ni,nj,nk,ngc      !< Grid dimensions.
   integer(I4P), intent(in)              :: nv                !< Number of q variables.
   integer(I4P), intent(in)              :: blocks_number     !< Number of current blocks.
   real(R8P),    intent(in)              :: dxyz(1:,1:)       !< Space steps.
   real(R8P),    intent(in)              :: f(1:,    &
                                              1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:,&
                                              1:)             !< Forcing distribution.
   real(R8P),    intent(inout)           :: q(1:,    &
                                              1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:,&
                                              1:)             !< Field variables.
   real(R8P),    intent(inout), optional :: dq(1:,    &
                                               1-ngc:,&
                                               1-ngc:,&
                                               1-ngc:,&
                                               1:)            !< Residuals.
   real(R8P),    intent(inout), optional :: dq_max            !< Maximum residual.
   integer(I4P), intent(in),    optional :: iterations_init   !< Smoothing iterations to initialize guess.
   integer(I4P), intent(in),    optional :: iterations_fine   !< Smoothing iterations for fine grid.
   integer(I4P), intent(in),    optional :: iterations_coarse !< Smoothing iterations for coarse grid.
   integer(I4P)                          :: iterations_       !< Smoothing iterations, local var.
   real(R8P)                             :: dq_max_           !< Maximum residual, local var.
   real(R8P)                             :: dx2,dy2,dz2       !< Square space steps.
   real(R8P)                             :: factor            !< Jacobi relaxation factor.
   real(R8P)                             :: omega             !< Overrelaxation parameter.
   real(R8P)                             :: q_old             !< Previous q.
   integer(I4P)                          :: i,j,k,b,v,iter    !< Counter.

   iterations_ = 1 ; if (present(iterations_fine)) iterations_ = iterations_fine
   omega = 1.8_R8P
   dq_max_ = 0._R8P
   do iter=1, iterations_
      call apply_bc_dirichlet(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, q=q)
      do b=1, blocks_number
         dx2 = dxyz(1,b)*dxyz(1,b)
         dy2 = dxyz(2,b)*dxyz(2,b)
         dz2 = dxyz(3,b)*dxyz(3,b)
         factor = 1._R8P / (2._R8P * (1._R8P/dx2 + 1._R8P/dy2 + 1._R8P/dz2))
         do k = 1, nk
         do j = 1, nj
         do i = 1, ni
            do v=1, nv
               q_old = q(v,i,j,k,b)
               q(v,i,j,k,b) = factor * ((q(v,i+1,j,  k  ,b) +  q(v,i-1,j,  k  ,b)) / dx2 + &
                                        (q(v,i,  j+1,k  ,b) +  q(v,i,  j-1,k  ,b)) / dy2 + &
                                        (q(v,i,  j,  k+1,b) +  q(v,i,  j,  k-1,b)) / dz2 - &
                                         f(v,i,  j,  k  ,b))
               q(v,i,j,k,b) = q_old + omega * (q(v,i,j,k,b) - q_old)
               dq_max_ = max(dq_max, abs(q(v,i,j,k,b) - q_old))
            enddo
         enddo
         enddo
         enddo
      enddo
   enddo
   if (present(dq_max)) dq_max = dq_max_
   endsubroutine compute_smoothing_sor

   subroutine compute_smoothing_sor_omp(ni, nj, nk, ngc, nv, blocks_number, dxyz, f, q, dq, dq_max, &
                                        iterations_init, iterations_fine, iterations_coarse)
   !< Compute smoothing by SOR, Successive Overrelaxation (parallel OpenMP) method: red-black 3D scheme.
   integer(I4P), intent(in)              :: ni,nj,nk,ngc      !< Grid dimensions.
   integer(I4P), intent(in)              :: nv                !< Number of q variables.
   integer(I4P), intent(in)              :: blocks_number     !< Number of current blocks.
   real(R8P),    intent(in)              :: dxyz(1:,1:)       !< Space steps.
   real(R8P),    intent(in)              :: f(1:,    &
                                              1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:,&
                                              1:)             !< Forcing distribution.
   real(R8P),    intent(inout)           :: q(1:,    &
                                              1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:,&
                                              1:)             !< Field variables.
   real(R8P),    intent(inout), optional :: dq(1:,    &
                                               1-ngc:,&
                                               1-ngc:,&
                                               1-ngc:,&
                                               1:)            !< Residuals.
   real(R8P),    intent(inout), optional :: dq_max            !< Maximum residual.
   integer(I4P), intent(in),    optional :: iterations_init   !< Smoothing iterations to initialize guess.
   integer(I4P), intent(in),    optional :: iterations_fine   !< Smoothing iterations for fine grid.
   integer(I4P), intent(in),    optional :: iterations_coarse !< Smoothing iterations for coarse grid.
   integer(I4P)                          :: iterations_       !< Smoothing iterations, local var.
   real(R8P)                             :: dq_max_           !< Maximum residual, local var.
   real(R8P)                             :: factor            !< Jacobi relaxation factor.
   real(R8P)                             :: dx2,dy2,dz2       !< Square space steps.
   real(R8P)                             :: omega             !< Overrelaxation parameter.
   integer(I4P)                          :: i,j,k,b,v,iter    !< Counter.
   integer(I4P)                          :: color             !< Color counter.
   real(R8P)                             :: residual          !< Residual.

   iterations_ = 1 ; if (present(iterations_fine)) iterations_ = iterations_fine
   omega = 1.8_R8P
   dq_max_ = 0._R8P
   do iter=1, iterations_
      call apply_bc_dirichlet(ni=ni, nj=nj, nk=nk, ngc=ngc, blocks_number=blocks_number, q=q)
      do b=1, blocks_number
         dx2 = dxyz(1,b)*dxyz(1,b)
         dy2 = dxyz(2,b)*dxyz(2,b)
         dz2 = dxyz(3,b)*dxyz(3,b)
         factor = 1._R8P / (2._R8P * (1._R8P/dx2 + 1._R8P/dy2 + 1._R8P/dz2))
         do color=0, 7
            !$omp parallel do firstprivate(b,factor,omega,color) private(residual) shared(q,dq) reduction(max:dq_max) collapse(3)
            do k = 1, nk
            do j = 1, nj
            do i = 1, ni
               if (mod(i+j+k,8)==color) then
                  do v=1, nv
                     residual = factor * ((q(v,i+1,j,  k  ,b) +  q(v,i-1,j,  k  ,b)) / dx2 + &
                                          (q(v,i,  j+1,k  ,b) +  q(v,i,  j-1,k  ,b)) / dy2 + &
                                          (q(v,i,  j,  k+1,b) +  q(v,i,  j,  k-1,b)) / dz2 - &
                                           f(v,i,  j,  k  ,b))
                     dq(v,i,j,k,b) = q(v,i,j,k,b) + omega * (residual - q(v,i,j,k,b))
                     dq_max_ = max(dq_max_, abs(residual - q(v,i,j,k,b)))
                  enddo
               endif
            enddo
            enddo
            enddo
            !$omp parallel do firstprivate(b,color) shared(q,dq) collapse(3)
            do k = 1, nk
            do j = 1, nj
            do i = 1, ni
               if (mod(i+j+k,8)==color) then
                  do v=1, nv
                     q(v,i,j,k,b) = dq(v,i,j,k,b)
                  enddo
               endif
            enddo
            enddo
            enddo
         enddo
      enddo
   enddo
   if (present(dq_max)) dq_max = dq_max_
   endsubroutine compute_smoothing_sor_omp

   subroutine update_smoothing_progress(progress_label, progress_counter, progress_total, progress_last_percent, &
                                        converged)
   character(len=*),   intent(in),    optional :: progress_label        !< Progress message label.
   integer(I4P),       intent(inout), optional :: progress_counter      !< Completed smoothing sweeps.
   integer(I4P),       intent(in),    optional :: progress_total        !< Planned smoothing sweeps.
   integer(I4P),       intent(inout), optional :: progress_last_percent !< Last printed progress percentage.
   logical,            intent(in),    optional :: converged             !< Force 100% when converged.
   integer(I4P)                                :: counter_              !< Local progress counter.
   integer(I4P)                                :: total_                !< Local planned sweeps.
   integer(I4P)                                :: progress              !< Current progress percentage.
   logical                                     :: converged_            !< Local convergence flag.

   if (.not. present(progress_label)) return
   if (.not. present(progress_counter)) return
   if (.not. present(progress_total)) return
   if (.not. present(progress_last_percent)) return

   converged_ = .false.
   if (present(converged)) converged_ = converged

   total_ = max(1_I4P, progress_total)
   counter_ = min(max(progress_counter + merge(0_I4P, 1_I4P, converged_), 0_I4P), total_)
   if (.not. converged_) progress_counter = counter_
   progress = int(100._R8P * real(counter_, R8P) / real(total_, R8P), kind=I4P)
   if (converged_) progress = 100_I4P
   progress = min(100_I4P, max(0_I4P, progress))

   if (progress > progress_last_percent) then
      call mpih%print_message(trim(progress_label)//' progress: '//trim(str(progress,.true.))//'% ('// &
                              trim(str(counter_,.true.))//'/'//trim(str(total_,.true.))//')')
      progress_last_percent = progress
   endif
   endsubroutine update_smoothing_progress
endmodule adam_flail_object
