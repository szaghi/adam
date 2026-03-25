!< ADAM, FLAIL Fortran Linear Algebra Interface Library class definition, CPU common backend.
module adam_flail_object
!< ADAM, FLAIL Fortran Linear Algebra Interface Library class definition, CPU common backend.

! ADAM singleton objects
use :: adam_global_mpih, only : mpih
! third party modules
use :: finer
use :: penf

implicit none
private
public :: flail_object
public :: compute_smoothing_interface
public :: compute_smoothing_gauss_seidel
public :: compute_smoothing_multigrid
public :: compute_smoothing_sor
public :: compute_smoothing_sor_omp

character(9),  parameter, public :: SMOOTHING_MULTIGRID   ='MULTIGRID'    !< Smoothing multigrid parameter.
character(12), parameter, public :: SMOOTHING_GAUSS_SEIDEL='GAUSS-SEIDEL' !< Smoothing Gauss-Seidel parameter.
character(3),  parameter, public :: SMOOTHING_SOR         ='SOR'          !< Smoothing SOR parameter.
character(7),  parameter, public :: SMOOTHING_SOR_OMP     ='SOR-OMP'      !< Smoothing SOR-OpenMP parameter.

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
   integer(I4P)                :: b             !< Counter.

   do b=1, blocks_number
      q(:,0,:,:,b) = 0._R8P ; q(:,ni+1,:   ,:   ,b) = 0._R8P
      q(:,:,0,:,b) = 0._R8P ; q(:,:   ,nj+1,:   ,b) = 0._R8P
      q(:,:,:,0,b) = 0._R8P ; q(:,:   ,:   ,nk+1,b) = 0._R8P
   enddo
   endsubroutine apply_bc_dirichlet

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

   subroutine compute_smoothing_gauss_seidel(ni, nj, nk, ngc, nv, blocks_number, dxyz, f, q, dq, dq_max, &
                                             iterations_init, iterations_fine, iterations_coarse)
   !< Compute smoothing by Gauss-Seidel method.
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
   real(R8P)                             :: q_old             !< Previous q.
   real(R8P)                             :: factor            !< Jacobi relaxation factor.
   integer(I4P)                          :: i,j,k,b,v,iter    !< Counter.

   iterations_ = 1 ; if (present(iterations_fine)) iterations_ = iterations_fine
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
   enddo
   if (present(dq_max)) dq_max = dq_max_
   endsubroutine compute_smoothing_gauss_seidel

   subroutine compute_smoothing_multigrid(ni, nj, nk, ngc, nv, blocks_number, dxyz, f, q, dq, dq_max, &
                                          iterations_init, iterations_fine, iterations_coarse)
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
                                       iterations_fine=iterations_init)
   call compute_laplacian_residuals(ni=ni,nj=nj,nk=nk,ngc=ngc,nv=nv,blocks_number=blocks_number,dxyz=dxyz,f=f,q=q,dq=dq)
   call compute_restriction(ngc=ngc,nic=ni/2,njc=nj/2,nkc=nk/2,nv=nv,blocks_number=blocks_number,fine=dq,coarse=f_c)
   do b=1, blocks_number
      q_c(:,:,:,:,b) = 0._R8P
   enddo
   call compute_smoothing_gauss_seidel(ni=ni/2,nj=nj/2,nk=nk/2,ngc=ngc,nv=nv,blocks_number=blocks_number,dxyz=dxyz*2_R8P,&
                                       f=f_c,q=q_c,iterations_fine=iterations_coarse)
   call compute_prolongation(ngc=ngc,nic=ni/2,njc=nj/2,nkc=nk/2,nv=nv,blocks_number=blocks_number,coarse=q_c,fine=dq)
   do b=1, blocks_number
      q(:,:,:,:,b) = q(:,:,:,:,b) + dq(:,:,:,:,b)
   enddo
   call compute_smoothing_gauss_seidel(ni=ni,nj=nj,nk=nk,ngc=ngc,nv=nv,blocks_number=blocks_number,dxyz=dxyz,f=f,q=q,dq_max=dq_max,&
                                       iterations_fine=iterations_fine)
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
endmodule adam_flail_object
