!< ADAM, MHD Initial Conditions class definition, CPU backend.
module adam_mhd_ic_object
!< ADAM, MHD Initial Conditions class definition, CPU backend.

use adam_field_object, only : field_object
use adam_mpih_object, only : mpih_object
use adam_mhd_physics_object, only : mhd_physics_object
use finer
use penf

implicit none
private
public :: IC_TYPE_RP
public :: IC_TYPE_IVORTEX
public :: IC_TYPE_UNIFORM
public :: mhd_ic_object

character(len=18), parameter :: INI_SECTION_NAME="initial_conditions" !< INI (config) file section name containing IC configs.
character(len=7),  parameter :: IC_TYPE_UNIFORM="uniform"             !< Uniform IC TYPE parameter.
character(len=17), parameter :: IC_TYPE_IVORTEX="isentropic-vortex"   !< Isentropic Vortex IC TYPE parameter.
character(len=15), parameter :: IC_TYPE_RP="riemann-problem"          !< Riemann Problem IC TYPE parameter.

type :: mhd_ic_object
   !< Initial Conditions class definition, CPU backend.
   type(mpih_object)         :: mpih                 !< MPI handler.
   character(:), allocatable :: ic_type              !< IC type.
   integer(I4P)              :: regions_number=1     !< Number of IC regions.
   real(R8P), allocatable    :: q(:,:)               !< Primitive variables (r,u,v,w,p), s fluid specie index at IC for each region.
   real(R8P), allocatable    :: emin(:,:), emax(:,:) !< IC regions bounding box.
   contains
      ! public methods
      procedure, pass(self) :: description            !< Return pretty-printed object description.
      procedure, pass(self) :: initialize             !< Initialize IC.
      procedure, pass(self) :: load_from_file         !< Load config from file.
      procedure, pass(self) :: set_initial_conditions !< Set initial conditions on mhd fields.
endtype mhd_ic_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(mhd_ic_object), intent(in) :: self             !< IC.
   character(len=:), allocatable      :: desc             !< Description.
   character(len=1), parameter        :: NL=new_line('a') !< New line character.
   integer(I4P)                       :: r, v             !< Counter.

   desc =       self%mpih%myrankstr//'IC main data'//NL
   desc = desc//self%mpih%myrankstr//'  regions number: '//trim(str(self%regions_number))
   do r=1, self%regions_number
      desc = desc//NL//self%mpih%myrankstr//'  region('//trim(str(r,.true.))//')'
      do v=1, size(self%q, dim=1)
         desc = desc//NL//self%mpih%myrankstr//'    q('//trim(str(v,.true.))//'): '//trim(str(self%q(v,r)))
      enddo
         desc = desc//NL//self%mpih%myrankstr//'    emin: '//trim(str(self%emin(:,r)))
         desc = desc//NL//self%mpih%myrankstr//'    emax: '//trim(str(self%emax(:,r)))
   enddo
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize the equation.
   class(mhd_ic_object), intent(inout) :: self            !< IC.
   type(file_ini),         intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'mhd_ic_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'mhd_ic_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(mhd_ic_object), intent(inout)        :: self            !< IC.
   type(file_ini),         intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                      :: go_on_fail_     !< Go on if load fails.
   character(:), allocatable                    :: sname           !< Section name.
   integer(I4P)                                 :: i               !< Counter.
   integer(I4P)                                 :: error           !< Error status.
   character(99)                                :: buff_char       !< Option character buffer.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='type', val=buff_char, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(type)')
   self%ic_type = trim(adjustl(buff_char))
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='regions_number', val=self%regions_number, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(regions_number)')

   if (self%regions_number>=1) then
      allocate(   self%q(1:6, 1:self%regions_number))
      allocate(self%emin(1:3, 1:self%regions_number))
      allocate(self%emax(1:3, 1:self%regions_number))
      do i=1, self%regions_number
         sname = INI_SECTION_NAME//'_region_'//trim(str(i,.true.))
         call file_parameters%get(section_name=sname, option_name='r', val=self%q(1,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(r)')
         call file_parameters%get(section_name=sname, option_name='u', val=self%q(2,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(u)')
         call file_parameters%get(section_name=sname, option_name='v', val=self%q(3,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(v)')
         call file_parameters%get(section_name=sname, option_name='w', val=self%q(4,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(w)')
         call file_parameters%get(section_name=sname, option_name='p', val=self%q(5,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(p)')
         call file_parameters%get(section_name=sname, option_name='s', val=self%q(6,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(s)')
         call file_parameters%get(section_name=sname, option_name='emin_x', val=self%emin(1,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(emin_x)')
         call file_parameters%get(section_name=sname, option_name='emin_y', val=self%emin(2,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(emin_y)')
         call file_parameters%get(section_name=sname, option_name='emin_z', val=self%emin(3,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(emin_z)')
         call file_parameters%get(section_name=sname, option_name='emax_x', val=self%emax(1,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(emax_x)')
         call file_parameters%get(section_name=sname, option_name='emax_y', val=self%emax(2,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(emax_y)')
         call file_parameters%get(section_name=sname, option_name='emax_z', val=self%emax(3,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(emax_z)')
      enddo
   endif
   endsubroutine load_from_file

   subroutine set_initial_conditions(self, physics, field)
   !< Set initial conditions on mhd fields.
   class(mhd_ic_object),     intent(in)    :: self                 !< IC.
   type(mhd_physics_object), intent(in)    :: physics              !< Fluids physiscs.
   type(field_object),         intent(inout) :: field                !< Field object.
   integer(I4P)                              :: b, i, j, k, s, ri    !< Counter.
   real(R8P)                                 :: rn                   !< Random number.
   real(R8P)                                 :: cv, R, g, delta, gm1 !< Fluid physics constants.
   real(R8P)                                 :: rho, u, v, w, p      !< Primitive variables.
   real(R8P)                                 :: xc, yc               !< Isentropic vortex center.
   real(R8P)                                 :: Rv, Rv2, r2, expr2   !< Isentropic vortex radius.
   real(R8P)                                 :: av                   !< Isentropic vortex speed of sound.
   real(R8P)                                 :: Mv, Mv2              !< Isentropic vortex Mach number.

   associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, &
             q=>field%q, x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell)
   select case(self%ic_type)
   case(IC_TYPE_UNIFORM) ! uniform, only one region (s=1); q(6,1) is the base level for random velocity perturbation
      cv = physics%eos(1)%cv
      R  = physics%eos(1)%R
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                                                                 rho = self%q(1,1)
                  call random_number(rn) ; rn = rn*self%q(6,1) ; u   = self%q(2,1)+rn
                  call random_number(rn) ; rn = rn*self%q(6,1) ; v   = self%q(3,1)+rn
                  call random_number(rn) ; rn = rn*self%q(6,1) ; w   = self%q(4,1)+rn
                                                                 p   = self%q(5,1)
                  q(1,i,j,k,b) = rho
                  q(2,i,j,k,b) = rho*u
                  q(3,i,j,k,b) = rho*v
                  q(4,i,j,k,b) = rho*w
                  q(5,i,j,k,b) = rho*(cv*p/(rho*R) + 0.5_R8P*(u**2+v**2+w**2))
               enddo
            enddo
         enddo
      enddo
   case(IC_TYPE_IVORTEX) ! isentropic vortex, only one region (s=1); q(6,1) is the vortex radius and emin(x,y) is the vortex center
      cv    = physics%eos(1)%cv
      R     = physics%eos(1)%R
      g     = physics%eos(1)%g
      delta = physics%eos(1)%delta
      gm1   = physics%eos(1)%gm1

      xc  = self%emin(1,1)
      yc  = self%emin(2,1)
      Rv  = self%q(6,1)
      Rv2 = Rv**2
      av  = sqrt(g*self%q(5,1)/self%q(1,1))
      Mv  = self%q(2,1)/av
      Mv2 = Mv**2
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  r2     = ((x_cell(i,b) - xc)**2 + (y_cell(j,b) - yc)**2)/Rv2
                  expr2  = exp((1._R8P - r2)*0.5_R8P)

                  rho    = self%q(1,1)*(1._R8P - delta*Mv2*expr2**2)**(1._R8P/gm1)
                  u      = Mv*(1._R8P  - (y_cell(j,b)-yc)/Rv*expr2)*av
                  v      = Mv*(          (x_cell(i,b)-xc)/Rv*expr2)*av
                  w      = 0._R8P
                  p      = self%q(5,1)*(rho/self%q(1,1))**g

                  q(1,i,j,k,b) = rho
                  q(2,i,j,k,b) = rho*u
                  q(3,i,j,k,b) = rho*v
                  q(4,i,j,k,b) = 0._R8P
                  q(5,i,j,k,b) = rho*(cv*p/(rho*R) + 0.5_R8P*(u**2+v**2+w**2))
               enddo
            enddo
         enddo
      enddo
   case(IC_TYPE_RP)      ! Riemann Problem like (uniform regions)
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do ri=1, self%regions_number
                     if ((x_cell(i,b) > self%emin(1,ri).and.x_cell(i,b) <= self%emax(1,ri)).and. &
                         (y_cell(j,b) > self%emin(2,ri).and.y_cell(j,b) <= self%emax(2,ri)).and. &
                         (z_cell(k,b) > self%emin(3,ri).and.z_cell(k,b) <= self%emax(3,ri))) then
                        q(1:5,i,j,k,b) = physics%eos(int(self%q(6,ri)))%primitive2conservative(primitive=self%q(1:5,ri))
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
   endselect
   endassociate
   endsubroutine set_initial_conditions
endmodule adam_mhd_ic_object
