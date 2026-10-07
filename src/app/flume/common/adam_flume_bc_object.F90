!< ADAM, FLUME boundary conditions class definition.
module adam_flume_bc_object
!< ADAM, FLUME boundary conditions class definition.
!<
!< `periodic` is the library `BC_PERIODIC`, so the tree builds true periodic neighbors and the ghost exchange fills
!< periodic ghosts across blocks and ranks (verified by `src/tests/amr/test_periodic_ghost`); the other kinds are
!< filled by the backends on the boundary crown maps.
!<
!< **Realm edges and corners.** A crown row beyond several realm faces (an edge or a corner of the realm) is filled by
!< the kind of one of those faces, `realm_edge_face`, applied to the donor `realm_edge_donor`: mirrored (wall) or
!< clamped (extrapolation) along that face's axis, so the donor lies beyond the other faces only and is a ghost
!< already filled. The backends fill the rows beyond one face first, then two, then three (issue #65 P0).
!<
!< **No-slip walls** (issue #65, D-M4-7). `wall-noslip` (adiabatic) and `wall-isothermal` mirror the interior cell like
!< `wall-inviscid` but reflect the whole velocity about the wall velocity, `u_g = 2 u_w - u` (`wall_u, wall_v, wall_w`,
!< tangential, default 0; a normal component is fatal), so the mean of the two cells is the wall velocity. The pressure
!< is mirrored. The adiabatic wall mirrors the temperature (and the density); the isothermal wall sets the ghost
!< temperature `T_g = 2 wall_temperature - T`, so the mean is the wall temperature, and the density from the mirrored
!< pressure, `rho_g = p / (R T_g)`: the ghost stays positive while `T < 2 wall_temperature`. The field follows the MHD
!< wall rule of `wall-inviscid` (normal component odd), psi is even. Both are second-order at the wall.

! ADAM classes, libraries, parameters
use :: adam_parameters,           only : BC_PERIODIC, BC_SEAM, FEC_TO_DELTA
! ADAM singleton objects
use :: adam_mpih_global,          only : mpih
! FLUME modules
use :: adam_flume_parameters,     only : IQ_BX, IQ_BZ, IQ_R, IQ_RE, IQ_RU, IQ_RW, MODEL_EULER, MODEL_MHD,            &
                                        MODEL_MHD_EGLM, MODEL_MHD_GLM, strip_control
use :: adam_flume_physics_object, only : flume_physics_object, primitive_state_to_conservative
! third party modules
use :: finer,                     only : file_ini
use :: penf,                      only : I4P, R8P, str

implicit none
private
public :: flume_bc_object
public :: BC_EXTRAPOLATION
public :: BC_INFLOW
public :: BC_WALL_INVISCID
public :: BC_WALL_NOSLIP
public :: BC_WALL_ISOTHERMAL
public :: wall_noslip_ghost
public :: BC_PERIODIC
public :: realm_edge_donor
public :: realm_edge_face

integer(I4P), parameter :: BC_EXTRAPOLATION = 1_I4P !< Zeroth-order extrapolation.
integer(I4P), parameter :: BC_INFLOW        = 2_I4P !< Prescribed state.
integer(I4P), parameter :: BC_WALL_INVISCID = 3_I4P !< Inviscid (slip) wall: mirror with normal momentum negated.
integer(I4P), parameter :: BC_WALL_NOSLIP = 4_I4P   !< No-slip adiabatic wall: velocity reflected about the wall's.
integer(I4P), parameter :: BC_WALL_ISOTHERMAL = 5_I4P !< No-slip isothermal wall: also the temperature set.

character(len=13), parameter :: BC_EXTRAPOLATION_STR="extrapolation"          !< Accepted spelling of BC_EXTRAPOLATION.
character(len=6),  parameter :: BC_INFLOW_STR="inflow"                        !< Accepted spelling of BC_INFLOW.
character(len=13), parameter :: BC_WALL_INVISCID_STR="wall-inviscid"          !< Accepted spelling of BC_WALL_INVISCID.
character(len=8),  parameter :: BC_PERIODIC_STR="periodic"                    !< Accepted spelling of BC_PERIODIC.
character(len=11), parameter :: BC_WALL_NOSLIP_STR="wall-noslip"              !< Accepted spelling of BC_WALL_NOSLIP.
character(len=15), parameter :: BC_WALL_ISOTHERMAL_STR="wall-isothermal"      !< Accepted spelling of BC_WALL_ISOTHERMAL.
character(len=6),  parameter :: WALL_VELOCITY_KEY(3)=['wall_u', 'wall_v', 'wall_w'] !< Wall velocity keys.
character(len=2),  parameter :: INFLOW_KEY(8)=['r ', 'u ', 'v ', 'w ', &
                                               'p ', 'bx', 'by', 'bz']       !< Inflow primitive state keys (MHD: all 8).
character(len=8),  parameter :: SECTION_NAME(6)=['bc_x_min', 'bc_x_max', &
                                                 'bc_y_min', 'bc_y_max', &
                                                 'bc_z_min', 'bc_z_max']      !< INI section names of the 6 faces.

type :: flume_bc_object
   !< FLUME boundary conditions class definition.
   integer(I4P)           :: bc_type(6)=0_I4P !< Boundary condition type of each face.
   real(R8P), allocatable :: q_inflow(:,:)    !< Conservative inflow state of each face [nv, 6].
   real(R8P), allocatable :: wall_sign(:,:)   !< Wall mirror sign of each variable per direction [nv, 3] (+1 or -1).
   real(R8P)              :: wall_velocity(3,6)=0._R8P !< No-slip walls: wall velocity of each face.
   real(R8P)              :: wall_temperature(6)=0._R8P !< Isothermal walls: wall temperature of each face.
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize boundary conditions.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype flume_bc_object

contains
   ! public methods
   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(flume_bc_object), intent(in) :: self             !< Boundary conditions.
   character(len=:), allocatable      :: desc             !< Description.
   character(len=1), parameter        :: NL=new_line('a') !< New line character.
   integer(I4P)                       :: f                !< Counter.

   desc = mpih%myrankstr//'Boundary conditions main data'
   do f=1, 6
      desc = desc//NL//mpih%myrankstr//'  '//SECTION_NAME(f)//': '//trim(str(self%bc_type(f)))
   enddo
   endfunction description

   subroutine initialize(self, file_parameters, physics)
   !< Initialize boundary conditions.
   class(flume_bc_object),     intent(inout) :: self            !< Boundary conditions.
   type(file_ini),             intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   type(flume_physics_object), intent(in)    :: physics         !< Physics (for the inflow state conversion).

   print '(A)', mpih%myrankstr//'flume_bc_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters, physics=physics)
   print '(A)', self%description()
   print '(A)', mpih%myrankstr//'flume_bc_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, physics)
   !< Load config from file; `type` is required on every face and an unknown value is fatal.
   class(flume_bc_object),     intent(inout) :: self            !< Boundary conditions.
   type(file_ini),             intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   type(flume_physics_object), intent(in)    :: physics         !< Physics (for the inflow state conversion).
   character(999)                            :: buff            !< Option value buffer.
   character(:), allocatable                 :: bc_str          !< BC type string.
   real(R8P)                                 :: prim(8)         !< Inflow primitive state (first nprim used).
   integer(I4P)                              :: nprim           !< Primitive keys number: 5 (Euler) or 8 (MHD).
   integer(I4P)                              :: error           !< Error status.
   integer(I4P)                              :: d, f, k         !< Counters.

   if (allocated(self%q_inflow)) deallocate(self%q_inflow)
   if (allocated(self%wall_sign)) deallocate(self%wall_sign)
   allocate(self%q_inflow(physics%nv,6), self%wall_sign(physics%nv,3))
   self%q_inflow = 0._R8P
   ! wall rule (issue #41, section 3.6): the mirror state negates the wall-normal momentum; MHD (perfectly conducting
   ! reflecting wall, D-12) also negates the wall-normal magnetic field, psi stays even
   self%wall_sign = 1._R8P
   select case(physics%model)
   case(MODEL_EULER)
      nprim = 5_I4P
      do d=1, 3
         self%wall_sign(IQ_RU+d-1,d) = -1._R8P
      enddo
   case(MODEL_MHD, MODEL_MHD_GLM, MODEL_MHD_EGLM)
      nprim = 8_I4P
      do d=1, 3
         self%wall_sign(IQ_RU+d-1,d) = -1._R8P
         self%wall_sign(IQ_BX+d-1,d) = -1._R8P
      enddo
   case default
      nprim = 0_I4P
      call mpih%error_stop(msg=': no boundary conditions for physical model "'//physics%physical_model//'"')
   endselect
   prim = 0._R8P
   do f=1, 6
      call file_parameters%get(section_name=SECTION_NAME(f), option_name='type', val=buff, error=error)
      if (error > 0) call mpih%error_stop(msg=': failed to load ['//SECTION_NAME(f)//'].(type)')
      bc_str = trim(adjustl(strip_control(buff)))
      select case(bc_str)
      case(BC_EXTRAPOLATION_STR)
         self%bc_type(f) = BC_EXTRAPOLATION
      case(BC_INFLOW_STR)
         self%bc_type(f) = BC_INFLOW
         do k=1, nprim
            call file_parameters%get(section_name=SECTION_NAME(f), option_name=trim(INFLOW_KEY(k)), val=prim(k), &
                                     error=error)
            if (error > 0) &
               call mpih%error_stop(msg=': failed to load ['//SECTION_NAME(f)//'].('//trim(INFLOW_KEY(k))//')')
         enddo
         call primitive_state_to_conservative(model=physics%model, gamma=physics%gamma, prim=prim, q=self%q_inflow(:,f))
      case(BC_WALL_INVISCID_STR)
         self%bc_type(f) = BC_WALL_INVISCID
      case(BC_WALL_NOSLIP_STR, BC_WALL_ISOTHERMAL_STR)
         self%bc_type(f) = BC_WALL_NOSLIP
         if (bc_str == BC_WALL_ISOTHERMAL_STR) self%bc_type(f) = BC_WALL_ISOTHERMAL
         ! the wall velocity: tangential, optional (a resting wall by default)
         do d=1, 3
            call file_parameters%get(section_name=SECTION_NAME(f), option_name=WALL_VELOCITY_KEY(d), &
                                     val=self%wall_velocity(d,f), error=error)
            if (error > 0) self%wall_velocity(d,f) = 0._R8P
         enddo
         if (self%wall_velocity((f+1)/2,f) /= 0._R8P) &
            call mpih%error_stop(msg=': ['//SECTION_NAME(f)//'].('//WALL_VELOCITY_KEY((f+1)/2)//') is the normal '// &
                                     'component of the wall velocity and must be 0 (the wall does not move through '// &
                                     'the domain)')
         if (self%bc_type(f) == BC_WALL_ISOTHERMAL) then
            call file_parameters%get(section_name=SECTION_NAME(f), option_name='wall_temperature', &
                                     val=self%wall_temperature(f), error=error)
            if (error > 0) call mpih%error_stop(msg=': failed to load ['//SECTION_NAME(f)//'].(wall_temperature)')
            if (.not.(self%wall_temperature(f) > 0._R8P)) &
               call mpih%error_stop(msg=': ['//SECTION_NAME(f)//'].(wall_temperature) must be positive')
         endif
      case(BC_PERIODIC_STR)
         self%bc_type(f) = BC_PERIODIC
      case default
         call mpih%error_stop(msg=': unknown ['//SECTION_NAME(f)//'].(type) "'//bc_str//'"; expected one of '// &
                                  BC_EXTRAPOLATION_STR//', '//BC_INFLOW_STR//', '//BC_WALL_INVISCID_STR//', '// &
                                  BC_WALL_NOSLIP_STR//', '//BC_WALL_ISOTHERMAL_STR//', '//BC_PERIODIC_STR)
      endselect
   enddo
   do f=1, 5, 2
      if ((self%bc_type(f) == BC_PERIODIC) .neqv. (self%bc_type(f+1) == BC_PERIODIC)) &
         call mpih%error_stop(msg=': ['//SECTION_NAME(f)//'] and ['//SECTION_NAME(f+1)//'] must be both periodic or neither')
   enddo
   endsubroutine load_from_file
   ! public procedures
   pure function realm_edge_face(fec, face_kind) result(face)
   !< Face whose condition fills a crown row beyond the realm faces of the edge or corner `fec`: the first inflow face
   !< among them (an inflow ghost holds the inflow state whatever else it lies beyond), else the first physical one,
   !< 0 when every one is a seam (the row lies outside the forest).
   integer(I4P), intent(in) :: fec          !< Tree boundary fec of the row (7..26).
   integer(I4P), intent(in) :: face_kind(6) !< Kind of each realm face (-x, +x, -y, +y, -z, +z), BC_SEAM for a seam.
   integer(I4P)             :: face         !< Face (1..6), 0 if none.
   integer(I4P)             :: d, f         !< Axis, face.
   !$acc routine seq
   !$omp declare target

   face = 0_I4P
   do d=1_I4P, 3_I4P
      if (FEC_TO_DELTA(d,fec) == 0_I4P) cycle
      f = 2_I4P * d - 1_I4P + (FEC_TO_DELTA(d,fec) + 1_I4P) / 2_I4P
      if (face_kind(f) == BC_INFLOW) then
         face = f
         return
      endif
   enddo
   do d=1_I4P, 3_I4P
      if (FEC_TO_DELTA(d,fec) == 0_I4P) cycle
      f = 2_I4P * d - 1_I4P + (FEC_TO_DELTA(d,fec) + 1_I4P) / 2_I4P
      if (face_kind(f) /= BC_SEAM) then
         face = f
         return
      endif
   enddo
   endfunction realm_edge_face

   pure subroutine realm_edge_donor(face, face_kind, ni, nj, nk, i, j, k, iref, jref, kref)
   !< Donor of the ghost `(i, j, k)` filled by the condition of `face`: the cell mirrored about the face (walls) or the
   !< first interior cell along its normal (extrapolation); the indexes along the other axes are kept.
   integer(I4P), intent(in)  :: face             !< Face (1..6: -x, +x, -y, +y, -z, +z).
   integer(I4P), intent(in)  :: face_kind        !< Kind of the face.
   integer(I4P), intent(in)  :: ni, nj, nk       !< Grid dimensions.
   integer(I4P), intent(in)  :: i, j, k          !< Ghost cell.
   integer(I4P), intent(out) :: iref, jref, kref !< Donor cell.
   logical                   :: mirror           !< Wall: mirror, else clamp.
   !$acc routine seq
   !$omp declare target

   mirror = any(face_kind == [BC_WALL_INVISCID, BC_WALL_NOSLIP, BC_WALL_ISOTHERMAL])
   iref = i ; jref = j ; kref = k
   select case(face)
   case(1_I4P)
      iref = merge(1_I4P - i, 1_I4P, mirror)
   case(2_I4P)
      iref = merge(2_I4P * ni + 1_I4P - i, ni, mirror)
   case(3_I4P)
      jref = merge(1_I4P - j, 1_I4P, mirror)
   case(4_I4P)
      jref = merge(2_I4P * nj + 1_I4P - j, nj, mirror)
   case(5_I4P)
      kref = merge(1_I4P - k, 1_I4P, mirror)
   case(6_I4P)
      kref = merge(2_I4P * nk + 1_I4P - k, nk, mirror)
   endselect
   endsubroutine realm_edge_donor
   pure subroutine wall_noslip_ghost(nv, d, gamma, R, psi_energy, isothermal, wall_velocity, wall_temperature, q, qg)
   !< Ghost state of a no-slip wall (normal axis `d`) from its mirrored interior state `q`: the velocity reflected about
   !< the wall velocity, `u_g = 2 u_w - u`; the pressure mirrored; the temperature mirrored (adiabatic) or set to
   !< `2 T_w - T` (isothermal, the density then `p / (R T_g)`); the normal field odd, psi even.
   !<
   !< Model-agnostic: the field components are `IQ_BX..min(IQ_BZ, nv)` (none for Euler) and psi the slots beyond them;
   !< `psi_energy` (1 with EGLM, 0 otherwise) says whether the total energy holds `psi^2 / 2`.
   integer(I4P), intent(in)  :: nv                  !< Conservative variables number.
   integer(I4P), intent(in)  :: d                   !< Wall-normal axis.
   real(R8P),    intent(in)  :: gamma               !< Specific heats ratio.
   real(R8P),    intent(in)  :: R                   !< Gas constant.
   real(R8P),    intent(in)  :: psi_energy          !< 1 if the total energy holds psi^2 / 2 (EGLM), else 0.
   logical,      intent(in)  :: isothermal          !< Isothermal (else adiabatic) wall.
   real(R8P),    intent(in)  :: wall_velocity(3)    !< Wall velocity.
   real(R8P),    intent(in)  :: wall_temperature    !< Wall temperature (isothermal).
   real(R8P),    intent(in)  :: q(nv)               !< Mirrored interior state.
   real(R8P),    intent(out) :: qg(nv)              !< Ghost state.
   real(R8P)                 :: u(3), ug(3)         !< Interior and ghost velocity.
   real(R8P)                 :: e_mag               !< Magnetic (and EGLM psi) energy, unchanged by the mirror.
   real(R8P)                 :: p                   !< Mirrored pressure.
   real(R8P)                 :: rho_g               !< Ghost density.
   integer(I4P)              :: v                   !< Counter.
   !$acc routine seq
   !$omp declare target

   u = q(IQ_RU:IQ_RW) / q(IQ_R)
   ug = 2._R8P * wall_velocity - u
   e_mag = 0._R8P
   do v=IQ_BX, min(IQ_BZ, nv)
      e_mag = e_mag + 0.5_R8P * q(v) * q(v)
      qg(v) = q(v)
      if (v - IQ_BX + 1_I4P == d) qg(v) = -q(v)
   enddo
   do v=IQ_BZ + 1_I4P, nv
      e_mag = e_mag + psi_energy * 0.5_R8P * q(v) * q(v)
      qg(v) = q(v)
   enddo
   p = (gamma - 1._R8P) * (q(IQ_RE) - 0.5_R8P * q(IQ_R) * dot_product(u, u) - e_mag)
   rho_g = q(IQ_R)
   if (isothermal) rho_g = p / (R * (2._R8P * wall_temperature - p / (q(IQ_R) * R)))
   qg(IQ_R) = rho_g
   qg(IQ_RU:IQ_RW) = rho_g * ug
   qg(IQ_RE) = p / (gamma - 1._R8P) + 0.5_R8P * rho_g * dot_product(ug, ug) + e_mag
   endsubroutine wall_noslip_ghost
endmodule adam_flume_bc_object
