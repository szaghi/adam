!< ADAM, forest manifest parser — reads forest.ini into a transient struct.
module adam_forest_manifest
!< ADAM, forest manifest parser — reads forest.ini into a transient struct.
!<
!< Phase D of issue #10 (design in issue #13). A **forest manifest** is the
!< small INI that lists the realms making up a forest and how they couple.
!< This module reads such a manifest and returns the data the driver needs
!< to allocate the realm array and the topology data the forest needs to
!< populate per-realm `maps%inter_realm_neighbors` after each realm has
!< been initialized.
!<
!< Manifest shape (see issue #13 for full spec):
!<
!<```
!<   [forest]
!<   realms_number = 2
!<
!<   [realm.1]
!<   ini = realm_1.ini       ; path relative to forest.ini's directory
!<
!<   [realm.2]
!<   ini = realm_2.ini
!<
!<   [forest.topology]
!<   inter_realm_faces_number = 1
!<
!<   [forest.topology.face_1]
!<   realm_a  = 1
!<   face_a   = +x
!<   realm_b  = 2
!<   face_b   = -x
!<   coupling = mirror       ; mirror | periodic | interpolate
!<```
!<
!< Single-INI detection: callers MUST check whether the file actually IS a
!< manifest before invoking this parser — a plain PRISM input.ini has no
!< [forest] section and `is_forest_manifest` returns .false. without
!< loading any state. This lets the driver short-circuit to the legacy
!< single-INI path for current rmf usage.
!<
!< Architectural invariants:
!<
!<   * **No singletons read.** This module is a pure data-transform on a
!<     file path; it does not touch mpih, grid, or any other program-scope
!<     state. That keeps the manifest parser callable BEFORE any of those
!<     singletons have been initialized.
!<   * **No realm-typed components.** The forest_manifest_t struct carries
!<     only intrinsic-typed components — realm paths as character strings,
!<     face/coupling codes as integers. The realm objects themselves are
!<     allocated by the driver from `realms_number`; the topology entries
!<     are translated to `inter_realm_neighbor_t` instances on the
!<     realm-side later (Phase D follow-up).

! ADAM classes, libraries, parameters
use :: adam_maps_object, only : FACE_X_MAX, FACE_X_MIN, FACE_Y_MAX, FACE_Y_MIN, FACE_Z_MAX, FACE_Z_MIN, &
                                 COUPLING_MIRROR, COUPLING_PERIODIC, COUPLING_INTERPOLATE
! third party modules
use :: finer, only : file_ini
use :: penf,  only : I4P, str

implicit none
private
public :: forest_manifest_t
public :: forest_face_pair_t
public :: is_forest_manifest
public :: read_forest_manifest

integer(I4P), parameter :: MAX_INI_PATH_LEN = 256_I4P !< Buffer length for realm INI paths.

type :: forest_face_pair_t
   !< One inter-realm face coupling as parsed from a [forest.topology.face_<n>] block.
   !<
   !< The struct is symmetric: both realms know about the same coupling. Later
   !< (during forest%initialize) each realm's maps%inter_realm_neighbors is
   !< populated with TWO inter_realm_neighbor_t entries per face pair — one
   !< from realm_a's perspective and one from realm_b's perspective.
   integer(I4P) :: realm_a  = 0_I4P            !< First realm in the pair (1-based).
   integer(I4P) :: face_a   = 0_I4P            !< Face code on realm_a (FACE_X_MAX..FACE_Z_MIN).
   integer(I4P) :: realm_b  = 0_I4P            !< Second realm in the pair.
   integer(I4P) :: face_b   = 0_I4P            !< Face code on realm_b.
   integer(I4P) :: coupling = COUPLING_MIRROR  !< Coupling kind (default MIRROR for the first use case).
endtype forest_face_pair_t

type :: forest_manifest_t
   !< Parsed contents of a forest.ini manifest.
   !<
   !< Transient data structure: built by `read_forest_manifest`, consumed
   !< once by the driver (to allocate `realm(:)`) and the forest (to wire
   !< inter-realm topology), then discarded. NOT stored as a long-lived
   !< component on `forest_object` (which is behavior-only by design).
   integer(I4P)                                :: realms_number = 0_I4P !< Number of realms in the forest.
   character(len=MAX_INI_PATH_LEN), allocatable :: realm_ini(:)         !< Per-realm INI path (1:realms_number).
   type(forest_face_pair_t),         allocatable :: face_pairs(:)        !< Inter-realm couplings; unallocated when none declared.
endtype forest_manifest_t

contains
   function is_forest_manifest(filename) result(yes)
   !< Return .true. iff the file is recognizable as a forest manifest.
   !<
   !< Detection rule: the file loads successfully via FiNeR AND contains a
   !< `[forest] realms_number = N` entry with N >= 1. A plain PRISM input.ini
   !< has no `[forest]` section and returns .false. — the driver then falls
   !< back to the legacy single-INI path.
   character(*), intent(in) :: filename !< Path to test.
   logical                  :: yes      !< Detection result.
   type(file_ini)           :: handle   !< Local FiNeR handle (discarded after probe).
   integer(I4P)             :: realms_number !< Probed value (ignored beyond presence check).
   integer(I4P)             :: error         !< FiNeR error code.

   yes = .false.
   call handle%initialize(filename=trim(filename))
   call handle%load
   call handle%get(section_name='forest', option_name='realms_number', val=realms_number, error=error)
   if (error == 0 .and. realms_number >= 1_I4P) yes = .true.
   call handle%free
   endfunction is_forest_manifest

   subroutine read_forest_manifest(filename, manifest)
   !< Parse a forest manifest INI into `manifest`.
   !<
   !< Caller is responsible for confirming the file IS a manifest via
   !< `is_forest_manifest` first; this routine error-stops on any malformed
   !< or missing required entry (so a caller that gets here trusts the
   !< parse to succeed).
   character(*),               intent(in)  :: filename !< Path to the manifest INI.
   type(forest_manifest_t),    intent(out) :: manifest !< Parsed result.
   type(file_ini)                          :: handle   !< FiNeR handle for the manifest.
   character(MAX_INI_PATH_LEN)             :: buff     !< Character buffer for path / kind reads.
   integer(I4P)                            :: error    !< FiNeR error code.
   integer(I4P)                            :: faces_number !< Count of [forest.topology.face_<n>] blocks.
   integer(I4P)                            :: is       !< Realm index counter.
   integer(I4P)                            :: f        !< Face-pair counter.
   character(:), allocatable               :: sname    !< Section-name builder.

   call handle%initialize(filename=trim(filename))
   call handle%load

   ! [forest] realms_number — required
   call handle%get(section_name='forest', option_name='realms_number', val=manifest%realms_number, error=error)
   if (error /= 0 .or. manifest%realms_number < 1_I4P) &
      error stop 'adam_forest_manifest: [forest] realms_number missing or < 1 in '//trim(filename)
   allocate(manifest%realm_ini(manifest%realms_number))

   ! [realm.<is>] ini — required, once per realm
   do is = 1_I4P, manifest%realms_number
      sname = 'realm.'//trim(str(is, .true.))
      call handle%get(section_name=sname, option_name='ini', val=buff, error=error)
      if (error /= 0) &
         error stop 'adam_forest_manifest: missing ['//sname//'] ini = ... in '//trim(filename)
      manifest%realm_ini(is) = trim(adjustl(buff))
   enddo

   ! [forest.topology] inter_realm_faces_number — optional; absent means no inter-realm coupling
   faces_number = 0_I4P
   call handle%get(section_name='forest.topology', option_name='inter_realm_faces_number', val=faces_number, error=error)
   if (error == 0 .and. faces_number > 0_I4P) then
      allocate(manifest%face_pairs(faces_number))
      do f = 1_I4P, faces_number
         sname = 'forest.topology.face_'//trim(str(f, .true.))
         call read_face_pair(handle=handle, section_name=sname, pair=manifest%face_pairs(f), context=filename)
      enddo
   endif

   call handle%free
   endsubroutine read_forest_manifest

   subroutine read_face_pair(handle, section_name, pair, context)
   !< Read one [forest.topology.face_<n>] section into a face-pair struct.
   type(file_ini),            intent(inout) :: handle       !< FiNeR handle.
   character(*),              intent(in)    :: section_name !< Section name to read from.
   type(forest_face_pair_t),  intent(out)   :: pair         !< Parsed face pair.
   character(*),              intent(in)    :: context      !< Caller's filename for error messages.
   character(MAX_INI_PATH_LEN)              :: buff         !< Character buffer.
   integer(I4P)                             :: error        !< FiNeR error code.

   call handle%get(section_name=section_name, option_name='realm_a', val=pair%realm_a, error=error)
   if (error /= 0) error stop 'adam_forest_manifest: ['//trim(section_name)//'] realm_a missing in '//trim(context)
   call handle%get(section_name=section_name, option_name='realm_b', val=pair%realm_b, error=error)
   if (error /= 0) error stop 'adam_forest_manifest: ['//trim(section_name)//'] realm_b missing in '//trim(context)
   call handle%get(section_name=section_name, option_name='face_a', val=buff, error=error)
   if (error /= 0) error stop 'adam_forest_manifest: ['//trim(section_name)//'] face_a missing in '//trim(context)
   pair%face_a = parse_face(trim(adjustl(buff)), section_name, context)
   call handle%get(section_name=section_name, option_name='face_b', val=buff, error=error)
   if (error /= 0) error stop 'adam_forest_manifest: ['//trim(section_name)//'] face_b missing in '//trim(context)
   pair%face_b = parse_face(trim(adjustl(buff)), section_name, context)
   ! coupling — optional, default MIRROR
   call handle%get(section_name=section_name, option_name='coupling', val=buff, error=error)
   if (error == 0) pair%coupling = parse_coupling(trim(adjustl(buff)), section_name, context)
   endsubroutine read_face_pair

   function parse_face(text, section_name, context) result(code)
   !< Translate "+x" / "-x" / "+y" / "-y" / "+z" / "-z" to the face code constant.
   character(*), intent(in) :: text         !< Face spec from the INI.
   character(*), intent(in) :: section_name !< Section name (for error message).
   character(*), intent(in) :: context      !< Manifest filename (for error message).
   integer(I4P)             :: code         !< Resolved face code.

   select case (text)
   case ('+x', '+X', 'x+', 'X+'); code = FACE_X_MAX
   case ('-x', '-X', 'x-', 'X-'); code = FACE_X_MIN
   case ('+y', '+Y', 'y+', 'Y+'); code = FACE_Y_MAX
   case ('-y', '-Y', 'y-', 'Y-'); code = FACE_Y_MIN
   case ('+z', '+Z', 'z+', 'Z+'); code = FACE_Z_MAX
   case ('-z', '-Z', 'z-', 'Z-'); code = FACE_Z_MIN
   case default
      error stop 'adam_forest_manifest: ['//trim(section_name)//'] unknown face "'//trim(text)//'" in '//trim(context)
   endselect
   endfunction parse_face

   function parse_coupling(text, section_name, context) result(code)
   !< Translate "mirror" / "periodic" / "interpolate" to the coupling code constant.
   character(*), intent(in) :: text         !< Coupling spec from the INI.
   character(*), intent(in) :: section_name !< Section name (for error message).
   character(*), intent(in) :: context      !< Manifest filename (for error message).
   integer(I4P)             :: code         !< Resolved coupling code.

   select case (text)
   case ('mirror',      'Mirror',      'MIRROR'     ); code = COUPLING_MIRROR
   case ('periodic',    'Periodic',    'PERIODIC'   ); code = COUPLING_PERIODIC
   case ('interpolate', 'Interpolate', 'INTERPOLATE'); code = COUPLING_INTERPOLATE
   case default
      error stop 'adam_forest_manifest: ['//trim(section_name)//'] unknown coupling "'//trim(text)//'" in '//trim(context)
   endselect
   endfunction parse_coupling
endmodule adam_forest_manifest
