!< ADAM, compressible MHD application solver (FLUME), FNL (OpenACC / OpenMP offload) backend.
program adam_flume_fnl
!< ADAM, compressible MHD application solver (FLUME), FNL (OpenACC / OpenMP offload) backend.
!<
!< Two driver paths, auto-detected from the input file: a plain FLUME `input.ini` runs one realm through
!< `forest%simulate`; a `forest.ini` manifest (a `[forest] realms_number = N` section) runs N realms through
!< `forest%simulate_from_manifest`. See `docs/concepts/forest.md` for the manifest schema.

! ADAM classes, libraries, parameters
use :: adam_forest_object,    only : forest_object
use :: adam_forest_manifest,  only : forest_manifest_t, is_forest_manifest, read_forest_manifest
use :: adam_realm_object,     only : realm_pointer
! FLUME modules
use :: adam_flume_fnl_object, only : flume_fnl_object

implicit none

type(flume_fnl_object), allocatable, target :: realm(:)          !< Realms: 1 for a single input, N for a manifest.
type(realm_pointer),    allocatable         :: realm_pointers(:) !< A pointer to each realm: what the forest takes (#96).
type(forest_object)                         :: forest            !< Orchestrator that drives the realm array.
type(forest_manifest_t)                     :: manifest          !< Parsed manifest (manifest path only).
integer                                     :: na                !< Number of command line arguments.
character(999)                              :: input_file_name   !< Input file name.

na = command_argument_count()
if (na == 0) then
   input_file_name = 'input.ini'
else
   call get_command_argument(1, input_file_name)
   input_file_name = trim(adjustl(input_file_name))
endif

if (is_forest_manifest(trim(input_file_name))) then
   call read_forest_manifest(filename=trim(input_file_name), manifest=manifest)
   allocate(realm(manifest%realms_number))
   call point_realms
   call forest%simulate_from_manifest(realm=realm_pointers, manifest=manifest)
else
   allocate(realm(1))
   call point_realms
   call forest%simulate(realm=realm_pointers, filename=trim(input_file_name))
endif
contains
   subroutine point_realms
   !< Point one realm_pointer at each realm. The indexing here is on the monomorphic array `realm(:)`, which nvfortran
   !< addresses correctly; the forest then never indexes a polymorphic array (issue #96).
   integer :: i !< Counter.

   allocate(realm_pointers(size(realm)))
   do i = 1, size(realm)
      realm_pointers(i)%ptr => realm(i)
   enddo
   endsubroutine point_realms
endprogram adam_flume_fnl
