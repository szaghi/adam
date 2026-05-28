!< ADAM, Maxwell application solver, CPU backend.
program adam_prism_cpu
!< ADAM, Maxwell application solver, CPU backend.
!<
!< Two driver paths, auto-detected from the input file:
!<
!<   * **Single-INI**: the input file is a plain PRISM `input.ini`. The
!<     driver allocates `realm(1)` and calls
!<     `forest%simulate(realm, filename=...)`. N=1 fast path; no inter-
!<     realm machinery is exercised.
!<
!<   * **Manifest**: the input file is a `forest.ini` manifest pointing
!<     at one PRISM INI per realm. The driver reads `realms_number` and
!<     per-realm INI paths, allocates `realm(realms_number)`, and calls
!<     `forest%simulate_from_manifest(realm, manifest)`.
!<
!< See `docs/guide/forest.md` for the manifest schema and the multi-realm
!< coupling-cadence semantics.
!<
!< Detection is done by `is_forest_manifest`: the file is a manifest iff
!< it contains a `[forest] realms_number = N` (N>=1) section. A plain
!< PRISM input has no `[forest]` section and falls through to the
!< single-INI path.

! ADAM common library — forest orchestrator + manifest parser
use :: adam_forest_object,    only : forest_object
use :: adam_forest_manifest,  only : forest_manifest_t, is_forest_manifest, read_forest_manifest
! PRISM modules
use :: adam_prism_cpu_object, only : prism_cpu_object

implicit none

type(prism_cpu_object), allocatable :: realm(:)         !< Realm array; size set by input shape (1 for single-INI, N for manifest).
type(forest_object)                 :: forest           !< Orchestrator that drives the realm array.
type(forest_manifest_t)             :: manifest         !< Parsed manifest (populated only in manifest path).
integer                             :: na               !< Number of command line arguments.
character(999)                      :: input_file_name  !< Input file name.

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
   call forest%simulate_from_manifest(realm=realm, manifest=manifest)
else
   allocate(realm(1))
   call forest%simulate(realm=realm, filename=trim(input_file_name))
endif
endprogram adam_prism_cpu
