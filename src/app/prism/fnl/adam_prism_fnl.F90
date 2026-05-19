!< ADAM, Maxwell application solver, GPU (FNL) backend.
program adam_prism_fnl
!< ADAM, Maxwell application solver, GPU (FNL) backend.
!<
!< Two driver paths, auto-detected from the input file:
!<
!<   * **Single-INI (legacy)**: the input file is a plain PRISM `input.ini`.
!<     The driver allocates `realm(1)` and calls
!<     `forest%simulate(realm, filename=...)`. Bit-identical to the
!<     Phase C single-realm rmf behaviour.
!<
!<   * **Manifest**: the input file is a `forest.ini` manifest pointing
!<     at one PRISM INI per realm (Phase D of issue #10, design in
!<     issue #13). The driver reads `realms_number` and per-realm INI
!<     paths, allocates `realm(realms_number)`, and calls
!<     `forest%simulate_from_manifest(realm, manifest)`.
!<
!< Detection is done by `is_forest_manifest`: the file is a manifest iff
!< it contains a `[forest] realms_number = N` (N>=1) section. A plain
!< PRISM input has no `[forest]` section and falls through to the
!< legacy path.
!<
!< The seven legacy `adam_*_global` shims (`grid`, `field`, ...)
!< continue to alias the LAST realm initialized inside
!< `forest%initialize` (`realm(N)`); for N=1 this is `realm(1)` and
!< matches Phase C behaviour. For N>1 the shims alias `realm(N)` —
!< that is the deferred C.3 gap documented in [#10] Phase C status;
!< the gap is fully load-bearing once N>1 actually runs and is the
!< next Phase D follow-up after the rmf-2realm regression case (D.4).

! ADAM common library — forest orchestrator + manifest parser
use :: adam_forest_object,    only : forest_object
use :: adam_forest_manifest,  only : forest_manifest_t, is_forest_manifest, read_forest_manifest
! PRISM modules
use :: adam_prism_fnl_object, only : prism_fnl_object

implicit none

type(prism_fnl_object), allocatable :: realm(:)         !< Realm array; size set by input shape (1 for single-INI, N for manifest).
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
endprogram adam_prism_fnl
