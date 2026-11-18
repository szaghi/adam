!< ADAM, the generic driver of an application: realms of the application's type, the forest, the run.
module adam_run_library
!< ADAM, the generic driver of an application: realms of the application's type, the forest, the run.
!<
!< Every application's main program did the same (issue #80 P3, §5(f)): read the input file name (default
!< `input.ini`), detect a forest manifest, allocate one realm (single input) or one per manifest realm, point the
!< forest at them and simulate. `run_forest` does it for any realm type, given a prototype of it:
!<```fortran
!< program adam_flume_cpu
!< use :: adam_run_library,      only : run_forest
!< use :: adam_flume_cpu_object, only : flume_cpu_object
!< implicit none
!< type(flume_cpu_object) :: prototype !< The realm type to run (only its type is used).
!<
!< call run_forest(prototype=prototype)
!< endprogram adam_flume_cpu
!<```

use :: adam_forest_manifest, only : forest_manifest_t, is_forest_manifest, read_forest_manifest
use :: adam_forest_object,   only : forest_object
use :: adam_realm_object,    only : realm_object, realm_pointer

implicit none
private
public :: run_forest

contains
   subroutine run_forest(prototype)
   !< Run the application whose realm type is the dynamic type of `prototype`: one realm for a single input file, one
   !< per realm of a forest manifest (the first command-line argument, default `input.ini`).
   !<
   !< Each realm is a scalar of its own, allocated through the forest's pointer with `mold=prototype` (default-initialized,
   !< the prototype's values are not copied): there is no array of `class(realm_object)`, which nvfortran indexes wrongly
   !< (issue #96), and, as the program-scope arrays they replace, the realms are never deallocated, so no finalization or
   !< device free runs after the forest has finalized MPI.
   class(realm_object), intent(in)  :: prototype         !< A realm of the type to run.
   type(realm_pointer), allocatable :: realm_pointers(:) !< The realms, one pointer each: what the forest takes.
   type(forest_object)              :: forest            !< Orchestrator that drives the realms.
   type(forest_manifest_t)          :: manifest          !< Parsed manifest (manifest path only).
   character(999)                   :: input_file_name   !< Input file name.
   integer                          :: realms_number     !< Realms: 1 for a single input, N for a manifest.
   integer                          :: i                 !< Counter.

   input_file_name = 'input.ini'
   if (command_argument_count() > 0) call get_command_argument(1, input_file_name)
   input_file_name = trim(adjustl(input_file_name))
   realms_number = 1
   if (is_forest_manifest(trim(input_file_name))) then
      call read_forest_manifest(filename=trim(input_file_name), manifest=manifest)
      realms_number = manifest%realms_number
   endif
   allocate(realm_pointers(realms_number))
   do i=1, realms_number
      allocate(realm_pointers(i)%ptr, mold=prototype)
   enddo
   if (is_forest_manifest(trim(input_file_name))) then
      call forest%simulate_from_manifest(realm=realm_pointers, manifest=manifest)
   else
      call forest%simulate(realm=realm_pointers, filename=trim(input_file_name))
   endif
   endsubroutine run_forest
endmodule adam_run_library
