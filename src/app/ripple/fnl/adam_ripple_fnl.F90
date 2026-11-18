!< RIPPLE, scalar linear advection on the ADAM library, FNL (OpenACC) backend.
program adam_ripple_fnl
!< RIPPLE, scalar linear advection on the ADAM library, FNL (OpenACC) backend.
!<
!< The library's generic driver (`adam_run_library%run_forest`, issue #80 P3) on this application's realm type: one
!< realm for a single input file, one per realm of a forest manifest (the first argument, default `input.ini`).

use :: adam_run_library, only : run_forest
use :: adam_ripple_fnl_object, only : ripple_fnl_object

implicit none

type(ripple_fnl_object) :: prototype !< The realm type to run (only its type is used).

call run_forest(prototype=prototype)
endprogram adam_ripple_fnl
