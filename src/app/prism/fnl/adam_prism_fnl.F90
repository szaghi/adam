!< ADAM, Maxwell application solver, GPU (FNL) backend.
program adam_prism_fnl
!< ADAM, Maxwell application solver, GPU (FNL) backend.
!<
!< The library's generic driver (`adam_run_library%run_forest`, issue #80 P3) on this application's realm type: one
!< realm for a single input file, one per realm of a forest manifest (the first argument, default `input.ini`).

use :: adam_run_library, only : run_forest
use :: adam_prism_fnl_object, only : prism_fnl_object

implicit none

type(prism_fnl_object) :: prototype !< The realm type to run (only its type is used).

call run_forest(prototype=prototype)
endprogram adam_prism_fnl
