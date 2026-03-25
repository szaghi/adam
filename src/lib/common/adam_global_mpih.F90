!< ADAM, global MPI handler singleton — single program-scope mpih_object instance.
module adam_global_mpih
!< ADAM, global MPI handler singleton — single program-scope mpih_object instance.

! ADAM classes, libraries, parameters
use :: adam_mpih_object, only : mpih_object

implicit none
private
public :: mpih

type(mpih_object), target :: mpih !< Program-scope MPI handler singleton.
endmodule adam_global_mpih
