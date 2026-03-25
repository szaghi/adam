!< ADAM, global grid singleton — single program-scope grid_object instance.
module adam_global_grid
!< ADAM, global grid singleton — single program-scope grid_object instance.

! ADAM classes, libraries, parameters
use :: adam_grid_object, only : grid_object

implicit none
private
public :: grid

type(grid_object), target :: grid !< Program-scope grid singleton.
endmodule adam_global_grid
