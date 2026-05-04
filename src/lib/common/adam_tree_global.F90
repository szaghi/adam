!< ADAM, global tree singleton — single program-scope tree_object instance.
module adam_tree_global
!< ADAM, global tree singleton — single program-scope tree_object instance.

! ADAM classes, libraries, parameters
use :: adam_tree_object, only : tree_object

implicit none
private
public :: tree

type(tree_object), target :: tree !< Program-scope tree singleton.
endmodule adam_tree_global
