!< ADAM, objects classes.
module adam_objects
!< ADAM, objects classes.

use adam_tree_node_object
use adam_tree_bucket_object
use adam_tree_object

implicit none
private
! public :: KEY_LEN
public :: destroy_tree_node
! public :: key_str, key_int
public :: tree_node_object, TO_BE_REFINED, TO_BE_DEREFINED
public :: tree_bucket_object, len
public :: iterator_interface
public :: tree_object
endmodule adam_objects
