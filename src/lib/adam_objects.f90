!< ADAM, objects classes.
module adam_objects
!< ADAM, objects classes.

use adam_dictionary_node_object
use adam_dictionary_object
use adam_hash_table_object
use adam_tree_topology

implicit none
private
public :: KEY_LEN
public :: destroy_dictionary_node
public :: key_str, key_int
public :: dictionary_node_object
public :: dictionary_object
public :: iterator_interface
public :: hash_table_object
public :: coordinates_to_morton, morton_to_coordinates
public :: child
public :: child_local
public :: first_at_level
public :: last_at_level
public :: level
public :: parent
public :: path
public :: siblings
endmodule adam_objects
