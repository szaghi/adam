!< ADAM, objects classes.
module adam_objects
!< ADAM, objects classes.

use adam_dictionary_node_object , only : KEY_LEN, destroy_dictionary_node, key_str, key_int, dictionary_node_object
use adam_dictionary_object, only : dictionary_object, iterator_interface
use adam_hash_table_object, only : hash_table_object

implicit none
private
public :: KEY_LEN
public :: destroy_dictionary_node
public :: key_str, key_int
public :: dictionary_node_object
public :: dictionary_object
public :: iterator_interface
public :: hash_table_object
endmodule adam_objects
