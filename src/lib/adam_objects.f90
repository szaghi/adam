!< ADAM, objects classes.
module adam_objects
!< ADAM, objects classes.

use adam_dictionary_node_object , only : KEY_LEN, destroy_dictionary_node, dictionary_node_object
use adam_dictionary_object, only : dictionary_object, iterator_interface

implicit none
private
public :: KEY_LEN
public :: destroy_dictionary_node
public :: dictionary_node_object
public :: dictionary_object
public :: iterator_interface
endmodule adam_objects
