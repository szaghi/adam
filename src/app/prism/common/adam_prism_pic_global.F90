!< ADAM, global PRISM PIC singleton — single program-scope prism_pic_object instance.
module adam_prism_pic_global
!< ADAM, global PRISM PIC singleton — single program-scope prism_pic_object instance.

! PRISM modules
use :: adam_prism_pic_object, only : prism_pic_object

implicit none
private
public :: pic

type(prism_pic_object), target :: pic !< Program-scope PRISM PIC singleton.
endmodule adam_prism_pic_global
