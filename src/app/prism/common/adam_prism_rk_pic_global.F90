!< ADAM, global PRISM RK-PIC singleton — single program-scope prism_rk_pic_object instance.
module adam_prism_rk_pic_global
!< ADAM, global PRISM RK-PIC singleton — single program-scope prism_rk_pic_object instance.

! PRISM modules
use :: adam_prism_rk_pic_object, only : prism_rk_pic_object

implicit none
private
public :: rk_pic

type(prism_rk_pic_object), target :: rk_pic !< Program-scope PRISM RK-PIC singleton.
endmodule adam_prism_rk_pic_global
