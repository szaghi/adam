!< ADAM, global PRISM leapfrog-PIC singleton — single program-scope prism_leapfrog_pic_object instance.
module adam_prism_leapfrog_pic_global
!< ADAM, global PRISM leapfrog-PIC singleton — single program-scope prism_leapfrog_pic_object instance.

! PRISM modules
use :: adam_prism_leapfrog_pic_object, only : prism_leapfrog_pic_object

implicit none
private
public :: leapfrog_pic

type(prism_leapfrog_pic_object), target :: leapfrog_pic !< Program-scope PRISM leapfrog-PIC singleton.
endmodule adam_prism_leapfrog_pic_global
