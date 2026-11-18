!< ADAM, FLUME time handler class definition.
module adam_flume_time_object
!< ADAM, FLUME time handler class definition.
!<
!< The library time handler as is (`adam_time_object`): step counter, time, time step and the required `[time]` keys
!< `it_max`, `time_max`, `CFL`; FLUME adds nothing (issue #80 P3).

use :: adam_time_object, only : time_object

implicit none
private
public :: flume_time_object

type, extends(time_object) :: flume_time_object
   !< FLUME time handler class definition.
endtype flume_time_object
endmodule adam_flume_time_object
