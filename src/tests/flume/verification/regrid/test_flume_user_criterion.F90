!< A registered refinement criterion for the test of `[amr_marker_N] mode = 6` (issue #110, P4).
module test_flume_user_criterion_box
!< The built-in primitive-box marker written as a registered criterion (a type-bound procedure needs a module).
use :: adam_amr_criterion_object, only : amr_criterion_object
use :: adam_parameters,           only : TO_BE_DEREFINED, TO_BE_REFINED, TO_NOT_TOUCH
use :: adam_realm_object,         only : realm_object
use :: penf,                      only : I4P, R8P

implicit none
private
public :: box_criterion

type, extends(amr_criterion_object) :: box_criterion
   !< The built-in primitive-box marker as a registered criterion.
   real(R8P)    :: box_emin(3)=[0.75_R8P, 0._R8P, 0._R8P] !< Box minimum corner.
   real(R8P)    :: box_emax(3)=[1._R8P, 1._R8P, 1._R8P]   !< Box maximum corner.
   integer(I4P) :: target_level=3_I4P                     !< Refine below this level inside the box.
   contains
      procedure, pass(self) :: vote
endtype box_criterion

contains
   function vote(self, realm, b) result(mark)
   !< The vote of the built-in primitive-box marker for block `b`.
   class(box_criterion), intent(inout) :: self        !< The criterion.
   class(realm_object),  intent(in)    :: realm       !< The realm.
   integer(I4P),         intent(in)    :: b           !< Local block index.
   integer(I4P)                        :: mark        !< The vote.
   real(R8P)                           :: centroid(3) !< Block centroid.
   integer(I4P)                        :: level       !< Block level.

   centroid = 0.5_R8P * (realm%adam%field%emin(:,b) + realm%adam%field%emax(:,b))
   level = realm%adam%tree%level(realm%adam%field%code(b))
   mark = TO_BE_DEREFINED
   if (all(centroid >= self%box_emin) .and. all(centroid <= self%box_emax)) then
      if (level < self%target_level) then
         mark = TO_BE_REFINED
      elseif (level == self%target_level) then
         mark = TO_NOT_TOUCH
      endif
   endif
   endfunction vote
endmodule test_flume_user_criterion_box

!< FLUME driver with a registered refinement criterion (issue #110, P4): the test of `[amr_marker_N] mode = 6`.
program test_flume_user_criterion
!< FLUME driver with a registered refinement criterion (issue #110, P4): the test of `[amr_marker_N] mode = 6`.
!<
!< The FLUME CPU driver (adam_flume_cpu.F90) plus one line per realm, `register_amr_criterion`: the way a user adds a
!< refinement criterion without editing the library or FLUME. The criterion, `box-criterion`, votes as the built-in
!< primitive-box marker does with the box x >= 0.75 and target level 3: refine inside the box below the target level,
!< keep at it, derefine elsewhere. A run whose marker is `mode = 6, name = box-criterion` must therefore be bitwise the
!< same run with the built-in box marker (`multirealm/check.sh` leg 10, on one realm and on a 2-realm forest).
use :: adam_forest_object,            only : forest_object
use :: adam_forest_manifest,          only : forest_manifest_t, is_forest_manifest, read_forest_manifest
use :: adam_realm_object,             only : realm_pointer
use :: adam_flume_cpu_object,         only : flume_cpu_object
use :: test_flume_user_criterion_box, only : box_criterion

implicit none

type(flume_cpu_object), allocatable, target :: realm(:)          !< Realms: 1 for a single input, N for a manifest.
type(realm_pointer),    allocatable         :: realm_pointers(:) !< A pointer to each realm.
type(forest_object)                         :: forest            !< Orchestrator that drives the realm array.
type(forest_manifest_t)                     :: manifest          !< Parsed manifest (manifest path only).
type(box_criterion),                 target :: box               !< The registered criterion.
integer                                     :: na                !< Number of command line arguments.
character(999)                              :: input_file_name   !< Input file name.
integer                                     :: i                 !< Counter.

na = command_argument_count()
if (na == 0) then
   input_file_name = 'input.ini'
else
   call get_command_argument(1, input_file_name)
   input_file_name = trim(adjustl(input_file_name))
endif
box%name = 'box-criterion'

if (is_forest_manifest(trim(input_file_name))) then
   call read_forest_manifest(filename=trim(input_file_name), manifest=manifest)
   allocate(realm(manifest%realms_number))
   call point_realms
   call forest%simulate_from_manifest(realm=realm_pointers, manifest=manifest)
else
   allocate(realm(1))
   call point_realms
   call forest%simulate(realm=realm_pointers, filename=trim(input_file_name))
endif
contains
   subroutine point_realms
   !< Point one realm_pointer at each realm, and register the criterion with every realm.
   allocate(realm_pointers(size(realm)))
   do i = 1, size(realm)
      realm_pointers(i)%ptr => realm(i)
      call realm(i)%register_amr_criterion(box)
   enddo
   endsubroutine point_realms

endprogram test_flume_user_criterion
