!< ADAM, PRISM absorbing-layer geometry helpers shared by fWLayer and PML.
module adam_prism_absorbing_layer_geometry

use :: penf, only : I4P, R8P

implicit none
private
public :: compute_absorbing_face_range

contains

   pure subroutine compute_absorbing_face_range(face, width, ni, nj, nk, domain_emin, domain_emax, block_emin, block_emax, &
                                                dxyz, ni_range, nj_range, nk_range, cells, last_center_distance, profile_extent)
   !< Return the active cell range of an absorbing layer on a block face.
   integer(I4P), intent(in)  :: face
   real(R8P),    intent(in)  :: width
   integer(I4P), intent(in)  :: ni, nj, nk
   real(R8P),    intent(in)  :: domain_emin(3)
   real(R8P),    intent(in)  :: domain_emax(3)
   real(R8P),    intent(in)  :: block_emin(3)
   real(R8P),    intent(in)  :: block_emax(3)
   real(R8P),    intent(in)  :: dxyz(3)
   integer(I4P), intent(out) :: ni_range(2)
   integer(I4P), intent(out) :: nj_range(2)
   integer(I4P), intent(out) :: nk_range(2)
   integer(I4P), intent(out) :: cells
   real(R8P),    intent(out) :: last_center_distance
   real(R8P),    intent(out) :: profile_extent
   real(R8P)                 :: ds
   real(R8P)                 :: distance0
   real(R8P)                 :: remaining

   ni_range = 0_I4P
   nj_range = 0_I4P
   nk_range = 0_I4P
   cells = 0_I4P
   last_center_distance = 0._R8P
   profile_extent = 0._R8P

   if (width <= 0._R8P) return

   select case (face)
   case (1_I4P)
      ds = dxyz(1)
      distance0 = block_emin(1) - domain_emin(1)
      remaining = width - distance0
      if (remaining <= 0._R8P) return
      cells = min(ni, ceiling(remaining / ds, kind=I4P))
      ni_range = [1_I4P, cells]
      nj_range = [1_I4P, nj]
      nk_range = [1_I4P, nk]
   case (2_I4P)
      ds = dxyz(1)
      distance0 = domain_emax(1) - block_emax(1)
      remaining = width - distance0
      if (remaining <= 0._R8P) return
      cells = min(ni, ceiling(remaining / ds, kind=I4P))
      ni_range = [ni - cells + 1_I4P, ni]
      nj_range = [1_I4P, nj]
      nk_range = [1_I4P, nk]
   case (3_I4P)
      ds = dxyz(2)
      distance0 = block_emin(2) - domain_emin(2)
      remaining = width - distance0
      if (remaining <= 0._R8P) return
      cells = min(nj, ceiling(remaining / ds, kind=I4P))
      ni_range = [1_I4P, ni]
      nj_range = [1_I4P, cells]
      nk_range = [1_I4P, nk]
   case (4_I4P)
      ds = dxyz(2)
      distance0 = domain_emax(2) - block_emax(2)
      remaining = width - distance0
      if (remaining <= 0._R8P) return
      cells = min(nj, ceiling(remaining / ds, kind=I4P))
      ni_range = [1_I4P, ni]
      nj_range = [nj - cells + 1_I4P, nj]
      nk_range = [1_I4P, nk]
   case (5_I4P)
      ds = dxyz(3)
      distance0 = block_emin(3) - domain_emin(3)
      remaining = width - distance0
      if (remaining <= 0._R8P) return
      cells = min(nk, ceiling(remaining / ds, kind=I4P))
      ni_range = [1_I4P, ni]
      nj_range = [1_I4P, nj]
      nk_range = [1_I4P, cells]
   case (6_I4P)
      ds = dxyz(3)
      distance0 = domain_emax(3) - block_emax(3)
      remaining = width - distance0
      if (remaining <= 0._R8P) return
      cells = min(nk, ceiling(remaining / ds, kind=I4P))
      ni_range = [1_I4P, ni]
      nj_range = [1_I4P, nj]
      nk_range = [nk - cells + 1_I4P, nk]
   endselect

   if (cells <= 0_I4P) return

   last_center_distance = distance0 + real(cells - 1_I4P, R8P) * ds
   profile_extent = distance0 + real(cells, R8P) * ds
   endsubroutine compute_absorbing_face_range

endmodule adam_prism_absorbing_layer_geometry
