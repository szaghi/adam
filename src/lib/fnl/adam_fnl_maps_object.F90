!< ADAM, maps class definition, FNL backend.
module adam_fnl_maps_object
!< ADAM, maps class definition, FNL backend.

! ADAM objects
use :: adam_maps_object,              only : maps_object
use :: adam_seam_interpolation_library, only : SEAM_FILL_INJECTION, SEAM_RESTRICT_MEAN
! ADAM FNL singleton objects
use :: adam_fnl_mpih_global, only : mpih_fnl
! third party modules
use :: fundal
use :: penf
! sdk modules
use :: mpi

implicit none
private
public :: maps_fnl_object

type :: maps_fnl_object
   !< Maps class, FNL backend.
   ! GPU data.
   ! The `=> null()` initializers are mandatory, not cosmetic: dev_assign_to_device
   ! does `if (associated(dst)) call dev_free(dst)` before allocating, and
   ! `associated()` on a pointer that was never nullified is undefined behaviour
   ! (nvfortran traps it as "Null pointer for tmp$r"). A default-initialised
   ! pointer component starts every instance in a defined, disassociated state.
   integer(I8P), pointer :: local_map_ghost_cell_gpu(:,:)     => null() !< Local map for ghost cells updating, cells order.
   integer(I4P), pointer :: local_map_ghost_cell_prol_gpu(:)  => null() !< The flag-4 (coarse->fine) rows of
                                                                        !< local_map_ghost_cell_gpu (#157 B2; none: null).
   real(R8P),    pointer :: seam_rho_gpu(:)                   => null() !< Clamped density per row of
                                                                        !< local_map_ghost_cell_gpu (its flag-4 rows): the
                                                                        !< seam fill guard's part once per ghost (#157 B2).
   integer(I8P), pointer :: comm_map_recv_ghost_cell_gpu(:,:) => null() !< Communication map, `fec` information, cell order.
   integer(I8P), pointer :: comm_map_send_ghost_cell_gpu(:,:) => null() !< Communication map, `fec` information, cell order.
   integer(I4P), pointer :: comm_map_send_density_rows_gpu(:) => null() !< The density rows of the flag-4 cells of
                                                                        !< comm_map_send_ghost_cell_gpu (#157 B3; none: null).
   real(R8P),    pointer :: seam_rho_send_gpu(:)              => null() !< Clamped density per row of
                                                                        !< comm_map_send_ghost_cell_gpu (at those rows).
   real(R8P),    pointer :: send_buffer_ghost_gpu(:)          => null() !< Send buffer of ghost cells.
   real(R8P),    pointer :: recv_buffer_ghost_gpu(:)          => null() !< Receive buffer of ghost cells.
   integer(I8P), pointer :: local_map_bc_crown_gpu(:,:,:)     => null() !< Local map for face BC ghost cells, "crown" order.
   ! Inter-realm seam ghost-fill — device-resident counterparts of the host
   ! seam maps/buffers on `maps_object` (see adam_maps_object.F90). The host
   ! arrays drive layout; these are sized/populated 1:1 by `copy_cpu_gpu`.
   ! Per-peer buffers shaped `(buf_len, n_peers)` so concurrent peers do not
   ! alias.
   integer(I4P), pointer :: seam_local_map_ghost_cell_gpu(:,:)  => null() !< Per-cell seam ghost map (sorted by peer_realm).
   real(R8P),    pointer :: seam_local_send_buf_gpu(:,:)        => null() !< Per-peer pack buffer, device-resident.
   real(R8P),    pointer :: seam_local_recv_buf_gpu(:,:)        => null() !< Per-peer unpack buffer, device-resident.
   ! Cross-rank seam rows (issue #40), device twins of `maps_object%seam_mpi_{send,recv}_cell`: the pack and unpack
   ! kernels gather and scatter the cross-rank seam cells through them.
   integer(I4P), pointer :: seam_mpi_send_cell_gpu(:,:) => null() !< Cross-rank send rows [rank, b, i, j, k].
   integer(I4P), pointer :: seam_mpi_recv_cell_gpu(:,:) => null() !< Cross-rank receive rows [rank, b, i, j, k].
   ! Intra-realm coarse->fine seam ghost fill (issue #21 N2 / #22 F3): host-side
   ! mirror of `maps_object%seam_ghost_fill`, the active fill regime consumed by
   ! the flag-4 branches of the ghost kernels (passed as a scalar kernel
   ! argument, not device-resident). Refreshed on every `copy_cpu_gpu`.
   integer(I4P)          :: seam_ghost_fill = SEAM_FILL_INJECTION             !< Active seam ghost-fill regime (host mirror).
   integer(I4P)          :: seam_density_var = 0_I4P                          !< Seam fill density guard (host mirror of
                                                                             !< `maps_object%seam_density_var`, #101).
   integer(I4P)          :: seam_restriction = SEAM_RESTRICT_MEAN             !< Fine->coarse seam restriction (host mirror of
                                                                              !< `maps_object%seam_restriction`, #68 C2).
   integer(I4P)          :: refine_ratio(3) = 2_I4P                           !< Per-axis 2:1 refinement factor (host mirror,
                                                                              !< issue #46): 1 on axes the tree never refines.
   contains
      procedure, pass(self) :: copy_cpu_gpu !< Copy data from (maps global singleton) CPU to (maps_fnl_object) GPU.
      procedure, pass(self) :: destroy      !< Free device data owned by the helper.
      procedure, pass(self) :: initialize   !< Initialize MPI handler data.
endtype maps_fnl_object

contains
   subroutine copy_cpu_gpu(self, maps, verbose)
   !< Copy data from the (realm-local) CPU `maps` to (maps_fnl_object) GPU.
   class(maps_fnl_object), intent(inout)        :: self     !< The maps.
   type(maps_object),      intent(in)           :: maps     !< Realm-local CPU maps; was the `maps` singleton.
   logical,                intent(in), optional :: verbose  !< Flag to activate verbose mode.
   logical                                      :: verbose_ !< Flag to activate verbose mode, local var.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call mpih_fnl%print_message('maps_fnl_object%copy_cpu_gpu start')
   ! Each source map is an `allocatable` component of the (realm-local) CPU
   ! `maps`, populated conditionally by maps_object's make_* routines.
   ! On a single MPI rank the inter-rank communication maps and buffers
   ! have nothing to build and are legitimately left unallocated. Passing
   ! an unallocated allocatable as the non-optional assumed-shape
   ! `intent(in)` `src` of dev_assign_to_device is illegal Fortran. Guard
   ! every copy with allocated(): an empty CPU map leaves the GPU pointer
   ! null, the correct device-side state for "nothing to exchange". A copy
   ! after a runtime regrid (issue #74) frees and nullifies the device map of
   ! a CPU map that is no longer allocated: the kernels size their loops from
   ! the device arrays, and a stale map would replay the old grid's rows.
   ! Intra-realm seam ghost fill (issue #22 F3): flag-4 map rows (coarse->fine
   ! seam interpolation, metadata column included in the verbatim map copy) are
   ! handled on device by the flag-4 branches of adam_fnl_field_kernels; the
   ! active regime is mirrored here for the kernel callers.
   self%seam_ghost_fill = maps%seam_ghost_fill
   self%refine_ratio    = maps%refine_ratio
   self%seam_restriction = maps%seam_restriction
   self%seam_density_var = maps%seam_density_var
   if (associated(self%local_map_ghost_cell_prol_gpu)) then
      call dev_free(self%local_map_ghost_cell_prol_gpu, mydev)
      nullify(self%local_map_ghost_cell_prol_gpu)
   endif
   if (associated(self%seam_rho_gpu)) then
      call dev_free(self%seam_rho_gpu, mydev)
      nullify(self%seam_rho_gpu)
   endif
   if (allocated(maps%local_map_ghost_cell)) then
      call dev_assign_to_device(dst=self%local_map_ghost_cell_gpu, src=maps%local_map_ghost_cell)
      ! the seam fill guard once per fine ghost (#157 B2): the flag-4 rows and a scratch value per row
      call make_prolongation_rows(map=maps%local_map_ghost_cell, cflag=9, prol_gpu=self%local_map_ghost_cell_prol_gpu, &
                                  rho_gpu=self%seam_rho_gpu)
      if (verbose_) call mpih_fnl%print_message('copy local_map_ghost_cell_gpu done ('// &
         trim(str(count(maps%local_map_ghost_cell(:,9) == 4_I8P)))//' seam flag-4 rows)')
   else
      if (associated(self%local_map_ghost_cell_gpu)) then
         call dev_free(self%local_map_ghost_cell_gpu, mydev)
         nullify(self%local_map_ghost_cell_gpu)
      endif
      if (verbose_) call mpih_fnl%print_message('skip local_map_ghost_cell_gpu (CPU map not allocated)')
   endif
   if (associated(self%comm_map_send_density_rows_gpu)) then
      call dev_free(self%comm_map_send_density_rows_gpu, mydev)
      nullify(self%comm_map_send_density_rows_gpu)
   endif
   if (associated(self%seam_rho_send_gpu)) then
      call dev_free(self%seam_rho_send_gpu, mydev)
      nullify(self%seam_rho_send_gpu)
   endif
   if (allocated(maps%comm_map_send_ghost_cell)) then
      call dev_assign_to_device(dst=self%comm_map_send_ghost_cell_gpu, src=maps%comm_map_send_ghost_cell)
      ! the seam fill guard once per fine ghost (#157 B3): the density row of each flag-4 cell (its variables are
      ! contiguous rows, in order) and a scratch value per row
      if (maps%seam_density_var > 0_I4P) &
         call make_prolongation_rows(map=maps%comm_map_send_ghost_cell, cflag=7, prol_gpu=self%comm_map_send_density_rows_gpu, &
                                     rho_gpu=self%seam_rho_send_gpu, cvar=5, var=maps%seam_density_var)
      if (verbose_) call mpih_fnl%print_message('copy comm_map_send_ghost_cell done ('// &
         trim(str(count(maps%comm_map_send_ghost_cell(:,7) == 4_I8P)))//' seam flag-4 rows)')
   else
      if (associated(self%comm_map_send_ghost_cell_gpu)) then
         call dev_free(self%comm_map_send_ghost_cell_gpu, mydev)
         nullify(self%comm_map_send_ghost_cell_gpu)
      endif
      if (verbose_) call mpih_fnl%print_message('skip comm_map_send_ghost_cell_gpu (CPU map not allocated)')
   endif
   if (allocated(maps%comm_map_recv_ghost_cell)) then
      call dev_assign_to_device(dst=self%comm_map_recv_ghost_cell_gpu, src=maps%comm_map_recv_ghost_cell)
      if (verbose_) call mpih_fnl%print_message('copy comm_map_recv_ghost_cell done')
   else
      if (associated(self%comm_map_recv_ghost_cell_gpu)) then
         call dev_free(self%comm_map_recv_ghost_cell_gpu, mydev)
         nullify(self%comm_map_recv_ghost_cell_gpu)
      endif
      if (verbose_) call mpih_fnl%print_message('skip comm_map_recv_ghost_cell_gpu (CPU map not allocated)')
   endif
   if (allocated(maps%send_buffer_ghost)) then
      call dev_assign_to_device(dst=self%send_buffer_ghost_gpu, src=maps%send_buffer_ghost)
      if (verbose_) call mpih_fnl%print_message('copy send_buffer_ghost done')
   else
      if (associated(self%send_buffer_ghost_gpu)) then
         call dev_free(self%send_buffer_ghost_gpu, mydev)
         nullify(self%send_buffer_ghost_gpu)
      endif
      if (verbose_) call mpih_fnl%print_message('skip send_buffer_ghost_gpu (CPU map not allocated)')
   endif
   if (allocated(maps%recv_buffer_ghost)) then
      call dev_assign_to_device(dst=self%recv_buffer_ghost_gpu, src=maps%recv_buffer_ghost)
      if (verbose_) call mpih_fnl%print_message('copy recv_buffer_ghost done')
   else
      if (associated(self%recv_buffer_ghost_gpu)) then
         call dev_free(self%recv_buffer_ghost_gpu, mydev)
         nullify(self%recv_buffer_ghost_gpu)
      endif
      if (verbose_) call mpih_fnl%print_message('skip recv_buffer_ghost_gpu (CPU map not allocated)')
   endif
   if (allocated(maps%local_map_bc_crown)) then
      call dev_assign_to_device(dst=self%local_map_bc_crown_gpu, src=maps%local_map_bc_crown)
      if (verbose_) call mpih_fnl%print_message('copy local_map_bc_crown done')
   else
      if (associated(self%local_map_bc_crown_gpu)) then
         call dev_free(self%local_map_bc_crown_gpu, mydev)
         nullify(self%local_map_bc_crown_gpu)
      endif
      if (verbose_) call mpih_fnl%print_message('skip local_map_bc_crown_gpu (CPU map not allocated)')
   endif
   ! Inter-realm seam ghost-fill — device-resident counterparts.
   if (allocated(maps%seam_local_map_ghost_cell)) then
      call dev_assign_to_device(dst=self%seam_local_map_ghost_cell_gpu, src=maps%seam_local_map_ghost_cell)
      if (verbose_) call mpih_fnl%print_message('copy seam_local_map_ghost_cell done')
   else
      if (associated(self%seam_local_map_ghost_cell_gpu)) then
         call dev_free(self%seam_local_map_ghost_cell_gpu, mydev)
         nullify(self%seam_local_map_ghost_cell_gpu)
      endif
      if (verbose_) call mpih_fnl%print_message('skip seam_local_map_ghost_cell_gpu (CPU map not allocated)')
   endif
   if (allocated(maps%seam_local_send_buf)) then
      call dev_assign_to_device(dst=self%seam_local_send_buf_gpu, src=maps%seam_local_send_buf)
      if (verbose_) call mpih_fnl%print_message('copy seam_local_send_buf done')
   else
      if (associated(self%seam_local_send_buf_gpu)) then
         call dev_free(self%seam_local_send_buf_gpu, mydev)
         nullify(self%seam_local_send_buf_gpu)
      endif
      if (verbose_) call mpih_fnl%print_message('skip seam_local_send_buf_gpu (CPU map not allocated)')
   endif
   if (allocated(maps%seam_local_recv_buf)) then
      call dev_assign_to_device(dst=self%seam_local_recv_buf_gpu, src=maps%seam_local_recv_buf)
      if (verbose_) call mpih_fnl%print_message('copy seam_local_recv_buf done')
   else
      if (associated(self%seam_local_recv_buf_gpu)) then
         call dev_free(self%seam_local_recv_buf_gpu, mydev)
         nullify(self%seam_local_recv_buf_gpu)
      endif
      if (verbose_) call mpih_fnl%print_message('skip seam_local_recv_buf_gpu (CPU map not allocated)')
   endif
   ! Cross-rank seam rows (issue #40): allocated only on ranks that have cross-rank seam rows.
   if (allocated(maps%seam_mpi_send_cell)) then
      call dev_assign_to_device(dst=self%seam_mpi_send_cell_gpu, src=maps%seam_mpi_send_cell)
      if (verbose_) call mpih_fnl%print_message('copy seam_mpi_send_cell done')
   elseif (associated(self%seam_mpi_send_cell_gpu)) then
      call dev_free(self%seam_mpi_send_cell_gpu, mydev)
      nullify(self%seam_mpi_send_cell_gpu)
   endif
   if (allocated(maps%seam_mpi_recv_cell)) then
      call dev_assign_to_device(dst=self%seam_mpi_recv_cell_gpu, src=maps%seam_mpi_recv_cell)
      if (verbose_) call mpih_fnl%print_message('copy seam_mpi_recv_cell done')
   elseif (associated(self%seam_mpi_recv_cell_gpu)) then
      call dev_free(self%seam_mpi_recv_cell_gpu, mydev)
      nullify(self%seam_mpi_recv_cell_gpu)
   endif
   if (verbose_) call mpih_fnl%print_message('maps_fnl_object%copy_cpu_gpu finish')
   endsubroutine copy_cpu_gpu

   subroutine destroy(self)
   !< Free device data owned by the maps FNL helper.
   class(maps_fnl_object), intent(inout) :: self !< The FNL helper.

   if (associated(self%local_map_ghost_cell_gpu)) then
      call dev_free(self%local_map_ghost_cell_gpu, mydev)
      nullify(self%local_map_ghost_cell_gpu)
   endif
   if (associated(self%seam_rho_gpu)) then
      call dev_free(self%seam_rho_gpu, mydev)
      nullify(self%seam_rho_gpu)
   endif
   if (associated(self%local_map_ghost_cell_prol_gpu)) then
      call dev_free(self%local_map_ghost_cell_prol_gpu, mydev)
      nullify(self%local_map_ghost_cell_prol_gpu)
   endif
   if (associated(self%comm_map_recv_ghost_cell_gpu)) then
      call dev_free(self%comm_map_recv_ghost_cell_gpu, mydev)
      nullify(self%comm_map_recv_ghost_cell_gpu)
   endif
   if (associated(self%comm_map_send_ghost_cell_gpu)) then
      call dev_free(self%comm_map_send_ghost_cell_gpu, mydev)
      nullify(self%comm_map_send_ghost_cell_gpu)
   endif
   if (associated(self%comm_map_send_density_rows_gpu)) then
      call dev_free(self%comm_map_send_density_rows_gpu, mydev)
      nullify(self%comm_map_send_density_rows_gpu)
   endif
   if (associated(self%seam_rho_send_gpu)) then
      call dev_free(self%seam_rho_send_gpu, mydev)
      nullify(self%seam_rho_send_gpu)
   endif
   if (associated(self%send_buffer_ghost_gpu)) then
      call dev_free(self%send_buffer_ghost_gpu, mydev)
      nullify(self%send_buffer_ghost_gpu)
   endif
   if (associated(self%recv_buffer_ghost_gpu)) then
      call dev_free(self%recv_buffer_ghost_gpu, mydev)
      nullify(self%recv_buffer_ghost_gpu)
   endif
   if (associated(self%local_map_bc_crown_gpu)) then
      call dev_free(self%local_map_bc_crown_gpu, mydev)
      nullify(self%local_map_bc_crown_gpu)
   endif
   if (associated(self%seam_local_map_ghost_cell_gpu)) then
      call dev_free(self%seam_local_map_ghost_cell_gpu, mydev)
      nullify(self%seam_local_map_ghost_cell_gpu)
   endif
   if (associated(self%seam_local_send_buf_gpu)) then
      call dev_free(self%seam_local_send_buf_gpu, mydev)
      nullify(self%seam_local_send_buf_gpu)
   endif
   if (associated(self%seam_local_recv_buf_gpu)) then
      call dev_free(self%seam_local_recv_buf_gpu, mydev)
      nullify(self%seam_local_recv_buf_gpu)
   endif
   if (associated(self%seam_mpi_send_cell_gpu)) then
      call dev_free(self%seam_mpi_send_cell_gpu, mydev)
      nullify(self%seam_mpi_send_cell_gpu)
   endif
   if (associated(self%seam_mpi_recv_cell_gpu)) then
      call dev_free(self%seam_mpi_recv_cell_gpu, mydev)
      nullify(self%seam_mpi_recv_cell_gpu)
   endif
   endsubroutine destroy

   subroutine make_prolongation_rows(map, cflag, prol_gpu, rho_gpu, cvar, var)
   !< The device list of the flag-4 (coarse->fine) rows of a ghost-cell map and a scratch value per row, for the seam
   !< fill guard computed once per fine ghost (#157 B2, B3); with `cvar` and `var` (a map with a row per variable), the
   !< flag-4 rows of variable `var` only. A map with no such row leaves both null. The caller frees any previous
   !< allocation.
   integer(I8P),          intent(in)           :: map(1:,1:)  !< Ghost-cell map (host).
   integer(I4P),          intent(in)           :: cflag       !< Column of the kind flag.
   integer(I4P), pointer, intent(inout)        :: prol_gpu(:) !< The flag-4 rows (device).
   real(R8P),    pointer, intent(inout)        :: rho_gpu(:)  !< Scratch, one value per row (device).
   integer(I4P),          intent(in), optional :: cvar        !< Column of the variable.
   integer(I4P),          intent(in), optional :: var         !< Variable.
   integer(I4P), allocatable                   :: rows(:)     !< The flag-4 rows (host).
   integer(I4P)                                :: r           !< Counter.

   if (present(cvar) .and. present(var)) then
      rows = pack([(r, r=1, size(map, dim=1))], map(:,cflag) == 4_I8P .and. map(:,cvar) == int(var, I8P))
   else
      rows = pack([(r, r=1, size(map, dim=1))], map(:,cflag) == 4_I8P)
   endif
   if (size(rows, dim=1) == 0) return
   call dev_assign_to_device(dst=prol_gpu, src=rows)
   call dev_assign_to_device(dst=rho_gpu, src=spread(0._R8P, dim=1, ncopies=size(map, dim=1)))
   endsubroutine make_prolongation_rows

   subroutine initialize(self, maps)
   !< Initialize maps from the (realm-local) CPU `maps`.
   !< Requires `mpih_fnl` (adam_fnl_mpih_global) initialized and `maps` populated before calling.
   class(maps_fnl_object), intent(inout) :: self !< The maps, FNL backend.
   type(maps_object),      intent(in)    :: maps !< Realm-local CPU maps; was the `maps` singleton.

   call mpih_fnl%print_message('maps_fnl_object%initialize start')
   call self%copy_cpu_gpu(maps=maps, verbose=.true.)
   call mpih_fnl%print_message('maps_fnl_object%initialize finish')
   endsubroutine initialize
endmodule adam_fnl_maps_object
