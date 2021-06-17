module memorysaver
   use iso_c_binding
   implicit none
   interface
      subroutine getmemory(memtot, memavail) bind(C, name="getmemory")
         import :: C_LONG
         integer(C_LONG), intent(in) :: memtot, memavail
      endsubroutine getmemory
   endinterface

   contains

   subroutine save_memory(it, rank)
   integer,           intent(in) :: it       !< Temporal iteration.
   integer, optional, intent(in) :: rank     !< Rank or integer to label filename.
   character(999)                :: fmem     !< File for memory saves.
   integer(C_LONG)               :: memtot   !< Total RAM memory of process.
   integer(C_LONG)               :: memavail !< Avail RAM memory of process.
   integer                       :: rank_    !< Rank or integer to label filename, local var.

   rank_ = 0 ; if(present(rank)) rank_ = rank
   call getmemory(memtot, memavail)
   fmem = "memory_usage_?????.dat"
   write(fmem(14:18), "(I5.5)") rank_
   open(unit=125, file=fmem, position="append")
   write(125,*) it, memavail, memtot
   close(125)
   endsubroutine save_memory
endmodule memorysaver
