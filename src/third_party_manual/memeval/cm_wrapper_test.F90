program getmemory_test
use getmemory_mod
integer(C_LONG) :: memtot, memavail
call getmemory(memtot, memavail)
print*,'memavail/memtot: ',memavail,'/',memtot
endprogram getmemory_test
