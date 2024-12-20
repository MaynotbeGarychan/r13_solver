!
!***  EASY dynamic memory management
!
!***  C STYLE INCLUDE because of the use of MINDEX
!
!***  ALLOCATE MEMORY:
!         call ez_mem_allocm(i,len,literal_string) 
!         call ez_mem_allocmchk(i,len,literal_string)
!
!       ez_mem_allocm puts a value of 0 in iezp(i) and mezp(i) to indicate
!                  that the memory allocation failed just as mem_allocm
!                  returns a 0 for the memory pointer
!
!       ez_mem_allocmchk terminates the code on failure.
!
!       where 
!         i is the equivalent of the integer used with addlen
!         len is the number of words of memory
!         literal_string is a 'description'
!
!         iezp(i) is the handle you are using
!         mezp(i) is the corresponding MINDEX pointer to i_mem and r_mem 
!
!       therefore iezp points to mezp through memh_ptr, i.e.,
!         mezp(i)=memh_ptr(iezp(i))
!
!***  REALLOCATE MEMORY:
!         call ez_mem_realloc(i,len,literal_string)
!
!***  FREE MEMORY:
!         automatic at end of run
!       or
!         call ez_mem_free(i)
!
!***  RESTARTS:
!         iezp is written to the restart file and
!         mezp is automatically updated during the restart
!
!         memory allocated using this method will automatically 
!         be handled during restart writes and reads
!
!***  IEZP AND MEZP DOCUMENTATION:
!         the equivalent of the ioshl comments in dynai
!         are located in routine ez_mem_allocmchk in atemp.F
!
!***  INITIALIZATION:
!         all arrays initialized to zero in ez_mem_init
!
      integer LENEZMEM  
      parameter (LENEZMEM=1000)
      integer iezp
      MINDEX  mezp
      common/iezpcom/iezp(LENEZMEM)
      common/mezpcom/mezp(LENEZMEM)
