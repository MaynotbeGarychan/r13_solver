!
!     dynamic memory allocation stuff.  This has to be a C style 
!     include because of the use of MINDEX
!
#ifdef NULL_IMEM
      integer i_mem(0:*)
      real r_mem(0:*)
      integer*4 i4_mem(0:*)
      integer*8 i8_mem(0:*)
      real*4 r4_mem(0:*)
      real*8 r8_mem(0:*)
      real*8 real8_mem(0:*)
      MINDEX d_mem,m_mem(0:*),mem8_blk
      common/dynmem/d_mem,mem8_blk(15)
      pointer (i_wew,i_mem)
      pointer (r_wew,r_mem)
      pointer (i8_wew,i8_mem)
      pointer (r8_wew,r8_mem)
      pointer (real8_wew,real8_mem)
      pointer (i4_wew,i4_mem)
      pointer (r4_wew,r4_mem)
      pointer (m_wew,m_mem)
      equivalence (i_wew,r_wew,i8_wew,r8_wew,real8_wew,i4_wew,r4_wew,   &
     &  m_wew,d_mem)
#else
      integer i_mem(0:1)
      real r_mem(0:1)
      integer*4 i4_mem(0:1)
      integer*8 i8_mem(0:1)
      real*4 r4_mem(0:1)
      real*8 r8_mem(0:1)
      real*8 real8_mem
      MINDEX m_mem(0:1)
      MINDEX mem8_blk(0:15)
      common/dynmem/real8_mem(0:15)
      equivalence (i_mem,r_mem,i8_mem,r8_mem,real8_mem,i4_mem,r4_mem,   &
     &  m_mem,mem8_blk)
#endif
      MINDEX mem_alloc, mems_alloc,mem_realloc, memh_ptr
      integer mem_length,memh_alloc,memh_realloc,memh_length
      integer memh_enum,memh_enums,memh_attach
!      integer mem_newmark,mem_getmark,memh_getmark
      MINDEX memh_detach,memh_ralloc
      MINDEX m_to_m8
      MINDEX m_to_mm
      MINDEX m8_to_m
      MINDEX m_to_m4
      MINDEX mptr_to_m
#ifndef PCWIN
      external mem_alloc, mems_alloc, mem_realloc, mem_length
      external memh_alloc, memh_realloc, memh_length, memh_ptr
      external memh_enum,memh_attach,memh_detach,memh_ralloc
      external memh_enums
      external m_to_m8
      external m_to_mm
      external m8_to_m
      external m_to_m4
      external mptr_to_m
#endif
      MINDEX dm_xtpz
      MINDEX dm_n45
      integer len_n45,len_n46
      common/dynmem1/                                                   &
     & dm_xtpz,                                                         &
     & dm_n45          
      common/dynmem1l/len_n45,len_n46
!
      MINDEX mem_allocm,mems_allocm,mem_reallocm,mem_loc
      integer memh_allocm,memh_reallocm
#ifndef PCWIN
      external mem_allocm,mems_allocm,mem_reallocm,mem_loc
      external memh_allocm,memh_reallocm
#endif
      MINDEX mem_allocmchk,mem_allocchk,mems_allocmchk
      MINDEX mems_allocchk,mem_reallocmchk,mem_reallocchk
      integer memh_allocmchk,memh_allocchk
#ifndef PCWIN
      external mem_allocmchk,mem_allocchk,mems_allocmchk
      external mems_allocchk,mem_reallocmchk,mem_reallocchk
      external memh_allocmchk,memh_allocchk
#endif
