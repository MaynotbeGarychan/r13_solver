# 1 "init_dyn21.F"
# 1 ".\define.inc" 1 


# 26
















# 44


# 48


# 53


# 57


# 62









# 73









# 85



# 90






























# 122

# 126




# 136



# 180


# 188


# 223


# 257


# 283






# 299






# 311

# 318









# 329




# 341










# 370


# 1 ".\adios.h" 1 
# 9



















# 373 ".\define.inc" 2 


# 379

# 382






# 408



# 417

# 424





# 2 "init_dyn21.F" 2 

      subroutine usercomm
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     to define threadprivate common (aux14loc) for the history variable
c     and any new threadprivate common blocks by the user
c
c     I/O unit 241~245 are reversed for user material usage
c
c     the maximum number of history variables (NHISVAR) is defined in
c     the include file (nhisparm.inc)
c     the user can increase the number of history variables by
c     editing NHISVAR in the include (nhisparm.inc)
c
c Note:
c  For SMP or MPP version, all threadprivate common blocks must be loaded
c  and initialized in the beginning of job
c
c  If the user adds any new threadprivate common blocks, the following
c  three lines must be added:
c   1. common/xxxloc/xxxzz(nlq) # define common block
c   2. c$omp threadprivate (/xxxloc/)  # declare it as threadprivate common
c   3. xxxzz(1)=0.0              # initailze the first variable
c
c   4. common /auxstr/ atores elemnent history variables in initialization routine.
c      The size of the array should be (7 + number of history variables)*nip.
c      User should adjust the size to hold element data.
c
c For example, to add a new threadprivate common/userloc/
c
c 123456789012345678901234567890123456789012345678901234567890
c      common/user1loc/user1zz1(nlq),user1zz2(nlq)
c c$omp threadprivate (/user1loc/)
c      user1zz1(1)=0.0
c
c---------------------------------------------------------------
# 48

c
      include 'nlqparm'
      include 'nhisparm.inc'
c
      parameter (NBEL8A=max(100,2*(NHISVUE+45)))
      parameter (NBEL7 =max(2601,
     . 1180+576*NXDOFUE+64*NXDOFUE*NXDOFUE+12*5))
# 58

      parameter (NBUX13=max(62,49+NHISVAR))

      common/auxstr/ ajunk(5000)
      common/aux14loc/ax14zz1(nlq,7+NHISVAR)
      common/bel8loc/bel8lz1(nlq,2*(NHISVUE+45))
      common/bel8aloc/bel8lz2(nlq,NBEL8A)
      common/bel7loc/bel7lz(nlq,NBEL7)
      common/bux13loc/bx13zz(nlq,NBUX13)

c
c$omp threadprivate (/bux13loc/)
c$omp threadprivate (/bel7loc/)
c$omp threadprivate (/bel8aloc/)
c$omp threadprivate (/bel8loc/)
c$omp threadprivate (/aux14loc/)
      ax14zz1(1,1)=0.0
      bel8lz1(1,1)=0.0
      bel8lz2(1,1)=0.0
      bel7lz(1,1)=0.0
      bx13zz(1,1)=0.0
      ajunk(1)   =0.0
c
      return
      end
      subroutine usercomm1
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
      parameter (NBEL8A=max(100,2*(NHISVUE+45)))
      parameter (NBEL7 =max(2601,
     . 1180+576*NXDOFUE+64*NXDOFUE*NXDOFUE+12*5))
      common/usercom/iusercomm(4)
      parameter (nconst_max=max(1474,(NHISVAR+26)*8+41+2*24))
# 103

      common/strhx/ideh,mxh,ixhx(8),xe(3,8),ve(3,8),strnh(6),
     . ressh,hisvar(nconst_max)

# 109

      ideh=0
c
c     define the length of 4 common blocks
      iusercomm(1)=nlq*(7+NHISVAR)
      iusercomm(2)=nlq*2*(NHISVUE+45)
      iusercomm(3)=nlq*NBEL8A
      iusercomm(4)=nlq*NBEL7
c
      return
      end
