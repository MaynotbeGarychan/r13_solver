      integer iptadj,ipta2,ipta3,ipta4,ipta5,ipta6,ipta7,ipta8,         &
     & ipta9,ipta10,ipta11,ipta12,ipta13
      real optsz,optszsv
      common/adjcom/iptadj,ipta2,ipta3,ipta4,ipta5,ipta6,ipta7,ipta8,   &
     & ipta9,ipta10,ipta11,ipta12,ipta13,optsz,optszsv
      integer conparmdim
      parameter (conparmdim =32)
!
! ipta13 points to a chunk of memory that is dimensioned
! (comparmdim,numsv), to hold whatever contact related
! parameters need to be added as capabilities are added.
!
!
!  Current storage allocations:
!
!  1- 6 Used for the Contraction Joint capability.
!  7- 7 Flag to activate local coordinate system for slave box
!  8-16 Coordinate transformation matrix for slave box
! 17-19 Transformed origin for slave box
! 20-20 Flag to activate local coordinate system for master box
! 21-29 Coordinate transformation matrix for master box
! 30-32 Transformed origin for master box
!
