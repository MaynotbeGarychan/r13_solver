      integer igmprt,ngmprt,massout,iprtms,iprtlc,iprtop,lcoutp,iglstat
      integer sw2_grpc
      real prtims,prtlst,prtnxt,gmfilts
      integer iascii_file
      parameter (iascii_file=100)
      common/ptimes/igmprt,ngmprt,massout(10),                          &
     & prtims(iascii_file),prtlst(iascii_file),                         &
     & iprtms(iascii_file),iprtlc(iascii_file),                         &
     & prtnxt(iascii_file),iprtop(iascii_file),                         &
     & lcoutp(45),gmfilts,iglstat,sw2_grpc
!
!     ascii output control data
!
!     prtims-fixed ascii output intervals
!      ascii files:
!           ( 1)-SECFORC cross section forces
!           ( 2)-RWFORC rigid wall forces
!           ( 3)-NODOUT nodal data
!           ( 4)-ELOUT  element data
!           ( 5)-GLSTAT global data
!           ( 6)-DEFORC discrete elements
!           ( 7)-MATSUM material energies
!           ( 8)-NCFORC nodal interface forces
!           ( 9)-RCFORC resultant interface forces
!           (10)-DEFGEO smug animator
!           (11)-SPCFORC spc reaction forces
!           (12)-SWFORC nodal constraint resultant forces
!           (13)-ABSTAT airbag statistics
!           (14)-AVSFLT avs database
!           (15)-NODFOR nodal force groups
!           (16)-BNDOUT output intervals for nodal boundary conditions
!           (17)-RBDOUT
!           (18)-GCEOUT
!           (19)-MPGS
!           (20)-MOVIE
!           (21)-SLEOUT
!           (22)-SBTOUT
!           (23)-JNTFORC
!           (24)-TRHIST
!           (25)-DCFAIL
!           (26)-SPHOUT
!           (27)-H3OUT
!           (28)-mcol
!           (29)-mcolmod
!           (30)-mcolkine
!           (31)-NODOUT
!           (32)-TPRINT 
!           (33)-
!           (34)-SSSTAT
!           (35)-MASSOUT
!           (36)-ELOUTDET
!           (37)-ATDOUT
!           (38)-BEPRES
!           (39)-BKSTAT
!           (40)-CURVOUT
!           (41)-SENSOR_CPM
!           (42)-ELOUT_SSD binout only
!           (43)-NODOUT_SSD binout only
!           (44)-ELOUT_SPCM
!           (45)-NODOUT_SPCM
!           (46)-ELOUT_PSD
!           (47)-NODOUT_PSD
!           (48)-DISBOUT
!           (49)-PLLYOUT
!           (50)-PBSTAT
!           (51)-BEARINGOUT
!           (52)-DEMASSFLOW
!           (53)-SPHMASSFLOW
!           (54)-SPHVICINITY
!           (63)-BEARING
!           (65)-PGSTAT
!           (66)-PGSNSR
!           (67)-PRTUBE
!           (68)-NODFOR_SSD
!           (69)-DESTAT global des data
!           (70)-PYRO
!           (75)-PBLAST sensor binout only
!           (78)-S-ALE sensor
!           (79)-S-ALE mmg switch
!           (80)-S-ALE tracer
!           (81)-S-ALE fsi database
!           (83)-D3MAX
!           (85)-DBFSI
!           (86)-DBSENSOR
!     prtlst-current ascii output time
!     prtnxt-next ascii output time
!     iprtms-output type  for ascii output
!     iprtlc-curve id foer ascii output interval
!     iprtop-option for ascii output time
!     lcoutp-misc. numbers
!     igmprt-flag for any ascii output 0-no 1-yes
!     ngmprt-number of ascii output files
!     iglstat = 0/1 flag, 1 if we are preparing to output glstat
!     sw2_grpc = 0/1 flag, 1 sw2 output was requested by grpc client
!
