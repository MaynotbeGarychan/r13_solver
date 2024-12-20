      integer                                                           &
     & iprelax,nxale,nyale,nzale,nalebn,nalegp,                         &
     & iale,ialetp,ialenc,ialead,newall,naltie,ialefst,jalefst,         &
     & jale00,jale01,jale02,jale03,jale04,jale05,jale06,jale07,         &
     & jale08,jale09,jale10,jale11,jale12,jale13,jale14,jale15,         &
     & jale16,jale17,jale18,jale19,jale20,jale21,jale22,jale23,         &
     & jale24,jale25,jale26,jale27,jale28,jale29,jale30,jale31,         &
     & jale32,jale33,jale34,jale35,jale36,jale37,jale38,jale39,         &
     & lenale,ntrace,nalesc,nintrf,jntrf,jalenp,nalenp,nvoid,           &
     & idotie,ltie1,ltie1b,ltie2,ltie3,ltie4,ltie5,ltie6,               &
     & ltie7,ltie8,ltie9,ltie10,ltie11,ltie12,ltie12e,ltie13,ltie13b,   &
     & ltie13c,ltie13d,ltie14,ltie15,ltie15b,ltie15c,ltie15d,           &
     & ltie15e,ltie15f,ltie15g,ltie16,                                  &
     & ltie17a,ltie17b,ltie18a,ltie18b,ltie19,ltie20,ltie21,ltie22,     &
     & ltie23,ltie24,ltie25,ltie26,ltie27,ltie28,ltie29,ltie30,         &
     & ltie31,ltie32,ltie33,ltie34,ltie35,ltie36,mxtquad,               &
     & nalebcsg,ialrst,nlagran,ibrfrm,iotflw,nofp,                      &
     & lvoid1,lvoid2,lvoid3,lvoid4,iframv,jnk2,jnk3,jnk4,               &
     & jnk5,jnk6,jnk7,jnk8,jnk9,jnk10,jnk11,jnk12,                      &
     & jnk13,jnk14,jnk15,iale_act,iale_mt41skpinc,jale,jalesies,        &
     & numrf,ialefp1,ialefp2,ialefp3,ialefp4,                           &
     & ialefp5,ialefp6,ialefp7,ialefp8,ialefp9,ialefp10,ialeprs1,       &
     & ialeprs2,                                                        &
     & ialeprs3,ialeprs4,ialeprs5,ialeprs6,ialeprs7,jalebkt1,jalebkt2,  &
     & jalebkt3,                                                        &
     & jalebkt4,jalebkt5,ialebc,ialebc1,ialehv,                         &
     & icplfq,ibktfq,ltie37,ltie38,ltie39,                              &
     & ltie40,ltie41,ltie42,ltie43,ltie44,ltie45,ltie46,ltie47,ltie48,  &
     & jprefbymmg,ialetie_nerod,ialetie_nerodprev,ialehydro_maxnode,    &
     & ialethk,ialesft,naleswmg,nalelp,lalelp,iale_edge,                &
     & ialepar,naltur,ltur1,ltur2,ltur3,nambhydro,iphydro,ninihydro,    &
     & iihydro,nale_his,nale30e,nale30f,ialemvbnd,iale_2d,iptadj2d,     &
     & iale_dtx,iale_mono,                                              &
     & nalefailswmg,jalefailswmg,jalefailswmg01,iadvec_mat8_tb,iale_1d, &
     & iptadj1d,naleebc,nalelpe,ialevap,jalevap,ntrale,jale12a,         &
     & ialectrl,noalectype2,ialtie_mem4ales3e,iale_volwgt,iale_nshare,  &
     & ialerez_nfac,ialectype2,ialephase,ialetie_nbyedge,ialetie_nedges,&
     & ialtie_maxelh
      real ale_chkr
      real                                                              &
     & alerz1,alerz2,alerz3,alerz4,alebeg,aleend,voidmin,voidmin2,      &
     & voidmas,salefac,alemin,alemax,aleming,alemaxg,                   &
     & alerrr,voidfac,xiautobc,                                         &
     & ecoupl,alerefp,dtmsale,ethcoupl,alefac9,alectrl
!
      common/alecom1/                                                   &
     & alerz1,alerz2,alerz3,alerz4,alebeg,aleend,voidmin,voidmin2,      & !  1
     & voidmas,salefac,alemin(3),alemax(3),                             & !  9
     & alerrr,iprelax,voidfac,nxale,nyale,nzale,nalebn,nalegp,          & ! 17
     & iale,ialetp,ialenc,ialead,newall,naltie,ialefst,jalefst,         & ! 25
     & jale00,jale01,jale02,jale03,jale04,jale05,jale06,jale07,         & ! 33
     & jale08,jale09,jale10,jale11,jale12,jale13,jale14,jale15,         & ! 41
     & jale16,jale17,jale18,jale19,jale20,jale21,jale22,jale23,         & ! 49
     & jale24,jale25,jale26,jale27,jale28,jale29,jale30,jale31,         & ! 57
     & jale32,jale33,jale34,jale35,jale36,jale37,jale38,jale39,         & ! 65
     & lenale,ntrace,nalesc,nintrf,jntrf,jalenp,nalenp,nvoid,           & ! 73
     & idotie,ntrale,jale12a,aleming(3),alemaxg(3)                        ! 81
      common/alecom2a/                                                  &
     & ltie1,ltie1b,ltie2,ltie3,ltie4,ltie5,ltie6,                      & !  1
     & ltie7,ltie8,ltie9,ltie10,ltie11,ltie12,ltie13,ltie13b,           & !  8
     & ltie13c,ltie13d,ltie14,ltie15,ltie15b,ltie15c,ltie15d,           & ! 16
     & ltie15e,ltie15f,ltie15g,ltie16,                                  & ! 23
     & ltie17a,ltie17b,ltie18a,ltie18b,ltie19,ltie20,ltie21,ltie22,     & ! 27
     & ltie23,ltie24,ltie25,ltie26,ltie27,ltie28,ltie29,ltie30,         & ! 35
     & ltie31,ltie32,ltie33,ltie34,ltie35,ltie36,mxtquad,xiautobc,      & ! 43
     & ialtie_mem4ales3e,ltie12e,ialetie_nbyedge,ialetie_nedges,        & ! 51 
     & ialtie_maxelh                                                      ! 55
      common/alecom3/                                                   &
     & nalebcsg,ialrst,nlagran,ibrfrm,iotflw,nofp(19),                  & !  1 
     & lvoid1,lvoid2,lvoid3,lvoid4,iframv,jnk2,jnk3,jnk4,               & ! 25
     & jnk5,jnk6,jnk7,jnk8,jnk9,jnk10,jnk11,jnk12,                      & ! 33
     & jnk13,jnk14,jnk15,iale_act,iale_mt41skpinc,jale(4),jalesies        ! 41
      common/alecom4/                                                   &
     & numrf,ecoupl,ialefp1,ialefp2,ialefp3,ialefp4,                    & !  1 
     & ialefp5,ialefp6,ialefp7,ialefp8,ialefp9,ialefp10,ialeprs1,       & !  7
     & ialeprs2,                                                        & ! 14 
     & ialeprs3,ialeprs4,ialeprs5,ialeprs6,ialeprs7,jalebkt1,jalebkt2,  & ! 15
     & jalebkt3,                                                        & ! 22
     & jalebkt4,jalebkt5,alerefp,dtmsale,ialebc,ialebc1,ialehv,ethcoupl,& ! 23
     & icplfq,ibktfq,ltie37,ltie38,ltie39,                              & ! 31
     & ltie40,ltie41,ltie42,ltie43,ltie44,ltie45,ltie46,ltie47,ltie48,  & ! 36
     & jprefbymmg,ialetie_nerod,ialetie_nerodprev,ialehydro_maxnode       ! 45
      common/alecom5/                                                   &
     & ialethk,ialesft,naleswmg,nalelp,lalelp,iale_edge,                & !  1
     & ialepar,naltur,ltur1,ltur2,ltur3,nambhydro,iphydro,ninihydro,    & !  7 
     & iihydro,nale_his,ale_chkr,nale30e,nale30f,ialemvbnd,iale_2d,     & ! 15
     & iptadj2d,iale_dtx,iale_mono,                                     & ! 22
     & nalefailswmg,jalefailswmg,jalefailswmg01,iadvec_mat8_tb,         & ! 25
     & iale_1d,iptadj1d,alefac9,alectrl(10),ialectrl(10),               & ! 29
     & naleebc,nalelpe,ialevap,jalevap,noalectype2,iale_volwgt,         & ! 52
     & iale_nshare,ialerez_nfac,ialectype2,ialephase                      ! 58
!     
      integer                                                           &
     & nalcpl,nsldcpl,iregcpl,                                          &
     & lcpl01,lcpl02,lcpl03,lcpl04,lcpl05,lcpl06,lcpl07,lcpl08,         &
     & lcpl09,lcpl10,lcpl11,lcpl12,lcpl13,lcpl14,lcpl15,lcpl16,         &
     & lcpl17,lcpl18,lcpl19,lcpl20,lcpl21,lcpl22,lcpl23,lcpl24,         &
     & lcpl25,lcpl26,lcpl27,lcpl28,lcpl29,lcpl30,lcpl31,lcpl32,         &
     & lcpl33,lcpl34,lcpl35,lcpl36,lcpl37,lcpl38,lcpl39,lcpl40,         &
     & lcpl41,lcpl42,lcpl43,lcpl44,lcpl45,lcpl46,lcpl47,lcpl48
      common/alecpl/                                                    &
     & nalcpl,nsldcpl,iregcpl,                                          &
     & lcpl01,lcpl02,lcpl03,lcpl04,lcpl05,lcpl06,lcpl07,lcpl08,         &
     & lcpl09,lcpl10,lcpl11,lcpl12,lcpl13,lcpl14,lcpl15,lcpl16,         &
     & lcpl17,lcpl18,lcpl19,lcpl20,lcpl21,lcpl22,lcpl23,lcpl24,         &
     & lcpl25,lcpl26,lcpl27,lcpl28,lcpl29,lcpl30,lcpl31,lcpl32,         &
     & lcpl33,lcpl34,lcpl35,lcpl36,lcpl37,lcpl38,lcpl39,lcpl40,         &
     & lcpl41,lcpl42,lcpl43,lcpl44,lcpl45,lcpl46,lcpl47,lcpl48
!
      integer                                                           &
     & nalercp,                                                         &
     & lebc00,lebc01,lebc02,lebc03,lebc04,lebc05,lebc06,lebc07
      common/aleebc/                                                    &
     & nalercp,                                                         &
     & lebc00,lebc01,lebc02,lebc03,lebc04,lebc05,lebc06,lebc07
!
      integer                                                           &
     & lrcp00,lrcp01,lrcp02,lrcp03,lrcp04,lrcp05,lrcp06,lrcp07
      common/alercp/                                                    &
     & lrcp00,lrcp01,lrcp02,lrcp03,lrcp04,lrcp05,lrcp06,lrcp07
!
! FSI and sensor
      integer        ialedb1,ialedb1b,ialedb2,ialedb3,ialedb4,ialedb5,  & ! 1
     & ialedb6,ialedb7,ialedb8,ialedb9,ialedb10,ialedb11,ialedb12,      & ! 7
     & ialedb13,ialedb14,ialedb15,ialedb16,ialedb17,ialedb17b,ialedb18, & ! 14
     & ialedb19,ialedb20,ialedb21,ialedb22,ialedb23,ialedb24,ialedb25,  & ! 21
     & ialeunit,ialesr1,ialesor,jaledbmom                                 ! 28
      common/aledata/ialedb1,ialedb1b,ialedb2,ialedb3,ialedb4,ialedb5,  & ! 1
     & ialedb6,ialedb7,ialedb8,ialedb9,ialedb10,ialedb11,ialedb12,      & ! 7
     & ialedb13,ialedb14,ialedb15,ialedb16,ialedb17,ialedb17b,ialedb18, & ! 14
     & ialedb19,ialedb20,ialedb21,ialedb22,ialedb23,ialedb24,ialedb25,  & ! 21
     & ialeunit,ialesr1,ialesor,jaledbmom                                 ! 28
!
      integer        iale_ups,nale_ups,ialeid
      real           aleups_st,aleups_tmas
      common/aleups/ iale_ups,aleups_st,aleups_tmas,nale_ups,ialeid(8)
!
      integer        ialeps1,ialeps2,ialeps3,ialeps4,ialeps5,ialeps6,   & ! 1
     & ialeps7,ialeps8,ialeps9,ialeps10                                   ! 7
      common/aleabag/ialeps1,ialeps2,ialeps3,ialeps4,ialeps5,ialeps6,   & ! 1
     & ialeps7,ialeps8,ialeps9,ialeps10                                   ! 7
!
      integer          iale_modeps, nale_modeps, iale_modeps_ptr,       &
     & mpp_modeps_ptr,mpp_modact,mpp_nmax,mpp_epsnmax
      common/alemodeps/iale_modeps, nale_modeps, iale_modeps_ptr(100),  &
     & mpp_modeps_ptr(97),mpp_modact,mpp_nmax,mpp_epsnmax
!
      integer           nale_fsiproj, iale_fsiproj_ptr,                 &
     & mpp_fsiproj_ptr
      integer ale_pjptr,mpp_npseg,mpp_npmax
      common/alefsiproj/nale_fsiproj, iale_fsiproj_ptr(100),            &
     & mpp_fsiproj_ptr(97),mpp_npseg,mpp_npmax,ale_pjptr
!
!  for the new ale advection and ctype 24
!  ADD ALE NEW COMMON BLOCK
!
      integer ialenew,iadvnew,ictype24,iale_segn,iale_sndn,iale_sndv,   &
     &        ifake12,iale_nofsifull
      common/alenew/ialenew,iadvnew,ictype24,iale_segn,iale_sndn,       &
     & iale_sndv,ifake12,iale_nofsifull
!
!  io buffer for ale forces file
      character*256   nameale
      real            aleint,alecut
      integer*8       ialeadd64
      integer         ioba9,id3ale,nifpale,ntotsg,ialeifp,inameale,     &
     &                ialeint_cid
      equivalence (inameale,nameale)
      common/aleiobuf/ioba9,id3ale,ialeadd64,aleint,alecut,             &
     & nifpale,ntotsg,ialeifp,inameale(64),ialeint_cid
!
!  io buffer for ale interface linking file
      character*256   name_alnk
      real            alnkint,alnkcut,fsilnk_t0,fsilnk_t1
      integer*8       ialnkadd64,ialnkadd64_sav
      integer         iobal9,id3alnk,nifpalnk,ntotsglnk,                &
     &                ialnkifp,iname_alnk,nzfsilnk,ngfsilnk,            &
     &                ialnkint_cid
      equivalence (iname_alnk,name_alnk)
      common/fsilnkiobuf/ iobal9,id3alnk,ialnkadd64,ialnkadd64_sav,     &
     & alnkint,alnkcut,nifpalnk,ntotsglnk,ialnkifp,nzfsilnk,ngfsilnk,   &
     & fsilnk_t0,fsilnk_t1,ialnkint_cid,                                &
     & iname_alnk(64)
!
      real           alest2,rlengmax
      common/sumtime/alest2,rlengmax
!
      integer nmtconx,nmcx,nmtconx1
      integer nmtconxo,nmcxo,nmtconx1o
      integer MAX_ALEG,MAX_ALEP
      parameter (MAX_ALEG=100,MAX_ALEP=1000)
      common/alenmt/                                                    &
     &   nmtconx(MAX_ALEG), nmcx(MAX_ALEG+1), nmtconx1(MAX_ALEP)
      common/alenmto/                                                   &
     &  nmtconxo(MAX_ALEG),nmcxo(MAX_ALEG+1),nmtconx1o(MAX_ALEP)
!
! storing MMG energy
      real ale_matsum
      common/aleenergy/ ale_matsum(10,MAX_ALEP)
!
! ALE fragmentation
      integer nalefrag,jalefrag
      common /ammgfrag/ nalefrag,jalefrag
!
! ALE burn switch mmg
      integer                  nalebrnswmg,ialebrnswmg_nvartmp,         &
     &                         jalebrnswmg1,jalebrnswmg2,jalebrnswmg3,  &
     &                         jalebrnswmg4,jalebrnswmg5,jalebrnswmg6,  &
     &                         jalebrnswmg7,jalebrnswmg8
      common/aleburnswitchmmg1/nalebrnswmg,ialebrnswmg_nvartmp,         &
     &                         jalebrnswmg1,jalebrnswmg2,jalebrnswmg3,  &
     &                         jalebrnswmg4,jalebrnswmg5,jalebrnswmg6,  &
     &                         jalebrnswmg7,jalebrnswmg8
!
!  xcell2d pointer  
      integer jale2dgmcsub,jale2dgvcsub,jale2dgmpsub
      common/alexcell/jale2dgmcsub,jale2dgvcsub,jale2dgmpsub     
!
!  
      integer icntlunit
      common /aleunit/ icntlunit
!
!  ale fsi 2 load node
      integer nalefsi2loadnode,jfsi2ld
      integer ifsi2ld_comm,ifsi2ld_npp,ifsi2ld_iam,ifsi2ld_iam0
      integer   ifsi2ld_iob,ifsi2ld_iobx,ifsi2ld_dir
      integer*8 ifsi2ld_add64,ifsi2ld_addx64
      character*80 fsi2ld_dir
      equivalence (ifsi2ld_dir,fsi2ld_dir)
      common/fsi2ldcom1/nalefsi2loadnode,jfsi2ld(9,20),ifsi2ld_dir(10)
      common/fsi2ldcom2/ifsi2ld_add64,ifsi2ld_addx64
      common/fsi2ldcom3/ifsi2ld_comm,ifsi2ld_npp,ifsi2ld_iam,           &
     &                  ifsi2ld_iob,ifsi2ld_iobx
!
! database ale mat
      integer ndatalemat,idatalemat_iob                                 &
     &,idatalemat_ptr1,idatalemat_ptr2,idatalemat_ptr3,idatalemat_ptr4
      common/datalemat1/ndatalemat,idatalemat_iob                       &
     &,idatalemat_ptr1,idatalemat_ptr2,idatalemat_ptr3,idatalemat_ptr4
      integer*8 idatalemat_addr,idatalemat_addw
      common/datalemat2/idatalemat_addr,idatalemat_addw
!
! database ale
      integer ndatale,idatale_iob,idatale_flg,idatale_mod3plot          &
     &,idatale_ptr1,idatale_ptr2,idatale_ptr3,idatale_ptr4,idatale_ptr5 &
     &,idatale_ptr6,idatale_ptr7,idatale_ptr8,idatale_ptr9,idatale_ptr10
      common/datale/ndatale,idatale_iob,idatale_mod3plot                &
     &,idatale_ptr1,idatale_ptr2,idatale_ptr3,idatale_ptr4,idatale_ptr5 &
     &,idatale_ptr6,idatale_ptr7,idatale_ptr8,idatale_ptr9,idatale_ptr10
      integer*8 idatale_addr,idatale_addw
      common/datale1/idatale_addr,idatale_addw
!
! database ale operation
      integer ndataleop,idataleop_nby,idataleop_flg,idataleop_iob       &
     &,idataleop_ptr1,idataleop_ptr2,idataleop_ptr3,idataleop_ptr4      &
     &,idataleop_ptr5,idataleop_ptr6,idataleop_ptr7,idataleop_ptr8
      common/dataleop/ndataleop,idataleop_iob                           &
     &,idataleop_ptr1,idataleop_ptr2,idataleop_ptr3,idataleop_ptr4      &
     &,idataleop_ptr5,idataleop_ptr6,idataleop_ptr7,idataleop_ptr8
      integer*8 idataleop_addr,idataleop_addw
      common/dataleop1/idataleop_addr,idataleop_addw
      parameter(idataleop_nby=10)
!
! ale switch mmg
      integer MAX_ALESWMMG,MAX_ALESWMMG1
      parameter(MAX_ALESWMMG=10,MAX_ALESWMMG1=4)
      integer naleswitchmmg,ialeswitchmmg_i,jaleswitchmmg
      common/aleswitchmmg1/naleswitchmmg,ialeswitchmmg_i,               &
     &jaleswitchmmg(MAX_ALESWMMG,MAX_ALESWMMG1)
!
! lagrangian porous element
      integer nlagpor,                                                  &
     &        jlagpor1,jlagpor2,jlagpor3,jlagpor4,jlagpor5
      common/lagpor1/nlagpor,                                           &
     &        jlagpor1,jlagpor2,jlagpor3,jlagpor4,jlagpor5
!
! initial ale eos
      integer ninieos,inieos_setele,                                    &
     &        jinieos1,jinieos2,jinieos3,jinieos4,jinieos5,             &
     &        jinieos6,jinieos7,jinieos8,jinieos9,jinieos10
      common/inieos1/ninieos,inieos_setele,                             &
     &        jinieos1,jinieos2,jinieos3,jinieos4,jinieos5,             &
     &        jinieos6,jinieos7,jinieos8,jinieos9,jinieos10
!
! ale smoothing
      integer jale09a,jale09b,jale09c,jale09nsnd,jale09nrcv,jale09dsnd, &
     &        jale09drcv,jale09isnd,jale09ircv,jale09xsnd,jale09xrcv,   &
     &        nalescmax
      common/alesmoothing/jale09a,jale09b,jale09c,jale09nsnd,jale09nrcv,&
     &                    jale09dsnd,jale09drcv,jale09isnd,jale09ircv,  &
     &                    jale09xsnd,jale09xrcv,nalescmax
!
! boundary ambient eos
      integer ialeaad
      common/aleaad/ialeaad(25)
      integer jaleaad1,jaleaad2,jaleaad3,jaleaad4,jaleaad5,             &
     &        jaleaad6,jaleaad7,jaleaad8,jaleaad9,jaleaad10
      common/aleaad1/                                                   &
     &        jaleaad1,jaleaad2,jaleaad3,jaleaad4,jaleaad5,             &
     &        jaleaad6,jaleaad7,jaleaad8,jaleaad9,jaleaad10
!
! boundary ambient
      integer iamboun_nbyheader
      parameter (iamboun_nbyheader=16)
      integer namboun,nalenrf,naleamb,                                  &
     &        jamboun1,jamboun2,jamboun3,jamboun4,jamboun5,             &
     &        jamboun6,jamboun7,jamboun8,jamboun9,jamboun10
      common/amboun1/namboun,nalenrf,naleamb,                           &
     &        jamboun1,jamboun2,jamboun3,jamboun4,jamboun5,             &
     &        jamboun6,jamboun7,jamboun8,jamboun9,jamboun10
!
! ale prescribed motion
      integer nalemvlc,                                                 &
     &        jalemvlc1,jalemvlc2,jalemvlc3,jalemvlc4,jalemvlc5,        &
     &        jalemvlc6,jalemvlc7,jalemvlc8,jalemvlc9,jalemvlc10
      common/alemvlc1/nalemvlc,                                         &
     &        jalemvlc1,jalemvlc2,jalemvlc3,jalemvlc4,jalemvlc5,        &
     &        jalemvlc6,jalemvlc7,jalemvlc8,jalemvlc9,jalemvlc10
!
! ale injection
      integer injctale_nbyheader
      parameter (injctale_nbyheader=44)
      integer njctale,injctale_nnd,                                     &
     &        injctale_nbyseg,injctale_nsegbynd,rnjctale_sizmax,        &
     &        jctale1,jctale2,jctale3,jctale4,jctale5,                  &
     &        jctale6,jctale7,jctale8,jctale9,jctale10
      common/njctale1/njctale,injctale_nnd,                             &
     &        injctale_nbyseg,injctale_nsegbynd,rnjctale_sizmax,        &
     &        jctale1,jctale2,jctale3,jctale4,jctale5,                  &
     &        jctale6,jctale7,jctale8,jctale9,jctale10
!
! control segments in ale coupling
      integer nalexcplsg,ialexcplsg_nclsgtot,                           &
     &        ialexcplsg_nndsnd,ialexcplsg_nsgrcvtot,                   &
     &        jalexcplsg1,jalexcplsg2,jalexcplsg3,jalexcplsg4,          &
     &        jalexcplsg5,jalexcplsg6,jalexcplsg7,jalexcplsg8,          &
     &        jalexcplsg9
      common/alexcplsg1/nalexcplsg,ialexcplsg_nclsgtot,                 &
     &        ialexcplsg_nndsnd,ialexcplsg_nsgrcvtot,                   &
     &        jalexcplsg1,jalexcplsg2,jalexcplsg3,jalexcplsg4,          &
     &        jalexcplsg5,jalexcplsg6,jalexcplsg7,jalexcplsg8,          &
     &        jalexcplsg9
!
! ale and sph mesh box and cylinder
      integer ialemeshcy,ialemeshbx,ialemesh1,ialemesh2
      real    uu1,ralemesh1,uu2,ralemesh2
      common/alemesh0/ialemeshcy,uu1(9),ialemesh1(10),ralemesh1(5)
      common/alemesh1/ialemeshbx,uu2(9),ialemesh2(5),ralemesh2(6)
!
! ale mesh interface
      integer ialemshint_nbyheader
      parameter (ialemshint_nbyheader=10)
      integer nalemshint,                                               &
     &        jalemshint1,jalemshint2,jalemshint3,jalemshint4,          &
     &        jalemshint5,jalemshint6,jalemshint7,jalemshint8
      common/alemshint1/nalemshint,                                     &
     &        jalemshint1,jalemshint2,jalemshint3,jalemshint4,          &
     &        jalemshint5,jalemshint6,jalemshint7,jalemshint8
!
! s-ale
      integer mmalemsh,mmaledecomp,mmalenhsv,numelh_sale,mmalealtdir
      common/mmale/mmalemsh,mmaledecomp,mmalenhsv,numelh_sale,          &
     &        mmalealtdir
      integer mmsh01,mmsh02,mmsh03,mmsh04,mmsh05,mmsh06,mmsh07,mmsh08,  &
     &        mmsh09,mmsh10,mmsh11,mmsh12,mmsh13,mmsh14,mmsh15,mmsh16
      common/mmmsh/                                                     &
     &        mmsh01,mmsh02,mmsh03,mmsh04,mmsh05,mmsh06,mmsh07,mmsh08,  &
     &        mmsh09,mmsh10,mmsh11,mmsh12,mmsh13,mmsh14,mmsh15,mmsh16
      integer mmale01,mmale02,mmale03,mmale04,mmale05,mmale06,mmale07,  &
     &        mmale08,mmale09,mmale10,mmale11,mmale12,mmale13,mmale14,  &
     &        mmale15,mmale16,mmale17,mmale18,mmale19,mmale20,mmale21,  &
     &        mmale22,mmale23,mmale24,mmale25,mmale26,mmale27,mmale28,  &
     &        mmale29,mmale30,mmale31,mmale32,mmale33,mmale34,mmale35,  &
     &        mmale36,mmale37,mmale38,mmale39,mmale40,mmale41,mmale42 
      common/mmalemem/                                                  &
     &        mmale01,mmale02,mmale03,mmale04,mmale05,mmale06,mmale07,  &
     &        mmale08,mmale09,mmale10,mmale11,mmale12,mmale13,mmale14,  &
     &        mmale15,mmale16,mmale17,mmale18,mmale19,mmale20,mmale21,  &
     &        mmale22,mmale23,mmale24,mmale25,mmale26,mmale27,mmale28,  &
     &        mmale29,mmale30,mmale31,mmale32,mmale33,mmale34,mmale35,  &
     &        mmale36,mmale37,mmale38,mmale39,mmale40,mmale41,mmale42
! s-ale volfilling
      integer mmalevolfil
      common/mmalevlf/mmalevolfil
      character*80 salevfcnm
      integer ilsdavfc
      common/mmaletmp/salevfcnm,ilsdavfc
! s-ale ambient
      integer mmaleamb,mmaleab01,mmaleab02,mmaleab03,mmaleab04,         &
     &  mmaleab05,mmaleab06,mmaleab07,mmaleab08,mmaleab09,mmaleab10,    &
     &  mmalenrf
      common/mmaleab/mmaleamb,mmaleab01,mmaleab02,mmaleab03,mmaleab04,  &
     &  mmaleab05,mmaleab06,mmaleab07,mmaleab08,mmaleab09,mmaleab10,    &
     &  mmalenrf
! s-ale mesh motion
      integer mmalemov,mmalemv01,mmalemv02,mmalemv03,mmalemv04,         &
     & mmalemv05,mmalemv06,mmalemv07,mmalemv08,mmalemv09,mmalemv10
      common/mmalemot/mmalemov,mmalemv01,mmalemv02,mmalemv03,mmalemv04, &
     & mmalemv05,mmalemv06,mmalemv07,mmalemv08,mmalemv09,mmalemv10
! s-ale mpp
      integer mmale_nelhsa,mmalesf_size,mmalesn_size,salecomm,saleiam,  &
     & mpp_nsnsa,mmalesfa_size,salenod_dch,saleface_dch,salemppt1
      common/mmalempp/mmale_nelhsa,mmalesf_size,mmalesn_size,salecomm,  &
     & saleiam,mpp_nsnsa,mmalesfa_size,salenod_dch,saleface_dch,        &
     & salemppt1
! s-ale misc
      integer mmalemshup,mmalemis
      common/mmalemisc/mmalemshup,mmalemis(9)
!
! cbis mpp commnunication put here tempoarily
      integer cbiscomm,cbisnpp,cbisiam,dbfcbis,dbcbiscomm,              &
     & dbcbisiam,dbcbispat1,dbcbispat2,ndbcbis,idbcbis,cbisnod_dch,     &
     & cbisnorcomm,cbisnor_dch,cbist3,cbist4
      common /cbismpp/cbiscomm,cbisnpp,cbisiam,dbfcbis,dbcbiscomm,      &
     & dbcbisiam,dbcbispat1,dbcbispat2,ndbcbis,idbcbis,cbisnod_dch,     &
     & cbisnorcomm,cbisnor_dch,cbist3,cbist4
