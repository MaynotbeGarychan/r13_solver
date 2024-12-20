!
! sph memory pointers
      integer  ldsph00,ldsph01,ldsph02,ldsph03,ldsph04,ldsph05,         &
     & ldsph06,ldsph07,ldsph08,ldsph09,ldsph10,ldsph11,ldsph12,         &
     & ldsph13,ldsph14,ldsph15,ldsph16,ldsph17,ldsph18,ldsph19,         &
     & ldsph20,ldsph21,ldsph22,ldsph23,ldsph24,ldsph25,ldsph26,         &
     & ldsph27,ldsph28,ldsph29,ldsph30,ldsph31,ldsph32,ldsph33,         &
     & ldsph34,ldsph35,ldsph36,ldsph37

      common/sph17/ldsph00,ldsph01,ldsph02,ldsph03,ldsph04,ldsph05,     &
     & ldsph06,ldsph07,ldsph08,ldsph09,ldsph10,ldsph11,ldsph12,         &
     & ldsph13,ldsph14,ldsph15,ldsph16,ldsph17,ldsph18,ldsph19,         &
     & ldsph20,ldsph21,ldsph22,ldsph23,ldsph24,ldsph25,ldsph26,         &
     & ldsph27,ldsph28,ldsph29,ldsph30,ldsph31,ldsph32,ldsph33,         &
     & ldsph34,ldsph35,ldsph36,ldsph37
!
! sph bk01
      real bgnsph,sphxl,sphxu,sphyl,sphyu,sphzl,sphzu,                  &
     & dtsphini,expdmsph,strtm,xvsph,fsphshw,                           &
     & sph_boxc(3),sph_gl(3,3),psphcoeff,tisphmax,tisphavg,tisphromax
      integer iactsph,nsph_extra,ispherd,isphcont,nzsphcpl,             &
     & nsphdecpl,nsph_react,nsph_reg,isphtmp(5),                        &
     & isphterode,isphavisc,isphshw,msph_extra,                         &
     & nsph_ntot,irf_sph,isph_osu,nsph_inj,isph_inj,isphthk,            &
     & isphstab,nsphidh,nzpesphcp,nsph_gfirst,nsph_lfirst,nzsphmfp,     &
     & nsph_av,ipsmoot,nsph_act,ntot_act,nadpsph,isphpfem,nzsphvicsens
      common/sphbk01/bgnsph,sphxl,sphxu,sphyl,sphyu,sphzl,sphzu,        &
     & dtsphini,expdmsph,strtm,xvsph,iactsph,nsph_extra,ispherd,        &
     & isphcont,nzsphcpl,nsph_react,isphtmp,isphterode,isphavisc,       &
     & isphshw,fsphshw,msph_extra,nsph_ntot,nsphdecpl,nsph_reg,         &
     & ipsmoot,sph_boxc,sph_gl,irf_sph,isph_osu,nsph_inj,isph_inj,      &
     & isphthk,isphstab,nsphidh,psphcoeff,nzpesphcp,nsph_gfirst,        &
     & nsph_lfirst,nzsphmfp,nsph_av,nsph_act,ntot_act,nadpsph,          &
     & isphpfem,tisphmax,tisphavg,tisphromax,nzsphvicsens
!
! sph bk03 adaptice blocks
! ksphsurf(7): unused scratch space
      real    adpsph
      integer nfstsph,nzadpsph,iadpsph,ifailcouple,nzsphsurf,           &
     & ihsphsurf,isphsurfc,isphsurfe,ksphsurf(7),ihsphav,nzsphobj
      common /sphbk03/ nfstsph,nzadpsph,adpsph(4),iadpsph(2),           &
     & ifailcouple,nzsphsurf,ihsphsurf,isphsurfc,isphsurfe,ksphsurf,    &
     & ihsphav,nzsphobj
!     data count per SPH OBJ       
      integer NCSTOBJ
      parameter (NCSTOBJ=7)
!
! sph bk04 sph shape functions for adp
      integer MAX_ADSPH,ADP_OFF
      parameter (MAX_ADSPH=8)
      parameter (ADP_OFF=MAX_ADSPH**3)
      real    sphh
      common /sphbk04/ sphh(8,3*ADP_OFF*MAX_ADSPH)
!
! sph extra variables for MPP
      integer mpp_sphflow,mpp_sphflw,mpp_dlen,mpp_ntsphh,mpp_ntsph,     &
     & mpp_ntot,mpp_sphg2lh,mpp_sphl2gh,mpp_sndmax,mpp_rcvmax,          &
     & mpp_sphtmp(12),memsize_sph
      common /mpp_sphbk01/ mpp_sphflow,mpp_sphflw,mpp_dlen,             &
     & mpp_ntsphh,mpp_ntsph,mpp_ntot,mpp_sphg2lh,mpp_sphl2gh,           &
     & mpp_sndmax,mpp_rcvmax,mpp_sphtmp,memsize_sph
! sph variable s for axisymmetry and 2.5 D option
      integer itwoh
      common /sphaxisym/itwoh
!
! sph extran memory handlers
      integer isphtempp,isphtbch,isphtbfh,isphbulkq,isphbqsz,           &
     & nsphtempp,hlistsmp,nlistsmp,hsph_smp,nsph_smp,isphtbcch,         &
     & isphtbrh,nsph_01,nblistsize
      common/sphtempt/isphtempp,isphtbch,isphtbfh,isphbulkq,isphbqsz,   &
     &     nsphtempp,hlistsmp,nlistsmp,hsph_smp,nsph_smp,isphtbcch,     &
     &     isphtbrh,nsph_01,nblistsize
! sph variables for INJECTION keyword
      integer isphinject,isphinjstart
      common/sphinject/isphinject,isphinjstart
! sph variables for STOCHASTIC keyword
      integer isphtable,isphhlist
      common/sphstochastic/isphtable,isphhlist
! sph variables for PERI coupling keyword
      integer ippesphcpl
      common/perisphcp/ippesphcpl
! sph variables for DEFINE_SPH_BOX and DEFINE_SPH_CYLINDER
      integer isphmeshcy,isphmeshbx
      common/sphmesh0/isphmeshcy,isphmeshbx
! sph variables for shifting stablizing keyword
      integer istabshift,isphtv,iptvsph
      common/shiftingstab/istabshift,isphtv,iptvsph
! sph memory buffer handler for IISPH, ITEM_IISPH=4 from sync_sph_cplfrcin
      integer ITEM_IISPH
      parameter (ITEM_IISPH=4)
      integer nh_iisphsync,ns_iisph,nr_iishp,ifstiisph,iisph_tmp
      common/mpp_membug_iisph/ nh_iisphsync,ns_iisph,nr_iishp,          &
     &      ifstiisph,iisph_tmp(16)
!sph to sph coupling      
      integer isph2sphcont
      common/sph2sphcont/isph2sphcont
!
! block for isph interface force file
!
      character*256   nameisph
      logical         lisphifp
      integer*8       iisphadd64,iisphadd64_save
      real            sphint,sphcut
      integer         iobisph9,id3isph,nifpisph,ntotsgisph,iisphifp,    &
     & inameisph,ifpisphnodes,mpp_isphsurfc,mpp_isphsurfe,isphifpvar,   &
     & iisphpb
      equivalence (inameisph,nameisph)
      common/isphiobuf/iobisph9,id3isph,iisphadd64,iisphadd64_save,     &
     & sphint,sphcut,nifpisph,ntotsgisph,iisphifp,inameisph(64),        &
     & lisphifp,ifpisphnodes,mpp_isphsurfc,mpp_isphsurfe,isphifpvar,    &
     &  iisphpb
!
