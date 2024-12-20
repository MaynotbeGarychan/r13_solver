c
c ale database data
      integer ivolf,ialedb,ialeproj,nalemms,ialemms,naleffi,
     . ndbale,npspart,npsnods,ialereg,ivolf1,nalegasm,ndbsor,
     . iale_gpname,ndbale_t,ndbsor_t,ialextrcomblk,ihcple
      integer ndale_min,ndale_max,ndale_tmp,ndale_noseq
      common/alemisc/ ivolf,ialedb,ialeproj,nalemms,ialemms,naleffi,
     . ndbale,npspart,npsnods,ialereg,ivolf1,nalegasm,ndbsor,
     . iale_gpname,ndbale_t,ndbsor_t,ialextrcomblk,ndale_min,ndale_max,
     . ndale_tmp(12),ndale_noseq,ihcple
c
c ale trim
      real    trimtime
      integer notrim
      common/aletrim/ trimtime(10),notrim
c
c ale nlf
      real    alefstr,alefrad
      integer ialefps,nalefhx
      common/alenlf/ ialefps,nalefhx,alefstr,alefrad
c
c alelco
      integer lcoupl,ik1,ik4,iregflg,nregflg,ncaletie
      common/alelco/lcoupl(20),ik1(40),ik4(40),iregflg,nregflg,ncaletie
c
c ale d3plot d3drlf 
      integer ialemat
      common/aled3plot/ialemat
