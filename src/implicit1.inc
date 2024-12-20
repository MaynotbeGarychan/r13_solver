!   ... implicit common ...
!
      character*80    sc_filename, moddyn_filename, modal_op2_filename,  &
     &                rv_filename
      character*10    sc_mass, sc_damp, sc_stiff, sc_inert
      common/bki01ch/ sc_filename, sc_mass, sc_damp, sc_stiff, sc_inert, &
     &                moddyn_filename(2), modal_op2_filename,            &
     &                rv_filename
!
      integer imauto,iteopt,lauto,mthsol,ilimit,maxref,icnvrg,           &
     & igdiv,nwebuf,neql,neqt,imterm,imphas,nbfgs,                       &
     & numupd,istif,itrlas,imerr,imdof,neqtgl,ilsmtd,lsdir,ptr_aspc_h,   &
     & impcntctl,ipmaxx,igpx,itermp,itermp1,h_lp,h_rd,ptr_rf_data_h,     &
     & moddyn_flag, moddyn_mdrange, moddyn_mdlen, moddyn_mdfreq,         &
     & moddyn_mdfrqlen, moddyn_aux,                                      &
     & mf2_ptrnum_f95, mf2_ptrsym_f95, h_rsvc, im_numnp, mpp_im_numnp,   &
     & h_frf
      common/bki01i/imauto,iteopt,lauto,mthsol,ilimit,maxref,icnvrg,     &
     & igdiv,nwebuf,neql,neqt,imterm,imphas,nbfgs,                       &
     & numupd,istif,itrlas,imerr,imdof,neqtgl,ilsmtd,lsdir,ptr_aspc_h,   &
     & impcntctl,ipmaxx,igpx,itermp,itermp1,h_lp,h_rd,ptr_rf_data_h,     &
     & moddyn_flag, moddyn_mdrange, moddyn_mdlen, moddyn_mdfreq,         &
     & moddyn_mdfrqlen, moddyn_aux(13),                                  &
     & mf2_ptrnum_f95, mf2_ptrsym_f95, h_rsvc, im_numnp, mpp_im_numnp,   &
     & h_frf
!
      real dtimp,dtimp0,timeim,dtmnim,dtmxim,cvtl,ectl,rctl,             &
     & tolls,dnorm2,dtprnt,dtplot,dtiter,dtrefm,ppmax,pmaxr,errord,umax, &
     & sumx,xrte,moddyn_valdmp,moddyn_output_time,ctmxs,htmxs,xkei,      &
     & tollls,pretim
      common/bki01r/dtimp,dtimp0,timeim,dtmnim,dtmxim,cvtl,ectl,rctl,    &
     & tolls,dnorm2,dtprnt(2),dtplot(2),dtiter(2),dtrefm(2),             &
     & ppmax,pmaxr,errord,umax,sumx,xrte,moddyn_valdmp,                  &
     & moddyn_output_time(2),ctmxs(8),htmxs(8),xkei(-2:2),tollls,pretim
!
      real ascntl,r78tags
      common/bki02r/ascntl(200),r78tags(20)
!
      logical lsensw,jlsmtd
      common/bki01l/lsensw(20),jlsmtd
!
      integer imip,isolvr,icwrb,elem_tbl,i78tags,moddyn_fast,numsmv
      common/bki02i/imip(300),isolvr(500),icwrb(50),elem_tbl(4,65),      &
     & i78tags(20),moddyn_fast,numsmv
