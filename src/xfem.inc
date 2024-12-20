c     Include file for fracture analysis with xfem
c     ictrxfem     0 not an xfem analysis
c                  2 2d/shell only
c                  3 3d only
c                  5 shell and 3d
      integer        ictrxfem,ictrxfem2,iformxfem,icohesive,ixfemgotip,
     &               ixfemerosion,ixfeminit,ixfemrepl,ndimxfem,
     &               ixfempropcr,ixfemcrklen,ixfemdam,ixfemjction
      common/xfemctr/ictrxfem,ictrxfem2,iformxfem,icohesive,ixfemgotip,
     &               ixfemerosion,ixfeminit,ixfemrepl,ndimxfem,
     &               ixfempropcr,ixfemcrklen,ixfemdam,ixfemjction
c
      integer       ixfemmass,ixfemndfc,ixfemsmorder,iformjump,icontip
      common/xfemfm/ixfemmass,ixfemndfc,ixfemsmorder,iformjump,icontip
c
c     cohesive material constants
c     a(ndtype)    mctype(*)
c     a(ndcohm)    cmcoh(48,*)
c     nmtcox       number of internal variables in each xfem
c     nmtcoh       number of internal variables in each cohesive zone
      integer       ndtype,ndcohm,nmtcox,nmtcoh
      common/cohmat/ndtype,ndcohm,nmtcox,nmtcoh
c       
c     element status
c     xf_inxfems   inxfems(*) local (xfem parts) eid to global eid
c     xf_ixfems    ixfems(*) for 2d/shell elements
c     xf_isbfem    isbfem(*) =1 boundary elements
c     xf_ckopen    ckopen(35,*) crack opening angle, location & nodal forces
c     xf_neighbors neighbors(4,*) direct front neighboring elements of xfem eid
c     xf_numbors   numbors(4,*) number of corner neighboring elements of xfem eid
c     xf_neighbor2 neighbor2(6,4,*) corner neighboring elements of xfem eid
c     xf_nbjction  nbjction(2,4,*) front neighbors from junction part
c     xf_numjction numjction(4,*) number of corner neighbors from junction
c     xf_nb2jction nb2jction(6,4,*) corner neighbors from junction part
c     xf_thicks    original shell thickness
c     xf_ifails    number of times fem fails but cannot become xfem
c     xf_rfails    ratio of failure value to smoothed value
      MINDEX xf_ixfems,xf_isbfem,xf_ckopen,xf_neighbors,xf_inxfems,
     &       xf_thicks,xf_ifails,xf_rfails,xf_numbors,xf_neighbor2,
     &       xf_nbjction,xf_numjction,xf_nb2jction
      common/xfemsts/xf_ixfems,xf_isbfem,xf_ckopen,xf_neighbors,
     &               xf_inxfems,xf_thicks,xf_ifails,xf_rfails,
     &               xf_numbors,xf_neighbor2,
     &               xf_nbjction,xf_numjction,xf_nb2jction
c
c     xfem part info
c     ngxfem       number of xfem groups (including triangle group for base 16)
c     npxfem       number of xfem parts
c     ntxfem       number of elements in xfem parts
c     xf_nelgxfem       group of eid
c     xf_mxexfem        mxe of the group
c     xf_mxpxfem        mxe of the part
c     xf_itegxfem       starting xfem eid of each group 
c     xf_nbpart         number of xfem parts connected to current part
c     xf_nbparts        part ids (among npxfem) of the neighboring xfem parts
c     xf_props          xfem part related data
      integer        ngxfem,ntxfem,npxfem
      MINDEX xf_nelgxfem,xf_mxexfem,xf_mxpxfem,xf_itegxfem,
     &       xf_nbpart,xf_nbparts,xf_props
      common/xfeminfo/ngxfem,ntxfem,npxfem
      common/xfemprts/xf_nelgxfem,xf_mxexfem,xf_mxpxfem,xf_itegxfem,
     &                xf_nbpart,xf_nbparts,xf_props
c
c     non-local strain 
c     xf_xeps(2*ntxfem) non-local strain:local strain
c     xf_nipxfem        nip for strain regularization
c     xf_nelxfem        eid in xfem parts for strain reg.
c     xf_dcjp           support
c     xf_xyzc           coordinates of element centers 
c     xf_ntipfem        if tip front element, pointed to xytipfem       
c     xf_xytipfem       local tip coordinate to reg. strain
c     xf_nipavg         nip for strain center
c     xf_nelavg         eid in xfem parts for strain center
      MINDEX xf_xeps,xf_nipxfem,xf_nelxfem,xf_dcjp,
     &       xf_xyzc,xf_ntipfem,xf_xytipfem,xf_nipavg,xf_nelavg
      common/xfemsmth/xf_xeps,xf_nipxfem,xf_nelxfem,xf_dcjp,
     &                xf_xyzc,xf_ntipfem,xf_xytipfem,xf_nipavg,xf_nelavg
      integer MXFEMIP
c     parameter (MXFEMIP=3500)
      parameter (MXFEMIP=1200)
c
c     crack description
c     ncracks      number of cracks
c     mcracks      allocated number of cracks
c     nxfems       number of xfem elements
c     mxfems       allocated number of xfem elements
c     npnods       number of phantom nodes
c     mpnods       allocated number of phantom nodes
c     nxact        number of active xfems
c     nactcks      number of active cracks (after merging)
      integer      ncracks,nxfems,npnods,mcracks,mxfems,mpnods,
     &             nxact,nprexfem,iphstat,nactcks
      real         crklength
      common/cracks/ncracks,nxfems,npnods,mcracks,mxfems,mpnods,
     &              nxact,nprexfem,iphstat,nactcks,crklength
c     ck_nxfems    number of xfems for each crack
c     ck_pairs     element pairs (orig & phantom) for cracks
c     ck_x         crack position in reference coordinates
c     ck_nside     element edges cut by crack
c     ck_ixesta    crack start element (xfem id)
c     ck_ixeend    crack end element (xfem id)
c     ck_cracks    crack id for active cracks
      MINDEX ck_nxfems,ck_pairs,ck_x,ck_ixesta,ck_ixeend,ck_cracks,
     &       ck_nside
      common/ckindex/ck_nxfems,ck_pairs,ck_x,ck_ixesta,ck_ixeend,
     &               ck_cracks,ck_nside
c
c     xfem database
c     xf_phnods    phantom nodes ids
c     xf_m         phantom node mass
c     xf_mr        phantom node mass (rotation)
c     xf_x         phantom node position
c     xf_disp      phantom node displacements
c     xf_v         phantom node velocities
c     xf_a         phantom node acceleration/rhs
c     xf_nbc       boundary conditions on phantom nodes
c     xf_nphsta    phamtom node status
c     xf_npcxfm    nodal constraint at crack tip
c     xf_cpcxfm    nodal constraint parameter at crack tip
c     xf_mast      added mass for phantom node
c     xf_masr      added mass for phantom node (rotation)
c     xf_dt        time step for xfem
      MINDEX xf_phnods,xf_m,xf_x,xf_disp,xf_v,xf_a,xf_nbc,xf_mr
      MINDEX xf_nphsta,xf_npcxfm,xf_cpcxfm,xf_mast,xf_masr,xf_dt
      common/xfnods/xf_phnods,xf_m,xf_x,xf_disp,xf_v,xf_a,xf_nbc,xf_mr,
     &  xf_nphsta,xf_npcxfm,xf_cpcxfm,xf_mast,xf_masr,xf_dt
c     xf_phelms    phantom element connectivity
c     xf_auxvec    memory for internal variables in phantom elements
c     xf_lavn      start location for internal variables in ph-elements
c     xf_auxcoh    memory for internal variables in cohesive interface
c     xf_lavcoh    start location for int. variables in cohesive zone
c     xf_aratio    area ratios of the two sub-region to whole element
c     xf_ipxfem    number of integration points for xfem (subdomain)
c     xf_ndxfem    nodal numbers in subdomains 
      MINDEX xf_phelms,xf_auxvec,xf_lavn,xf_auxcoh,xf_lavcoh,
     &       xf_aratio,xf_ipxfem,xf_ndxfem
      common/xfelms/xf_phelms,xf_auxvec,xf_lavn,xf_auxcoh,xf_lavcoh,
     &              xf_aratio,xf_ipxfem,xf_ndxfem
c
      integer        ntncog,ncoir
      common/coeleir/ntncog,ncoir(20)
      integer nxfemco,nxfem2
      common/xfemtype2/nxfemco,nxfem2(20)
#ifdef MPP
      integer nteshmppin,nteshmppout,nphshmpp
      integer neinmppxfem,neoutmppxfem
      common/xfemmpp1/nteshmppin,nteshmppout,nphshmpp,
     &       neinmppxfem,neoutmppxfem
      MINDEX mppxfemprocin,mppxfemneshin,xyzshmppxfem,mppxfemphsh,
     &       mppxfemprocout,mppxfemneshout,mppxfemshne,mppxfembfem,
     &       xytipmppxfem,crkangmppxfem,xshmppxfem,fvalmppxfem
      common/xfemmpp2/mppxfemprocin,mppxfemneshin,xyzshmppxfem,
     &       mppxfemprocout,mppxfemneshout,mppxfemshne,mppxfemphsh,
     &       mppxfembfem,xytipmppxfem,crkangmppxfem,xshmppxfem,
     &       fvalmppxfem
      MINDEX mppxfemnein,mppxfemneinproc,xyzc0mppxfem,
     &       mppxfemneout,mppxfemneoutproc,xyztipmppxfem,
     &       xyztip2mppxfem,mppcounter,strnmppxfem,xyzcmppxfem
      common/xfemmpp3/mppxfemnein,mppxfemneinproc,xyzc0mppxfem,
     &       mppxfemneout,mppxfemneoutproc,xyztipmppxfem,
     &       xyztip2mppxfem,mppcounter,strnmppxfem,xyzcmppxfem
      MINDEX ctipmppxfem,cangmppxfemout,cangmppxfemin,ncelmmppxfem
      common/xfemmpp4/ctipmppxfem,cangmppxfemout,cangmppxfemin,
     &       ncelmmppxfem
#endif
