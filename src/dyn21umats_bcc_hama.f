#include "define.inc"
#include "define2.inc"
      subroutine umatBccHama(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1   temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
!============================================================
! Declaration of constitutive variables
!------------------------------------------------------------
      implicit none
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
c     UMAT variables
      double precision cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*)
      double precision cma(*),qmat(3,3)
      integer nnpcrv(*)
      integer ::nhsv=600
      double precision dt1
      character*5 etype
      logical failel,reject
      INTEGER idele
c     
      integer umatType
c     IO
      integer oriType
c     Define the crystal
      integer,parameter:: typeCry=0
      integer,parameter:: numSys=24
c     Intermedia variables
      integer i,j,k,l
      double precision, parameter :: pi=acos(-1.0d0)
      double precision idm(3,3)
c     Kinetic model variables
      double precision f(3,3),f_n1(3,3),df_n1(3,3),f_n1_inv(3,3)
      double precision f_det
      double precision L_n1(3,3)
      double precision D(3,3),Dv(6)
      double precision W(3,3),Wv(3)
      double precision We(3,3),Wev(3)
      double precision We1
      double precision exp_we(3,3)
c     Kinematic variables
      double precision s11(3,numSys),m11(3,numSys)
      double precision s11_n1(3,numSys),m11_n1(3,numSys)
      double precision r(3,3),RL(6,6),r_n1(3,3)
      double precision Pa(3,3,numSys),eschmid(6,numSys)
      double precision Wa(3,3,numSys),wschmid(3,numSys)
      double precision euler(3),euler_n1(3)
      double precision Dp(6),Wp(3)
c     Slip system constitutive model variables
      double precision tau(numSys),dgamma(numSys)
      double precision gamma_slip(numSys)
      double precision dgamma_tol,gamma_n1
      double precision g_crss(numSys)
      double precision dgamma0,mval,dgamma_lim
c     Stress Update model variables
      double precision ym,pr,bk,sm
      double precision L_ela(6,6),L_ela_cry(6,6)
      double precision dsig_0(6),sig_jau(6),sig_r(6),sig_n1(6)
c     Hardening model variables
      integer hardType
      double precision g0,gs,h0,hs,q
c     Solid element variables
      double precision g,g2,gc,q1,q3,davg,p,deti,c22i,c23i,fac
      double precision temper,elsiz,epsp,capa,tt
!============================================================
! Declaration of non-constitutive variables
!------------------------------------------------------------
c     Stress state vaiables
      double precision sig_m,sig_eq,st
c     Strain variables
      double precision peeq
!============================================================
! Declaration of utilized functions
!------------------------------------------------------------
      double precision calSigEq,calSigMean,calSigTri,calPeeq
!============================================================
! Obtain variables from materials constants
!------------------------------------------------------------
      call      umatBccHamaGetMc(cm,ym,pr,bk,sm,
     1       umatType,dgamma0,mval,dgamma_lim,
     2       oriType,euler,
     3       hardType,g0,gs,h0,hs,q)
!============================================================
! Initial step: ncrycle = 0
!------------------------------------------------------------   
      if (.not.failel) then
      if(ncycle==0) then
c     Init crystal orientation, slip system vectors
            call initCrystal(oriType,typeCry,euler,
     1           numSys,r,s11,m11)
c     Initialize the hsv list
            call umatBccHamaInitHsv(g0,r,s11,m11,numSys,nhsv,
     1      hsv,sig)
      else
!============================================================
! Calculation begins: ncycle > 0
!------------------------------------------------------------ 
c     Obtain defromation gradient from hsv
      call getDispGradfromHsv(hsv,nhsv,f,f_n1)
c     Obtain CRSS from hsv
      call umatBccHamaGetHsv(hsv,numSys,g_crss,r,s11,m11,
     1       gamma_n1,gamma_slip)

!============================================================
! Kinetic model
!------------------------------------------------------------
      call calDeforSpinRateByDispGrad(f,f_n1,dt1,Dv,Wv)

!============================================================
! Constitutive model for slip at slip system
!------------------------------------------------------------ 
      call calSchmidTensor(s11,m11,numSys,
     1     Pa,Wa,eschmid,wschmid)
      call calSigRss(sig(1:6),eschmid,numSys,tau)
c      call calSlipRateVp(tau,g_crss,mval,dgamma0,
c     1   dgamma_lim,numSys,dgamma)
      call calSlipRateVpHeatAct(tau,g_crss,tau0110,tau0112,
     1   dgk0,dgamma0,pval,qval,tval,kb,numSys,dgamma)
      call updateCss(dgamma,numSys,dt1,
     1     dgamma_tol,gamma_slip,gamma_n1)
c     Project the slip deformation into macro deformation and spin
c     as plastic corrector
      call calStrainRateBySlip(dgamma,eschmid,numSys,Dp)
      call calSpinRateBySlip(dgamma,wschmid,numSys,Wp)

c   Update dislocation
      

!============================================================
! Update cauchy stress
!------------------------------------------------------------
c     Fourth-order elastic tensor at crystal coordinate
      call calElasIsoTensorCrystalGlobal(ym,pr,r,L_ela)
c     Update the Cauchy stress tensor
      call mat33Det(f,f_det)
      call updateSigJaum(sig,eschmid,wschmid,dgamma,Wv,
     1     Dv,L_ela,f_det,numSys,dt1,sig_n1)

!============================================================
! Rotation model
!------------------------------------------------------------
c     calculate the rotation rate for further rotation
      call updateOrientation(Wv,Wp,s11,m11,r,numSys,dt1,
     1  s11_n1,m11_n1,r_n1)
c     Extract the euler angle from the tranformation matrix
      call calEulerbyTransMat(r_n1,euler_n1)

!============================================================
! Hardening model
!------------------------------------------------------------
      call updateCrssFcc(g0,gs,h0,hs,q,gamma_n1,dgamma,
     1           dt1,g_crss)
!============================================================
! Calculation of non-constitutive variables
!------------------------------------------------------------
c     Calculate stress state
      sig_m=calSigMean(sig(1:6))
      sig_eq=calSigEq(sig(1:6))
      st=calSigTri(sig_m,sig_eq)
c     Calculate strain state
      peeq=calPeeq(Dp(1:6),hsv(208),dt1)
!============================================================
! Give constitutive variables to hsv
!------------------------------------------------------------
      call umatBccHamaUpdateHsv(numSys,sig_n1,f_n1,g_crss,
     1       r_n1,s11_n1,m11_n1,gamma_slip,gamma_n1,
     2       hsv,sig)
!============================================================
! Give non-constitutive variables to hsv
!------------------------------------------------------------
      call umatBccHamaUpdateVar(euler_n1,st,peeq,hsv) 
      endif
      endif
!============================================================
! End of cpfem
!------------------------------------------------------------
      end subroutine umatBccHama

      subroutine umatBccHamaGetMc(cm,ym,pr,bk,sm,
     1       umatType,dgamma0,mval,dgamma_lim,
     2       oriType,euler,
     3       hardType,g0,gs,h0,hs,q)
      implicit none
      double precision cm(*)
      double precision ym,pr,bk,sm
      double precision umatType,dgamma0,mval,dgamma_lim
      double precision oriType,euler(3)
      double precision hardType,g0,gs,h0,hs,q
c     cm(1 ~ 8)   Constitutive parameters
      ym=cm(1)
      pr=cm(2)
      bk=cm(3)
      sm=cm(4)       ! Shear modulus
c     cm(9 ~ 16) basic crystal plasticity model
      umatType=cm(9)
      dgamma0=cm(10)
      mval=cm(11)
      dgamma_lim=cm(12)
c     cm(17 ~ 24) orientation information
      oriType=cm(17)
      euler=cm(18:20)
c     cm(25 ~ 32) hardening 
      hardType=cm(25)
      g0=cm(26)
      gs=cm(27)
      h0=cm(28)
      hs=cm(29)
      q=cm(30)
      end subroutine umatBccHamaGetMc

      subroutine umatBccHamaInitHsv(g0,r,s11,m11,numSys,nhsv,hsv,sig)
      implicit none
      integer nhsv,numSys
      integer i,l,k
      double precision g0,hsv(nhsv),r(3,3),sig(6)
      double precision s11(3,numSys),m11(3,numSys)
c     Initialize the hsv list
            do l=1,nhsv
                  hsv(l)=0.
            enddo
c     Diagonal part of deformation gradient
            hsv(1)=1.
            hsv(5)=1.
            hsv(9)=1.
c     CRSS
            do l=1,numSys
                  hsv(9+l)=g0
            enddo
c     Transformation matrix for orientation
            hsv(22)=r(1,1)
            hsv(23)=r(2,1)
            hsv(24)=r(3,1)
            hsv(25)=r(1,2)
            hsv(26)=r(2,2)
            hsv(27)=r(3,2)
            hsv(28)=r(1,3)
            hsv(29)=r(2,3)
            hsv(30)=r(3,3)
c     Slip system vectors
            do l=1,numSys
                  k=(l-1)*3
                  hsv(31+k)=s11(1,l)
                  hsv(32+k)=s11(2,l)
                  hsv(33+k)=s11(3,l)
                  hsv(67+k)=m11(1,l)
                  hsv(68+k)=m11(2,l)
                  hsv(69+k)=m11(3,l)
            enddo
c     Cauchy stress tensor
            do i=1,6
                  sig(i)=0.
            enddo
      end subroutine umatBccHamaInitHsv

C     Hsv List
c     f(3,3)            <- hsv(1 ~ 9)
c     g_crss(12)        <- hsv(10 ~ 21)
c     r(3,3)            <- hsv(22 ~ 30)
c     s11(3,12)         <- hsv(31 ~ 66)
c     m11(3,12)         <- hsv(67 ~ 102)
c     gamma_slip(12)    <- hsv(103 ~ 114)
c     gamma_n1          <- hsv(115)
      subroutine umatBccHamaGetHsv(hsv,numSys,g_crss,r,s11,m11,
     1       gamma_n1,gamma_slip)
      implicit none
      integer numSys
      integer l,k
      double precision hsv(*)
      double precision g_crss(numSys),r(3,3)
      double precision s11(3,numSys),m11(3,numSys)
      double precision gamma_n1
      double precision gamma_slip(numSys)
c     Obtain CRSS from hsv
      do l=1,numSys
            g_crss(l)=hsv(9+l)
      enddo
c     Obtain slip volume from hsv
c     Obtain orientation info
      r(1,1)=hsv(22)
      r(2,1)=hsv(23)
      r(3,1)=hsv(24)
      r(1,2)=hsv(25)
      r(2,2)=hsv(26)
      r(3,2)=hsv(27)
      r(1,3)=hsv(28)
      r(2,3)=hsv(29)
      r(3,3)=hsv(30)
      do l=1,12
            k=(l-1)*3
            s11(1,l)=hsv(31+k)
            s11(2,l)=hsv(32+k)
            s11(3,l)=hsv(33+k)
            m11(1,l)=hsv(67+k)
            m11(2,l)=hsv(68+k)
            m11(3,l)=hsv(69+k)
      enddo
c     Obtain slip volume from hsv
      gamma_n1=hsv(115)
      do l=1,numSys
            gamma_slip(l)=hsv(102+l)
      enddo
      end subroutine umatBccHamaGetHsv
      
      subroutine umatBccHamaUpdateHsv(numSys,sig_n1,f_n1,g_crss,
     1       r_n1,s11_n1,m11_n1,gamma_slip,gamma_n1,
     2       hsv,sig)
      implicit none
      integer l,k
      integer numSys
      double precision hsv(*)
      double precision sig(6),sig_n1(6)
      double precision f_n1(3,3),g_crss(numSys),r_n1(3,3)
      double precision s11_n1(3,numSys),m11_n1(3,numSys)
      double precision gamma_slip(numSys),gamma_n1
      sig(1)=sig_n1(1)
      sig(2)=sig_n1(2)
      sig(3)=sig_n1(3)
      sig(4)=sig_n1(4)
      sig(5)=sig_n1(5)
      sig(6)=sig_n1(6)
c
      hsv(1)=f_n1(1,1)
      hsv(2)=f_n1(2,1)
      hsv(3)=f_n1(3,1)
      hsv(4)=f_n1(1,2)
      hsv(5)=f_n1(2,2)
      hsv(6)=f_n1(3,2)
      hsv(7)=f_n1(1,3)
      hsv(8)=f_n1(2,3)
      hsv(9)=f_n1(3,3)
c     
      do l=1,numSys
            hsv(9+l)=g_crss(l)
      enddo
c     
      hsv(22)=r_n1(1,1)
      hsv(23)=r_n1(2,1)
      hsv(24)=r_n1(3,1)
      hsv(25)=r_n1(1,2)
      hsv(26)=r_n1(2,2)
      hsv(27)=r_n1(3,2)
      hsv(28)=r_n1(1,3)
      hsv(29)=r_n1(2,3)
      hsv(30)=r_n1(3,3)
c     
      do l=1,numSys
            k=(l-1)*3
            hsv(31+k)=s11_n1(1,l)
            hsv(32+k)=s11_n1(2,l)
            hsv(33+k)=s11_n1(3,l)
            hsv(67+k)=m11_n1(1,l)
            hsv(68+k)=m11_n1(2,l)
            hsv(69+k)=m11_n1(3,l)
      enddo
c
      do l=1,numSys
            hsv(102+l)=gamma_slip(l)
      enddo
      hsv(115)=gamma_n1
      end subroutine umatBccHamaUpdateHsv

      subroutine umatBccHamaUpdateVar(euler_n1,st,peeq,hsv)
      implicit none
      double precision hsv(*)
      double precision euler_n1(3),st,peeq
      hsv(201)=euler_n1(1)
      hsv(202)=euler_n1(2)
      hsv(203)=euler_n1(3)
      hsv(204)=st
      hsv(205)=peeq
      end subroutine umatBccHamaUpdateVar