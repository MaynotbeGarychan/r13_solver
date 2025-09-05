#include "define.inc"
#include "define2.inc"
      subroutine umat_bcc(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1   temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
!============================================================
! Declaration of constitutive variables
!------------------------------------------------------------
! Note:
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
      integer umat_type
c     IO
      integer ori_type
c     Define the crystal
      integer,parameter:: num_ss=48
      integer,parameter:: cry_type=1
c     Intermedia variables
      integer i,j,k,l
      double precision, parameter :: pi=acos(-1.0d0)
      double precision idm(3,3)
c     Kinetic model variables
      double precision f(3,3),f_n1(3,3),df_n1(3,3),f_n1_inv(3,3)
      double precision f_det
      double precision fe(3,3),fe_n1(3,3)
      double precision fp(3,3),fp_inv(3,3)
      double precision L_n1(3,3)
      double precision D(3,3),Dv(6)
      double precision W(3,3),Wv(3)
      double precision We(3,3),Wev(3)
      double precision We1
      double precision exp_we(3,3)
c     Kinematic variables
      double precision s11(3,num_ss),m11(3,num_ss)
      double precision s11e(3,num_ss),m11e(3,num_ss)
      double precision s11_n1(3,num_ss),m11_n1(3,num_ss)
      double precision r(3,3),RL(6,6),r_n1(3,3)
      double precision Pa(3,3,num_ss),eschmid(6,num_ss)
      double precision Wa(3,3,num_ss),wschmid(3,num_ss)
      double precision euler(3)
      double precision phi1,lphi,phi2
      double precision Dp(6),Wp(3)
c     Slip system constitutive model variables
      double precision tau(num_ss),dgamma(num_ss)
      double precision tau_abs(num_ss)
      double precision gamma_slip(num_ss)
      double precision dgamma_tol,gamma_n1
      double precision g_crss(num_ss)
      double precision dgamma_0,mval,dgamma_lim
c     Stress Update model variables
      double precision ym,pr,bk,sm
      double precision L_ela(6,6),L_ela_cry(6,6)
      double precision dsig_0(6),sig_jau(6),sig_r(6),sig_n1(6)
c     Hardening model variable
      integer hard_type
      double precision g0,gs,h0,hs,q
c     Solid element variables
      double precision g,g2,gc,q1,q3,davg,p,deti,c22i,c23i,fac
      double precision temper,elsiz,epsp,capa,tt
!============================================================
! Declaration of non-constitutive variables
!------------------------------------------------------------
! Note:
!------------------------------------------------------------

!============================================================
! Declaration of utilized functions
!------------------------------------------------------------
! Note:
!------------------------------------------------------------
      integer array_maxidx,ss2sp_fcc
      double precision calSigEq,calSigMean,calSigTri,cal_lode
      double precision calPeeq
!============================================================
! Obtain variables from materials constants
!------------------------------------------------------------
! Note:
!------------------------------------------------------------
!------------------------------------------------------------
c     cm(1 ~ 8)   Constitutive parameters
      ym=cm(1)
      pr=cm(2)
      bk=cm(3)
      sm=cm(4)       ! Shear modulus
c     cm(9 ~ 16) basic crystal plasticity model
      umat_type=cm(9)
      dgamma_0=cm(10)
      mval=cm(11)
      dgamma_lim=cm(12)
c     cm(17 ~ 24) orientation information
      ori_type=cm(17)
      euler=cm(18:20)
c     cm(25 ~ 32) hardening 
      hard_type=cm(25)
      g0=cm(26)
      gs=cm(27)
      h0=cm(28)
      hs=cm(29)
      q=cm(30)
!============================================================
! Initialize Hsv list
!------------------------------------------------------------
! Note:
!------------------------------------------------------------     
      if (.not.failel) then
      if(ncycle==0) then
            call INITCRY(ori_type,cry_type,euler,
     1           num_ss,r,s11,m11)
c     Cauchy stress tensor
            do i=1,6
                  sig(i)=0.
            enddo
c     Initialize the hsv list
            do i=1,nhsv
                  hsv(i)=0.
            enddo
c     Diagonal part of deformation gradient
            hsv(1)=1.
            hsv(5)=1.
            hsv(9)=1.
c     CRSS
            do l=1,num_ss
                  hsv(9+l)=g0
            enddo
c     
            hsv(58)=r(1,1)
            hsv(59)=r(2,1)
            hsv(60)=r(3,1)
            hsv(61)=r(1,2)
            hsv(62)=r(2,2)
            hsv(63)=r(3,2)
            hsv(64)=r(1,3)
            hsv(65)=r(2,3)
            hsv(66)=r(3,3)
            do l=1,num_ss
                  k=(l-1)*3
                  hsv(67+k)=s11(1,l)
                  hsv(68+k)=s11(2,l)
                  hsv(69+k)=s11(3,l)
                  hsv(211+k)=m11(1,l)
                  hsv(212+k)=m11(2,l)
                  hsv(213+k)=m11(3,l)
            enddo
c     Cauchy stress tensor
            do i=1,6
                  sig(i)=0.
            enddo
c     elastic deformation gradient tensor
            ! hsv(404)=1.
            ! hsv(408)=1.
            ! hsv(412)=1.
      else
!============================================================
! Calculation begins
!------------------------------------------------------------
! Note:
!------------------------------------------------------------   
c     Obtain defromation gradient from hsv
      f(1,1)=hsv(1)
      f(2,1)=hsv(2)
      f(3,1)=hsv(3)
      f(1,2)=hsv(4)
      f(2,2)=hsv(5)
      f(3,2)=hsv(6)
      f(1,3)=hsv(7)
      f(2,3)=hsv(8)
      f(3,3)=hsv(9)
      f_n1(1,1)=hsv(nhsv+1)
      f_n1(2,1)=hsv(nhsv+2)
      f_n1(3,1)=hsv(nhsv+3)
      f_n1(1,2)=hsv(nhsv+4)
      f_n1(2,2)=hsv(nhsv+5)
      f_n1(3,2)=hsv(nhsv+6)
      f_n1(1,3)=hsv(nhsv+7)
      f_n1(2,3)=hsv(nhsv+8)
      f_n1(3,3)=hsv(nhsv+9)
c     Obtain CRSS from hsv
      do l=1,num_ss
            g_crss(l)=hsv(9+l)
      enddo
c     
      r(1,1)=hsv(58)
      r(2,1)=hsv(59)
      r(3,1)=hsv(60)
      r(1,2)=hsv(61)
      r(2,2)=hsv(62)
      r(3,2)=hsv(63)
      r(1,3)=hsv(64)
      r(2,3)=hsv(65)
      r(3,3)=hsv(66)
      do l=1,num_ss
            k=(l-1)*3
            s11(1,l)=hsv(67+k)
            s11(2,l)=hsv(68+k)
            s11(3,l)=hsv(69+k)
            m11(1,l)=hsv(211+k)
            m11(2,l)=hsv(212+k)
            m11(3,l)=hsv(213+k)
      enddo
c     Obtain slip volume from hsv
      gamma_n1=hsv(403)
      do l=1,num_ss
            gamma_slip(l)=hsv(354+l)
      enddo
c     elastic deformation gradient tensor
      ! fe(1,1)=hsv(404)
      ! fe(2,1)=hsv(405)
      ! fe(3,1)=hsv(406)
      ! fe(1,2)=hsv(407)
      ! fe(2,2)=hsv(408)
      ! fe(3,2)=hsv(409)
      ! fe(1,3)=hsv(410)
      ! fe(2,3)=hsv(411)
      ! fe(3,3)=hsv(412)

!============================================================
! Kinetic model
!------------------------------------------------------------
! Note:
!------------------------------------------------------------   
c     calculate deformation gradient rate
      call deformation_gradient_rate(f,f_n1,dt1,df_n1)
c     calculate the inverse of deformation gradient
      call matInverse(f_n1,3,f_n1_inv)
c     calculate velocity gradient
      call matInnProd(df_n1,f_n1_inv,3,3,3,L_n1)
c     decompose velocity gradient
      call velocity_gradient_decompose(L_n1,Dv,Wv)

!============================================================
! Constitutive model for slip at slip system
!------------------------------------------------------------
! Note:
!------------------------------------------------------------  
c      call ss_vec_to_configuration(fe,s11,m11,num_ss,
c     1   s11e,m11e)
      call calSchmidTensor(s11,m11,num_ss,
     1     Pa,Wa,eschmid,wschmid)
      call calSigRss(sig(1:6),eschmid,num_ss,tau)
      call calSlipRate(tau,g_crss,mval,dgamma_0,
     1   dgamma_lim,num_ss,dgamma)
      call updateCss(dgamma,num_ss,dt1,
     1     dgamma_tol,gamma_slip,gamma_n1)
c     Project the slip deformation into macro deformation and spin
      call calStrainRateBySlip(dgamma,eschmid,num_ss,Dp)
      call calSpinRateBySlip(dgamma,wschmid,num_ss,Wp)

!============================================================
! Update cauchy stress
!------------------------------------------------------------
! Note:
!------------------------------------------------------------ 
c     Fourth-order elastic tensor at crystal coordinate
      call calElasIsoTensor(ym,pr,L_ela_cry)
c     Transform elastic tensor to material coordinate
      call calTransMatFourthOrd(r,RL)
      call tranElasTensorLocal2Global(L_ela_cry,RL,L_ela)
c     Update the Cauchy stress tensor
      call mat33Det(f,f_det)
      call updateSigJaum(sig,eschmid,wschmid,dgamma,Wv,
     1     Dv,L_ela,f_det,num_ss,dt1,sig_n1)

!============================================================
! Calculate elastic deformation gradient
!------------------------------------------------------------
! Note:
!------------------------------------------------------------ 
c      call deformation_gradient_plastic(dgamma,s11,m11,
c     1     num_ss,fp)
c      call matInverse(fp,3,fp_inv)
c      call matInnProd(f_n1,fp_inv,3,3,3,fe_n1)

!============================================================
! Rotation model
!------------------------------------------------------------
! Note:
!------------------------------------------------------------ 
c     calculate the rotation rate for further rotation
      call cal_We_We1(Wv,Wp,We,We1)
      call matIdentity(3,idm)
      if(We1==0.) then
            do j=1,3
                  do i=1,3
                        exp_we(i,j)=idm(i,j)
                  enddo
            enddo
      else
            do j=1,3
                  do i=1,3
                  exp_we(i,j)=idm(i,j)+(sin(We1*dt1)/We1)*We(i,j)
                  do k=1,3
                        exp_we(i,j)=exp_we(i,j)
     1                     +((1-cos(We1*dt1))/(We1**2))
     2                     *We(i,k)*We(k,j)
                  enddo
                  enddo
            enddo
      endif
c     Rotate the tranformation matrix
      call rot_tran_mat(r,exp_we,r_n1)
c     Rotate the slip systsme vectors
      call rot_slip_vec(s11,m11,exp_we,num_ss,
     1     s11_n1,m11_n1)
c     Extract the euler angle from the tranformation matrix
      call calEulerbyTransMat(r_n1,pi,
     1     phi1,lphi,phi2)

!============================================================
! Hardening model
!------------------------------------------------------------
! Note:
!------------------------------------------------------------ 
      call hardening_bcc(g0,gs,h0,hs,q,gamma_n1,dgamma,
     1           dt1,g_crss)
!============================================================
! Calculation of non-constitutive variables
!------------------------------------------------------------
! Note:
!------------------------------------------------------------ 
c     Calculate stress state

!============================================================
! Give constitutive variables to hsv
!------------------------------------------------------------
! Note:
!------------------------------------------------------------ 

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
      do l=1,num_ss
            hsv(9+l)=g_crss(l)
      enddo
c     
      hsv(58)=r_n1(1,1)
      hsv(59)=r_n1(2,1)
      hsv(60)=r_n1(3,1)
      hsv(61)=r_n1(1,2)
      hsv(62)=r_n1(2,2)
      hsv(63)=r_n1(3,2)
      hsv(64)=r_n1(1,3)
      hsv(65)=r_n1(2,3)
      hsv(66)=r_n1(3,3)
c     
      do l=1,num_ss
            k=(l-1)*3
            hsv(67+k)=s11_n1(1,l)
            hsv(68+k)=s11_n1(2,l)
            hsv(69+k)=s11_n1(3,l)
            hsv(211+k)=m11_n1(1,l)
            hsv(212+k)=m11_n1(2,l)
            hsv(213+k)=m11_n1(3,l)
      enddo
c
      do l=1,num_ss
            hsv(354+l)=gamma_slip(l)
      enddo
      hsv(403)=gamma_n1
c     
      ! hsv(404)=fe_n1(1,1)
      ! hsv(405)=fe_n1(2,1)
      ! hsv(406)=fe_n1(3,1)
      ! hsv(407)=fe_n1(1,2)
      ! hsv(408)=fe_n1(2,2)
      ! hsv(409)=fe_n1(3,2)
      ! hsv(410)=fe_n1(1,3)
      ! hsv(411)=fe_n1(2,3)
      ! hsv(412)=fe_n1(3,3)
!============================================================
! Give non-constitutive variables to hsv
!------------------------------------------------------------
! Note:
!------------------------------------------------------------ 
c     phi1,lphi,phi2    <- hsv(201 ~ 203)
      hsv(501)=phi1
      hsv(502)=lphi
      hsv(503)=phi2
      endif
      endif
!============================================================
! End of cpfem
!------------------------------------------------------------
! Note:
!------------------------------------------------------------ 
      end subroutine umat_bcc