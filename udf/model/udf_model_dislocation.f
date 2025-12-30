      subroutine updateDsl(rho_r,numSys,dt1,rho)
      !============================================================
      ! Update of the dislocation density
      !------------------------------------------------------------
      ! input:
      ! rho_r(numSys)   - Dislocation glide rate
      ! numSys          - Number of slip system
      ! dt1             - Time step
      ! output:
      ! rho(numSys)     - Dislocation density of slip system
      !------------------------------------------------------------
      implicit none
      include 'define_cp.inc'
      integer l
      integer numSys
      double precision dt1
      double precision rho_r(maxSys)
      double precision rho(maxSys)
      double precision val
      do l=1,numSys
          val= rho_r(l)*dt1
          rho(l)=rho(l)+val
          if (rho(l) .lt. 0.d0) then
              rho(l)=0.d0 ! safeguard
              print *, 'Warning: negative dislocation density'
              print *, 'corrected to zero on slip system.'
          endif
      enddo
      end subroutine updateDsl
      
      subroutine calDslLineVec(m11,s11,numSys,l11)
      !============================================================
      ! calculate the dislocation line vector
      !------------------------------------------------------------
      ! input: 
      ! m11(3,numSys)   - slip direction vector
      ! s11(3,numSys)   - slip plane vector
      ! numSys      - number of slip system
      ! output:
      ! l11(3,numSys)   - dislocation line vector
      !============================================================
        implicit none
        include 'define_cp.inc'
        integer numSys
        integer l
        double precision m11(3,maxSys),s11(3,maxSys)
        double precision l11(3,maxSys)
        do l=1,numSys
            call vector_cross_product(m11(:,l),
     1   s11(:,l),l11(:,l))
        enddo
      end subroutine calDslLineVec

      subroutine calDslEvolKocks(ka,kb,bv,dgamma,numSys,dt1,
     1       rho,rho_r,rho_tol)
        !============================================================
        ! Evolution model for dislocation
        ! Laws for work-hardening and low-temperature creep.
        !  J. Eng. Mater. Tech., ASMEH 98, 76., Kocks
        !------------------------------------------------------------
        ! input: 
        ! ka, kb          - Fitting constants
        ! dgamma(numSys)  - Slip rate
        ! numSys          - Number of slip system
        ! rho(numSys)     - Dislocation density of slip system
        ! output:
        ! rho(numSys)     - Dislocation density of slip system
        ! rho_r(numSys)   - Rate ofDislocation density of slip system
        ! rho_tol         - Total dislocation (for checking)
        !------------------------------------------------------------
        implicit none
        integer l
        integer numSys
        double precision rho_tol
        double precision rho(numSys)
        double precision rho_r(numSys)
        double precision dgamma(numSys)
        double precision ka,kb,bv
        double precision dt1
        double precision val
        rho_tol=0.
        do l=1,numSys
              rho_tol=rho_tol+rho(l)
        enddo
        do l=1,numSys
              val=sqrt(rho_tol)/ka-kb*rho(l)
              rho_r(l)=abs(dgamma(l))/bv*val
              rho(l)=rho(l)+rho_r(l)*dt1
        enddo
        end subroutine calDslEvolKocks

        subroutine calDslRateHama(yc,km,rho,typeCry,numSys,
     1      dgamma,rho_r)
      !============================================================
      ! Calculate the dislocation glide rate
      ! for Hama sensei's model
      !------------------------------------------------------------
      ! input:
      ! dsl_L           - Dislocation mean free path
      ! yc              - Annihilation coefficient
      ! rho(numSys)     - Dislocation density of slip system
      ! dgamma(numSys)  - Slip rate at current step
      ! numSys          - Number of slip system
      ! output:
      ! rho_r(numSys)   - Dislocation glide rate
      !---------------------------------------------------------
      implicit none
      include 'define_cp.inc'
      integer l,j
      integer numSys,typeCry
      double precision yc,km
      double precision rho(maxSys)
      double precision dgamma(maxSys)
      double precision rho_r(maxSys)
      double precision dsl_L(maxSys)
      double precision val
      double precision matInteract(maxSys,maxSys)

      call getDslInteractionMatrix(typeCry,numSys,matInteract)

      do l=1,numSys
          val=1.0d0
          do j=1,numSys
            if (j .ne. l) then
              val=val+matInteract(l,j)*rho(j)
            endif
          enddo
          dsl_L(l)=1/sqrt(val)*km
      enddo

      do l=1,numSys
        val=(1/dsl_L(l))-2*yc*rho(l)
        rho_r(l)=(1/CST_BV)*val*abs(dgamma(l))
      enddo

      end subroutine calDslRateHama

      subroutine getDslInteractionMatrix(typeCry,numSys,matInteract)
      !============================================================
      ! Get dislocation (latent) hardening interaction matrix
      ! For BCC {110}<111> slip systems (12 systems)
      !
      ! matInteract(i,j) = h( ih(i,j) )
      !
      !============================================================
      implicit none
      include 'define_cp.inc'
      integer typeCry
      integer numSys
      double precision matInteract(maxSys,maxSys)
      integer i, j
      double precision h(17)
      integer ih(12,12)
!============================================================
! 1. Latent hardening parameters h(k)
!    (from the table you provided)
!============================================================
      data h / 
     & 0.1d0,    ! h(1)
     & 0.1d0,    ! h(2)
     & 0.45d0,   ! h(3)
     & 5.5d0,    ! h(4)
     & 0.4d0,    ! h(5)
     & 0.6d0,    ! h(6)
     & 1.5d0,    ! h(7)
     & 0.1d0,    ! h(8)
     & 0.01d0,   ! h(9)
     & 3.5d0,    ! h(10)
     & 1.0d0,    ! h(11)
     & 1.0d0,    ! h(12)
     & 0.1d0,    ! h(13)
     & 0.08d0,   ! h(14)
     & 0.05d0,   ! h(15)
     & 0.1d0,    ! h(16)
     & 1.0d0     ! h(17)
     & /

!============================================================
! 2. Interaction index matrix ih(i,j)
!    Corresponds exactly to the {110}<111> table
!============================================================
            data ih /
!        1  2  3  4  5  6  7  8  9 10 11 12
     &  1, 2, 3, 3, 5, 4, 5, 6, 6, 5, 5, 4,   ! 1
     &  2, 1, 3, 3, 6, 5, 4, 5, 5, 4, 6, 5,   ! 2
     &  3, 3, 1, 2, 4, 5, 6, 5, 4, 5, 5, 6,   ! 3
     &  3, 3, 2, 1, 5, 6, 5, 4, 5, 6, 4, 5,   ! 4
     &  5, 6, 4, 5, 1, 2, 3, 3, 4, 5, 6, 5,   ! 5
     &  4, 5, 5, 6, 2, 1, 3, 3, 5, 6, 5, 4,   ! 6
     &  5, 4, 6, 5, 3, 3, 1, 2, 5, 4, 5, 6,   ! 7
     &  6, 5, 5, 4, 3, 3, 2, 1, 6, 5, 4, 5,   ! 8
     &  6, 5, 4, 5, 4, 5, 5, 6, 1, 2, 3, 3,   ! 9
     &  5, 4, 5, 6, 5, 6, 4, 5, 2, 1, 3, 3,   ! 10
     &  5, 6, 5, 4, 6, 5, 5, 4, 3, 3, 1, 2,   ! 11
     &  4, 5, 6, 5, 5, 4, 6, 5, 3, 3, 2, 1    ! 12
     & /

!============================================================
! 3. Build interaction matrix
!============================================================
      if (numSys .ne. 12) then
          print *, 'Error: BCC {110}<111> requires numSys = 12'
          stop
      endif

      do i = 1, numSys
          do j = 1, numSys
              matInteract(i,j) = h( ih(i,j) )
          end do
      end do

      return
      end subroutine getDslInteractionMatrix

      subroutine calDslRcvyRateKohenert(rho, temp, mu,
     &     kappa1, kappa2, numSys, drho)
      !============================================================
      ! Calculate the dislocation recovery rate, Kohenert
      !------------------------------------------------------------
      ! input:
      ! rho        - dislocation density
      ! T          - temperature
      ! mu         - shear modulus
      ! kappa1     - material parameter κ1
      ! kappa2     - material parameter κ2
      ! numSys     - number of slip systems
      ! output:
      ! drho       - recovery rate dρ/dt
      !------------------------------------------------------------
      implicit none
      include 'define_cp.inc'
      double precision rho(maxSys)   ! dislocation density
      double precision temp          ! temperature
      double precision mu         ! shear modulus
      double precision diffCoef   ! diffusion coefficient
      double precision kappa1     ! material parameter κ1
      double precision kappa2     ! material parameter κ2
      double precision drho(maxSys)    ! recovery rate dρ/dt
      double precision expo, arg
      integer numSys
      integer l
!------------------------------------------------------------
! Safety check
!------------------------------------------------------------
      call calDslDiffCoeffArrhen(temp,diffCoef)
      do l=1,numSys
      if (rho(l) .le. 0.d0 .or. temp .le. 0.d0) then
         drho(l) = 0.d0
         continue
      end if
!------------------------------------------------------------
! Exponential argument
!------------------------------------------------------------
      arg = kappa2*mu*(CST_BV**4)*sqrt(rho(l))/(CST_BZ*temp)
! Optional: avoid overflow
      arg = min(arg, 80.d0)
      expo = exp(arg) - 1.d0
      drho(l)=(kappa1*diffCoef/CST_BV)*(rho(l)**1.5d0)*expo
      enddo
      end subroutine calDslRcvyRateKohenert

      subroutine calDslDiffCoeffArrhen(temp,diffCoef)
      !============================================================
      ! Calculate the diffusion coefficient
      ! Arrhenius type
      !------------------------------------------------------------
      ! input:
      ! temp          - Temperature
      ! output:
      ! diffCoef      - Diffusion coefficient
      !------------------------------------------------------------
        implicit none
        include 'define_cp.inc'
        double precision temp
        double precision diffCoef
        diffCoef = CST_DIFF0*exp(-CST_QFeCr/(CST_RGAS*temp))
      end subroutine calDslDiffCoeffArrhen


