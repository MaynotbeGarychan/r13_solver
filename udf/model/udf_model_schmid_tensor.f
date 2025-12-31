      subroutine calSchmidTensor(s11e,m11e,numSys,
     1     sfmate,wfmate,sfve,wfve)
      !============================================================
      ! Calculate Schmid tensor its Vogit notation
      ! based on slip systems vectors in reference scheme
      !------------------------------------------------------------
      ! input: 
      ! m11e(3,numSys),s11e(3,numSys) - local slip system vectors
      ! numSys        - num of slip systems
      ! output: 
      ! sfmate(3,3,numSys),wfmate(3,3,numSys) - Schmid tensor
      ! sfve(6,numSys),wfve(3,numSys) - Schmid tensor in Vog
      !============================================================
      implicit none
      include 'define_cp.inc'
      integer l,i,j
      integer numSys
      double precision m11e(3,maxSys),s11e(3,maxSys)
      double precision sfmate(3,3,maxSys),wfmate(3,3,maxSys)
      double precision sfve(6,maxSys),wfve(3,maxSys)
      do l=1,numSys
            do i=1,3
                  do j=1,3
                      sfmate(i,j,l)=(s11e(i,l)*m11e(j,l)+
     1                        s11e(j,l)*m11e(i,l))/2.
                      wfmate(i,j,l)=(s11e(i,l)*m11e(j,l)-
     1                        s11e(j,l)*m11e(i,l))/2.
                  enddo
              enddo
            sfve(1,l)=sfmate(1,1,l)
            sfve(2,l)=sfmate(2,2,l)
            sfve(3,l)=sfmate(3,3,l)
            sfve(4,l)=sfmate(1,2,l)*2.
            sfve(5,l)=sfmate(2,3,l)*2.
            sfve(6,l)=sfmate(3,1,l)*2.
            wfve(1,l)=wfmate(1,2,l)
            wfve(2,l)=wfmate(2,3,l)
            wfve(3,l)=wfmate(3,1,l)
      enddo

      end subroutine calSchmidTensor

      subroutine calNonSchmidTensor(s11e, m11e, numSys, nsfve)
      !============================================================
      ! Calculate Non-Schmid tensor η for each slip system
      ! based on slip systems vectors in reference scheme
      !------------------------------------------------------------
      ! input:
      ! m11e(3,numSys),s11e(3,numSys) - local slip system vectors
      ! numSys        - num of slip systems
      ! output:
      ! nschmid(3,3,numSys) - Non-Schmid tensor for each slip sys
      !============================================================
      implicit none
      include 'define_cp.inc'
      integer numSys
      integer i, j, l
      double precision s11e(3,maxSys), m11e(3,maxSys)
      double precision nsfmate(3,3,maxSys) 
      double precision nsfve(6,maxSys)
      double precision z(3)
! =================== build tensor for each slip system ===================
      do l = 1, numSys

! ---- compute z = s x m
         z(1) = s11e(2,l)*m11e(3,l) - s11e(3,l)*m11e(2,l)
         z(2) = s11e(3,l)*m11e(1,l) - s11e(1,l)*m11e(3,l)
         z(3) = s11e(1,l)*m11e(2,l) - s11e(2,l)*m11e(1,l)

! ---- assemble η tensor
         do i = 1, 3
            do j = 1, 3
               nsfmate(i,j,l) =
     &           CST_NSFSS * s11e(i,l)*s11e(j,l)
     &         + CST_NSFMM * m11e(i,l)*m11e(j,l)
     &         + CST_NSFZZ * z(i)*z(j)
     &         + CST_NSFSZ * ( s11e(i,l)*z(j) + z(i)*s11e(j,l) )
     &         + CST_NSFMZ * ( m11e(i,l)*z(j) + z(i)*m11e(j,l) )
            end do
         end do
! ---- convert to Vogit notation
         nsfve(1,l) = nsfmate(1,1,l)
         nsfve(2,l) = nsfmate(2,2,l)
         nsfve(3,l) = nsfmate(3,3,l)
         nsfve(4,l) = nsfmate(1,2,l)*2.0d0
         nsfve(5,l) = nsfmate(2,3,l)*2.0d0
         nsfve(6,l) = nsfmate(3,1,l)*2.0d0
      end do

      end subroutine calNonSchmidTensor


      SUBROUTINE CALSF(SLPDIR,SLPNOR,NUMSS,SLPDEF)
      !============================================================
      ! Calculate Schmid factor (Vogit)
      ! based on local slip systems vectors
      !------------------------------------------------------------
      ! input: 
      ! NUMSS - number of slip system
      ! SLPDIR(3,NUMSS),SLPNOR(3,NUMSS)-glob slip direciton, normal
      ! output: 
      ! SLPDEF(6,NUMSS)  - Vogit notation of schmid factor
      !------------------------------------------------------------
      INTEGER NUMSS
      INTEGER J
      DOUBLE PRECISION SLPDIR(3,NUMSS),SLPNOR(3,NUMSS)
      DOUBLE PRECISION SLPDEF(6,NUMSS)
      
      DO J=1,NUMSS
            SLPDEF(1,J)=SLPDIR(1,J)*SLPNOR(1,J)
            SLPDEF(2,J)=SLPDIR(2,J)*SLPNOR(2,J)
            SLPDEF(3,J)=SLPDIR(3,J)*SLPNOR(3,J)
            SLPDEF(4,J)=SLPDIR(1,J)*SLPNOR(2,J)+SLPDIR(2,J)*SLPNOR(1,J)
            SLPDEF(5,J)=SLPDIR(1,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(1,J)
            SLPDEF(6,J)=SLPDIR(2,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(2,J)
      END DO

      END SUBROUTINE CALSF

      SUBROUTINE CALSFSPIN(SLPDIR,SLPNOR,NUMSS,SLPSPN)
      !============================================================
      ! Calculate Schmid factor for spin (Vogit)
      ! based on local slip systems vectors
      !------------------------------------------------------------
      ! input: 
      ! NUMSS - number of slip system
      ! SLPDIR(3,NUMSS),SLPNOR(3,NUMSS)-glob slip direciton, normal
      ! output: 
      ! SLPSPN(6,NUMSS)  - Vogit notation of schmid factor
      !------------------------------------------------------------

      INTEGER NUMSS
      INTEGER J
      DOUBLE PRECISION SLPDIR(3,NUMSS),SLPNOR(3,NUMSS)
      DOUBLE PRECISION SLPSPN(3,NUMSS)

      DO J=1,NUMSS
            SLPSPN(1,J)=0.5*(SLPDIR(1,J)*SLPNOR(2,J)-
     2                       SLPDIR(2,J)*SLPNOR(1,J))
            SLPSPN(2,J)=0.5*(SLPDIR(3,J)*SLPNOR(1,J)-
     2                       SLPDIR(1,J)*SLPNOR(3,J))
            SLPSPN(3,J)=0.5*(SLPDIR(2,J)*SLPNOR(3,J)-
     2                       SLPDIR(3,J)*SLPNOR(2,J))
            END DO

      END SUBROUTINE CALSFSPIN







