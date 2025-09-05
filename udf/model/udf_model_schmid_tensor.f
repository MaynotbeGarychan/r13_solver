      subroutine calSchmidTensor(s11e,m11e,num_ss,
     1     sfmate,wfmate,sfve,wfve)
      !============================================================
      ! Calculate Schmid tensor its Vogit notation
      ! based on slip systems vectors in reference scheme
      !------------------------------------------------------------
      ! input: 
      ! m11e(3,num_ss),s11e(3,num_ss) - local slip system vectors
      ! num_ss        - num of slip systems
      ! output: 
      ! sfmate(3,3,num_ss),wfmate(3,3,num_ss) - Schmid tensor
      ! sfve(6,num_ss),wfve(3,num_ss) - Schmid tensor in Vog
      !============================================================
      implicit none
      integer l,i,j
      integer num_ss
      double precision m11e(3,num_ss),s11e(3,num_ss)
      double precision sfmate(3,3,num_ss),wfmate(3,3,num_ss)
      double precision sfve(6,num_ss),wfve(3,num_ss)
      do l=1,num_ss
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

      subroutine normal_tensor(m11,num_sp,num_ss,nfmat,nfv)
      !============================================================
      ! Calculate non-Schmid tensor Na its vogit 
      ! based on local slip systems vectors
      !------------------------------------------------------------
      ! input: 
      ! m11(3,num_ss) - local slip plane vectors
      ! num_ss        - num of slip systems
      ! num_sp        - num of slip plane
      ! output: 
      ! nfmat(3,3,num_sp) - Nomral tensor
      ! nfv(6,num_sp) - Normal tensor based on vogit
      !------------------------------------------------------------
      implicit none
      integer i,j,l,k
      integer num_ss,num_sp
      double precision m11(3,num_ss)
      double precision nfmat(3,3,num_sp)
      double precision nfv(6,num_sp)

      do l=1,num_sp
            do i=1,3
               do j=1,3
                     k=3*l
                     nfmat(i,j,l)=m11(i,k)*m11(j,k)
               enddo
            enddo

            nfv(1,l)=nfmat(1,1,l)
            nfv(2,l)=nfmat(2,2,l)
            nfv(3,l)=nfmat(3,3,l)
            nfv(4,l)=nfmat(1,2,l)*2.
            nfv(5,l)=nfmat(2,3,l)*2.
            nfv(6,l)=nfmat(3,1,l)*2.
      enddo

      end subroutine normal_tensor

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







