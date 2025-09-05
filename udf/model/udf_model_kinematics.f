      subroutine calDeforSpinRateByDispGrad(f,f_n1,dt1,Dv,Wv)
      !============================================================
      ! Main program to calculate the deformation and spin rate
      ! tensor from displacement gradient
      !------------------------------------------------------------
      ! input: 
      ! f(3,3) - deformation gradient tensor at current step
      ! f_n1(3,3) - deformation gradient tensor at next step
      ! dt1       - time step
      ! output: 
      ! Dv(6)     - strain rate tensor in Vogit notation
      ! Wv(3)     - spin rate tensor in Vogit notation
      !============================================================
      implicit none
      double precision f(3,3),f_n1(3,3),df_n1(3,3),f_n1_inv(3,3)
      double precision L_n1(3,3),Dv(6),Wv(6)
      double precision dt1
c     calculate the deformation gradieent rate from disp grad
      call deformation_gradient_rate(f,f_n1,dt1,df_n1)
c     calculate the inverse of deformation gradient
      call matInverse(f_n1,3,f_n1_inv)
c     calculate velocity gradient
      call matInnProd(df_n1,f_n1_inv,3,3,3,L_n1)
c     decompose velocity gradient
      call velocity_gradient_decompose(L_n1,Dv,Wv)
      end subroutine calDeforSpinRateByDispGrad

      subroutine getDispGradfromHsv(hsv,nhsv,f,f_n1)
      implicit none
      integer nhsv
      double precision f(3,3),f_n1(3,3)
      double precision hsv(*)
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
      end subroutine getDispGradfromHsv
      
      subroutine deformation_gradient_rate(f,f_n1,dt1,df_n1)
      !============================================================
      ! calculate the deformation gradient
      !------------------------------------------------------------
      ! input: 
      ! f(3,3) - deformation gradient tensor at current step
      ! f_n1(3,3) - deformation gradient tensor at next step
      ! dt1       - time step
      ! output: 
      ! df_n1(3,3)  - rate deformation gradient tensor at curr step
      !============================================================
      implicit none
      double precision f(3,3),f_n1(3,3)
      double precision dt1
      double precision df_n1(3,3)
      integer i,j

      do i=1,3
            do j=1,3
                  df_n1(i,j)=(f_n1(i,j)-f(i,j))/dt1
            enddo
      enddo

      end subroutine deformation_gradient_rate
      
      subroutine velocity_gradient_decompose(l_n1,Dv,Wv)
      !============================================================
      ! decompose the velocity gradient into D,W
      !------------------------------------------------------------
      ! input: 
      ! l_n1(3,3) - velocity gradient tensor
      ! output: 
      ! Dv(6)     - strain rate tensor in Vogit notation
      ! Wv(3)     - spin rate tensor in Vogit notation
      !============================================================
      implicit none
      double precision l_n1(3,3)
      double precision Dv(6),Wv(3)
      integer i

      do i=1,3
            Dv(i)=l_n1(i,i)
      enddo
      Dv(4)=l_n1(1,2)+l_n1(2,1)
      Dv(5)=l_n1(2,3)+l_n1(3,2)
      Dv(6)=l_n1(3,1)+l_n1(1,3)

      Wv(1)=(l_n1(1,2)-l_n1(2,1))/2.
      Wv(2)=(l_n1(2,3)-l_n1(3,2))/2.
      Wv(3)=(l_n1(3,1)-l_n1(1,3))/2.

      end subroutine velocity_gradient_decompose

      subroutine deformation_gradient_plastic(dgamma,s11,m11,num_ss,
     1     fp)
      !============================================================
      ! calculate the deformation gradient plastic part
      !------------------------------------------------------------
      ! input: 
      ! dgamma(num_ss) - slip rate at current step
      ! s11(3,num_ss)  - slip system vectors
      ! m11(3,num_ss)  - slip plane normal vectors
      ! num_ss
      ! output: 
      ! fp(3,3)     - deformation gradient plastic part
      !============================================================
      implicit none
      integer l,i,j
      integer num_ss
      double precision dgamma(num_ss)
      double precision s11(3,num_ss)
      double precision m11(3,num_ss)
      double precision fp(3,3)

      call matIdentity(3,fp)
      do l=1,num_ss
            do i=1,3
                  do j=1,3
                        fp(i,j)=fp(i,j)+s11(i,l)*m11(j,l)*dgamma(l)
                  enddo
            enddo
      enddo

      end subroutine deformation_gradient_plastic

      subroutine calStrainRateBySlip(dgamma,eschmid,num_ss,
     1     Dpv)
      !============================================================
      ! Calculate strain tensor, plastic part
      !------------------------------------------------------------
      ! input: 
      ! dgamma(num_ss)      - slip rate at current step
      ! eschmid(6,num_ss)   - schmid tensor in vogit for each sys
      ! num_ss              - number of slip system
      ! output: 
      ! Dpv(6)              - strain tensor in vogit
      !============================================================
      implicit none
      integer l,i
      integer num_ss
      double precision dgamma(num_ss)
      double precision eschmid(6,num_ss)
      double precision Dpv(6)

      do i=1,6
            Dpv(i)=0.
            do l=1,num_ss
                  Dpv(i)=Dpv(i)+eschmid(i,l)*dgamma(l)
            enddo
      enddo

      end subroutine calStrainRateBySlip

      subroutine calSpinRateBySlip(dgamma,wschmid,num_ss,
     1     Wpv)
      !============================================================
      ! Calculate spin tensor, plastic part
      !------------------------------------------------------------
      ! input: 
      ! dgamma(num_ss)      - slip rate at current step
      ! eschmid(6,num_ss)   - schmid tensor in vogit for each sys
      ! num_ss              - number of slip system
      ! output: 
      ! Dpv(6)              - strain tensor in vogit
      !============================================================
      implicit none
      integer l,i
      double precision dgamma(num_ss)
      double precision wschmid(3,num_ss)
      integer num_ss
      double precision Wpv(3)

      do i=1,3
            Wpv(i)=0.
            do l=1,num_ss
                  Wpv(i)=Wpv(i)+wschmid(i,l)*dgamma(l)
            enddo
      enddo

      end subroutine calSpinRateBySlip