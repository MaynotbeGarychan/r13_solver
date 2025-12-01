      function ss2sp_fcc(ss) result(sp)
        !============================================================
        ! Return the index of corresponding slip plane for a specific
        ! slip system
        !------------------------------------------------------------
        ! input: 
        ! ss       - index of the slip system
        ! output:
        ! sp       - index of the slip plane
        !------------------------------------------------------------
        implicit none
        integer ss
        integer sp

        sp=(ss-1)/3+1

      end function ss2sp_fcc

      subroutine ROTMATBYEULER(phi1,lphi,phi2,r)
        !============================================================
        ! Calculate the tranformation matrix 
        ! to tranform the slip system vectors to local based on 
        ! crystal orientations (euler angles)
        !------------------------------------------------------------
        ! input: 
        ! phi1,lphi,phi2  - euler angle of the crystal (element)
        ! output:
        ! r(3,3)          - tranformation matrix
        !------------------------------------------------------------
        
        implicit none
        double precision phi1,lphi,phi2
        double precision c1,s1,c2,s2,c3,s3
        double precision r(3,3)

        c1=cos(phi1)
        s1=sin(phi1)
        c2=cos(lphi)
        s2=sin(lphi)
        c3=cos(phi2)
        s3=sin(phi2)
        
        r(1,1)=c1*c3
     1        -s1*c2*s3
        r(1,2)=-c1*s3
     1        -s1*c2*c3
        r(1,3)=s1*s2
        r(2,1)=s1*c3
     1         +c1*c2*s3
        r(2,2)=-s1*s3
     1        +c1*c2*c3
        r(2,3)=-c1*s2
        r(3,1)=s2*s3
        r(3,2)=s2*c3
        r(3,3)=c2

      end subroutine ROTMATBYEULER

      subroutine calEulerbyTransMat(r,euler)
         !============================================================
         ! Calculate the tranformation matrix 
         ! to tranform the slip system vectors to local based on 
         ! crystal orientations (euler angles)
         !------------------------------------------------------------
         ! input: 
         ! r(3,3)          - tranformation matrix
         ! pi              - constant of pi
         ! output:
         ! euler(3)  - euler angle of the crystal (element
         !------------------------------------------------------------
         implicit none
         double precision r(3,3)
         double precision euler(3)
         double precision, parameter :: pi=acos(-1.0d0)

         euler(2)=acos(r(3,3))
         if(euler(2)==0.) then
            euler(1)=atan2(-r(1,2),r(1,1))
            euler(2)=0.
         else if(euler(2)==pi) then
            euler(1)=atan2(r(1,2),r(1,1))
            euler(2)=0.
         else
            euler(1)=atan2(r(1,3),-r(2,3))
            euler(2)=atan2(r(3,1),r(3,2))
         endif
         if(euler(1) < 0.) then
            euler(1)=euler(1)+2.*pi
         endif
         if(euler(2) < 0.) then
            euler(2)=euler(2)+2.*pi
         endif
      end subroutine calEulerbyTransMat

      subroutine calTransMatFourthOrd(r,RL)
        !============================================================
        ! Calculate the tranformation matrix for four-order
        ! tensor
        !------------------------------------------------------------
        ! input: 
        ! r(3,3)     - tranformation matrix
        ! output: 
        ! RL(6,6)    - tranformation matrix for four-order
        !------------------------------------------------------------

         implicit none
         integer i,j
         double precision r(3,3)
         double precision RL(6,6)

         do i=1,3
            do j=1,3
               RL(i,j)=r(i,j)**2
            enddo
               RL(i,4)=2.*r(i,1)*r(i,2)
               RL(i,5)=2.*r(i,2)*r(i,3)
               RL(i,6)=2.*r(i,3)*r(i,1)
         enddo
         do j=1,3
            RL(4,j)=r(1,j)*r(2,j)
            RL(5,j)=r(2,j)*r(3,j)
            RL(6,j)=r(3,j)*r(1,j)
         enddo
         RL(4,4)=r(1,1)*r(2,2)+r(2,1)*r(1,2)
         RL(4,5)=r(1,2)*r(2,3)+r(2,2)*r(1,3)
         RL(4,6)=r(1,1)*r(2,3)+r(2,1)*r(1,3)
         RL(5,4)=r(2,1)*r(3,2)+r(3,1)*r(2,2)
         RL(5,5)=r(2,2)*r(3,3)+r(3,2)*r(2,3)
         RL(5,6)=r(2,1)*r(3,3)+r(3,1)*r(2,3)
         RL(6,4)=r(1,1)*r(3,2)+r(3,1)*r(1,2)
         RL(6,5)=r(1,2)*r(3,3)+r(3,2)*r(1,3)
         RL(6,6)=r(1,1)*r(3,3)+r(3,1)*r(1,3)

      end subroutine calTransMatFourthOrd

      subroutine trans_fourth_order_tensor(mat,rl,n,mat_trans)
      !============================================================
      ! transform the fourth order tensor
      !  mat_trans = rl * mat * rlt
      !------------------------------------------------------------
      ! input: 
      ! mat(n,n)   - initial matrix
      ! rl(n,n)    - transformation matrix
      ! n          - num of col and row
      ! output: 
      ! mat_trans(n,n)  - matrix after transformation
      !------------------------------------------------------------
         implicit none
         integer n
         double precision mat(n,n)
         double precision rl(n,n),rlt(n,n)
         double precision mat_trans(n,n)
         double precision temp(n,n)

         call matInnProd(rl,mat,n,n,n,temp)
         call matTranspose(rl,n,n,rlt)
         call matInnProd(temp,rlt,n,n,n,mat_trans)

      end subroutine trans_fourth_order_tensor

      subroutine trans_matrix_by_axis_angle(axis,angle,r)
      !============================================================
      ! Calculate the tranformation matrix 
      ! for a rotation repsect to an axis for a angle
      !------------------------------------------------------------
      ! input: 
      ! axis(3)  - rotation axis
      ! angle    - total angle for rotation (in radians)
      ! output:
      ! r(3,3)          - tranformation matrix
      !------------------------------------------------------------
         implicit none
         double precision axis(3)
         double precision angle
         double precision r(3,3)
         double precision c,s
         double precision ux,uy,uz

         call vector_normalize(axis,3,axis)
         ux = axis(1)
         uy = axis(2)
         uz = axis(3)
         c=cos(angle)
         s=sin(angle)

         r(1,1) = c + ux*ux*(1.0d0 - c)
         r(1,2) = ux*uy*(1.0d0 - c) - uz*s
         r(1,3) = ux*uz*(1.0d0 - c) + uy*s

         r(2,1) = uy*ux*(1.0d0 - c) + uz*s
         r(2,2) = c + uy*uy*(1.0d0 - c)
         r(2,3) = uy*uz*(1.0d0 - c) - ux*s

         r(3,1) = uz*ux*(1.0d0 - c) - uy*s
         r(3,2) = uz*uy*(1.0d0 - c) + ux*s
         r(3,3) = c + uz*uz*(1.0d0 - c)

      end subroutine trans_matrix_by_axis_angle

      subroutine trans_matrix_2d(angle,r)
      !============================================================
      ! Calculate the tranformation matrix in 2d
      ! counter clockwise is positive
      !------------------------------------------------------------
      ! input: 
      ! angle    - total angle for rotation (in radians)
      ! output:
      ! r(2,2)          - tranformation matrix
      !------------------------------------------------------------
         double precision angle
         double precision r(2,2)
         double precision c,s

         c=cos(angle)
         s=sin(angle)

         r(1,1)=c
         r(1,2)=-s
         r(2,1)=s
         r(2,2)=c

      end subroutine trans_matrix_2d

