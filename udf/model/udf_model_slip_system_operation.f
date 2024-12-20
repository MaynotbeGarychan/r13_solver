      subroutine ss_vec_to_global(s1,m1,r,num_ss,s11,m11)
      !============================================================
      ! Transform slip system from global to local by tranformation
      ! matrix
      !------------------------------------------------------------
      ! input: 
      ! m1(3,num_ss),s1(3,num_ss)  - global slip systems vectors
      ! r(3,3)        - tranformation matrix
      ! num_ss        - num of slip systems
      ! output:
      ! m11(3,num_ss),s11(3,num_ss) - local slip system vectors
      !------------------------------------------------------------
      implicit none
      integer l,i
      integer num_ss
      double precision m1(3,num_ss),s1(3,num_ss)
      double precision m11(3,num_ss),s11(3,num_ss)
      double precision r(3,3)

      do l=1,num_ss
            do i=1,3
                  s11(i,l)=0.
                  m11(i,l)=0.
                  s11(i,l)=s11(i,l)+r(i,1)*s1(1,l)+
     1                 r(i,2)*s1(2,l)+r(i,3)*s1(3,l)
                  m11(i,l)=m11(i,l)+r(i,1)*m1(1,l)+
     1                 r(i,2)*m1(2,l)+r(i,3)*m1(3,l)
            enddo
      enddo

      end subroutine ss_vec_to_global