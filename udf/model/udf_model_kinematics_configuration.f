      subroutine ss_vec_to_configuration(fe,s11,m11,num_ss,
     1   s11e,m11e)
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
      double precision fe(3,3),fe_inv(3,3)
      double precision s11(3,num_ss),m11(3,num_ss)
      integer num_ss
      double precision s11e(3,num_ss),m11e(3,num_ss)
      integer l,i

      call matrix_inverse(fe,3,fe_inv)
      do l=1,num_ss
        do i=1,3
            s11e(i,l)=dot_product(fe(i,:),s11(:,l))
            m11e(i,l)=dot_product(fe_inv(i,:),m11(:,l))
        enddo
      enddo

      do l=1,num_ss
        call vector_normalize(s11e(:,l),3,s11e(:,l))
        call vector_normalize(m11e(:,l),3,m11e(:,l))
      enddo
      
      end subroutine ss_vec_to_configuration