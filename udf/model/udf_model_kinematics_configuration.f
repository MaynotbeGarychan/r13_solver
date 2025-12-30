      subroutine ss_vec_to_configuration(fe,s11,m11,numSys,
     1   s11e,m11e)
      !============================================================
      ! Calculate spin tensor, plastic part
      !------------------------------------------------------------
      ! input: 
      ! dgamma(numSys)      - slip rate at current step
      ! eschmid(6,numSys)   - schmid tensor in vogit for each sys
      ! numSys              - number of slip system
      ! output: 
      ! Dpv(6)              - strain tensor in vogit
      !============================================================
      implicit none
      double precision fe(3,3),fe_inv(3,3)
      double precision s11(3,numSys),m11(3,numSys)
      integer numSys
      double precision s11e(3,numSys),m11e(3,numSys)
      integer l,i

      call matInverse(fe,3,fe_inv)
      do l=1,numSys
        do i=1,3
            s11e(i,l)=dot_product(fe(i,:),s11(:,l))
            m11e(i,l)=dot_product(fe_inv(i,:),m11(:,l))
        enddo
      enddo

      do l=1,numSys
        call vecNormalize(s11e(:,l),3,s11e(:,l))
        call vecNormalize(m11e(:,l),3,m11e(:,l))
      enddo
      
      end subroutine ss_vec_to_configuration