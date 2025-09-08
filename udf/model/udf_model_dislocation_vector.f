      subroutine edge_dsl_line_vec(m11,s11,numSys,l11)
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
        integer numSys
        integer l
        double precision m11(3,numSys),s11(3,numSys)
        double precision l11(3,numSys)
        
        do l=1,numSys
            call vector_cross_product(m11(:,l),
     1   s11(:,l),l11(:,l))
        enddo

      end subroutine edge_dsl_line_vec