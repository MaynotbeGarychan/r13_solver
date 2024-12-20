      subroutine edge_dsl_line_vec(m11,s11,num_ss,l11)
      !============================================================
      ! calculate the dislocation line vector
      !------------------------------------------------------------
      ! input: 
      ! m11(3,num_ss)   - slip direction vector
      ! s11(3,num_ss)   - slip plane vector
      ! num_ss      - number of slip system
      ! output:
      ! l11(3,num_ss)   - dislocation line vector
      !============================================================
        implicit none
        integer num_ss
        integer l
        double precision m11(3,num_ss),s11(3,num_ss)
        double precision l11(3,num_ss)
        
        do l=1,num_ss
            call vector_cross_product(m11(:,l),
     1   s11(:,l),l11(:,l))
        enddo

      end subroutine edge_dsl_line_vec