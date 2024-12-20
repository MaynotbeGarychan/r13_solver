      function dge_to_rad(deg) result(rad)
        !============================================================
        ! convert degree to radian
        !------------------------------------------------------------
        ! input: 
        ! deg     - degrees
        ! output:
        ! rad     - radian
        !------------------------------------------------------------
        implicit none
        double precision deg
        double precision rad
        double precision, parameter :: pi=acos(-1.0d0)

        rad=deg*pi/180

      end function 

      subroutine deg_to_rad_arr(arr,num)

        implicit none
        integer num
        double precision arr(num)
        integer i
        double precision dge_to_rad

        do i=1,num
          arr(i)=dge_to_rad(arr(i))
        enddo
        
      end subroutine deg_to_rad_arr