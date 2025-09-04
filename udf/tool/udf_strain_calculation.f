       function calPeeq(Dp,peeq,dt1) result(peeq_n1)
        !============================================================
        ! Calculate the peeq based on increament
        !------------------------------------------------------------
        ! input: 
        ! Dp(6)   - Plastic strain tensor in Voigt
        ! peeq    - PEEQ at current step
        ! dt1     - Time increment
        ! output:
        ! peeq_n1 - PEEQ after increament
        !------------------------------------------------------------

        implicit none
        double precision Dp(6)
        double precision peeq
        double precision dpeeq
        double precision dt1
        double precision peeq_n1

        dpeeq=(SQRT(2.)/3.)*SQRT((Dp(1)-Dp(2))
     1            **2+(Dp(2)-Dp(3))**2+(Dp(3)-Dp(1))**2+
     2            1.5*(Dp(4)**2+Dp(5)**2+Dp(6)**2))
        peeq_n1=peeq+dpeeq*dt1

       end function calPeeq

       function calPeeqr(Dp) result(peeqr)
        !============================================================
        ! Calculate the peeqr
        !------------------------------------------------------------
        ! input: 
        ! Dp(6)   - Plastic strain tensor in Voigt
        ! output:
        ! peeqr - rate of peeq
        !------------------------------------------------------------

         implicit none
         double precision Dp(6)
         double precision peeqr

         peeqr=(SQRT(2.)/3.)*SQRT((Dp(1)-Dp(2))
     1            **2+(Dp(2)-Dp(3))**2+(Dp(3)-Dp(1))**2+
     2            1.5*(Dp(4)**2+Dp(5)**2+Dp(6)**2))

       end function calPeeqr
         


