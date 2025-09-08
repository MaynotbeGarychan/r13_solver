      subroutine calSlipRate(tau,g_crss,mval,dgamma_0,
     1   dgamma_lim,numSys,dgamma)
        !============================================================
        ! A power-law style visco-plasticity slip model 
        !------------------------------------------------------------
        ! input: 
        ! tau(numSys)     - RSS at each slip system
        ! g_crss(numSys)  - CRSS at each slip system
        ! mval            - Sensitivity of slip rate
        ! dgamma_0        - Reference slip rate
        ! dgamma_lim      - Limit of slip rate
        ! numSys          - Num of slip system
        ! output:
        ! dgamma(numSys)  - Slip rate at current step
        !------------------------------------------------------------
        implicit none
        integer l
        integer numSys
        double precision dgamma_0
        double precision tau(numSys)
        double precision g_crss(numSys)
        double precision mval
        double precision dgamma_lim
        double precision dgamma(numSys)

        do l=1,numSys
            dgamma(l)=dgamma_0*(tau(l)/g_crss(l))
     1         *abs(tau(l)/g_crss(l))**((1/mval)-1)
        enddo

        do l=1,numSys
            if(dgamma(l) .gt. dgamma_lim) then
                  dgamma(l)=dgamma_lim
            elseif(dgamma(l) .lt. -dgamma_lim) then
                  dgamma(l)=-dgamma_lim
            endif
        enddo


      end subroutine calSlipRate


      subroutine updateCss(dgamma,numSys,dt1,
     1      dgamma_tol,gamma_slip,gamma_n1)
        !============================================================
        ! Update of the slip deformation volume
        !------------------------------------------------------------
        ! input: 
        ! dgamma(numSys)  - Slip rate at current step
        ! numSys          - Num of slip system
        ! dt1             - Time step
        ! gamma_slip(numSys)    - Cumulative slip deformation at sysm
        ! output:
        ! dgamma_tol      - Slip deformation of all sysm at curr step
        ! gamma_slip(numSys)    - Cumulative slip deformation at sysm
        ! gamma_n1        - Cumulative slip deformation of all sysm
        !------------------------------------------------------------
        implicit none
        integer l
        integer numSys
        double precision dt1
        double precision dgamma(numSys)
        double precision gamma_slip(numSys)
        double precision dgamma_tol,gamma_n1
        double precision val

        val=0
        dgamma_tol=0
        do l=1,numSys
            val=abs(dgamma(l))*dt1
            gamma_slip(l)=gamma_slip(l)+val
            dgamma_tol=dgamma_tol+val
            gamma_n1=gamma_n1+val
        enddo

      end subroutine updateCss

