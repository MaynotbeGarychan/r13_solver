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
        include 'define_cp.inc'
        integer l
        integer numSys
        double precision dgamma_0
        double precision tau(maxSys)
        double precision g_crss(maxSys)
        double precision mval
        double precision dgamma_lim
        double precision dgamma(maxSys)

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
      
      subroutine calSlipRateHeatAct(tau,g_crss,tau0110,tau0112,
     1   dgk0,dgamma_0,pval,qval,tval,kb,numSys,dgamma)

        implicit none
        include 'define_cp.inc'
        integer i,numSys
        double precision tau(maxSys),g_crss(maxSys)
        double precision tauEff(maxSys),dgk(maxSys)
        double precision dgamma(maxSys)
        double precision tau0110,tau0112
        double precision dgk0,pval,qval,tval,dgamma_0
        double precision kb

        do i=1,numSys
            tauEff(i)=abs(tau(i))-g_crss(i)
        enddo

c       1~12
        do i=1,12
            dgk(i)=dgk0*(1.-(tauEff(i)/tau0110)**pval)**qval
        enddo
c       13~24
        do i=13,24
            dgk(i)=dgk0*(1.-(tauEff(i)/tau0112)**pval)**qval
        enddo

        do i=1,numSys
            if(tauEff(i).gt.0)then
                dgamma(i)=dgamma_0*exp(0.-dgk(i)/(kb*tval))
            else
                dgamma(i)=0.
            endif
        enddo

        end subroutine

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
        include 'define_cp.inc'
        integer l
        integer numSys
        double precision dt1
        double precision dgamma(maxSys)
        double precision gamma_slip(maxSys)
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

