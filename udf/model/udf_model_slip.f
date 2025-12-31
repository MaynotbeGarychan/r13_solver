      subroutine calSlipRateVp(tau,g_crss,mval,dgamma_0,
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

      end subroutine calSlipRateVp

      subroutine calSlipRateHeatAct(tau,g_crss,dgamma_0,tau_0,
     1    Delta_Gk0,p,q,T,dgamma_lim,numSys,dgamma)
    !============================================================
    ! A visco-plasticity slip model based on the provided formula
    !------------------------------------------------------------
    ! input:
    ! tau(numSys)     - Stress at each slip system
    ! g_crss(numSys)  - Critical resolved shear stress at each slip system
    ! dgamma_0         - Reference slip rate
    ! tau_0           - Reference stress
    ! Delta_Gk0       - Reference change in Gibbs free energy
    ! p, q            - Model parameters
    ! T               - Temperature
    ! numSys          - Number of slip systems
    ! output:
    ! dgamma(numSys)  - Slip rate at the current step
    !------------------------------------------------------------
        implicit none
        include 'define_cp.inc'
        integer l
        integer numSys
        double precision tau(maxSys)
        double precision g_crss(maxSys)
        double precision dgamma_0,tau_0,Delta_Gk0,p,q,T
        double precision dgamma_lim
        double precision dgamma(maxSys)
        double precision tau_eff, delta_Gk

        do l = 1, numSys
            tau_eff=abs(tau(l))-g_crss(l)
            if (tau_eff.le.0) then
                dgamma(l)=0.d0
            elseif (tau_eff.ge.tau_0) then
                dgamma(l)=dgamma_0*sign(1.d0, tau(l))
                print *, 'Warning: Runs into athermal regime'
                print *, 'but handled as thermal activation.'
                print *, 'Slip system', l,' tau_eff >= tau_0'
            else
                delta_Gk=Delta_Gk0*(1.0-(tau_eff/tau_0)**p)**q
                dgamma(l)=dgamma_0*exp(-delta_Gk/(CST_BZ*T))
     1            *sign(1.d0,tau(l))
            endif
        end do

        do l=1,numSys
            if(dgamma(l) .gt. dgamma_lim) then
                  dgamma(l)=dgamma_lim
            elseif(dgamma(l) .lt. -dgamma_lim) then
                  dgamma(l)=-dgamma_lim
            endif
        enddo

      end subroutine calSlipRateHeatAct


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

