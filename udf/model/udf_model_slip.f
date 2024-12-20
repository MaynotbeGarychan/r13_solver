      subroutine cal_slip_pl(tau,g_crss,mval,dgamma_0,
     1   dgamma_lim,num_ss,dgamma)
        !============================================================
        ! A power-law style visco-plasticity slip model 
        !------------------------------------------------------------
        ! input: 
        ! tau(num_ss)     - RSS at each slip system
        ! g_crss(num_ss)  - CRSS at each slip system
        ! mval            - Sensitivity of slip rate
        ! dgamma_0        - Reference slip rate
        ! dgamma_lim      - Limit of slip rate
        ! num_ss          - Num of slip system
        ! output:
        ! dgamma(num_ss)  - Slip rate at current step
        !------------------------------------------------------------
        implicit none
        integer l
        integer num_ss
        double precision dgamma_0
        double precision tau(num_ss)
        double precision g_crss(num_ss)
        double precision mval
        double precision dgamma_lim
        double precision dgamma(num_ss)

        do l=1,num_ss
            dgamma(l)=dgamma_0*(tau(l)/g_crss(l))
     1         *abs(tau(l)/g_crss(l))**((1/mval)-1)
        enddo

        do l=1,num_ss
            if(dgamma(l) .gt. dgamma_lim) then
                  dgamma(l)=dgamma_lim
            elseif(dgamma(l) .lt. -dgamma_lim) then
                  dgamma(l)=-dgamma_lim
            endif
        enddo


      end subroutine cal_slip_pl


      subroutine update_ccs(dgamma,num_ss,dt1,
     1      dgamma_tol,gamma_slip,gamma_n1)
        !============================================================
        ! Update of the slip deformation volume
        !------------------------------------------------------------
        ! input: 
        ! dgamma(num_ss)  - Slip rate at current step
        ! num_ss          - Num of slip system
        ! dt1             - Time step
        ! gamma_slip(num_ss)    - Cumulative slip deformation at sysm
        ! output:
        ! dgamma_tol      - Slip deformation of all sysm at curr step
        ! gamma_slip(num_ss)    - Cumulative slip deformation at sysm
        ! gamma_n1        - Cumulative slip deformation of all sysm
        !------------------------------------------------------------
        implicit none
        integer l
        integer num_ss
        double precision dt1
        double precision dgamma(num_ss)
        double precision gamma_slip(num_ss)
        double precision dgamma_tol,gamma_n1
        double precision val

        val=0
        dgamma_tol=0
        do l=1,num_ss
            val=abs(dgamma(l))*dt1
            gamma_slip(l)=gamma_slip(l)+val
            dgamma_tol=dgamma_tol+val
            gamma_n1=gamma_n1+val
        enddo

      end subroutine update_ccs

!       subroutine cal_slip_thm_act(tau,g_crss,rho_ssdm,vec_b,eng_act,k_bolz
!      1   num_ss,dgamma)
!         !============================================================
!         ! A power-law style visco-plasticity slip model 
!         !------------------------------------------------------------
!         ! input: 
!         ! tau(num_ss)     - RSS at each slip system
!         ! g_crss(num_ss)  - CRSS at each slip system
!         ! num_ss          - Num of slip system
!         ! output:
!         ! dgamma(num_ss)  - Slip rate at current step
!         !------------------------------------------------------------
!         implicit none
!         integer l
!         integer num_ss
!         double precision dgamma_0
!         double precision tau(num_ss)
!         double precision g_crss(num_ss)
!         double precision mval
!         double precision dgamma_lim
!         double precision dgamma(num_ss)

!         do l=1,12
!             dgamma(l)=dgamma_0*(tau(l)/g_crss(l))
!      1         *abs(tau(l)/g_crss(l))**((1/mval)-1)
!         enddo

!         do l=1,12
!             if(dgamma(l) .gt. dgamma_lim) then
!                   dgamma(l)=dgamma_lim
!             elseif(dgamma(l) .lt. -dgamma_lim) then
!                   dgamma(l)=-dgamma_lim
!             endif
!         enddo

!       end subroutine cal_slip_thm_act