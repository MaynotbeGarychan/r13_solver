        subroutine cal_crss_by_dsl_Lee(ga,sm,bv,rho,m11,s11,l11,num_ss,
     1       g_crss)
        !============================================================
        ! Calculate the current crss by SSD
        ! A dislocation density-based single crystal constitutive 
        ! equation, IJP, Lee, 2010
        !------------------------------------------------------------
        ! input: 
        ! ga              - A=0.4
        ! sm              - shear modulus
        ! bv              - Burger vectors
        ! rho             - SSD dislocation
        ! m11(3,num_ss)   - Slip plane vectors
        ! s11(3,num_ss)   - Slip system vectors
        ! num_ss          - Number of slip system
        ! output:
        ! g_crss(num_ss)  - CRSS
        !------------------------------------------------------------
        implicit none
        integer l,k
        integer num_ss
        double precision rho(num_ss)
        double precision g_crss(num_ss)
        double precision m11(3,num_ss),s11(3,num_ss)
        double precision l11(3,num_ss)
        double precision ga,sm,bv
        double precision val

        do l=1,num_ss
                val=0.
                do k=1,num_ss
                        val=val+rho(k)*
     1            abs(dot_product(m11(:,l),l11(:,k)))
                enddo
                g_crss(l)=ga*sm*bv*sqrt(val)
        enddo

        end subroutine cal_crss_by_dsl_Lee


        subroutine dsl_evolution_kocks(ka,kb,bv,dgamma,num_ss,dt1,
     1       rho,rho_r,rho_tol)
        !============================================================
        ! Evolution model for dislocation
        ! Laws for work-hardening and low-temperature creep.
        !  J. Eng. Mater. Tech., ASMEH 98, 76., Kocks
        !------------------------------------------------------------
        ! input: 
        ! ka, kb          - Fitting constants
        ! dgamma(num_ss)  - Slip rate
        ! num_ss          - Number of slip system
        ! rho(num_ss)     - Dislocation density of slip system
        ! output:
        ! rho(num_ss)     - Dislocation density of slip system
        ! rho_r(num_ss)   - Rate ofDislocation density of slip system
        ! rho_tol         - Total dislocation (for checking)
        !------------------------------------------------------------
        implicit none
        integer l
        integer num_ss
        double precision rho_tol
        double precision rho(num_ss)
        double precision rho_r(num_ss)
        double precision dgamma(num_ss)
        double precision ka,kb,bv
        double precision dt1
        double precision val

        rho_tol=0.
        do l=1,num_ss
              rho_tol=rho_tol+rho(l)
        enddo
  
        do l=1,num_ss
              val=sqrt(rho_tol)/ka-kb*rho(l)
              rho_r(l)=abs(dgamma(l))/bv*val
              rho(l)=rho(l)+rho_r(l)*dt1
        enddo

        end subroutine dsl_evolution_kocks