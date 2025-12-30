      subroutine updateCrss(dcrss,dt1,numSys,crss)
        !============================================================
        ! Update the crss by its rate
        !------------------------------------------------------------
        ! input: 
        ! dcrss(numSys)   - crss rate
        ! dt1             - time step
        ! crss(numSys)    - crss for current step
        ! output:
        ! crss(numSys)    - crss for next step
        !------------------------------------------------------------
        implicit none
        include 'define_cp.inc'
        integer l
        integer numSys
        double precision dcrss(maxSys)
        double precision dt1
        double precision crss(maxSys)

        do l=1,numSys
            crss(l)=crss(l)+dcrss(l)*dt1
        enddo
        end subroutine updateCrss
      
      subroutine updateCrssFcc(g0,gs,h0,hs,q,gamma_n1,dgamma,dt1,
     1      g_crss)
        !============================================================
        ! Hardening model specifically for aluminum alloy
        ! Ref: Kadi?? Computational Material Science
        !------------------------------------------------------------
        ! input: 
        ! g0    - initial crss
        ! gs    - saturated crss
        ! h0    - initial hardening rate
        ! hs    - saturated hardening rate
        ! q     - coeffecient for co-planar and conjugate hardening
        ! gamma_n1  - cumulative shear strain at t=n+1
        ! dgamma(12)- cumulative shear strain at each slip system
        ! dt1   - time step
        ! g_crss(12)    - crss for current step
        ! output:
        ! g_crss(12)    - crss for next step
        !------------------------------------------------------------
        implicit none
        integer l,k
        double precision g0,gs,h0,hs,q,gamma_n1
        double precision dgamma(12)
        double precision dt1
        double precision term1,term2
        double precision h_n1
        double precision g_crss(12)
        
c   calculate hardening rate
        term1=(h0-hs)*gamma_n1/(gs-g0)
        term2=2.*exp(-term1)/(1.+exp(-2.*term1))
        h_n1=hs+(h0-hs)*term2**2
c   renew CRSS
        do l=1,12
            do k=1,3
                if(l .le. 3) then
                        g_crss(l)=g_crss(l)
     1                          +h_n1*abs(dgamma(k))*dt1
                else
                        g_crss(l)=g_crss(l)
     1                          +q*h_n1*abs(dgamma(k))*dt1
                endif
            enddo                       
            do k=4,6
                if(l .ge. 4 .and. l .le. 6) then
                        g_crss(l)=g_crss(l)
     1                          +h_n1*abs(dgamma(k))*dt1
                else
                        g_crss(l)=g_crss(l)
     1                          +q*h_n1*abs(dgamma(k))*dt1
                endif
            enddo
            do k=7,9
                if(l .ge. 7 .and. l .le. 9) then
                        g_crss(l)=g_crss(l)
     1                          +h_n1*abs(dgamma(k))*dt1
                else
                        g_crss(l)=g_crss(l)
     1                          +q*h_n1*abs(dgamma(k))*dt1
                endif
            enddo
            do k=10,12
                if(l .ge. 10 .and. l .le. 12) then
                        g_crss(l)=g_crss(l)
     1                          +h_n1*abs(dgamma(k))*dt1
                else
                        g_crss(l)=g_crss(l)
     1                          +q*h_n1*abs(dgamma(k))*dt1
                endif
            enddo
        enddo
        end subroutine updateCrssFcc

        subroutine hardening_bcc(g0,gs,h0,hs,q,gamma_n1,dgamma,dt1,
     1      g_crss)
        !============================================================
        ! Hardening model specifically for aluminum alloy
        ! Ref: Kadi?? Computational Material Science 
        ! for BCC, consider all 48 slip system
        !------------------------------------------------------------
        ! input: 
        ! g0    - initial crss
        ! gs    - saturated crss
        ! h0    - initial hardening rate
        ! hs    - saturated hardening rate
        ! q     - coeffecient for co-planar and conjugate hardening
        ! gamma_n1  - cumulative shear strain at t=n+1
        ! dgamma(48)- cumulative shear strain at each slip system
        ! dt1   - time step
        ! g_crss(48)    - crss for current step
        ! output:
        ! g_crss(48)    - crss for next step
        !------------------------------------------------------------
        implicit none
        integer l,k
        integer i,j
        double precision g0,gs,h0,hs,q,gamma_n1
        double precision dgamma(48)
        double precision dt1
        double precision term1,term2
        double precision h_n1
        double precision g_crss(48)
        double precision h_mat(48,48)

c       calculate hardening rate
        term1=(h0-hs)*gamma_n1/(gs-g0)
        term2=2.*exp(-term1)/(1.+exp(-2.*term1))
        h_n1=hs+(h0-hs)*term2**2
c       formulate a large matrix
        do l=1,48
                do k=1,48
                        h_mat(l,k)=q*h_n1
                enddo
        enddo
c       change self slip system and {110}<111> co planar
        do l=1,48
                h_mat(l,l)=h_n1
        enddo
        do i=1,6
                l=2*i
                h_mat(l,l-1)=h_n1
                h_mat(l-1,l)=h_n1
        enddo
c       compute hardening
        do l=1,48
                do k=1,48
                        g_crss(l)=g_crss(l)+
     1       h_mat(l,k)*abs(dgamma(k))*dt1
                enddo
        enddo

        end subroutine hardening_bcc

        subroutine updateCrssDslLee(ga,sm,bv,rho,m11,s11,l11,numSys,
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
        ! m11(3,numSys)   - Slip plane vectors
        ! s11(3,numSys)   - Slip system vectors
        ! numSys          - Number of slip system
        ! output:
        ! g_crss(numSys)  - CRSS
        !------------------------------------------------------------
        implicit none
        integer l,k
        integer numSys
        double precision rho(numSys)
        double precision g_crss(numSys)
        double precision m11(3,numSys),s11(3,numSys)
        double precision l11(3,numSys)
        double precision ga,sm,bv
        double precision val

        do l=1,numSys
                val=0.
                do k=1,numSys
                        val=val+rho(k)*
     1            abs(dot_product(m11(:,l),l11(:,k)))
                enddo
                g_crss(l)=ga*sm*bv*sqrt(val)
        enddo

        end subroutine updateCrssDslLee


        subroutine calCrssRateDslHama(rho,dgamma,alpha,mu,cst_k,yc,
     1         typeCry,numSys,dcrss)
        !============================================================
        ! Calculate the crss rate for Hama sensei's model
        !------------------------------------------------------------
        ! input:
        ! alpha           - interaction coeffecient
        ! mu              - shear modulus
        ! cst_k           - constant k
        ! yc              - coeffecient for dynamic recovery
        ! typeCry        - crystal type
        ! numSys          - number of slip system
        ! rho(numSys)    - dislocation density of slip system
        ! dgamma(numSys) - slip rate
        ! output:
        ! dcrss(numSys)  - crss rate
        !------------------------------------------------------------
        implicit none
        include 'define_cp.inc'

        integer l,k
        integer numSys,typeCry
        double precision alpha,mu,cst_k,yc
        double precision dcrss(maxSys)
        double precision dgamma(maxSys)
        double precision h(maxSys,maxSys)
        double precision matInteract(maxSys,maxSys)
        double precision rho(maxSys)
        double precision sumAll, sumExcel

c       initialize h
        do l=1,numSys
           do k=1,numSys
              h(l,k)=0.0d0
           enddo
        enddo

        call getDslInteractionMatrix(typeCry,numSys,matInteract)

c       calculation of h matrix
        do l=1,numSys
           sumAll   = 0.0d0
           sumExcel= 0.0d0

           do k=1,numSys
              sumAll = sumAll + matInteract(l,k)*rho(k)
              if(k.ne.l) then
                 sumExcel = sumExcel + matInteract(l,k)*rho(k)
              endif
           enddo

c          numerical safeguard
           sumAll   = max(sumAll  ,1.0d-20)
           sumExcel= max(sumExcel,1.0d-20)
           do k=1,numSys
              h(l,k) = (alpha*mu/2.0d0)*matInteract(l,k)
     &               * sumAll**(-0.5d0)
     &               * ( (1.0d0/cst_k)*(sumExcel**0.5d0)
     &               - 2.0d0*yc*rho(k) )
           enddo
        enddo

c       calculation of crss rate
        do l=1,numSys
           dcrss(l)=0.0d0
           do k=1,numSys
              dcrss(l)=dcrss(l)+h(l,k)*abs(dgamma(k))
           enddo
        enddo

        return
        end


