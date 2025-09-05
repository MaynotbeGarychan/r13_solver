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

        ! if(ncycle.eq.5) then
                ! open(17,file='hsv.txt',status='old')
                ! do i=1,48
                !         do j=1,48
                !                 write(17,*) h_mat(i,j)
                !         enddo
                ! enddo
                ! close(17)
        !   endif

        end subroutine hardening_bcc