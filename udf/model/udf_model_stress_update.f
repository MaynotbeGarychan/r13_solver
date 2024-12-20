        subroutine update_stress_by_jau(sig,eschmid,wschmid,dgamma,
     1    Wv,Dv,L_ela,f_det,num_ss,dt1,sig_n1)
       !============================================================
       ! update the stress by jaumann method
       !------------------------------------------------------------
       ! input: 
       ! sig(6)           - cauchy stress tensor in vogit notation
       ! Wv(3)            - spin rate tensor in vogit notation
       ! sig_r(6)         - rate of cauchy stress by rotation
       ! L_ela(6,6)       - elastic tensor
       ! dsig_0(6)        - viscosoplastic stress rate
       ! Dv(6)            - strain rate tensor in vogit notation
       ! dgamma(num_ss)   - slip rate
       ! Ra(6,num_ss)  - Mapping tensor
       ! eschmid(6,num_ss)     - Schmid tensor in vogit notation
       ! wschmid(6,num_ss)     - Schmid tensor in vogit notation
       ! f_det            - determinant of F tensor
       ! num_ss           - num of slip systems
       ! output:
       ! sig_n1(6)   - cauchy stress tensor in next step
       !============================================================
       implicit none
       integer i
       integer num_ss
       double precision sig(6)
       double precision dt1
       double precision sig_n1(6)
       double precision sig_jau(6)
       double precision sig_r(6)
       double precision Wv(3)
       double precision Dv(6)
       double precision dgamma(num_ss)
       double precision eschmid(6,num_ss)
       double precision wschmid(3,num_ss)
       double precision L_ela(6,6)
       double precision dsig_0(6)
       double precision Ra(6,num_ss)
       double precision f_det

       call vp_stress_rate_map_tensor(sig,L_ela,eschmid,
     1   wschmid,f_det,num_ss,Ra)
       call vp_stress_rate(dgamma,Ra,num_ss,dsig_0)
       call cauchy_stress_jau_rate(sig,L_ela,dsig_0,Dv,
     1   f_det,sig_jau)
       call cauchy_stress_rot_rate(sig,Wv,sig_r)
       
        do i=1,6
            sig_n1(i)=sig(i)+(sig_jau(i)+sig_r(i))*dt1
        enddo

        end subroutine update_stress_by_jau

        subroutine vp_stress_rate_map_tensor(sig,L_ela,eschmid,
     1   wschmid,f_det,num_ss,Ra)
        !============================================================
        ! Calculate the mapping tensor for calculating 
        ! viscosoplastic stress rate
        !------------------------------------------------------------
        ! input: 
        ! sig(6)           - cauchy stress tensor in vogit notation
        ! L_ela(6,6)       - elastic tensor
        ! eschmid(6,num_ss)     - Schmid tensor in vogit notation
        ! wschmid(6,num_ss)     - Schmid tensor in vogit notation
        ! num_ss           - num of slip systems
        ! f_det            - determinant of F tensor
        ! output:
        ! Ra(6,num_ss)  - Mapping tensor
        !============================================================
        implicit none
        double precision sig(6)
        double precision L_ela(6,6)
        double precision eschmid(6,num_ss)
        double precision wschmid(3,num_ss)
        integer num_ss
        double precision Ra(6,num_ss)
        double precision f_det
        integer i,j,l

        do l=1,num_ss
            do i=1,6
                Ra(i,l)=0.
                do j=1,6
                    Ra(i,l)=Ra(i,l)+L_ela(i,j)*eschmid(j,l)/f_det
c                    Ra(i,l)=Ra(i,l)+L_ela(i,j)*eschmid(j,l)
                enddo
            enddo
        enddo

        do l=1,num_ss
            Ra(1,l)=Ra(1,l)
     1             +2*(wschmid(1,l)*sig(4)-wschmid(3,l)*sig(6))
            Ra(2,l)=Ra(2,l)
     1             +2*(wschmid(2,l)*sig(5)-wschmid(1,l)*sig(4))
            Ra(3,l)=Ra(3,l)
     1             +2*(wschmid(3,l)*sig(6)-wschmid(2,l)*sig(5))
            Ra(4,l)=Ra(4,l)
     1             +wschmid(1,l)*sig(2)+wschmid(2,l)*sig(6)
     2             -wschmid(1,l)*sig(1)-wschmid(3,l)*sig(5)
            Ra(5,l)=Ra(5,l)
     1             +wschmid(2,l)*sig(3)+wschmid(3,l)*sig(4)
     2             -wschmid(2,l)*sig(2)-wschmid(1,l)*sig(6)
            Ra(6,l)=Ra(6,l)
     1             +wschmid(3,l)*sig(1)+wschmid(1,l)*sig(5)
     2             -wschmid(3,l)*sig(3)-wschmid(2,l)*sig(4)
        enddo

        end subroutine vp_stress_rate_map_tensor
        
        subroutine vp_stress_rate(dgamma,Ra,num_ss,dsig_0)
        !============================================================
        ! Calculate the viscosoplastic stress rate
        !------------------------------------------------------------
        ! input: 
        ! dgamma(num_ss)   - slip rate
        ! Ra(6,num_ss)  - Mapping tensor
        ! num_ss           - num of slip systems
        ! output:
        ! dsig_0(6)        - viscosoplastic stress rate
        !============================================================
        implicit none
        double precision Ra(6,num_ss)
        integer num_ss
        double precision dgamma(num_ss)
        double precision dsig_0(6)
        integer i,l

        do i=1,6
            dsig_0(i)=0.
            do l=1,num_ss
                dsig_0(i)=dsig_0(i)+Ra(i,l)*dgamma(l)
            enddo
        enddo

        end subroutine vp_stress_rate

        subroutine cauchy_stress_jau_rate(sig,L_ela,dsig_0,Dv,f_det,
     1   sig_jau)
        !============================================================
        ! Calculate the jaumann rate of cauchy stress
        !------------------------------------------------------------
        ! input: 
        ! sig(6)           - cauchy stress tensor in vogit notation
        ! L_ela(6,6)       - elastic tensor
        ! dsig_0(6)        - viscosoplastic stress rate
        ! Dv(6)            - strain rate tensor in vogit notation
        ! f_det            - determinant of F tensor
        ! output:
        ! sig_jau(6)       - jaumann rate of cauchy stress
        !============================================================
        implicit none
        double precision sig(6)
        double precision sig_jau(6)
        double precision L_ela(6,6)
        double precision dsig_0(6)
        double precision Dv(6)
        double precision f_det
        integer i,j
    
        do i=1,6
            sig_jau(i)=-dsig_0(i)-sig(i)*(Dv(1)+Dv(2)+Dv(3))
            do j=1,6
                sig_jau(i)=sig_jau(i)+L_ela(i,j)*Dv(j)/f_det
c                sig_jau(i)=sig_jau(i)+L_ela(i,j)*Dv(j)
            enddo
        enddo
        
        end subroutine cauchy_stress_jau_rate

        subroutine cauchy_stress_rot_rate(sig,Wv,sig_r)
        !============================================================
        ! Calculate the rate of cauchy stress by rotation
        !------------------------------------------------------------
        ! input: 
        ! sig(6)           - cauchy stress tensor in vogit notation
        ! Wv(3)            - spin rate tensor in vogit notation
        ! output:
        ! sig_r(6)         - rate of cauchy stress by rotation
        !============================================================
        implicit none
        double precision sig(6)
        double precision Wv(3)
        double precision sig_r(6)

        sig_r(1)=2.*(Wv(1)*sig(4)-Wv(3)*sig(6))
        sig_r(2)=2.*(Wv(2)*sig(5)-Wv(1)*sig(4))
        sig_r(3)=2.*(Wv(3)*sig(6)-Wv(2)*sig(5)) 
        sig_r(4)=Wv(1)*sig(2)+Wv(2)*sig(6)
     1       -Wv(1)*sig(1)-Wv(3)*sig(5)
        sig_r(5)=Wv(2)*sig(3)+Wv(3)*sig(4)
     1       -Wv(2)*sig(2)-Wv(1)*sig(6)
        sig_r(6)=Wv(3)*sig(1)+Wv(1)*sig(5)
     1       -Wv(3)*sig(3)-Wv(2)*sig(4)

        end subroutine cauchy_stress_rot_rate