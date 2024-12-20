      subroutine CALRSS(sig,sf_v,num_ss,tau)
        !============================================================
        ! Calculate the resolved shear stress array
        !------------------------------------------------------------
        ! input: 
        ! sig(6)          - array of the stress tensor
        ! sf_v(6,num_ss)  - schmid tensor in vogit notation for each
        !                   slip system
        ! num_ss          - num of slip system
        ! output:
        ! tau             - array of the resolved shear stress
        !------------------------------------------------------------
        implicit none
        integer i
        integer num_ss
        double precision sig(6)
        double precision sf_v(6,num_ss)
        double precision tau(num_ss)

        do i=1,num_ss
            tau(i)=dot_product(sf_v(:,i),sig)
        enddo

      end subroutine CALRSS

      function most_active_ss(tau,num_ss) result(idx)
        !============================================================
        ! Calculate the most activated slip system by tau
        !------------------------------------------------------------
        ! input: 
        ! tau(num_ss)     - array of resolved shear stress
        ! num_ss          - num of slip plane
        ! output:
        ! idx             - index of the most activated slip system
        !------------------------------------------------------------
        implicit none
        double precision tau(num_ss)
        double precision arr(num_ss)
        integer num_ss
        integer idx
        integer l
        integer array_maxidx

        do l=1,num_ss
          arr(l)=abs(tau(l))
        enddo

        idx=array_maxidx(tau,num_ss)

      end function most_active_ss

      subroutine cal_rns(sig,nf_v,num_sp,rns)
        !============================================================
        ! Calculate the resolved normal stress array
        !------------------------------------------------------------
        ! input: 
        ! sig(6)          - array of the stress tensor
        ! nf_v(6,num_ss)  - normal schmid tensor in vogit notation 
        !                   for each slip plane
        ! num_ss          - num of slip plane
        ! output:
        ! rns             - array of the resolved normal stress
        !------------------------------------------------------------
        implicit none
        integer i
        integer num_sp
        double precision sig(6)
        double precision nf_v(6,num_sp)
        double precision rns(num_sp)

        do i=1,num_sp
            rns(i)=dot_product(nf_v(:,i),sig)
        enddo

      end subroutine cal_rns

      subroutine cal_rr(num_sp,num_ss,rns,tau,rr)
        !============================================================
        ! Calculate the resolved ratio array
        !------------------------------------------------------------
        ! input: 
        ! rns(num_sp)     - array of the resolved normal stress
        ! tau(num_ss)     - array of the resolved shear stress 
        ! num_sp          - num of slip plane
        ! num_ss          - num of slip system
        ! output:
        ! rr(num_ss)      - array of the resolved ratio
        !------------------------------------------------------------
        implicit none
        integer i,j,l
        integer num_sp,num_ss
        double precision rns(num_sp)
        double precision tau(num_ss)
        double precision rr(num_ss)

        do i=1,4
          do j=1,3
            l=3*(i-1)+j
            rr(l)=rns(i)/abs(tau(l))
          enddo
        enddo

      end subroutine

    !   function cal_rr_ss_ut(sf,nf) 
    !  1   result(rr)
    !     implicit none
    !     double precision sf(6)
    !     double precision nf(6)
    !     double precision sig(6)
    !     double precision rr
    !     double precision tau,rns
        
    !     data sig/0.0,1.0,0.0,0.0,0.0,0.0/

    !     tau=dot_product(sf(:),sig)
    !     rns=dot_product(nf(:),sig)
    !     rr=rns/tau

    !   end function cal_rr_ss_ut

      function cal_sig_eq(sig) result(sig_eq)
        !============================================================
        ! Calculate the equivalent stress
        !------------------------------------------------------------
        ! input: 
        ! sig    - Cauchy stress tensor
        ! output:
        ! sig_eq - Equivalent stress
        !------------------------------------------------------------
        implicit none
        double precision sig(6)
        double precision sig_eq
        
        sig_eq=SQRT(0.5*((sig(1)-sig(2))**2+
     1     (sig(2)-sig(3))**2+(sig(3)-sig(1))
     2            **2+(sig(4)**2+sig(5)**2+sig(6)**2)*6.))
        
      end function cal_sig_eq

      function cal_sig_m(sig) result(sig_m)
        !============================================================
        ! Calculate the mean hydrostatic stress
        !------------------------------------------------------------
        ! input: 
        ! sig(6)    - Cauchy stress tensor
        ! output:
        ! sig_m - mean hydrostatic stress
        !------------------------------------------------------------
        implicit none
        double precision sig(6)
        double precision sig_m

        sig_m=(sig(1)+sig(2)+sig(3))/3

       end function cal_sig_m

       subroutine cal_sig_dev(sig,sig_d)
        !============================================================
        ! Calculate the deviatoric stress
        !------------------------------------------------------------
        ! input: 
        ! sig(6)    - Cauchy stress tensor
        ! output:
        ! sig_d(6) - Deviatoric stress
        !------------------------------------------------------------
        implicit none
        double precision sig(6)
        double precision sig_d(6)
        double precision sig_m
        double precision cal_sig_m
        integer i

        sig_m=cal_sig_m
        do i=1,6
          sig_d(i)=sig(i)
        enddo
        do i=1,3
          sig_d(i)=sig_d(i)-sig_m
        enddo

       end subroutine cal_sig_dev

       function cal_st(sig_m,sig_eq) result(st)
        !============================================================
        ! Calculate the mean hydrostatic stress
        !------------------------------------------------------------
        ! input: 
        ! sig_m  - mean hydrostatic stress
        ! sig_eq - Equivalent stress
        ! output:
        ! st     - Stress triaxiality
        !------------------------------------------------------------
        implicit none
        double precision sig_m
        double precision sig_eq
        double precision st

        st=sig_m/sig_eq

       end function cal_st

       function cal_lode(sig,sig_m,sig_eq) result(lode)
        !============================================================
        ! Calculate the lode parameters
        !------------------------------------------------------------
        ! input: 
        ! sig(6) - Cauchy stress tensor
        ! sig_m  - mean hydrostatic stress
        ! sig_eq - Equivalent stress
        ! output:
        ! lode     - Lode parameters
        !------------------------------------------------------------

        implicit none
        double precision sig(6)
        double precision sig_m
        double precision sig_eq
        double precision lode

        lode=-(27./2.)*((sig(1)-sig_m)*(sig(2)-
     1            sig_m)*(sig(3)-sig_m))/sig_eq

       end function cal_lode


