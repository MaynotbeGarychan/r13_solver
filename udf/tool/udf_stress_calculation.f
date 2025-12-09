      subroutine calSigRss(sig,sf_v,numSys,tau)
        !============================================================
        ! Calculate the resolved shear stress array
        !------------------------------------------------------------
        ! input: 
        ! sig(6)          - array of the stress tensor
        ! sf_v(6,numSys)  - schmid tensor in vogit notation for each
        !                   slip system
        ! numSys          - num of slip system
        ! output:
        ! tau             - array of the resolved shear stress
        !------------------------------------------------------------
        implicit none
        include './model/define_cp.inc'
        integer i
        integer numSys
        double precision sig(6)
        double precision sf_v(6,maxSys)
        double precision tau(maxSys)

        do i=1,numSys
            tau(i)=dot_product(sf_v(:,i),sig)
        enddo

      end subroutine calSigRss

      function calSigEq(sig) result(sig_eq)
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
        
      end function calSigEq

      function calSigMean(sig) result(sig_m)
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

       end function calSigMean

       subroutine calSigDev(sig,sig_d)
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
        double precision calSigMean
        integer i

        sig_m=calSigMean
        do i=1,6
          sig_d(i)=sig(i)
        enddo
        do i=1,3
          sig_d(i)=sig_d(i)-sig_m
        enddo

       end subroutine calSigDev

       function calSigTri(sig_m,sig_eq) result(st)
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

       end function calSigTri

       function calSigLode(sig,sig_m,sig_eq) result(lode)
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

       end function calSigLode


