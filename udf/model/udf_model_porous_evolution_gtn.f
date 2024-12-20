      subroutine porous_evl_gtn_growth(vf,Dpv,vfr)
        !============================================================
        ! Calculate the void growth for rate
        ! (1-f)tr(Dp)
        !------------------------------------------------------------
        ! input: 
        ! vf      - void fraction
        ! Dpv(6)  - plastic strain rate tensor in vogit
        ! output:
        ! vfr     - void fraction rate
        !------------------------------------------------------------
        implicit none
        double precision vf
        double precision Dpv(6)
        double precision vfr
        integer i

        vfr=(1.-vf)*(Dpv(1)+Dpv(2)+Dpv(3))

      end subroutine porous_evl_gtn_growth

      subroutine porous_evl_gtn_nucleation(peeq,peeqr,
     1   fn,peeqn,sn,vfr)
    !============================================================
    ! Calculate the void nucleation for rate
    ! An*peeqr
    !------------------------------------------------------------
    ! input: 
    ! peeq      - accumulative plastic strain PEEQ
    ! peeqr     - rate of PEEQ at current step
    ! fn,peeqn,sn  - three variables to control the distribution
    ! output:
    ! vfr     - void fraction rate
    !------------------------------------------------------------
      implicit none
      double precision peeq
      double precision peeqr
      double precision fn,sn,peeqn
      double precision vfr
      double precision distribution_nml
      double precision val

      val=distribution_nml(fn,peeqn,sn,peeq)
      vfr=val*peeqr

      end subroutine porous_evl_gtn_nucleation

      subroutine porous_evl_gtn_shear(vf,sig_d,sig_eq,lode,Dpv,
     1   vfr)
    !============================================================
    ! Calculate the shearing deformation for void growth
    ! kw*f*w(sig)*(sij*peeq_r_ij)/sig_eq
    ! Ref: Tension-torsion fracture experiments - Part II
    ! Zhenyu Xue, et. al., IJSS
    !------------------------------------------------------------
    ! input: 
    ! vf      - void fraction
    ! sig_d   - deviatoric stress tensor
    ! sig_eq  - equivalent stress
    ! lode    - lode parameter
    ! Dpv(6)  - plastic strain rate tensor in vogit notation
    ! output:
    ! vfr     - void fraction rate
    !------------------------------------------------------------
        implicit none
        double precision vf
        double precision sig_d(6)
        double precision sig_eq
        double precision lode
        double precision Dpv(6)
        double precision vfr
        double precision val
        double precision a,b
        double precision kw
        integer i
c       sij*epsilon_r_ij
        val=0.
        do i=1,3
            val=val+sig_d(i)*Dpv(i)
        enddo
        do i=4,6
            val=val+sig_d(i)*Dpv(i)*0.5
        enddo
c       w(sigma)
        a=lode**2
        b=(((a-1)**2))/((3+a)**3)*27.
c       vfr
        kw=1.0
        vfr=kw*vf*b*val/sig_eq

      end subroutine porous_evl_gtn_shear

      subroutine porous_evl_rice_tracey(fv,dpeeq,st,fv_n1)
    !============================================================
    ! Updateting void growth using rice tracey model
    !------------------------------------------------------------
    ! input: 
    ! fv      - void fraction
    ! dpeeq   - peeq rate
    ! st      - stress triaxiality
    ! output:
    ! fv_n1     - void fraction in next step
    !------------------------------------------------------------
        implicit none
        double precision fv,dpeeq,st,fv_n1
        if((st.lt.20.0).and.(st.gt.-20.0))then
                fv_n1=fv+0.283*dpeeq*exp(1.5*st)
        else
                fv_n1=fv
        endif
        end subroutine