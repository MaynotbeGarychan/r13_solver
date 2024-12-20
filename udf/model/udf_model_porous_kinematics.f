      subroutine defomration_gradient_porous(fv,fv_0,fvmat)
        !============================================================
        ! construct the Fv deformation gradient
        !------------------------------------------------------------
        ! input: 
        ! fv,fv_0       - current and initial void fraction
        ! output:
        ! fvmat(3,3)      - Fv deformation gradient
        !============================================================
        implicit none
        double precision fv,fv_0
        double precision fvmat(3,3)
        double precision coeff
        
        coeff=((1.-fv_0)/(1.-fv))**(1./3.)
        call matrix_identity(3,fvmat)
        call matrix_multipy_coeff(fvmat,coeff,3,3,fvmat)

      end subroutine defomration_gradient_porous

      subroutine strain_rate_tensor_porous(fv,fvr,Dvv)
        !============================================================
        ! construct the Dv strain rate tensor
        !------------------------------------------------------------
        ! input: 
        ! fv,fvr       - void fracrion, void fraction rate
        ! output:
        ! Dvv(3)      - volumetric strain rate tensor
        !============================================================
        implicit none
        double precision fv,fvr
        double precision Dvv(3)
        double precision val
        integer i

        val=fvr/(3.*(1.-fv))
        Dvv(1)=val
        Dvv(2)=val
        Dvv(3)=val

      end subroutine strain_rate_tensor_porous
