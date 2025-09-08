      subroutine calElasIsoTensorCrystalGlobal(ym,pr,r,L_ela)
        !============================================================
        ! Calculate the global elasticity tensor of the crystal
        ! based on its transformation matrix
        !------------------------------------------------------------
        ! input: 
        ! ym,pr    - young's modulus, possion's ratio
        ! r(3,3)   - transformation matrix representing orientation
        ! output: 
        ! L_ela(6,6) - elastic tensor at the global coordinates
        !------------------------------------------------------------
      implicit none
      double precision ym,pr
      double precision r(3,3),RL(6,6)
      double precision L_ela(6,6),L_ela_cry(6,6)
c     calculate elastic tensor
      call calElasIsoTensor(ym,pr,L_ela_cry)
c     Transform elastic tensor to material coordinate
      call calTransMatFourthOrd(r,RL)
      call tranElasTensorLocal2Global(L_ela_cry,RL,L_ela)
      end subroutine calElasIsoTensorCrystalGlobal
      
      subroutine tranElasTensorLocal2Global(ela_ten_cry,rl,ela_ten_glb)
        !============================================================
        ! Transform the elastic tensor from crystal to global coor
        !------------------------------------------------------------
        ! input: 
        ! ela_ten_cry(6,6)  - elastic tensor in crystal coodinate
        ! output: 
        ! rl(6,6)  - Transformation matrix for four-order tensor
        ! ela_ten_glb(6,6)  - elastic tensor in global coodinate
        !------------------------------------------------------------
        implicit none
        double precision ela_ten_cry(6,6),ela_ten_glb(6,6)
        double precision rl(6,6)
        integer,parameter:: n=6

        call trans_fourth_order_tensor(ela_ten_cry,rl,n,ela_ten_glb)

      end subroutine tranElasTensorLocal2Global
      
      subroutine calElasIsoTensor(ym,pr,ela_ten)
        !============================================================
        ! Calculate the four-order elastic tensor
        !------------------------------------------------------------
        ! input: 
        ! ym     - young's modulus
        ! pr     - poisson's ratio
        ! output: 
        ! ela_ten(6,6)    - four-order elastic tensor
        !------------------------------------------------------------

        implicit none
        integer i,j
        double precision ym, pr
        double precision ela_ten(6,6)
        double precision ld,sm

        ld=pr*ym/((1+pr)*(1-2*pr))
        sm=ym/(2*(1+pr))

        ! a=pr*ym/((1+pr)*(1-2*pr))
        ! b=(1-pr)*ym/((1+pr)*(1-2*pr))
        ! c=ym/(2*(1+pr))

        do j=1,6
            do i=1,6
                ela_ten(i,j)=0.
            enddo
        enddo
        do j=1,3
            do i=1,3
                ela_ten(i,j)=ld
            enddo
        enddo
        do i=1,3
            ela_ten(i,i)=ela_ten(i,i)+2.*sm
        enddo
        do i=4,6
            ela_ten(i,i)=sm
        enddo

      end subroutine calElasIsoTensor

      subroutine compliance_tensor_otp(ym,pr,cmp_ten)
        !============================================================
        ! construct the fourth-order compliance tensor for
        ! orthotropic materials
        !------------------------------------------------------------
        ! input: 
        ! ym     - young's modulus
        ! pr     - poisson's ratio
        ! output: 
        ! cmp_ten(6,6)    - four-order compliance tensor
        !------------------------------------------------------------
        implicit none
        double precision ym,pr
        double precision cmp_ten(6,6)
        double precision a,b,c
        integer i,j

        a=1./ym
        b=2.*(1.+pr)*a
        c=-pr*a

        do i=1,6
            do j=1,6
                cmp_ten(i,j)=0.0
            enddo
        enddo

        do i=1,3
            cmp_ten(i,i)=a
        enddo
        cmp_ten(1,2)=c
        cmp_ten(1,3)=c
        cmp_ten(2,3)=c
        cmp_ten(2,1)=c
        cmp_ten(3,1)=c
        cmp_ten(3,2)=c
        do i=4,6
            cmp_ten(i,i)=b
        enddo

      end subroutine

      subroutine elastic_tensor_by_ec_iso(ec11,ec12,ec44,ela_ten)
        !============================================================
        ! Calculate the four-order elastic tensor
        !------------------------------------------------------------
        ! input: 
        ! ec11,ec12,ec44     - elastic constants for the tensor
        ! output: 
        ! ela_ten(6,6)    - four-order elastic tensor
        !------------------------------------------------------------
        implicit none
        integer i,j
        double precision ec11,ec12,ec44
        double precision ela_ten(6,6)

        do i=1,6
            do j=1,6
                ela_ten(i,j)=0.
            enddo
        enddo

        do i=1,3
            ela_ten(i,i)=ec11
        enddo

        ela_ten(1,2)=ec12
        ela_ten(1,3)=ec12
        ela_ten(2,1)=ec12
        ela_ten(2,3)=ec12
        ela_ten(3,1)=ec12
        ela_ten(3,2)=ec12

        do i=4,6
            ela_ten(i,i)=ec44
        enddo

      end subroutine elastic_tensor_by_ec_iso

      function sm_by_ym_pr(ym,pr) result(sm)
        !============================================================
        ! Calculate the shear modulus by youngs modulus and possion
        ! ratio
        !------------------------------------------------------------
        ! input: 
        ! ym,pr     - young's modulus and poisson ratio
        ! output: 
        ! sm        - shear modulus
        !------------------------------------------------------------
        implicit none
        double precision ym,pr
        double precision sm

        sm=0.5*ym/(1+pr)

      end function