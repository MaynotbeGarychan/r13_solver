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