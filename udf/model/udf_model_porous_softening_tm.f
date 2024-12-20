      subroutine tm_geometric_tensor(esbsmat,vf,geomat)
      !============================================================
      ! Construct the geometric tensor for tanaka mori method
      ! considering a isotropic two phase: matrix and sphere void
      ! geomat = [vf(I-S)-1 + (1-vf)I]-1
      !------------------------------------------------------------
      ! input: 
      ! esbsmat(6,6) - corresponding eshelby tensor
      ! vf           - fraction of void
      ! output:
      ! geomat(6,6)  - geometric matrix for tm method
      !============================================================
            implicit none
            double precision esbsmat(6,6)
            integer i,j
            double precision vf,mf
            double precision geomat(6,6)

            mf=1.0-vf

            do i=1,6
                  do j=1,6
                        geomat(i,j)=0.0
                  enddo
            enddo
            do i=1,6
                  geomat(i,i)=1.0
            enddo
            do i=1,6
                  do j=1,6
                        geomat(i,j)=geomat(i,j)-esbsmat(i,j)
                  enddo
            enddo

            call matrix_inverse(geomat,6,geomat)
            do i=1,6
                  do j=1,6
                        geomat(i,j)=geomat(i,j)*vf
                  enddo
            enddo

            do i=1,6
                  geomat(i,i)=geomat(i,i)+mf
            enddo

            call matrix_inverse(geomat,6,geomat)

      end subroutine tm_geometric_tensor

      subroutine tm_porous_sm_sphere(sm,vf,pr,sm_d)
      !============================================================
      ! Calculate the damaged shear modulus using TM method
      ! considering a isotropic two phase: matrix and sphere void
      !------------------------------------------------------------
      ! input: 
      ! sm         - shear modulus
      ! vf         - fraction of void
      ! pr         - poisson ratio
      ! output:
      ! sm_d       - damaged shear modulus
      !============================================================
            implicit none
            double precision sm,vf,pr
            double precision sm_d
            double precision b
            double precision coeff

            b=(8.-10.*pr)/(15.-15.*pr)
            coeff=1.-vf/(1.-b*(1.-vf))
            sm_d=sm*coeff

      end subroutine tm_porous_sm_sphere