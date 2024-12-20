      subroutine esbsmat_iso_sphere(pr,esbsmat)
    !============================================================
    ! calculate the eshelby tensor for iso sphere inclusion
    !------------------------------------------------------------
    ! input: 
    ! pr            - poisson ratio
    ! output:
    ! esbsmat(6,6)      - eshelby tensor
    !============================================================
        implicit none
        double precision pr
        double precision esbsmat(6,6)
        double precision a,b
        integer i,j

        a=(5.*pr-1.)/(15.-15.*pr)
        b=(4.-5.*pr)/(15.-15.*pr)
        do i=1,6
            do j=1,6
                esbsmat(i,j)=0.
            enddo
        enddo

c       ij,kl
        do i=1,3
            do j=1,3
                esbsmat(i,j)=esbsmat(i,j)+a
            enddo
        enddo
c       ik,jl
        do i=1,6
            esbsmat(i,i)=esbsmat(i,i)+b
        enddo
c       il,jk
        do i=1,3
            esbsmat(i,i)=esbsmat(i,i)+b
        enddo

      end subroutine esbsmat_iso_sphere


      subroutine esbsmat_elliptic_local(pr,r,f,e,esbmat)
    !============================================================
    ! calculate the eshelby tensor for elliptic inclusion
    ! ref: Evaluation of the Eshelby solution for the ellipsoidal
    !      inclusion ellipsoidal inclusion and heterogeneity,
    !       chunfang Meng, et. al., computers & geosciences,
    !       2012
    !------------------------------------------------------------
    ! input: 
    ! pr            - poisson ratio
    ! r(3)          - sorted array for radius of elliptic void
    ! f,e           - incomplete elliptic integral
    ! output:
    ! esbsmat(6,6)      - eshelby tensor
    !============================================================
        implicit none
        double precision pr
        double precision f,e
        double precision esbmat(6,6)
        double precision r(3)
        double precision v1,v2,v3
        double precision ifir(3),isec(3,3)
        integer i,j
        double precision,parameter:: pi=acos(-1.0d0)

c     Ifir
      v1=4*pi*r(1)*r(2)*r(3)
      v2=(r(1)**2-r(2)**2)*((r(1)**2-r(3)**2)**0.5)
      v3=f-e
      ifir(1)=v1/v2*v3
      v2=(r(2)**2-r(3)**2)*((r(1)**2-r(3)**2)**0.5)
      v3=(r(2)*((r(1)**2-r(3)**2)**0.5)/(r(1)*r(3)))-e
      ifir(3)=v1/v2*v3
      ifir(2)=4*pi-ifir(1)-ifir(3)
c     isec
      isec(1,2)=(ifir(2)-ifir(1))/(r(1)**2-r(2)**2)
      isec(2,1)=isec(1,2)
      isec(2,3)=(ifir(3)-ifir(2))/(r(2)**2-r(3)**2)
      isec(3,2)=isec(2,3)
      isec(3,1)=(ifir(3)-ifir(1))/(r(1)**2-r(3)**2)
      isec(1,3)=isec(3,1)
      isec(1,1)=(4*pi/(r(1)**2)-isec(1,2)-isec(1,3))/3
      isec(2,2)=(4*pi/(r(2)**2)-isec(2,3)-isec(2,1))/3
      isec(3,3)=(4*pi/(r(3)**2)-isec(3,1)-isec(3,2))/3
c     esbmat
      do i=1,6
            do j=1,6
                  esbmat(i,j)=0.
            enddo
      enddo
      v1=8.*pi*(1.-pr)
      v2=1./v1
      v3=(1.-2.*pr)/v1
      do i=1,3
            do j=1,3
                  if(i.eq.j) then
                        esbmat(i,j)=3*v2*(r(i)**2)*isec(i,i)+v3*ifir(i)
                  else
                        esbmat(i,j)=v2*(r(j)**2)*isec(i,j)-v3*ifir(i)
                  endif
            enddo
      enddo
      v1=16.*pi*(1.-pr)
      v2=1./v1
      v3=(1.-2.*pr)/v1
      esbmat(4,4)=v2*(r(1)**2+r(2)**2)*isec(1,2)+v3*(ifir(1)+ifir(2))
      esbmat(5,5)=v2*(r(2)**2+r(3)**2)*isec(2,3)+v3*(ifir(2)+ifir(3))
      esbmat(6,6)=v2*(r(3)**2+r(1)**2)*isec(3,1)+v3*(ifir(3)+ifir(1))

      end subroutine esbsmat_elliptic_local

      subroutine esbmat_to_crystal(esbmat,rv,rc,esbmat_cry)
    !============================================================
    ! transform the eshelby tensor from void coordinate into
    ! crystal coordinate sys.
    ! be careful with the definition of orientation
    !  cry to spe : rc       spe to void : rv
    ! therefore, void to cry: rv-1*rc-1
    !------------------------------------------------------------
    ! input: 
    ! esbmat(6,6)  - initial eshelby tensor at local
    ! rv(3,3)      - trans matrix to define void orientation
    ! rc(3,3)      - trans matrix to define crystal orientation
    ! output:
    ! esbmat_cry(6,6)    - eshelby tensor in crystal coor.
    !============================================================
        implicit none
        double precision esbmat(6,6)
        double precision rv(3,3),rc(3,3)
        double precision rvi(3,3),rci(3,3)
        double precision esbmat_cry(6,6)
        double precision r(3,3)
        double precision rl(6,6)

        call matrix_inverse(rv,3,rvi)
        call matrix_inverse(rc,3,rci)
        call matrix_inner_product(rvi,rci,3,3,3,r)
        call ROTMAT4ORD(r,rl)
        call trans_fourth_order_tensor(esbmat,rl,6,esbmat_cry)

      end subroutine esbmat_to_crystal


