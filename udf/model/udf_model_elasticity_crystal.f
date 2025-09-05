c slip system fcc
cdata m/-1.00,1.00,1.00, -1.00,1.00,1.00, -1.00,1.00,1.00,
c1         1.00,1.00,1.00, 1.00,1.00,1.00, 1.00,1.00,1.00,
c2        -1.00,-1.00,1.00, -1.00,-1.00,1.00, -1.00,-1.00,1.00,
c3        1.00,-1.00,1.00, 1.00,-1.00,1.00, 1.00,-1.00,1.00/
c       data s/0.00,-1.00,1.00, 1.00,0.00,1.00, 1.00,1.00,0.00,
c1        0.00,-1.00,1.00, -1.00,0.00,1.00, -1.00,1.00,0.00,
c2        0.00,1.00,1.00, 1.00,0.00,1.00, -1.00,1.00,0.00,
c3        0.00,1.00,1.00, -1.00,0.00,1.00, 1.00,1.00,0.00/
c cauchy stress array
c     sig(1)=local x  stress
c     sig(2)=local y  stress
c     sig(3)=local z  stress
c     sig(4)=local xy stress
c     sig(5)=local yz stress
c     sig(6)=local zx stress
      subroutine shear_modulus_ss_safe(ela_ten_cry,sm_arr)

        implicit none
        double precision ela_ten_cry(6,6)
        double precision ela_ten_rot(6,6)
        double precision sm_arr(12)
        double precision axis(3)
        double precision angle
        double precision r(3,3)
        double precision rl(6,6)
        double precision dge_to_rad
        integer i,j

c       0,-1,1  & 0,1,1
c       construct transformation matrix
        axis(1)=1.0
        axis(2)=0.0
        axis(3)=0.0
        angle=45
        angle=dge_to_rad(angle)
        call trans_matrix_by_axis_angle(axis,angle,r)
        call calTransMatFourthOrd(r,rl)
        call tranElasTensorLocal2Global(ela_ten_cry,rl,ela_ten_rot)
        sm_arr(1)=ela_ten_rot(4,4)
        sm_arr(4)=ela_ten_rot(4,4)
        sm_arr(7)=ela_ten_rot(6,6)
        sm_arr(10)=ela_ten_rot(6,6)

!         open(18,file='chk_sm.txt',status='replace')
!         write(18,*) 'r(3,3)'
!         do i=1,3
!                 write(18,*) r(i,1),r(i,2),r(i,3)
!         enddo
!         write(18,*) 'rl(6,6)'
!         do i=1,6
!                 write(18,*) rl(i,1),rl(i,2),rl(i,3),
!      1       rl(i,4),rl(i,5),rl(i,6)
!         enddo
!         close(18)

c       1,0,1  & -1,0,1
        axis(1)=0.0
        axis(2)=1.0
        axis(3)=0.0
        angle=45
        angle=dge_to_rad(angle)
        call trans_matrix_by_axis_angle(axis,angle,r)
        call calTransMatFourthOrd(r,rl)
        call tranElasTensorLocal2Global(ela_ten_cry,rl,ela_ten_rot)
        sm_arr(2)=ela_ten_rot(4,4)
        sm_arr(8)=ela_ten_rot(4,4)
        sm_arr(5)=ela_ten_rot(5,5)
        sm_arr(11)=ela_ten_rot(5,5)

c       1,1,0  &  -1,1,0
        axis(1)=0.0
        axis(2)=0.0
        axis(3)=1.0
        angle=45
        angle=dge_to_rad(angle)
        call trans_matrix_by_axis_angle(axis,angle,r)
        call calTransMatFourthOrd(r,rl)
        call tranElasTensorLocal2Global(ela_ten_cry,rl,ela_ten_rot)
        sm_arr(3)=ela_ten_rot(5,5)
        sm_arr(12)=ela_ten_rot(5,5)
        sm_arr(6)=ela_ten_rot(6,6)
        sm_arr(9)=ela_ten_rot(6,6)

      end subroutine shear_modulus_ss_safe

      subroutine shear_modulus_ss_fast(ela_ten_cry,sm_arr)

        implicit none
        double precision ela_ten_cry(6,6)
        double precision sm_arr(12)
        double precision axis(3)
        double precision angle
        double precision dge_to_rad
        double precision r(3,3),rl(6,6)
        double precision shear_modulus_ss_fast_func
        integer smid
        double precision val

        angle=45
        angle=dge_to_rad(angle)
c       0,-1,1  & 0,1,1
        axis(1)=1.0
        axis(2)=0.0
        axis(3)=0.0
        call trans_matrix_by_axis_angle(axis,angle,r)
        call calTransMatFourthOrd(r,rl)
        smid=4
        val=shear_modulus_ss_fast_func(ela_ten_cry,rl,smid)
        sm_arr(1)=val
        sm_arr(4)=val
        smid=6
        val=shear_modulus_ss_fast_func(ela_ten_cry,rl,smid)
        sm_arr(7)=val
        sm_arr(10)=val
c       1,0,1  & -1,0,1
        axis(1)=0.0
        axis(2)=1.0
        axis(3)=0.0
        call trans_matrix_by_axis_angle(axis,angle,r)
        call calTransMatFourthOrd(r,rl)
        smid=4
        val=shear_modulus_ss_fast_func(ela_ten_cry,rl,smid)
        sm_arr(2)=val
        sm_arr(8)=val
        smid=5
        val=shear_modulus_ss_fast_func(ela_ten_cry,rl,smid)
        sm_arr(5)=val
        sm_arr(11)=val
c       1,1,0  &  -1,1,0
        axis(1)=0.0
        axis(2)=0.0
        axis(3)=1.0
        call trans_matrix_by_axis_angle(axis,angle,r)
        call calTransMatFourthOrd(r,rl)
        smid=5
        val=shear_modulus_ss_fast_func(ela_ten_cry,rl,smid)
        sm_arr(3)=val
        sm_arr(12)=val
        smid=6
        val=shear_modulus_ss_fast_func(ela_ten_cry,rl,smid)
        sm_arr(6)=val
        sm_arr(9)=val

      end subroutine shear_modulus_ss_fast

      function shear_modulus_ss_fast_func(ela_ten_cry,
     1       rl,smid) result(val)

        implicit none
        double precision axis(3)
        integer smid
        double precision val
        double precision ela_ten_cry(6,6)
        double precision vec(3)
        double precision rl(6,6)
        integer i,j,a

        do i=1,3
                a=i+3
                vec(i)=dot_product(rl(smid,4:6),ela_ten_cry(4:6,a))
        enddo
        val=dot_product(vec(:),rl(smid,4:6))
      end function shear_modulus_ss_fast_func