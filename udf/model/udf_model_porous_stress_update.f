      subroutine updateSigJaum_porous(sig,eschmid,wschmid,dgamma,
     1    Wv,Dv,L_ela,f_det,num_ss,dt1,Dvv,sig_n1)
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
      ! Dvv(3)           - volumetric strain tensor in vogit
      ! output:
      ! sig_n1(6)   - cauchy stress tensor in next step
      !============================================================
      implicit none
      integer i,j
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
      double precision Dvv(3)

      call vp_stress_rate_map_tensor(sig,L_ela,eschmid,
     1   wschmid,f_det,num_ss,Ra)
      call vp_stress_rate(dgamma,Ra,num_ss,dsig_0)
      call cauchy_stress_jau_rate(sig,L_ela,dsig_0,Dv,
     1   f_det,sig_jau)
      call cauchy_stress_rot_rate(sig,Wv,sig_r)

      do i=1,6
        do j=1,3
            sig_jau(i)=sig_jau(i)-L_ela(i,j)*Dvv(j)/f_det
        enddo
      enddo
      
       do i=1,6
           sig_n1(i)=sig(i)+(sig_jau(i)+sig_r(i))*dt1
       enddo

       end subroutine updateSigJaum_porous