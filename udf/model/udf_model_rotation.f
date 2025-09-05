       SUBROUTINE updateOrientation(Wv,Wp,s11,m11,r,num_ss,dt1,
     1  s11_n1,m11_n1,r_n1)
        !============================================================
        ! Major subroutine to update the orientation
        !------------------------------------------------------------
        ! input: 
        ! Wv(3),Wp(3)    - Spin rate tensor
        ! s11(3,num_ss)  - Unit vector slip plane at this step
        ! m11(3,num_ss)  - Unit vector slip direction at this step
        ! r(3,3)         - Trans. matrix for orientation this step
        ! num_ss         - Number of the slip systems
        ! dt1            - Time step length
        ! output:
        ! s11_n1(3,num_ss) - Unit vector slip plane at next step
        ! m11_n1(3,num_ss) - Unit vector slip direction at next step
        ! r_n1(3,3)        - Trans. matrix for orientation next step
        !------------------------------------------------------------
        implicit none
        integer num_ss
        double precision dt1
        double precision Wv(3),Wp(3)
        double precision s11(3,num_ss),m11(3,num_ss)
        double precision s11_n1(3,num_ss),m11_n1(3,num_ss)
        double precision r(3,3),r_n1(3,3),exp_we(3,3)
        call calWeExp(Wv,Wp,dt1,exp_we)
        call rot_tran_mat(r,exp_we,r_n1)
        call rot_slip_vec(s11,m11,exp_we,num_ss,s11_n1,m11_n1)
        END SUBROUTINE updateOrientation

       SUBROUTINE calWeExp(Wv,Wp,dt1,exp_we)
        !============================================================
        ! Calculate spin exponent
        !------------------------------------------------------------
        ! input: 
        ! s11(3,num_ss)  - Unit vector slip plane at this step
        ! m11(3,num_ss)  - Unit vector slip direction at this step
        ! exp_we         - Matrix for rotation calculation (expW*dt)
        ! num_ss         - Number of the slip systems
        ! output:
        ! s11_n1(3,num_ss) - Unit vector slip plane at next step
        ! m11_n1(3,num_ss) - Unit vector slip direction at next step
        !------------------------------------------------------------
        implicit none
        integer i,j,k
        double precision Wv(3),Wp(3),Wev(3)
        double precision We(3,3),idm(3,3),exp_we(3,3)
        double precision We1,dt1
        do i=1,3
            Wev(i)=Wv(i)-Wp(i)
        enddo
        We(1,1)=0.0
        We(1,2)=Wev(1)
        We(1,3)=-Wev(3)
        We(2,1)=-Wev(1)
        We(2,2)=0.0
        We(2,3)=Wev(2)
        We(3,1)=Wev(3)
        We(3,2)=-Wev(2)
        We(3,3)=0.0
        We1=sqrt(We(1,2)*We(1,2)+We(2,3)*We(2,3)+We(3,1)*We(3,1))
        if(We1==0.) then
            call matIdentity(3,exp_we)
        else
            call matIdentity(3,idm)
            do j=1,3
                  do i=1,3
                  exp_we(i,j)=idm(i,j)+(sin(We1*dt1)/We1)*We(i,j)
                  do k=1,3
                        exp_we(i,j)=exp_we(i,j)
     1                     +((1-cos(We1*dt1))/(We1**2))
     2                     *We(i,k)*We(k,j)
                  enddo
                  enddo
            enddo
        endif
      END SUBROUTINE calWeExp
       
       
       subroutine rot_slip_vec(s11,m11,exp_we,num_ss,
     1   s11_n1,m11_n1)
        !============================================================
        ! Rotate the slip system vectors for next step
        !------------------------------------------------------------
        ! input: 
        ! s11(3,num_ss)  - Unit vector slip plane at this step
        ! m11(3,num_ss)  - Unit vector slip direction at this step
        ! exp_we         - Matrix for rotation calculation (expW*dt)
        ! num_ss         - Number of the slip systems
        ! output:
        ! s11_n1(3,num_ss) - Unit vector slip plane at next step
        ! m11_n1(3,num_ss) - Unit vector slip direction at next step
        !------------------------------------------------------------

        implicit none
        integer i,k,l
        integer num_ss
        double precision s11(3,num_ss)
        double precision m11(3,num_ss)
        double precision exp_we(3,3)
        double precision s11_n1(3,num_ss)
        double precision m11_n1(3,num_ss)

        do l=1,num_ss
            do i=1,3
                s11_n1(i,l)=0.
                m11_n1(i,l)=0.
                do k=1,3
                    s11_n1(i,l)=s11_n1(i,l)
     1                  +exp_we(i,k)*s11(k,l)
                    m11_n1(i,l)=m11_n1(i,l)
     1                  +exp_we(i,k)*m11(k,l)
                enddo
            enddo
        enddo

        end subroutine rot_slip_vec

        subroutine rot_tran_mat(r,exp_we,r_n1)
        !============================================================
        ! Rotate the tranformation matrix (Local->Material)
        !------------------------------------------------------------
        ! input: 
        ! r(3,3)  - Transformation matrix from local to material 
        !           at this step
        ! exp_we  - Matrix for rotation calculation (expW*dt)
        ! output:
        ! r_n1(3,3)  - Transformation matrix from local to material 
        !               at next step
        !------------------------------------------------------------
            
            implicit none
            integer i,j,k
            double precision r(3,3)
            double precision r_n1(3,3)
            double precision exp_we(3,3)

            do j=1,3
                do i=1,3
                        r_n1(i,j)=0.
                        do k=1,3
                            r_n1(i,j)=r_n1(i,j)
     1                         +exp_we(i,k)*r(k,j)
                        enddo
                enddo
          enddo

        end subroutine rot_tran_mat

        subroutine cal_We_We1(Wv,Wp,We,We1)
        !============================================================
        ! Calculate spin rate
        !------------------------------------------------------------
        ! input: 
        ! Wv  - 
        ! We1  - 
        ! output:
        ! We(3,3)  - 
        ! We1      - 
        !------------------------------------------------------------
            implicit none
            integer i
            double precision Wv(3),Wp(3),Wev(3)
            double precision We(3,3)
            double precision We1

            ! Wev(:)=0.0
            do i=1,3
                Wev(i)=Wv(i)-Wp(i)
            enddo

            We(1,1)=0.0
            We(1,2)=Wev(1)
            We(1,3)=-Wev(3)
            We(2,1)=-Wev(1)
            We(2,2)=0.0
            We(2,3)=Wev(2)
            We(3,1)=Wev(3)
            We(3,2)=-Wev(2)
            We(3,3)=0.0
            We1=sqrt(We(1,2)*We(1,2)+We(2,3)*We(2,3)+We(3,1)*We(3,1))
            
        end subroutine cal_We_We1


