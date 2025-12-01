      subroutine initCrystal(io_type,cryType,ori_arr,numSys,
     1       r,s11,m11)
        !============================================================
        ! Initialize the crystal orientation, and calculate the
        ! normalized slip systems vectors in materials coordinates.
        !------------------------------------------------------------
        ! input: 
        ! io_type        - approach to init the ori, 0-cm, 1-csv
        ! cryType       - type of crystal, 0-fcc, 1-bcc
        ! ori_arr(3)     - array (phi1,lphi,phi2) in radian
        ! output:
        ! m1(3,numSys),s1(3,numSys)  - normalized sp, ss vectors
        ! m11(3,numSys),s11(3,numSys)- normalized sp, ss vectors in 
        !                              local coordinates
        !------------------------------------------------------------
c       declare variables
        implicit none
        integer i,j,l
        integer io_type,cryType
        integer numSys
        double precision m1(3,numSys),s1(3,numSys) ! normalized vec
        double precision m11(3,numSys),s11(3,numSys) ! normalized local vec
        double precision ori_arr(3)
        double precision phi1,lphi,phi2
        double precision r(3,3) ! tranformation matrix
c       declare functions

c       init slip system vectors
        if(cryType.eq.0)then
                call getSlipSysVecFcc(m1,s1)
        elseif(cryType.eq.1)then
                call getSlipSysVecBcc(m1,s1)
        elseif(cryType.eq.2)then
                call getSlipSysVecBcc_24(m1,s1)
        endif
c       Obtain orientation from ori array
        phi1=ori_arr(1)
        lphi=ori_arr(2)
        phi2=ori_arr(3)
c       

c       calculate the tranformation matrix
        call ROTMATBYEULER(phi1,lphi,phi2,r)

c       transform the ss,sp vector to materials coordinates
        call calSlipSys2Global(s1,m1,r,numSys,s11,m11)

        end subroutine initCrystal

!         subroutine read_angle_fg(idele,phi1,lphi,phi2)
!         !============================================================
!         ! Read the euler angle from csv file, for fine grain model
!         !------------------------------------------------------------
!         ! input: 
!         ! output:
!         ! phi1,lphi,phi2  - euler angle of the crystal (element)
!         !------------------------------------------------------------
!         implicit none
!         integer i
!         integer idele
!         double precision phi1,lphi,phi2
!         double precision euler(1000,3)
!         character*50 fname
!         character*5 fidxchar
!         integer fidx,ridx

! c       obtain the posistion file-row index for the orientation input
!         ridx=mod(idele,1000)+1
!         fidx=idele/1000+1
! c       read data
!         print *,idele
!         write(fidxchar,"(I4)") fidx ! character index to integer
!         fname='./Euler_angle/set_'//trim(adjustL(fidxchar))
!      1                  //'.csv'
!         open(18,file=fname,status='old')
!         do i=1,1000
!                 read(18,*) euler(i,1),euler(i,2),euler(i,3)
!         enddo
!         close(18)
! c       give the data
!         phi1=euler(ridx,1)
!         lphi=euler(ridx,2)
!         phi2=euler(ridx,3)

!         end subroutine read_angle_fg
        
