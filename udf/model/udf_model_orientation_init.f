      subroutine initCrystal(io_type,typeCry,ori_arr,numSys,
     1       r,s11,m11)
        !============================================================
        ! Initialize the crystal orientation, and calculate the
        ! normalized slip systems vectors in materials coordinates.
        !------------------------------------------------------------
        ! input: 
        ! io_type        - approach to init the ori, 0-cm, 1-csv
        ! typeCry       - type of crystal, 0-fcc, 1-bcc
        ! ori_arr(3)     - array (phi1,lphi,phi2) in radian
        ! output:
        ! m1(3,numSys),s1(3,numSys)  - normalized sp, ss vectors
        ! m11(3,numSys),s11(3,numSys)- normalized sp, ss vectors in 
        !                              local coordinates
        !------------------------------------------------------------
c       declare variables
        implicit none
        include 'define_cp.inc'
        integer i,j,l
        integer io_type,typeCry
        integer numSys
        double precision m1(3,maxSys),s1(3,maxSys) ! normalized vec
        double precision m11(3,maxSys),s11(3,maxSys) ! normalized local vec
        double precision ori_arr(3)
        double precision phi1,lphi,phi2
        double precision r(3,3) ! tranformation matrix
c       declare functions

c       init slip system vectors
        if(typeCry.eq.0)then
                call getSlipSysVecFcc(m1,s1)
        elseif(typeCry.eq.1)then
                call getSlipSysVecBcc(m1,s1)
        elseif(typeCry.eq.2)then
                call getSlipSysVecBcc12(m1,s1,numSys)
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