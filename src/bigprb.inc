      logical large,mlarg,longs,lsdain,bigids,bigidb,bigidp,bigidq,i10s
      common/bigprb/large,mlarg,longs,lsdain,bigids,bigidb,bigidp,      &
     & bigidq,i10s
!
!     lsdain -replace structure input with lsda database
!     longs  -long structure input
!     i10s   -i10  structure input
!     bigids -use bigid format
!     bigidb -have big id (>214748364) for bigid option
!     bigidp -have big id (>214748364) for printing
!     bigidq -have big id (>214748364) for bigid option and special d3plot
!
      character*12 fmt16i85,fmt10i85,fmt16i5,fmt10i8,fmt8i10,fmt8e10,   &
     & fmt5e16,fmti10e10,fmt10xe10,fmt2xi8a,fmt2xi11a
      common/inpfmt/fmt16i85,fmt10i85,fmt16i5,fmt10i8,fmt8i10,fmt8e10,  &
     & fmt5e16,fmti10e10,fmt10xe10,fmt2xi8a,fmt2xi11a
!
!                 long          mlarge       88 format
!     fmt16i85    (10i20)       (10i8)       (16i5)
!     fmt10i85    (10i20)       (10i8)       (10i5)
!     fmt16i5     (8i20/8i20)   (16i5)       (16i5)
!     fmt10i8     (8i20)        (10i8)       (10i8)
!     fmt8i10     (4i20)        (8i10)       (8i10)
!     fmt8e10     (4e20.0)      (8e10.0)     (8e10.0)
!     fmt5e16     (4e20.0)      (5e16.0)     (8e10.0)
!     fmti10e10   (i20,7e20.0)  (i10,7e10.0) (i10,7e10.0)
!     fmt10xe10   (20x,7e20.0)  (10x,7e10.0) (10x,7e10.0)
!     fmt2xi8a    (2x,i20,a)    (2x,i8,a)    (2x,i8,a)
!     fmt2xi11a   (2x,i20,a)    (2x,i11,a)   (2x,i11,a)
!
!     output format for structural file
!
      character*12 sfmt16i85,sfmt10i85,sfmt16i5,sfmt10i8,sfmt8i10,      &
     &  sfmt8e10,sfmt5e16
      character*28 sfmti10e10,sfmta10e10
      common/outfmt/sfmt16i85,sfmt10i85,sfmt16i5,sfmt10i8,sfmt8i10,     &
     & sfmt8e10,sfmti10e10,sfmta10e10,sfmt5e16
!
!     output      mlarge        mlarge           88 format
!     sfmt16i85   (10i20)       (10i8)           (16i5)
!     sfmt10i85   (10i20)       (10i8)           (10i5)
!     sfmt16i5    (8i20/8i20)   (16i5)           (16i5)
!     sfmt10i8    (8e20.13)     (10i8)           (10i8)
!     sfmt8i10    (8e20.13)     (8i10)           (8i10)
!     sfmt8e10    (8e20.13)     (8e10.3)         (8e10.3)
!     sfmt5e16    (8e20.13)     (5e16.9)         (5e16.9)
!     sfmti10e10  (i20,7e20.13) (i10,7e10.3)     (i10,7e10.3)
!     sfmta10e10  (a20,7e20.13) (a10,7e10.3)     (a10,7e10.3)
!
!     output format for dynain file
!
      character*2 dfmtchar
      character*21 dfmt10i8,dfmt8i10,dfmti8e16i,                        &
     & dfmti10e10,dfmti10e10i,dfmti10e16i,dfmt2i8e16i,                  &
     & dfmt2i10e16i,dfmt8e10,dfmt5e16,dfmt2i10e10i,dfmt6ia10i,          &
     & dfmt4e16i,dfmt5i10a,dfmt2i3e10
      common/koutfmt/dfmtchar,dfmt10i8,dfmt8i10,dfmti8e16i,             &
     & dfmti10e10,dfmti10e10i,dfmti10e16i,dfmt2i8e16i,                  &
     & dfmt2i10e16i,dfmt8e10,dfmt5e16,dfmt2i10e10i,dfmt6ia10i,          &
     & dfmt4e16i,dfmt5i10a,dfmt2i3e10
!
!     output       small               i10                 long 
!     dfmtchar     '  '                ' %'                ' +'
!     dfmt10i8     (10i8)              (10i10)             (10i10)
!     dfmt8i10     (8i10)              (8i10)              (8i20)
!     dfmti8e16i   (i8,1p3e16.9,2i8)   (i10,1p3e16.9,2i10) (i20,1p3e20.13,2i20)
!     dfmti10e10   (i10,1p7e10.3)      (i10,1p7e10.3)      (i20,1p7e20.13)
!     dfmti10e10i  (i10,1p3e10.3,i10   (i10,1p3e10.3,i10)  (i20,1p3e20.13,i20)
!     dfmti10e16i  (i10,1p3e16.9,2i10) (i10,1p3e16.9,2i10) (i20,1p3e20.13,2i20)
!     dfmti10e10i  (i10,1p4e10.3,2i10) (i10,1p4e10.3,2i10) (i20,1p4e20.13,2i20)
!     dfmt2i10e16i (2i10,1p4e10.3,2i10) (2i10,1p4e10.3,2i10) (2i20,1p4e20.13,2i20)
!     dfmt8e10     (1p8e10.3)          (1p8e10.3)          (1p8e20.13)
!     dfmt5e16     (1p5e16.9)          (1p5e16.9)          (1p5e20.13)
!     dfmt2i10e10i (2i10,1p4e10.3,2i10) (2i10,1p4e10.3,2i10) (2i20,1p4e20.13,2i20)
!     dfmt6ia10i   (6i10,a10,i10)      (6i10,a10,i10)      (6i20,a20,i20)
!     dfmt4e16i    (1p4e16.8,i20)      (1p4e16.8,i20)      (1p4e20.13,i20)
!     dfmt5i10a    (5i10,9x,a)         (5i10,9x,a)         (5i20,19x,a)
!
      integer nbigp,nbigp1
!     nbigp - number of integer words for big IDs
!     nbigp1- nbigp-1
      common/bigid8/nbigp,nbigp1
