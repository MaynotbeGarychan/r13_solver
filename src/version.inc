      character*10 git_date,git_time,git_commit
      parameter (git_date     = '06/30/2022')
      parameter (git_time     = '15:57:06')
      parameter (git_commit   = 'e41832f8f0')
      character*32 git_branch
      parameter (git_branch   = 'releases/R13.1.1')
      character*10 git_short_branch
      parameter (git_short_branch = 'R13.1.1')
      integer git_branch_count
      parameter (git_branch_count = 6)
      integer git_master_count
      parameter (git_master_count = 83656)
      character*32 git_build_id
      parameter (git_build_id   = 'R13.1.1-6-ge41832f8f0')
      character*10 build_date
      parameter (build_date   = '07/25/2022')
      character*64 build_host
      parameter (build_host   =
     &'')
      character*64 build_src_path
      parameter (build_src_path=
     &'F:\git\ls971r1311b180d6u')
      character*64 build_path
      parameter (build_path=
     &'F:\git\ls971r1311b180d6u\smpd6480')
      character*1024 build_cmd
      parameter (build_cmd=
     &'i -j12 --jobserver-auth=3,13 -- ANSYS=ans222+ LSDYNA=lsdyna RELEA!
     &SE=no LSDAIN=yes MUMPS=yes DUALCESE=no BATTERY=yes MAKE_LAPACK=no !
     &MAKE_BLAS=no USA=no MADCL=no OPTION=sse FCC=180 MSC=no FFTW=no BLA!
     &SLIB=mkli8r8 SPRNG_BUILD=no SPRNG=no SPRNG_PROD=no LAPACK_PROD=no !
     &BLAS_PROD=no EXTENDED=yes PFEM=yes HAVE_CPP=yes HAVE_CXX=yes HAVES!
     &PV7=yes DOUBLE=yes USE_SPOOLESV7=yes ANSYSLIC=ans211 ANSYSOUT=no M!
     &AKEDIR=g:/make/femster/X64 SYSTEM=win64 BASEDIR=../')
