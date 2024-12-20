#ifndef HAVE_DEFINE_INC
#define HAVE_DEFINE_INC
#if 0
!
! MAXNUMMAT below is the number of materials in arrays like nconst0(*),nconst(*), etc.
!
! The purpose of this is to easily expand these arrays in the future, by increasing
! the number below you will get compile errors to the places you need to make changes.
!
! Please use
!
! #include "bk19.inc"
! #include "bk19a.inc"
! #include "legal.inc"
!
! whenever any common block in these include files are needed.
!
! Note that if changing number of digits in MAXNUMMAT, then bk19.inc/bk19a.inc/legal.inc needs
! to be modified to keep the ampersand at column 73 of continuation lines.
!
! MAXNUMMAT increased to 2000 for *MODULE_USE.
! Make sure that LMATFIXD in umatss.inc is the same as MAXNUMMAT
! since both are used interchangeably in the code. We should replace
! either, but I wait with this for now since using pre-processor-includes is ugly.
!
#endif
#define MAXNUMMAT 2000

#define C_NODIV
#define C_NOMOVE
#define C_OVERLAP
#define C_NOEXTENDREORDER
#define C_NOOPT
#define C_NOARRYXX
#define C_COMPRESS
#define HAVE_CHARACTER_POINTER
#define C_LARGE_NLQ 10000
#define C_LISTVEC
#define C_SMPDEFAULT none
#define C_VEC_NONTEMP

#ifdef  C_SYMBOL
#undef  C_SYMBOL
#endif

#ifdef  ALPHA
#define MINDEX integer*8
#endif

#ifdef  IA64
#undef  HAVE_CHARACTER_POINTER
#define MINDEX integer*8
#endif

#ifdef OPEN64
#undef  HAVE_CHARACTER_POINTER
#endif

#ifdef  OSX
#undef  HAVE_CHARACTER_POINTER
#define MINDEX integer*8
#endif

#ifdef  OPENMP
#ifndef C_SYMBOL
#define C_SYMBOL
#define C_IVDEP c
#define C_REPEAT c
#define C_VDOPT c
#define C_SCALAR c
#ifdef  IBMAIX
#define C_TASKCOMMON(x) cibm* threadlocal /x/ ; save /x/
#else
#define C_TASKCOMMON(x) c$omp threadprivate (/x/)
#endif
#define C_BLANK
#define C_BEGIN_PARALLEL c$omp parallel default(C_SMPDEFAULT)
#define C_END_PARALLEL c$omp end parallel
#define C_DO_PARALLEL c$omp do
#define C_DO_END_NOWAIT c$omp end do nowait
#define C_DOALL_PARALLEL c$omp paralleldo default(C_SMPDEFAULT)
#ifdef HYBRID
#define C_ENDDO_NOWAIT c$omp end do nowait
#define C_ENDDO_PARALLEL c$omp end do
#else
#define C_ENDDO_PARALLEL c
#endif
#ifdef FJPP
#define C_DYNAMIC c$omp& schedule(static,1)
#else
#define C_DYNAMIC c$omp& schedule(dynamic,1)
#endif
#define C_STATIC c$omp& schedule(static,1)
#define C_GUIDED(x) c$omp& schedule(guided,x)
#define C_CHUNKSIZE(x) c$omp& chunksize(x)
#define C_IF_PARALLEL c$omp& if
#define C_SHARED1 c$omp& shared
#define C_SHARED2 c$omp& shared
#define C_PRIVATE c$omp& private
#define C_BARRIER c$omp barrier
#define C_BEGIN_GUARD c$omp critical
#define C_END_GUARD c$omp end critical
#define C_ORDERED c$omp& ordered
#define C_BEGIN_ORDERED c$omp ordered
#define C_END_ORDERED c$omp end ordered
#define C_DOACROSS c$omp paralleldo default(C_SMPDEFAULT)
#define C_SHARE1 c$omp& shared
#define C_SHARE2 c$omp& shared
#define C_LOCAL c$omp& private
#define C_REDUCTION(x,y) c$omp& reduction(x:y)
#define C_MP_SETLOCK c$omp critical
#define C_MP_UNSETLOCK c$omp end critical
#define C_FLUSH c$omp flush
#if !(defined(ALPHA)||defined(OPEN64))
#define C_ATOMIC c$omp atomic
#endif
#endif
#endif

#ifndef C_FLUSH
#define C_FLUSH c
#endif
#ifndef C_ATOMIC
#define C_ATOMIC c
#define C_END_ATOMIC c
#else
#define C_END_ATOMIC c$omp end atomic
#endif

#ifdef  IA64
#define C_SYMBOL
#if defined A_VECTOR
#undef C_IVDEP
#define C_IVDEP !dir$ivdep
#endif
#endif
#endif

#ifdef  NEC
#define C_SYMBOL
#undef  HAVE_CHARACTER_POINTER
#define A_VECTOR
#ifdef  MPP
#define G_VECTOR
#endif
#undef  C_IVDEP
#define C_IVDEP !cdir nodep
#undef  C_SVDEP
#define C_SVDEP !cdir nodep
#undef  C_SCALAR
#define C_SCALAR !cdir novector
#undef  C_NOARRYXX
#define C_NOARRYXX !cdir gthreorder
#undef  C_COMPRESS
#define C_COMPRESS !cdir compress
#undef  C_LISTVEC
#define C_LISTVEC !cdir listvec
#undef  C_ALTCODE
#define C_ALTCODE !cdir altcode=loopcnt
#undef  C_SHORTLOOP
#define C_SHORTLOOP !cdir shortloop
#undef  C_NOMOVEIF
#define C_NOMOVEIF !cdir nomove
#undef  C_NOLOOPCHG
#define C_NOLOOPCHG !cdir noloopchg
#undef  C_UNROLL
#define C_UNROLL(x) !cdir unroll=x
#undef  C_OVERLAP
#define C_OVERLAP *option* -O overlap
#undef  C_NODIV
#define C_NODIV *option* -O nodiv
#undef  C_NOMOVE
#define C_NOMOVE *option* -O nomove
#undef  C_NOEXTENDREORDER
#define C_NOEXTENDREORDER *option* -O noextendreorder
#undef  C_NOOPT
#define C_NOOPT *option* -NO
#undef  C_LOOPCNT
#define C_LOOPCNT *option* -pvctl loopcnt=10000
#endif

#if defined NEC || defined IA64
#ifdef FTRACE
#define C_BLANK
#define C_FTRACE_BEGIN(x) C_BLANK C_BLANK C_BLANK C_BLANK C_BLANK C_BLANK CALL FTRACE_REGION_BEGIN('x')
#define C_FTRACE_END(x) C_BLANK C_BLANK C_BLANK C_BLANK C_BLANK C_BLANK CALL FTRACE_REGION_END('x')
#endif
#endif

#ifdef  NOVECTOR
#define C_SYMBOL
#define C_IVDEP c
#define C_REPEAT c
#define C_VDOPT c
#define C_SCALAR c
#define C_TASKCOMMON c /x/
#define C_BLANK c
#define C_BEGIN_PARALLEL c
#define C_END_PARALLEL c
#define C_DO_PARALLEL c
#define C_DO_END_NOWAIT c
#define C_ENDDO_PARALLEL c
#define C_DOALL_PARALLEL_DYNAMIC c
#define C_DOALL_PARALLEL c
#define C_CHUNKSIZE(x) c
#define C_IF_PARALLEL c
#define C_SHARED1 c
#define C_SHARED2 c
#define C_PRIVATE c
#define C_BARRIER c
#define C_BEGIN_GUARD c
#define C_END_GUARD c
#define C_ORDERED c
#define C_BEGIN_ORDERED c
#define C_END_ORDERED c
#define C_DOACROSS c
#define C_SHARE1 c
#define C_SHARE2 c
#define C_LOCAL c
#define C_REDUCTION(x,y) c reduction x,y
#define C_MP_SETLOCK c
#define C_MP_UNSETLOCK c
#endif

#ifndef C_SYMBOL
#define C_IVDEP c
#define C_REPEAT c
#define C_VDOPT c
#define C_SCALAR c
#define C_TASKCOMMON(x) c
#define C_BLANK
#define C_BEGIN_PARALLEL c
#define C_END_PARALLEL c
#define C_DO_PARALLEL c
#define C_DO_END_NOWAIT c
#define C_ENDDO_PARALLEL c
#define C_DOALL_PARALLEL_DYNAMIC c
#define C_DOALL_PARALLEL c
#define C_CHUNKSIZE(x) c
#define C_IF_PARALLEL c
#define C_SHARED1 c
#define C_SHARED2 c
#define C_PRIVATE c
#define C_BARRIER c
#define C_BEGIN_GUARD c
#define C_END_GUARD c
#define C_ORDERED c
#define C_BEGIN_ORDERED c
#define C_END_ORDERED c
#define C_DOACROSS c
#define C_SHARE1 c
#define C_SHARE2 c
#define C_LOCAL c
#define C_REDUCTION(x,y) c
#define C_MP_SETLOCK c
#define C_MP_UNSETLOCK c
#endif

#ifdef  FJPP
#define C_SVDEP c
#undef  A_VECTOR
#undef  G_VECTOR
#ifdef  COLOR_SORT
#define C_SVDEP !ocl simd, norecurrence
#define A_VECTOR
#ifdef  MPP
#define G_VECTOR
#endif
#endif
#undef  C_NOOPT
#define C_NOOPT !options -O1
#undef  C_SMPDEFAULT
#define C_SMPDEFAULT none
#define C_PREFETCH(xx) !OCL PREFETCH_READ(xx(j+64))
#define IEOR xor
#define IAND and
#undef C_IVDEP
#define C_IVDEP !ocl simd, norecurrence
#undef C_VECTOR
#define C_VECTOR !ocl simd
#undef C_SCALAR
#define C_SCALAR !ocl nosimd
#endif

#ifndef MINDEX
#define MINDEX integer*8
#endif

#if defined I4R8 || defined AUTODOUBLE
#if 0
#define RECLEN 65536
#define RECSTR 32768
#define ZRECSTR 512
#else
#define RECLEN 8192
#define RECSTR 4096
#define ZRECSTR 2048
#endif
#else
#define RECLEN 8192
#define RECSTR 1024
#define ZRECSTR 512
#endif

#if defined AUTODOUBLE
#define ZREAL complex*16
#define REALX real*8
#define REAL_SZ 2
#define ZREAL_SZ 4
#define MPI_FLT MPI_DOUBLE_PRECISION
#else
#if defined I4R8
#define ZREAL complex*16
#define REALX real*8
#define REAL_SZ 2
#define ZREAL_SZ 4
#define MPI_FLT MPI_DOUBLE_PRECISION
#else
#define ZREAL complex
#define REALX real
#define REAL_SZ 1
#define ZREAL_SZ 2
#define MPI_FLT MPI_REAL
#endif
#endif

#ifdef BIGID
#define INTEGER8 integer*8
#else
#define INTEGER8 integer
#endif

#if defined AUTODOUBLE
#define FLOAT(x) dble(x)
#define DBLE(x) x
#define SNGL(x) x
#define REAL4 real
#define REAL8 real
#define IPRECS 1
#define INTEGER4 integer*4
#else
#define FLOAT(x) float(x)
#define DBLE(x) dble(x)
#define SNGL(x) sngl(x)
#define REAL4 real*4
#define REAL8 real*8
#define IPRECS 2
#define INTEGER4 integer
#endif

#if defined IFORT || defined PGI || defined IA64 || defined NEC || defined FJPP
#if defined C_SCALAR
#undef C_SCALAR
#endif
#if defined PGI
#define C_SCALAR   cpgi$l novector
#else
#if defined NEC
#define C_SCALAR   !cdir novector
#else
#if defined FJPP
#define C_SCALAR   !ocl nosimd
#else
#if defined IFORT || defined IA64
#define C_SCALAR   !dir$ novector
#endif
#endif
#endif
#endif
#endif

#include "adios.h"
#define NCPARM 256

#if defined(ALPHA) || defined(OPEN64)
# define OLD_SYSTEM 1
# define POINTER_CONTIGUOUS pointer
# define MEMFLG setmemflag
#else
# if defined(PGI165)
#  define POINTER_CONTIGUOUS pointer
# else
#  define POINTER_CONTIGUOUS pointer,contiguous
# endif
# define MEMFLG memstat
#endif

#ifdef FFTW
# ifdef AUTODOUBLE
#  ifdef __GNUC__
#   define FFT(name) dfftw_/**/name
#  else
#   define FFT(name) dfftw_##name
#  endif
# else
#  ifdef __GNUC__
#   define FFT(name) sfftw_/**/name
#  else
#   define FFT(name) sfftw_##name
#  endif
# endif
# ifndef OLD_SYSTEM
#   define FFTW_NEW 1
#   if !defined(MPP) && !defined(WINX64)
#     define FFTW_OPENMP 1
#   endif
# endif
#endif

#ifndef ISNAN
#if defined PGI || defined FJPP || defined NEC
# ifdef AUTODOUBLE
#  define ISNAN isnand
# else
#  define ISNAN isnanf
# endif
#else
#if defined IBM
# ifdef AUTODOUBLE
#  define ISNAN _isnand
# else
#  define ISNAN _isnanr
# endif
#else
# define ISNAN isnan
#endif
#endif
#endif
