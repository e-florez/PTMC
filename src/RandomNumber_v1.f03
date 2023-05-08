!----------------------------------------------
!      Written by Edison Florez using the
!      Fortran 03 version of the language
!
!            edisonffh@gmail.com
!----------------------------------------------

REAL(8) FUNCTION RandomNumber(Seed)
!
!  Purpose:
!     "Minimal" random number generator of Park, S.K., and Miller, K.W. 1988,
!     Communications of the ACM, vol. 31, pp. 1192–1201, similar to
!     Press, W.H. et al., Numerical Recipes in Fortran 77: The Art of Scientific
!     Computing, Cambridge University Press. 2 ED. (1992) Ch 7, pp 266-277.
!     Note that this particular random number generator has a period of  
!     ca. 2.1E+09, which is presumably not much ...?!
!
!     Set or reset "Seed" to any integer value (except the unlikely value MASK)
!     to initialize the sequence; "Seed" must not be altered between calls for 
!     successive deviates in a sequence.
!
!  Record of the revisions:
!       Date            Programmer              Description of change
!       ========        ==============          ======================================
!       Sept-2017       E. Florez               Original code 
!
IMPLICIT NONE
!-Declare LOCAL parameter & Data Dictionary
  INTEGER(8),PARAMETER :: ia  =     16807
  INTEGER(8),PARAMETER :: im  =2147483647
  REAL(8),PARAMETER    :: am  =   1.00/im
  INTEGER(8),PARAMETER :: iq  =    127773
  INTEGER(8),PARAMETER :: ir  =      2836
  INTEGER(8),PARAMETER :: mask= 123459876
  INTEGER(8) :: k
  INTEGER(8) :: Seed

!-Random number geberator
  Seed=IEOR(Seed,mask)
  k=Seed/iq
  Seed=ia*(Seed-k*iq)-ir*k
  IF(Seed.LT.0) Seed=Seed+im

  RandomNumber=am*Seed
  Seed=IEOR(Seed,mask)

!  WRITE(120,'(I20)') Seed

  RETURN
END FUNCTION RandomNumber
