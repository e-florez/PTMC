!----------------------------------------------
!      Written by Edison Florez using the
!      Fortran 03 version of the language
!
!            edisonffh@gmail.com
!----------------------------------------------

INTEGER(8) FUNCTION ClockSeed()
!
!  Purpose:
!    Function to generate a SEED to use into the random number generator 
!    function. This function uses the argument "VALUES" from the Fortran 03
!    subroutine DATE_AND_TIME(DATE, TIME, ZONE, VALUES)
!
!  Record of the revisions:
!       Date            Programmer              Description of change
!       ========        ==============          ======================================
!       Sept-2017       E. Florez              	The seed to generate the random number is the concatenation of VALUE(5), VALUE(6),
!						VALUE(7) and VALUE(8); i.e. hours//minutes//seconds//milliseconds, where
!
!                                               VALUE(1):The year
!                                               VALUE(2):The month
!                                               VALUE(3):The day of the month
!                                               VALUE(4):Time difference with UTC in minutes
!                                               VALUE(5):The hour of the day
!                                               VALUE(6):The minutes of the hour
!                                               VALUE(7):The seconds of the minute
!                                               VALUE(8):The milliseconds of the second
IMPLICIT NONE
!-Declare LOCAL parameter & Data Dictionary
  INTEGER(8),DIMENSION(8) :: Values
  CHARACTER(len=40),DIMENSION(8) :: Values_String
  CHARACTER(len=40) :: Values_string2
!  INTEGER(8) :: ClockSeed

!-Getting the seed

  CALL DATE_AND_TIME(VALUES=Values)

  WRITE(Values_String,"(I5)") Values                               !conver INT to CHARACTER array
  Values_String2=TRIM(ADJUSTL(Values_String(8)))//TRIM(ADJUSTL(Values_String(6)))//&
                 TRIM(ADJUSTL(Values_String(7)))//TRIM(ADJUSTL(Values_String(5)))
  READ(values_string2,"(I35)") ClockSeed                            !convert character (after concatenation) to int

  RETURN
END FUNCTION ClockSeed
