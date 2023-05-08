!----------------------------------------------
!      Function to Assign the Atomic Mass 
!
!      Written by Edison Florez using the
!      Fortran 03 version of the language
!
!            edisonffh@gmail.com
!----------------------------------------------

REAL(8) FUNCTION AtomicMass( Zi )

!
!  Purpose:
!    Fuction to assign the atomic mass according the atomic number Z 
!    Taken from https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
!               https://www.science.co.il/elements/#Notes
!
!  Record of the revisions:
!       Date            Programmer              Description of change
!       ========        ==============          ======================================
!       Feb-2018	E. Florez		v1: Original code
!
!
IMPLICIT NONE

!-Declare CALLING parameter & Data Dictionary
  INTEGER(4),INTENT(IN) :: Zi

!--------------------------------------------------------------------------------------------
!-Input Programs

  SELECT CASE (Zi)
   CASE (1)
    AtomicMass = 1.00784 
   CASE (2)
    AtomicMass = 4.002602
   CASE (3) 
    AtomicMass = 6.938
   CASE (4)
    AtomicMass = 9.0121831
   CASE (5)
    AtomicMass = 10.806
   CASE (6)
    AtomicMass = 12.0116
   CASE (7)
    AtomicMass = 14.00728
   CASE (8)
    AtomicMass = 15.99977
   CASE (9)
    AtomicMass = 18.998403163 
   CASE (10)
    AtomicMass = 20.1797
   CASE (11)
    AtomicMass = 22.98976928
   CASE (12)
    AtomicMass = 24.307
   CASE (13)
    AtomicMass =  26.9815385
   CASE (14)
    AtomicMass = 28.086
   CASE (15)
    AtomicMass = 30.973761998
   CASE (16)
    AtomicMass = 32.076
   CASE (17)
    AtomicMass = 35.457
   CASE (18)
    AtomicMass = 39.948
   CASE (19)
    AtomicMass = 39.0983
   CASE (20)
    AtomicMass = 40.078
   CASE (21)
    AtomicMass = 44.955908
   CASE (22)
    AtomicMass = 47.867
   CASE (23)
    AtomicMass = 50.9415
   CASE (24)
    AtomicMass = 51.9961
   CASE (25)
    AtomicMass = 54.938044
   CASE (26)
    AtomicMass = 55.845
   CASE (27)
    AtomicMass = 58.933194
   CASE (28)
    AtomicMass = 58.6934
   CASE (29)
    AtomicMass = 63.546
   CASE (30)
    AtomicMass = 65.38
   CASE (31)
    AtomicMass = 69.723 
   CASE (32)
    AtomicMass = 72.630 
   CASE (33)
    AtomicMass = 74.921595 
   CASE (34)
    AtomicMass = 78.971 
   CASE (35)
    AtomicMass = 79.907 
   CASE (36)
    AtomicMass = 83.798
   CASE (37)
    AtomicMass = 85.4678
   CASE (38)
    AtomicMass = 87.62
   CASE (39)
    AtomicMass = 88.90584
   CASE (40)
    AtomicMass = 91.224
   CASE (41)
    AtomicMass = 92.90637 
   CASE (42)
    AtomicMass = 95.95
   CASE (43)
    AtomicMass = 98
   CASE (44)
    AtomicMass = 101.07
   CASE (45)
    AtomicMass = 102.90550
   CASE (46)
    AtomicMass = 106.42
   CASE (47)
    AtomicMass = 107.8682
   CASE (48)
    AtomicMass = 112.414
   CASE (49)
    AtomicMass = 114.818
   CASE (50)
    AtomicMass = 118.710
   CASE (51)
    AtomicMass = 121.760
   CASE (52)
    AtomicMass = 127.60
   CASE (53)
    AtomicMass = 126.90447
   CASE (54)
    AtomicMass = 131.293
   CASE (55)
    AtomicMass = 132.90545196
   CASE (56)
    AtomicMass = 137.327
   CASE (57)
    AtomicMass = 138.90547
   CASE (58)
    AtomicMass = 140.116
   CASE (59)
    AtomicMass = 140.90766
   CASE (60)
    AtomicMass = 144.242
   CASE (61)
    AtomicMass = 145 
   CASE (62)
    AtomicMass = 150.36
   CASE (63)
    AtomicMass = 151.964
   CASE (64)
    AtomicMass = 157.25
   CASE (65)
    AtomicMass = 158.92535
   CASE (66)
    AtomicMass = 162.500
   CASE (67)
    AtomicMass = 164.93033
   CASE (68)
    AtomicMass = 167.259
   CASE (69)
    AtomicMass = 168.93422
   CASE (70)
    AtomicMass = 173.054
   CASE (71)
    AtomicMass = 174.9668
   CASE (72)
    AtomicMass = 178.49
   CASE (73)
    AtomicMass = 180.94788
   CASE (74)
    AtomicMass = 183.84
   CASE (75)
    AtomicMass = 186.207
   CASE (76)
    AtomicMass = 190.23
   CASE (77)
    AtomicMass = 192.217
   CASE (78)
    AtomicMass = 195.084
   CASE (79)
    AtomicMass = 196.966569
   CASE (80)
    AtomicMass = 200.592
   CASE (81)
    AtomicMass = 204.385
   CASE (82)
    AtomicMass = 207.2
   CASE (83)
    AtomicMass = 208.98040
   CASE (84)
    AtomicMass = 209
   CASE (85)
    AtomicMass = 210
   CASE (86)
    AtomicMass = 222
   CASE (87)
    AtomicMass = 223
   CASE (88)
    AtomicMass = 226
   CASE (89)
    AtomicMass = 227
   CASE (90)
    AtomicMass = 232.0377
   CASE (91)
    AtomicMass = 231.03588
   CASE (92)
    AtomicMass = 238.02891
   CASE (93)
    AtomicMass = 237
   CASE (94)
    AtomicMass = 244
   CASE (95)
    AtomicMass = 243
   CASE (96)
    AtomicMass = 247
   CASE (97)
    AtomicMass = 247
   CASE (98)
    AtomicMass = 251
   CASE (99)
    AtomicMass = 252
   CASE (100)
    AtomicMass = 257
   CASE (101)
    AtomicMass = 258
   CASE (102)
    AtomicMass = 259
   CASE (103)
    AtomicMass = 262
   CASE (104)
    AtomicMass = 261
   CASE (105)
    AtomicMass = 262
   CASE (106)
    AtomicMass = 266
   CASE (107)
    AtomicMass = 264
   CASE (108)
    AtomicMass = 277
   CASE (109)
    AtomicMass = 268
   CASE DEFAULT 
    WRITE(10,'(A)')""
    WRITE(10,'(A)')"****************************************"
    WRITE(10,'(A,I3,A)')"ERROR: Atomic number ",Zi," unknown"
    WRITE(10,'(A)')"****************************************"
    WRITE(10,'(A)')""
    STOP
  END SELECT 

  RETURN
END FUNCTION AtomicMass
