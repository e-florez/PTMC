!----------------------------------------------
!    Molecular Displacement of CM & Rotation
!          for a Selectec Fragment
!
!      Written by Edison Florez using the
!      Fortran 03 version of the language
!
!            edisonffh@gmail.com
!----------------------------------------------

SUBROUTINE RotateCluster( r,NMolecules,MaxRotation,Seed,NTrajectories,t,MiddlePoint,&   !IN
                          r_rot                                                      )  !OUT

!
!  Purpose:
!    Subroutine to rotate a whole cluster 
!
!  Record of the revisions:
!       Date            Programmer              Description of change
!       ========        ==============          ======================================
!       Dec-2018	E. Florez		v1:  Original code
!       xxx-2018	E. Florez		v2:
!
IMPLICIT NONE

!-Declare CALLING parameter & Data Dictionary
  INTEGER(4),INTENT(IN) :: NMolecules                                         !Number of Molecules 
  INTEGER(4),INTENT(IN) :: NTrajectories,t                                    !Number of temperarures
  REAL(8),INTENT(IN) :: MaxRotation                                           !Displacement could be adjusted in the fly
  REAL(8),DIMENSION(0:NTrajectories,NMolecules,1,3),INTENT(IN) :: r           !Molecular Coordinades x,y,z; r(Molec,Atom,x), means rank=3
  INTEGER(8),INTENT(IN) :: Seed                                               !The Seed for the Random Number generator function based on clock system
  REAL(8),DIMENSION(3),INTENT(IN) :: MiddlePoint                              !coordines of the center of SBC to rotate with respect of it
  REAL(8),DIMENSION(NMolecules,3),INTENT(OUT) :: r_rot                        !New Molecular Coordinades x,y,z for selected molecule 
                                                                            
!-Declare LOCAL parameter & Data Dictionary                                 
  INTEGER(8) :: i,k                                                           !The indices mening are: i=Molecule, j=AtomLabel, k=Coordinates
  REAL(8) :: a,b,c                                                            !Euler angles
  REAL(8),DIMENSION(3) :: Rotation                                            !Random rotation according Euler angles
  REAL(8),DIMENSION(:,:),ALLOCATABLE :: r_ppal                                !Coordinades x,y,z of principal axes
  REAL(8) :: RandomNumber                                                     !Random Number generator function

!--------------------------------------------------------------------------------------------
!-Allocating variables
  ALLOCATE(r_ppal(NMolecules,3))

!--------------------------------------------------------------------------------------------
!-Euler rotations:
!
!  George Arfken Hans Weber Frank E. Harris. Mathematical Methods for Physicists:
!  A Comprehensive Guide. Academic Press. 7th Edition. (2012) pag 139, and
!  
!  WOLFRAMALPHA:https://www.wolframalpha.com/input/?i=%7B%7BCos%5Bc%5D,+Sin%5Bc%5D,
!               +0%7D,+%7B-Sin%5Bc%5D,+Cos%5Bc%5D,+0%7D,+%7B0,+0,+1%7D%7D+.+%7B%7B
!               Cos%5Bb%5D,+0,+-Sin%5Bb%5D%7D,+%7B0,+1,+0%7D,+%7BSin%5Bb%5D,+0,+Cos
!               %5Bb%5D%7D%7D+.+%7B%7BCos%5Ba%5D,+Sin%5Ba%5D,+0%7D,+%7B-Sin%5Ba%5D,
!               +Cos%5Ba%5D,+0%7D,+%7B0,+0,+1%7D%7D
!--------------------------------------------------------------------------------------------




  !computing principal axis to rotate
  DO k=1,3
     Rotation(k)=(2.0*RandomNumber(Seed)-1.0)*MaxRotation    
     DO i=1,NMolecules !Findind th principal axes for rotation
        r_ppal(i,k)=r(t,i,1,k)-MiddlePoint(k)
     ENDDO
  ENDDO

  !Rotating the whole cluster 
  a=Rotation(1) ; b=Rotation(2) ; c=Rotation(3)
  DO i=1,NMolecules
     r_rot(i,1)=r_ppal(i,1)*( COS(a)*COS(b)*COS(c) - SIN(a)*SIN(c) )+&    !Rotating principal axes
                r_ppal(i,2)*(-COS(a)*COS(b)*SIN(c) - SIN(a)*COS(c) )+&
                r_ppal(i,3)*( COS(a)*SIN(b)                        )
     r_rot(i,1)=r_rot(i,1)+MiddlePoint(1)                                        !Final coordinates (x, k=1) after move and rotate
     
     r_rot(i,2)=r_ppal(i,1)*( SIN(a)*COS(b)*COS(c) + COS(a)*SIN(c) )+&
                r_ppal(i,2)*(-SIN(a)*COS(b)*SIN(c) + COS(a)*COS(c) )+&
                r_ppal(i,3)*( SIN(a)*SIN(b)                        )
     r_rot(i,2)=r_rot(i,2)+MiddlePoint(2)
     
     r_rot(i,3)=r_ppal(i,1)*(-SIN(b)*COS(c)                        )+&
                r_ppal(i,2)*( SIN(b)*SIN(c)                        )+&
                r_ppal(i,3)*( COS(b)                               )
     r_rot(i,3)=r_rot(i,3)+MiddlePoint(3)
  ENDDO

RETURN
END SUBROUTINE RotateCluster
