!----------------------------------------------
!    Molecular Displacement of CM & Rotation
!          for a Selectec Fragment
!
!      Written by Edison Florez using the
!      Fortran 03 version of the language
!
!            edisonffh@gmail.com
!----------------------------------------------

SUBROUTINE DisplaceAtoms( r,NMolecules,maxDisplacement,SelectedMolecule, &   !IN 
                          Seed,radiusSBC2,MiddlePoint,NTrajectories,t,   &   !IN
                          r_move                                          )  !OUT

!
!  Purpose:
!    Subroutine to dispace atoms randomly
!
!  Record of the revisions:
!       Date            Programmer              Description of change
!       ========        ==============          ======================================
!       Nov-2018	E. Florez		v1: Original code
!       xxx-2018	E. Florez		v2:
!
IMPLICIT NONE

!-Declare CALLING parameter & Data Dictionary
  INTEGER(4),INTENT(IN) :: SelectedMolecule                                   !Molecule to Displace & Rotate. It can be iqual to sm or SM
  INTEGER(4),INTENT(IN) :: NMolecules                                         !Number of Molecules 
  INTEGER(4),INTENT(IN) :: NTrajectories,t                                    !Number of temperarures
  REAL(8),INTENT(IN) :: radiusSBC2                                            !radius of SBD squared
  REAL(8),DIMENSION(3),INTENT(IN) :: MiddlePoint
  REAL(8),INTENT(IN) :: maxDisplacement                                       !Displacement could be adjusted in the fly
  INTEGER(8),INTENT(IN) :: Seed                                               !The Seed for the Random Number generator function based on clock system
  REAL(8),DIMENSION(0:NTrajectories,NMolecules,1,3),INTENT(IN) :: r           !Molecular Coordinades x,y,z; r(Molec,Atom,x), means rank=3

  REAL(8),DIMENSION(3),INTENT(OUT) :: r_move                                  !New Molecular Coordinades x,y,z for selected molecule 
                                                                            
!-Declare LOCAL parameter & Data Dictionary                                 
  INTEGER(8) :: i,j,k,sm 
  REAL(8),DIMENSION(3) :: Displacement                                        !Random Displacement in x,y,z
  INTEGER(8) :: DisplaCounter,Counter                                         !Counter to know how many steps are necessaries to avoid atoms overlapping
  REAL(8) r_cut                                                               !Length to define the boundary conditions
  REAL(8),DIMENSION(:),ALLOCATABLE :: dist1,Rmin_covalent                     !Euclidean distance
  REAL(8) :: RandomNumber                                                     !Random Number generator function
  CHARACTER(len=10) :: DisplaCounter_char                                     !Printing Left-Justified Values, converting integer to character
  INTEGER(8) :: Overlapping                                                 

  INTEGER(8) :: fileunit                                                      !designating a UNIT number to write out 
  CHARACTER(len=50) :: comment_line                                           !Comment line for XYZ file

  INTEGER(4) :: Max_Atoms
  INTEGER(4),DIMENSION(NMolecules) :: NAtoms                       !Number of Atoms 
  INTEGER(4),DIMENSION(NMolecules,1) :: ZNucleus                   !Atomic Numbers, per atom in each molecule. Zn(Molec,Atom), means rank=2
  REAL(8) :: radiusSBC

!--------------------------------------------------------------------------------------------
!-Allocating variables
  ALLOCATE(dist1(NMolecules))                
  ALLOCATE(Rmin_covalent(NMolecules))       

!--------------------------------------------------------------------------------------------
!-Body of the procedure  
  DisplaCounter=0
  DO                                                                         !Avoiding atoms overlapping
     IF ( DisplaCounter < 100000 ) THEN                                      !1E5 is the max number of movements allowed to find a good displacement

        !Abbreviation
        sm=SelectedMolecule

        DO k=1,3
           r_move(k)=r(t,sm,1,k)
        ENDDO

        !-Displacement of CM inside boundary condition (spheric or cubic)
        Counter=0
        DO 
           DO k=1,3                                                              !Defining a new Rcm and atomic coord (r) to move in term of the old ones
!              Displacement(k)=(2.0*RandomNumber(Seed)-1.0)*maxDisplacement
              Displacement(k)=(RandomNumber(Seed)-0.5)*maxDisplacement           !according to Florent's code
              r_move(k)=r_move(k)+Displacement(k)
           ENDDO
          
           !Confined inside a Sphere, where the center (h,j,k), so (x-h)^2 + (y-j)^2 + (z-k)^2 = r^2 is the general equation.
           !h=MiddlePoint(1) ; j=MiddlePoint(2) ; k=MiddlePoint(3) & r=radiusSBC
           r_cut=( (r_move(1)-MiddlePoint(1))**2 + (r_move(2)-MiddlePoint(2))**2 + &                !(x-h)^2 + (y-j)^2 + (z-k)^2
                   (r_move(3)-MiddlePoint(3))**2  )

           IF ( r_cut < radiusSBC2 ) THEN                           !r_cut is not SQRT and radiusSBC2 is squared
              !Exiting from the main Do cycle
              EXIT
           ELSE
              !Recovering the old coordinates
              DO k=1,3
                 r_move(k)=r_move(k)-Displacement(k)
              ENDDO
              Counter=Counter+1
              IF ( Counter > 100000 ) THEN
                 WRITE(*,*)''
                 WRITE(*,*)''
                 WRITE(*,"(A)")' *** Warning! ***'
                 WRITE(DisplaCounter_char,"(I6)") Counter-1
                 WRITE(*,"(3A)")'More than ',TRIM(ADJUSTL(DisplaCounter_char)),&
                 ' movements of CM (OUT OF SBC) and no good configuration was found'
                 WRITE(*,*)''


                 OPEN(UNIT=3000,FILE="last_configuration.xyz",STATUS='REPLACE')
                 
                 comment_line='Last configuration, ATOMS OUT OF SBC'
                 fileunit=3000

                 DO i=1,NMolecules
                    ZNucleus(i,1)=10
                    NAtoms(i)=1 
                 ENDDO

                 Max_Atoms=1  
                 radiusSBC=SQRT(radiusSBC2)
                 
                 CALL DrawMolecularConfiguration( ZNucleus,radiusSBC,NMolecules,NAtoms, &
                                                  Max_Atoms,MiddlePoint,NTrajectories,r,&
                                                  t,fileunit,comment_line                 )
                 CLOSE(3000)

                 STOP
              ENDIF
           ENDIF
        ENDDO
 
        !Calculating distances
        Overlapping=0                                             !Initialization
        DO i=1,NMolecules

           IF ( i == sm ) CYCLE                                   !Avoiding self-interaction

           !Distance between sm and the others          
           dist1(i)=SQRT(                               &         !dist1_ab=sqrt[
                          ( r_move(1) - r(t,i,1,1) )**2+&         !              (xa-xb)^2+
                          ( r_move(2) - r(t,i,1,2) )**2+&         !              (ya-yb)^2+
                          ( r_move(3) - r(t,i,1,3) )**2   )       !              (za-zb)^2 ]

           Rmin_covalent(i)=0.0                                   !The minimum distance is teh covalent radius for HH, it is 0.32 Angstrom
           IF( dist1(i) < Rmin_covalent(i) )  Overlapping=Overlapping+1   !Overlapping could be TRUE (1) or FALSE (0)
        ENDDO

    !Exiting from the main Do cycle
    IF ( Overlapping == 0 ) EXIT !Avoiding atoms overlapping 

     ELSE  
        WRITE(*,*)''
        WRITE(*,*)''
        WRITE(*,"(A)")' *** Warning! ***'
        WRITE(DisplaCounter_char,"(I6)") DisplaCounter-1
        WRITE(*,"(3A)")'More than ',TRIM(ADJUSTL(DisplaCounter_char)),&
                        'movements of CM (OVERLAPPING) and no good configuration was found'
        WRITE(*,*)''

        OPEN(UNIT=3000,FILE="last_configuration.xyz",STATUS='REPLACE')

        comment_line='Last configuration, OVERLAPPING ATOMS'
        fileunit=3000

        DO i=1,NMolecules
           ZNucleus(i,1)=10
           NAtoms(i)=1
        ENDDO

        Max_Atoms=1
        radiusSBC=SQRT(radiusSBC2)

        CALL DrawMolecularConfiguration( ZNucleus,radiusSBC,NMolecules,NAtoms, &
                                         Max_Atoms,MiddlePoint,NTrajectories,r,&
                                         t,fileunit,comment_line                 )
        CLOSE(3000)


        STOP
     ENDIF

     DisplaCounter=DisplaCounter+1

  ENDDO !Atoms overlapping

RETURN
END SUBROUTINE DisplaceAtoms
