!----------------------------------------------
!        Drawing Molecular configuration
!
!      Written by Edison Florez using the
!      Fortran 03 version of the language
!
!            edisonffh@gmail.com
!----------------------------------------------
    
SUBROUTINE DrawMolecularConfiguration( ZNucleus,radiusSBC,NMolecules,NAtoms, &
                                       Max_Atoms,MiddlePoint,NTrajectories,r,&
                                       t,fileunit,comment_line                )
!
!  Purpose:
!    Subroutine has been writing to rewrite  molecular configuration
!    in XYZ format, i.e. a file like
!
!    <number of atoms>
!    comment line
!    <Atomic Number> <X> <Y> <Z>
!
!  Record of the revisions:
!       Date            Programmer              Description of change
!       ========        ==============          ======================================
!       Sept-2017	E. Florez               v0: It did not exist
!       Sept-2017	E. Florez		v1: build file "configuration_i.xyz" for each i step in Monte-Carlo cycle
!       Sept-2017	E. Florez              	v2: build file "all_configuration.xyz", it append every configuration into a single file
!       Sept-2017	E. Florez              	v3: The name of the XYZ file was changed by "accepted_configuration.xyz"
!       Sept-2017	E. Florez              	v4: Allocate array
!       Oct-2017	E. Florez              	v5: The middlepoint and box vertix are uniques, so they are only calculated for the initial configuration
!       Nov-2017	E. Florez              	v6: Input coordinates for rare gases (Xenon i.g.) are in atomic unit, so the output is in Angstrom
!       Nov-2017	E. Florez              	v7: The rotation of principal axes is only for fragments with more than one atom 
!       Feb-2018	E. Florez              	v8: Drawing a 'sphere' to visualize SBC
!       Jan-2019	E. Florez              	v9: Writing out the coordinates in a XYZ file for any UNIT
!
IMPLICIT NONE

!-Declare CALLING parameter & Data Dictionary
  INTEGER(4),INTENT(IN) :: NMolecules                                         !Number of Molecules 
  INTEGER(4),DIMENSION(NMolecules),INTENT(IN) :: NAtoms                       !Number of Atoms 
  INTEGER(4),INTENT(IN) :: Max_Atoms                                          !Maximun number of atoms to declare array dimension
  INTEGER(8),INTENT(IN) :: fileunit                                           !Designating a UNIT number to write out 
  CHARACTER(len=50),INTENT(IN) :: comment_line                                !Comment line 
  INTEGER(4),INTENT(IN) :: NTrajectories,t                                    !Number of temperarures
  INTEGER(4),DIMENSION(NMolecules,Max_Atoms),INTENT(IN) :: ZNucleus           !Atomic Numbers, per atom in each molecule. Zn(Molec,Atom), means rank=2
  REAL(8),DIMENSION(3),INTENT(IN) :: MiddlePoint                                        
  REAL(8),INTENT(IN) :: radiusSBC                                             !Cube's Length 
  REAL(8),DIMENSION(0:NTrajectories,NMolecules,Max_Atoms,3),INTENT(IN) :: r   !Molecular Coordinades x,y,z; r(temp,Molec,Atom,x), means rank=4
                                                                             
!-Declare LOCAL parameter & Data Dictionary                                  
  INTEGER(4) :: i,j                                                           !The indices mening are: i=MoleculeLabel, j=AtomLabel, k=Coordinates
  INTEGER(4) :: TotalAtoms                                                    !Total number of atoms in all system
  REAL(8),DIMENSION(3) :: d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,&                !coordinates of the sphere (SBC)
                          d12,d13,d14,d15,d16,d17,d18
  REAL(8) :: sin45                                                            !sin(45deg)=cos(45deg), so here we always use sin45. It used to calculate the SBC
  CHARACTER(len=*),PARAMETER :: FMT1="(I3,3(3X,F15.10))",&                     !Formats to print coordinates
                                FMT2="(A3,3(3X,F15.10))"

!--------------------------------------------------------------------------------------------
!-Body of the Subroutine

  TotalAtoms=0
  DO i=1,NMolecules                                                       !The number of atoms in the whole molecular conglomerate
     TotalAtoms=TotalAtoms+NAtoms(i)
  ENDDO 

  WRITE(fileunit,"(10X,I6)") TotalAtoms 
  WRITE(fileunit,"(A)") comment_line

  DO i=1,NMolecules
!     DO j=1,NAtoms(i)
        WRITE(fileunit,FMT1) ZNucleus(i,1),r(t,i,1,1),r(t,i,1,2),r(t,i,1,3)
!        WRITE(fileunit,FMT1) ZNucleus(i,j),r(t,i,j,1)*0.529177249,&
!                        r(t,i,j,2)*0.529177249,r(t,i,j,3)*0.529177249                  !1 bohr = 0.529177249 angstrom
!     ENDDO
  ENDDO

!d_i are the end points of a diameter for the sphere (x-h)**2+(y-j)**2+(z-k)**2=r**2, where h:MiddlePoint(x), j:MiddlePoint(y) and k:MiddlePoint(z)
!The x axis is horizontal, z is vertical and y is out. The points d_i (i=1,2,...,8) are in plane zx+j, the others are in the plane xy+k.
  sin45=1.0/SQRT(2.0) !sin(45deg)=cos(45deg), so here we always use sin45
 
  d1(1) =MiddlePoint(1)+radiusSBC       ; d1(2) =MiddlePoint(2)                 ; d1(3) =MiddlePoint(3)
  d2(1) =MiddlePoint(1)-radiusSBC       ; d2(2) =MiddlePoint(2)                 ; d2(3) =MiddlePoint(3)
  d3(1) =MiddlePoint(1)                 ; d3(2) =MiddlePoint(2)                 ; d3(3) =MiddlePoint(3)-radiusSBC
  d4(1) =MiddlePoint(1)                 ; d4(2) =MiddlePoint(2)                 ; d4(3) =MiddlePoint(3)+radiusSBC
  d5(1) =MiddlePoint(1)+radiusSBC*sin45 ; d5(2) =MiddlePoint(2)                 ; d5(3) =MiddlePoint(3)+radiusSBC*sin45
  d6(1) =MiddlePoint(1)-radiusSBC*sin45 ; d6(2) =MiddlePoint(2)                 ; d6(3) =MiddlePoint(3)-radiusSBC*sin45
  d7(1) =MiddlePoint(1)+radiusSBC*sin45 ; d7(2) =MiddlePoint(2)                 ; d7(3) =MiddlePoint(3)-radiusSBC*sin45
  d8(1) =MiddlePoint(1)-radiusSBC*sin45 ; d8(2) =MiddlePoint(2)                 ; d8(3) =MiddlePoint(3)+radiusSBC*sin45
  d9(1) =MiddlePoint(1)                 ; d9(2) =MiddlePoint(2)+radiusSBC       ; d9(3) =MiddlePoint(3)
  d10(1)=MiddlePoint(1)                 ; d10(2)=MiddlePoint(2)-radiusSBC       ; d10(3)=MiddlePoint(3)
  d11(1)=MiddlePoint(1)+radiusSBC*sin45 ; d11(2)=MiddlePoint(2)+radiusSBC*sin45 ; d11(3)=MiddlePoint(3)
  d12(1)=MiddlePoint(1)-radiusSBC*sin45 ; d12(2)=MiddlePoint(2)-radiusSBC*sin45 ; d12(3)=MiddlePoint(3)
  d13(1)=MiddlePoint(1)-radiusSBC*sin45 ; d13(2)=MiddlePoint(2)+radiusSBC*sin45 ; d13(3)=MiddlePoint(3)
  d14(1)=MiddlePoint(1)+radiusSBC*sin45 ; d14(2)=MiddlePoint(2)-radiusSBC*sin45 ; d14(3)=MiddlePoint(3)
  d15(1)=MiddlePoint(1)                 ; d15(2)=MiddlePoint(2)+radiusSBC*sin45 ; d15(3)=MiddlePoint(3)+radiusSBC*sin45
  d16(1)=MiddlePoint(1)                 ; d16(2)=MiddlePoint(2)-radiusSBC*sin45 ; d16(3)=MiddlePoint(3)-radiusSBC*sin45
  d17(1)=MiddlePoint(1)                 ; d17(2)=MiddlePoint(2)+radiusSBC*sin45 ; d17(3)=MiddlePoint(3)-radiusSBC*sin45
  d18(1)=MiddlePoint(1)                 ; d18(2)=MiddlePoint(2)-radiusSBC*sin45 ; d18(3)=MiddlePoint(3)+radiusSBC*sin45

  !coordinates of the sphere (SBC)
  WRITE(fileunit,FMT2)'X00 ',MiddlePoint
  WRITE(fileunit,FMT2)'X01 ',d1
  WRITE(fileunit,FMT2)'X02 ',d2
  WRITE(fileunit,FMT2)'X03 ',d3
  WRITE(fileunit,FMT2)'X04 ',d4
  WRITE(fileunit,FMT2)'X05 ',d5
  WRITE(fileunit,FMT2)'X06 ',d6
  WRITE(fileunit,FMT2)'X07 ',d7
  WRITE(fileunit,FMT2)'X08 ',d8
  WRITE(fileunit,FMT2)'X09 ',d9
  WRITE(fileunit,FMT2)'X10 ',d10
  WRITE(fileunit,FMT2)'X11 ',d11
  WRITE(fileunit,FMT2)'X12 ',d12
  WRITE(fileunit,FMT2)'X13 ',d13
  WRITE(fileunit,FMT2)'X14 ',d14
  WRITE(fileunit,FMT2)'X15 ',d15
  WRITE(fileunit,FMT2)'X16 ',d16
  WRITE(fileunit,FMT2)'X17 ',d17
  WRITE(fileunit,FMT2)'X18 ',d18

  RETURN
END SUBROUTINE DrawMolecularConfiguration
