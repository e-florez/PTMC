!----------------------------------------------
!      Written by Edison Florez using the
!      Fortran 03 version of the language
!
!            edisonffh@gmail.com
!----------------------------------------------

PROGRAM atomic_PTMC
!
!  Abstract:
!    This is a Parallel Tempering Monte Carlo (PTMC) a fortran program
!    to calculate themodinamical properties (canonical ensemble) displacing 
!    and rotating a random molecule, as a fixed set, of a molecular 
!    conglomerate, according the Markov chain. Each link in the Markov 
!    chain has a energy (Ei), so the program changes from Ei to Ej and 
!    the configuration j is accepted whether Ej<Ei or the Metropolis 
!    criterion is met, otherwise configuration j is rejected
!
!  Record of the revisions:
!       Date		Programmer		Description of change
!	========	==============   	======================================
!	Nov-2018	E. Florez		v1: Original code
!	Mar-2019	E. Florez		All the version is controlled by Git --> git@bitbucket.org:edisonflorez/atomic_ptmc_checkpoints.git
!
!  FORTRAN UNITS:
!        unit            FILE                    Comment
!        ========        ===============         ======================================
!  before
!        100             input FILE
!        300             initial config
!        200             restarting INPUT
!        400             energy function
!        500+t           config saved            saving configs after certain number of MC cycles (one per trajectory t)
!  MC
!        1000            status
!        1200            lowest energy config    
!  after
!        2000            histE.data              general info for the multihistogram analysis 
!        2100+t          histE                   energy histogram information (one per trajectory t)
!        2300+t          RDFhist                 global radial distribution analysis (one per trajectory t)
!        2400+t          RDFhist_x               X: radial distribution analysis (one per trajectory t)
!        2600+t          RDFhist_y               Y: radial distribution analysis (one per trajectory t)
!        2800+t          RDFhist_z               Z: radial distribution analysis (one per trajectory t)
!        10000+i         outout FILE             10000 is the global OUTPUT and 10000+i are checkpoints' output FILE (where i>0)
!        20000+i         restarting OUTPUT
!
!  Displ
!        3000            last config             last configuration after failing down (atoms overlapping or moved out of SBC)
!
! ** units '9000-9999' for TESTING purposes  
!
!----------------------------------------------

USE OMP_LIB

IMPLICIT NONE

!-Gobal Variables: Declaration & Data Dictionary
  INTEGER(8) :: SelectedMolecule,i,j,k,n                             !The indices mening are: i=MoleculeLabel, j=AtomLabel, k=Coordinates, n=dummy
  INTEGER(8) :: MCcycles,MCstep,Nequilibration                       !Number of MC cycles and step, "MCstep" means a single step.
  REAL(8),PARAMETER:: PI=3.141592653589793238462643D0
  INTEGER(4) :: NMolecules                                           !Number of Molecules
  INTEGER(4),ALLOCATABLE :: NAtoms(:)                                !Number of Atoms per molecule array(NMolecules)
  INTEGER(4) :: Max_Atoms                                            !Maximun number of atoms to declare array dimension
  INTEGER(4) :: TotalAtoms                                           !Total number of atoms in all system
  CHARACTER(len=20),DIMENSION(:),ALLOCATABLE :: MoleculeName         !Name of Atom/Molecule fragments
  INTEGER(4),DIMENSION(:,:),ALLOCATABLE :: ZNucleus                  !Atomic Number, per atom in each molecule. Zn(Molec,Atom). INT kind 8 index start from 0
  REAL(8) :: radiusSBC,radiusSBC2                                    !Radius for Spherical Boundary Conditions SQUARED 
  REAL(8),DIMENSION(:,:,:),ALLOCATABLE :: ro                         !Molecular Coordinades x,y,z for the starting config. r(Molec,Atom,x), means rank=3
  REAL(8),DIMENSION(:,:,:,:),ALLOCATABLE :: rf,r                     !Molecular Coordinades x,y,z for every temperature. r(Tempe,Molec,Atom,x), means rank=4
  REAL(8),DIMENSION(:,:,:,:),ALLOCATABLE :: rex                      !rex is a temporary variable to exchange trajectories, so rex=rf(j); rf(j)=rf(j+1); rf(j+1)=rex
  REAL(8),DIMENSION(3) :: r_move                                     !Molecular Coordinades x,y,z for Selected Atom 
  REAL(8),DIMENSION(:,:),ALLOCATABLE :: r_rot                        !Rotated coordinades x,y,z for the whole cluster
  REAL(8) :: MaxDisplacementRcm,MaxAngle                             !Max Displacement of Mass Center & Max Angle of Rotation
  REAL(8),DIMENSION(:),ALLOCATABLE :: maxDisplacement,MaxRotation    !Displacement of the CM and ROt of the whole cluster is a fuction of T and it is adjustable in the fly
  REAL(8),DIMENSION(3) :: MiddlePoint                                !Coordinates of the CM and the Middle point of the box
  REAL(8),DIMENSION(:,:),ALLOCATABLE :: Rcm                          !Coordinates of the CM, Rcm(Molec,x) Rcm(Molec,y) Rcm(Molec,z)
  REAL(8),DIMENSION(3) :: RcmMinimum,RcmMaximum                      !Maximum & Minimum according Center of Masses to calculate Meiddle point
  REAL(8),DIMENSION(:),ALLOCATABLE :: TotalWeight                    !Sum of all Atomic Masses Numbers for each specific fragment
  REAL(8),PARAMETER :: kB=1.38065040E-23                             !kB is 1.38065040E-23 J/K, 0.0083144621(75) kJ/(mol⋅K) or 3.16499D-6 Hartrees
  INTEGER(8) :: NTrajectories,t                                      !Number of Temperatures
  REAL(8) :: Tf,To                                                   !Maximun and minimum temperatures respectively
  INTEGER(8) :: TotalAccepted_Mov,TotalAccepted2_Mov,&               !Total accepted and exchanged configuration 
                TotalAccepted_Rot,TotalAccepted2_Rot,TotalExchanged
  REAL(8) :: acceptance,Adjust                                       !Metropolis acceptance criterion and adjusted step size
  INTEGER(8),DIMENSION(:),ALLOCATABLE :: Mov_accepted_eq,&
                                         Rot_accepted_eq,&           !Number of confuguration accepted per trajectory
                                         Mov_accepted,Rot_accepted   !After equilibration
  INTEGER(8),DIMENSION(:),ALLOCATABLE :: adj_Mov_accepted,&          !Acceptance ration must be in [0.4,0.6] we adjust it on the fly
                                         adj_Rot_accepted
  INTEGER(8),DIMENSION(:),ALLOCATABLE :: adj_Mov_accepted_eq,&       !before eq
                                         adj_Rot_accepted_eq
  CHARACTER(len=50) :: EnergyName                                    !Name of the Potential Energy procedure 
  CHARACTER(len=500) :: EnergyFullPath                                    !Full path Name of the Potential Energy procedure 
  integer :: last_slash_position
!  REAL(8) :: EnergyFunction                                          !The procedure to calculate the Energy
  REAL(8),DIMENSION(:),ALLOCATABLE :: Energy_Tot,Energy_Tot_fin      !The total energy for the configuration
  REAL(8) :: DeltaEnergy                                             !Delta Energy is "Energy_Tot_ini(t) - Energy_Tot_fin(t)"
  REAL(8),DIMENSION(:,:,:),ALLOCATABLE :: Energy_ab_ini,Energy_ab_fin
  REAL(8),DIMENSION(:,:),ALLOCATABLE :: Energy_2a_ini,Energy_2a_fin
  INTEGER(8) :: Seed2,Seed,ClockSeed                                       !The Seed for the Random Number generator function based on clock system
  INTEGER(8) :: Rotate                                               !rotate the whole cluster, Yes: >0 / No: <=0
  INTEGER(8) :: input_error,ierror                                   !OPEN "input": Status flag from I/O statements. And, checking input type 
  INTEGER(8) :: FILEunit                                             !designating a UNIT number to write out 
  CHARACTER(len=100) :: comment_line                                  !Comment line for XYZ FILE
  CHARACTER(len=100) :: energy_tot_char                               !converting real to character
  REAL(8) :: RandomNumber                                            !Random Number generator function
  REAL(8) :: sin45                                                   !sin(45deg)=cos(45deg), so here we always use sin45. It used to calculate the SBC
  REAL(8),DIMENSION(3) :: d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,&       !coordinates of the sphere (SBC)
                          d12,d13,d14,d15,d16,d17,d18
  REAL(8),DIMENSION(:),ALLOCATABLE :: Temperature,beta               !Temperature(t) is the temperature to evaluate and beta is iqual to 1/(kB*T)
  REAL(8),DIMENSION(:),ALLOCATABLE :: TotalEnergy,TotalEnergy2,&
                                      AverageEnergy,AverageEnergy2,&
                                      Energy_Average,Energy_Average2,&
                                      Energy_Average3,Cv,dCv
  REAL(8),DIMENSION(:),ALLOCATABLE :: MaxEnergy,MinEnergy            !Local Energy Max & Min (per temperature) for the  histogram analysis
  REAL(8) :: MaxEnergyTotal,MinEnergyTotal                           !Global Energy Max & Min
  REAL(8),PARAMETER :: JtokJmol=6.022140857E23/1.0E3                 !From J to kJ/mol, NA=6.022140857E23
  REAL(8),PARAMETER :: EhtoJ=4.35974465E-18                          !from Hartree to J, 1 hartree = 4.359744650E-18 J
  INTEGER(8) :: nBins,hit                                            !number of bins for the energy range and 'hit' is teh index for the frequency (histogram)
  INTEGER(8) :: nConfigELT,nConfigEGT                                !Counting the energies OUTSIDE the range (Emax-Emin) for the histogram
  REAL(8) :: DeltaEHistogram                                         !It is the 'Bin width'
  INTEGER(8),DIMENSION(:,:),ALLOCATABLE :: nhits                     !frequency counter. 1st index is for the temperature and the 2nd one is for 'hit'
  INTEGER(8),DIMENSION(:),ALLOCATABLE :: Exchange_tot,Exchange_acc   !counting the total and accepted exchange
  INTEGER(8) :: ExchangeTrajectory                                   !index to exchange adjacent trajectories
  REAL(8) :: DeltaBeta,ExchangeAccept,Energy_ex                      !DeltaBet=beta(i)-beta(i+1) and ExchangeAcc=exp( DeltaBeta*DeltaEnergy ). Energy_ex is a temporal var
  REAL(8) :: Energy_2a_ini_ex,Energy_ab_ini_ex 
  CHARACTER(len=40) :: MCcycles_char,MCstep_char,&                   !Printing Left-Justified Values, converting integer to character
                       Max_Atoms_char,NMolecules_char,&
                       t_char,Nequilibration_char
  CHARACTER(len=20),DIMENSION(:),ALLOCATABLE :: NAtoms_char 
  REAL(8) :: AtomicMass                                              !Atomic Mass Number
  CHARACTER(len=50) :: input_name                                    !Name of the input FILE using a max of 40 character
  CHARACTER(len=500) :: path                                    !Name of the input FILE using a max of 40 character
  CHARACTER(len=:),ALLOCATABLE :: input,output                       !Name of the input & output FILEs
  LOGICAL,DIMENSION(:),ALLOCATABLE :: accept                         !Is it this config accepted? .TRUE. or .FALSE.
  INTEGER(8) :: SaveConfig                                           !Saves configuration every certain number of MC cycles
  INTEGER(4) :: days,hours,minutes                                   !Running Time
  REAL(8) :: StartTime,EndTime,RunTime,seconds,&                     !Running time
             PartialStartRunTime,PartialEndRunTime
  CHARACTER(10) :: today,now
  INTEGER(8) :: Molecule_a,Molecule_b,read_error                     !indices for molecule A and B
  INTEGER(8) :: indexAngle,indexDistance,NumDistances,NumAngles
  INTEGER(8) :: indexAngle_ini,indexDistance_ini
  INTEGER(8) :: indexAngle_fin,indexDistance_fin
  REAL(8) :: MinDistance,SpacingDistances,MinAngle,SpacingAngles
  REAL(8) :: dist,theta
  REAL(8) :: dist_ini,theta_ini,dist_fin,theta_fin
  REAL(8),DIMENSION(:,:),ALLOCATABLE :: Energy_LookUpTable
  !RDF histogram
  INTEGER(8) :: indexRDF,j_a,j_b,BinsRDF                              !The indices menans for molecule A and B. the index j is for atoms inside A or B
  REAL(8) :: MinR,MaxR,DeltaR 
  REAL(8) :: distx,disty,distz
  INTEGER(8),DIMENSION(:,:),ALLOCATABLE :: RDFHistogram                 !Radial Distribution Function. The number of bins for the RDF hist is 1000
  INTEGER(8),DIMENSION(:,:),ALLOCATABLE :: RDFHistogramx
  INTEGER(8),DIMENSION(:,:),ALLOCATABLE :: RDFHistogramy
  INTEGER(8),DIMENSION(:,:),ALLOCATABLE :: RDFHistogramz
  CHARACTER(len=*),PARAMETER :: FMT1="(I3,3(3X,F15.10))",&           !Formats to print coordinates
                                FMT2="(A3,3(3X,F10.6))",&
                                FMT3="(I3,3(3X,F30.26))"
  ! === Restart ===
  !Checkpoint, Restarting option
  INTEGER(8) :: check,check2                            !counter and total number of point respectively
  INTEGER(8) :: Ncheckpoints_eq,Ncheckpoints_mc                            !counter and total number of point respectively
  INTEGER(8) :: check_split_eq,check_split_mc,&
                check_MCstep,check_MCcycles
  CHARACTER(len=40) :: check_char,Ncheckpoints_mc_char,&
                       Ncheckpoints_eq_char,check_MCcycles_char
  !checking FILE/directory exists and can be opened
  LOGICAL :: FILExist,direxist
 
  !input
  CHARACTER(len=10) :: restart
  CHARACTER(len=50) :: restart_FILE
  !read
  CHARACTER(len=200) :: string         
  INTEGER(8) :: MCcycles_restart,Nequilibration_restart,&
                current_restart,Seed_restart,Ncheckpoints_mc_restart,&
                Ncheckpoints_eq_restart,check_restart_eq,check_restart_mc
  REAL(8) :: eq_percent_restart,sampling_percent_restart
  

  
  !=== Parallelization == openMP
  REAL(8) :: wtime_start,wtime_end,part_wtime_start,part_wtime_end
  INTEGER(8) :: threads,aval_threads

!--------------------------------------------------------------------------------------------
!-Start time
  CALL DATE_AND_TIME(today,now)

  CALL CPU_TIME(StartTime)
  PartialStartRunTime = StartTime

  wtime_start = 0.0
  wtime_end = 0.0
  threads = 1
  aval_threads = 1

  !$ wtime_start = omp_get_wtime()
  !$ part_wtime_start = wtime_start
  !$ threads = omp_get_max_threads()
  !$ aval_threads = omp_get_num_procs()

!--------------------------------------------------------------------------------------------
!-Read input
  !Get the path of the executable
  CALL GET_COMMAND_ARGUMENT(0,path)
  last_slash_position = index(path, '/', back=.true.)

    if (last_slash_position > 0) then
        EnergyFullPath = path(1:last_slash_position)
    else
        EnergyFullPath = './'
    endif

  !Asign the value to input based on "input" FILE
  CALL GET_COMMAND_ARGUMENT(1,input_name)

  !Automatic allocation on assignment
  input=TRIM(ADJUSTL(input_name))

  !Input FILE has the extension '.in', so we remove its, i.e. the last 3 characters
  output=input(1:LEN(input)-3)

  OPEN (UNIT=100,FILE=input,STATUS='OLD',ACTION='READ',IOSTAT=input_error)   !Was the OPEN successful?
  IF ( input_error == 0 ) THEN
     READ(100,*,iostat=ierror) MCcycles,Ncheckpoints_mc
                IF ( MCcycles*Ncheckpoints_mc <= 0 .OR. ierror /= 0 ) THEN
                   WRITE(*,*)
                   WRITE(*,*)' *** ERROR ***'
                   WRITE(*,*)' Number of Cycles and Checkpoints for MC sampling'
                   WRITE(*,*)' must be positive integers'
                   WRITE(*,*)
                   STOP
                ENDIF

     READ(100,*,iostat=ierror) Nequilibration,Ncheckpoints_eq  
                IF ( Nequilibration*Ncheckpoints_eq <= 0 .OR. ierror /= 0 ) THEN
                   WRITE(*,*)
                   WRITE(*,*)' *** ERROR ***'
                   WRITE(*,*)' Number of Cycles and Checkpoints for Equilibration '
                   WRITE(*,*)' must be positive integers'
                   WRITE(*,*)
                   STOP
                ENDIF

     READ(100,*,iostat=ierror) To,Tf,NTrajectories
                IF ( NTrajectories < 2 .OR. ierror /= 0 ) THEN
                   WRITE(*,*) 
                   WRITE(*,*)' *** ERROR ***'
                   WRITE(*,*)' Number of Trajectories must be a positive '
                   WRITE(*,*)' integer lager than two' 
                   WRITE(*,*) 
                   STOP
                ENDIF

     READ(100,*,iostat=ierror) Seed                                               !The Seed for the Random Number generator
                IF ( ierror /= 0 ) THEN
                   WRITE(*,*)
                   WRITE(*,*)' *** ERROR ***'
                   WRITE(*,*)' Number of Cycles and Checkpoints for Equilibration '
                   WRITE(*,*)' must be positive integer, larger than zero'
                   WRITE(*,*)
                   STOP
                ENDIF
                IF ( Seed <= 0 ) Seed=ClockSeed()                  !it is based on the clock system

     READ(100,*) EnergyName
     EnergyFullPath=TRIM(ADJUSTL(EnergyFullPath))//"/EnergyDataBase/"//TRIM(ADJUSTL(EnergyName))

     READ(100,*,iostat=ierror) restart,restart_FILE
                IF ( restart == 'yes' ) THEN
                   OPEN (UNIT=200,FILE=TRIM(ADJUSTL(restart_FILE)),&
                         STATUS='OLD',ACTION='READ',IOSTAT=input_error)   !Was the OPEN successful?

                   IF ( ierror /= 0 .OR. input_error /= 0) THEN
                      WRITE(*,*)
                      WRITE(*,*)' *** ERROR ***'
                      WRITE(*,*)' Check restart settings'
                      WRITE(*,*)' Restart:',restart
                      WRITE(*,*)' File NOT FOUND: ',restart_FILE
                      WRITE(*,*)
                      STOP
                   ENDIF
                ENDIF
     READ(100,*,iostat=ierror) radiusSBC
                IF ( ierror /= 0 .OR. radiusSBC <= 0.0 ) THEN
                   WRITE(*,*)
                   WRITE(*,*)' *** ERROR ***'
                   WRITE(*,*)' Sphere radius for Boundary Conditions' 
                   WRITE(*,*)' must be positive float, larger than zero'
                   WRITE(*,*)
                   STOP
                ENDIF
                radiusSBC2=radiusSBC**2                            !Radius for SBC squared
     READ(100,*,iostat=ierror) NMolecules
                IF ( ierror /= 0 .OR. NMolecules <= 0.0 ) THEN
                   WRITE(*,*)
                   WRITE(*,*)' *** ERROR ***'
                   WRITE(*,*)' Number of Molecules must be a positive integer larger than zero'
                   WRITE(*,*)
                   STOP
                ENDIF
     READ(100,*,iostat=ierror) SaveConfig                                         !Save configs? Yes: SaveConfig>0 or No: SaveConfig<=0
                IF ( ierror /= 0 ) THEN
                   WRITE(*,*)
                   WRITE(*,*)' *** ERROR ***'
                   WRITE(*,*)' Number of cycles to saves configuration must be an integer'
                   WRITE(*,*)' Yes: >0 or No: <=0 '
                   WRITE(*,*)
                   STOP
                ENDIF
                IF ( SaveConfig < 0 ) SaveConfig = 0
  
     READ(100,*,iostat=ierror) MaxDisplacementRcm,MaxAngle
                IF ( ierror /= 0 ) THEN
                   WRITE(*,*)
                   WRITE(*,*)' *** ERROR ***'
                   WRITE(*,*)' Max displacement of the Centre of Mass [Angstrom] and'
                   WRITE(*,*)' Max. angle of Rotation [degrees] must be floats'
                   WRITE(*,*)
                   STOP
                ENDIF
                 !to rotate randomly the whole cluster
                 IF ( MaxAngle > 1.0 ) THEN
                    Rotate = 1   !rotate
                 ELSE
                    Rotate = 0   !NOT rotate
                 ENDIF
                !To convert degrees to radians. Trigonometric Fortran functions use radians      
                MaxAngle=MaxAngle*(PI/180.0)
     READ(100,*) 
   
     !Allocating variables to use immediately
     ALLOCATE(MoleculeName(NMolecules))
     ALLOCATE(NAtoms(NMolecules))
                   
     !It's only to account the max number of atoms in a fragment/molecule  
     Max_Atoms=0                          
     DO i=1,NMolecules
        READ(100,*) MoleculeName(i),NAtoms(i)
        Max_Atoms=MAX(Max_Atoms,NAtoms(i)) 
        DO j=1,NAtoms(i)
           READ(100,*) !ZNucleus(i,j),ro(i,j,1),ro(i,j,2),ro(i,j,3)      !They are not reading now 'cause we need to allocate them, so we need Max_Atom
        ENDDO
     ENDDO
   

     !Go back to initial point. Rewind FILE, only, to read ZNucleus and coordinates
     REWIND(UNIT=100)
   
     READ(100,*) !MCcycles                            !Already read
     READ(100,*) !Nequilibration,Ncheckpoints_eq      !Already read
     READ(100,*) !To,Tf,NTrajectories                 !Already read
     READ(100,*) !Seed                                !Already read 
     READ(100,*) !EnergyName                          !Already read
     READ(100,*) !restart                             !Already read
     READ(100,*) !radiusSBC                           !Already read
     READ(100,*) !NMolecules                          !Already read
     READ(100,*) !SaveConfig                          !Already read
     READ(100,*) !MaxDisplacementRcm,MaxAngle         !Already read
     READ(100,*) 
   
     !Allocating variables to use immediately
     ALLOCATE(ZNucleus(NMolecules,Max_Atoms))
     ALLOCATE(ro(NMolecules,Max_Atoms,3))            !ro(Molec,Atom,x) is just for the starting configuration
   
     DO i=1,NMolecules
        READ(100,*) !MoleculeName(i),NAtoms(i)        !Already read
        DO j=1,NAtoms(i)
           READ(100,*) ZNucleus(i,j),ro(i,j,1),ro(i,j,2),ro(i,j,3)  !starting configuration (ro)
        ENDDO
     ENDDO
   
  ELSE
     !Else FILE open failed. Tell user.
     WRITE (*,"(A)") 
     WRITE (*,"(A)")   ' *** ERROR *** '
     WRITE (*,"(2A)")' INPUT FILE open failed, check: ',input
     WRITE (*,"(A)") 
     STOP
  END IF
   
  CLOSE (UNIT=100)

!-The number of bins for the histogram by default for all systems. It is define here for allocate a nhits 
  nBins=500
  BinsRDF=100


!--------------------------------------------------------------------------------------------
!Renaming and Saving old result FILEs into a directory called 'result_yyyymmdd-hhmmss'
  
  !execution date and time: yyyymmdd-hhmmss
  string = today(1:4)//today(5:6)//today(7:8)//"-"//now(1:2)//now(3:4)//now(5:6)

  INQUIRE (FILE="./radial_distribution/",EXIST=direxist)
  IF ( direxist ) CALL SYSTEM ("mkdir -p result_"//string// &
                               " && mv -fu radial_distribution/ result_"//string )

  INQUIRE (FILE="./configurations/",EXIST=direxist)
  IF ( direxist ) CALL SYSTEM ("mkdir -p result_"//string// &
                               " && mv -fu configurations/ result_"//string )

  INQUIRE (FILE="./histograms/",EXIST=direxist)
  IF ( direxist ) CALL SYSTEM ("mkdir -p result_"//string// &
                               " && mv -fu histograms/ result_"//string )

  INQUIRE (FILE="./checkpoint/",EXIST=direxist)
  IF ( direxist ) CALL SYSTEM ("mkdir -p result_"//string// &
                               " && mv -fu checkpoint/ result_"//string )

  INQUIRE (FILE=output//'.out',EXIST=filexist)
  IF ( filexist ) CALL SYSTEM ("mkdir -p result_"//string// &
                               " && mv -fu "//output//".out"//" result_"//string )

  INQUIRE (FILE="status_"//output,EXIST=filexist)
  IF ( filexist ) CALL SYSTEM ("mkdir -p result_"//string// &
                               " && mv -fu status_"//output//" result_"//string )

  INQUIRE (FILE="initial_configuration.xyz",EXIST=filexist)
  IF ( filexist ) CALL SYSTEM ("mkdir -p result_"//string// &
                               " && mv -fu initial_configuration.xyz result_"//string )

  INQUIRE (FILE="configuration_1.xyz",EXIST=filexist)
  IF ( filexist ) CALL SYSTEM ("mkdir -p result_"//string// &
                               " && mv -fu configuration_*.xyz result_"//string )

  INQUIRE (FILE="lowest_Energy.xyz",EXIST=filexist)
  IF ( filexist ) CALL SYSTEM ("mkdir -p result_"//string// &
                               " && mv -fu lowest_Energy.xyz result_"//string )

  INQUIRE (FILE="histE.data",EXIST=filexist)
  IF ( filexist ) CALL SYSTEM ("mkdir -p result_"//string// &
                               " && mv -fu histE.* result_"//string )

  INQUIRE (FILE="RDFhist.1",EXIST=filexist)
  IF ( filexist ) CALL SYSTEM ("mkdir -p result_"//string// &
                               " && mv -fu RDFhist* result_"//string )


!--------------------------------------------------------------------------------------------
!-Allocating variables. Some of them were allocated previously
  ALLOCATE(NAtoms_char(NMolecules))
  ALLOCATE(Energy_ab_ini(0:NTrajectories,NMolecules,NMolecules))       !Pairs energy for the selected unit BEFORE it is moved randomly
  ALLOCATE(Energy_ab_fin(NTrajectories,NMolecules,NMolecules))         !Pairs energy for the selected unit AFTER it is moved randomly
  ALLOCATE(Energy_2a_ini(0:NTrajectories,NMolecules))                 
  ALLOCATE(Energy_2a_fin(NTrajectories,NMolecules))                 
  ALLOCATE(Energy_Tot(0:NTrajectories))                                !Total Energy, 0 means the startig configuration
  ALLOCATE(Energy_Tot_fin(0:NTrajectories))                            !Total Energy after rotation, 0 means the startig configuration
  ALLOCATE(Mov_accepted_eq(NTrajectories))
  ALLOCATE(Rot_accepted_eq(NTrajectories))
  ALLOCATE(Mov_accepted(NTrajectories))              !after eq
  ALLOCATE(Rot_accepted(NTrajectories))              !after eq
  ALLOCATE(adj_Mov_accepted(NTrajectories))              !after eq
  ALLOCATE(adj_Rot_accepted(NTrajectories))              !after eq
  ALLOCATE(adj_Mov_accepted_eq(NTrajectories))
  ALLOCATE(adj_Rot_accepted_eq(NTrajectories))
  ALLOCATE(beta(NTrajectories))
  ALLOCATE(maxDisplacement(NTrajectories))
  ALLOCATE(MaxRotation(NTrajectories))
  ALLOCATE(Temperature(NTrajectories))
  ALLOCATE(TotalEnergy(NTrajectories))
  ALLOCATE(TotalEnergy2(NTrajectories))
  ALLOCATE(AverageEnergy(NTrajectories))
  ALLOCATE(AverageEnergy2(NTrajectories))
  ALLOCATE(Cv(NTrajectories))
  ALLOCATE(Energy_Average(NTrajectories))
  ALLOCATE(Energy_Average2(NTrajectories))
  ALLOCATE(Energy_Average3(NTrajectories))
  ALLOCATE(dCv(NTrajectories))
  ALLOCATE(nhits(NTrajectories,0:nBins))                  !frequency counter for the histogram. 1st index is for the temperature and the 2nd one is for 'hit'
  ALLOCATE(Exchange_tot(NTrajectories))
  ALLOCATE(Exchange_acc(NTrajectories))
  ALLOCATE(MaxEnergy(NTrajectories))
  ALLOCATE(MinEnergy(NTrajectories))
  ALLOCATE(rf(0:NTrajectories,NMolecules,Max_Atoms,3))  !coordinates for accepted configuration for every temp.
  ALLOCATE(r(0:NTrajectories,NMolecules,Max_Atoms,3))   !coordinates for any configuration accepted or not for every temp
  ALLOCATE(rex(0:NTrajectories,NMolecules,Max_Atoms,3)) !Parallel Tempering: rex is a temporary variable for any configuration accepted or not
  ALLOCATE(Rcm(NMolecules,3))
  ALLOCATE(r_rot(NMolecules,3))
  ALLOCATE(TotalWeight(NMolecules))
  ALLOCATE(accept(NTrajectories))                       !Is it this config accepted? .TRUE. or .FALSE.
  ALLOCATE(RDFHistogram(NTrajectories,0:BinsRDF))           !Radial Distribution Function. The number of bins for the RDF hist is 1000
  ALLOCATE(RDFHistogramx(NTrajectories,0:BinsRDF))           !Radial Distribution Function. The number of bins for the RDF hist is 1000
  ALLOCATE(RDFHistogramy(NTrajectories,0:BinsRDF))           !Radial Distribution Function. The number of bins for the RDF hist is 1000
  ALLOCATE(RDFHistogramz(NTrajectories,0:BinsRDF))           !Radial Distribution Function. The number of bins for the RDF hist is 1000


!--------------------------------------------------------------------------------------------

  !-Calculate the Center of Mass for each Molecule
  DO i=1,NMolecules
     TotalWeight(i)=0.0                                                   !Initialization
     DO k=1,3                                                             !Loop for each coordinate x(k=1), y(k=2), z(k=3)
        Rcm(i,k)=0.0
        IF ( NAtoms(i) == 1 ) THEN                                        !For one atom fragments Rcm is not claculated
           Rcm(i,k)=ro(i,1,k)
        ELSE
           DO j=1,NAtoms(i)                                               !Loop for every Atom inside a Molecule(i)
              TotalWeight(i)=TotalWeight(i)+AtomicMass(ZNucleus(i,j))
              Rcm(i,k)=Rcm(i,k)+ro(i,j,k)*AtomicMass(ZNucleus(i,j))
           ENDDO
           Rcm(i,k)=Rcm(i,k)/TotalWeight(i)
        ENDIF
     ENDDO
  ENDDO
  
  !-Calculating the origin of the coordenates system based on the CMs
  DO k=1,3
     RcmMinimum(k)=Rcm(1,1)                                               !Initialization
     RcmMaximum(k)=Rcm(1,1)
     DO i=1,NMolecules
        RcmMinimum(k)=MIN(RcmMinimum(k),Rcm(i,k))
        RcmMaximum(k)=MAX(RcmMaximum(k),Rcm(i,k))
     ENDDO
     MiddlePoint(k)=0.5*(RcmMaximum(k)+RcmMinimum(k))
  ENDDO

  !-Radial Distribution Funtion
  MaxR= 2.0*radiusSBC
  MinR= 0.0
  DeltaR= (MaxR - MinR)/BinsRDF                     !The number of bins for the RDF hist is 1000
       
  !-Temperature Grating. The system is slowly "HEATED", the temperature goes from To to Tf, where Tf > To
  IF ( Tf > To .AND. Tf*To > 0.0 ) THEN
     DO t=1,NTrajectories
        Temperature(t)=To*(Tf/To)**(REAL(t-1)/REAL(NTrajectories-1))            !generate temperature list geometric
        !Temperature(t)=To + (t-1)*( (Tf-To)/(NTrajectories-1))                   !generate temperature list arithmetic
        beta(t)=1.0/(kB*Temperature(t))
     ENDDO
  ELSE
     WRITE(*,*)''
     WRITE(*,'(A)')         ' *** ERROR ***'
     WRITE(*,'(A,T20,F5.2)')' T initial: ', To
     WRITE(*,'(A,T20,F5.2)')' T final: ', Tf
     WRITE(*,'(A)')         ' ----'
     WRITE(*,'(A)')         ' T final is not larger than T initial!'
     WRITE(*,'(A)')         ' The system is slowly "HEATED", so the temperature'
     WRITE(*,'(A)')         ' goes from To to Tf, where Tf > To (both larger than zero)'
     WRITE(*,'(A)')         ' Please CHECK the input FILE'
     WRITE(*,*)''
     STOP
  ENDIF

  !Saving accepted configurations after a certain number of MC cycles (SaveConfig)
  IF ( SaveConfig > 0 ) THEN          !Save configs? Yes: SaveConfig>0 or No: SaveConfig<=0, SaveConfig defines a cycle to save the configurations         
     DO t=1,NTrajectories
        !molecular configuration in XYZ format
        WRITE(t_char,"(I10)") t
        OPEN(UNIT=500+t,FILE="configuration_"//TRIM(ADJUSTL(t_char))//".xyz",STATUS='UNKNOWN')
     ENDDO
  ENDIF

  !Reading the data from the inter and extrapolation (Look-Up table)
  OPEN (UNIT=400,FILE=EnergyFullPath,STATUS='OLD',ACTION='READ',IOSTAT=read_error)   !Was the OPEN successful?
  IF ( read_error == 0 ) THEN
     READ(400,*) NumDistances,MinDistance,SpacingDistances     !number of distances, minimum distance, spacing of distances
     READ(400,*) NumAngles,MinAngle,SpacingAngles              !number of angles, minimum angle, spacing of angles

     MinAngle=MinAngle*PI/180.0
     SpacingAngles=SpacingAngles*PI/180.0

     ALLOCATE(Energy_LookUpTable(0:NumAngles,0:NumDistances))

     DO IndexAngle=0,NumAngles-1
        DO IndexDistance=0,NumDistances-1
           READ(400,*) Energy_LookUpTable(IndexAngle,IndexDistance)
           Energy_LookUpTable(IndexAngle,IndexDistance)=Energy_LookUpTable(IndexAngle,IndexDistance)*EhtoJ
        ENDDO
     ENDDO
  ELSE
     !Else FILE open failed. Tell user.
     WRITE(*,*)''
     WRITE(*,*)'****************************************************'
     WRITE(*,"(A)")     ' Error! File '
     WRITE(*,"(10X,A)") '"'//EnergyFullPath//'"'
     WRITE(*,"(A)")     ' not found. Cannot evaluate the Energy'
     WRITE(*,*)'****************************************************'
     WRITE(*,*)''
     STOP
  END IF

  CLOSE (UNIT=400)

!--------------------------------------------------------------------------------------------
!Restarting Option

  !- Initialization
  MCcycles_restart = 0
  Ncheckpoints_mc_restart = 0
  Nequilibration_restart = 0
  Ncheckpoints_eq_restart = 0
  current_restart = 0
  eq_percent_restart = 0.0
  sampling_percent_restart = 0.0

  IF ( restart == 'yes' ) THEN
     IF ( input_error == 0 ) THEN

        READ(200,*) !comment line
        READ(200,*) !comment line
        READ(200,*) !comment line
        READ(200,*) !comment line
        READ(200,*) !comment line
        READ(200,*) !comment line

        READ(200,'(A)') string
        READ(string(30:),*,iostat=ierror) MCcycles_restart
                   IF ( MCcycles_restart <= 0 .OR. ierror /= 0 ) THEN
                      WRITE(*,*)
                      WRITE(*,*)' *** ERROR ***'
                      WRITE(*,*)' The total number of MC cycles'
                      WRITE(*,*)' must be positive integers'
                      WRITE(*,*) 'MC cycles/sampling: ',MCcycles_restart
                      WRITE(*,*)
                      WRITE(*,*)' Check restart FILE at line 7'
                      WRITE(*,*)' File: ',restart_FILE
                      WRITE(*,*)
                      STOP
                   ENDIF

        READ(200,'(A)') string
        READ(string(30:),*,iostat=ierror) Ncheckpoints_mc_restart
                   IF ( Ncheckpoints_mc_restart <= 0 .OR. ierror /= 0 ) THEN
                      WRITE(*,*)
                      WRITE(*,*)' *** ERROR ***'
                      WRITE(*,*)' Checkpoints for sampling must be a positive integer'
                      WRITE(*,*)' Checkpoint for sampling: ',Ncheckpoints_mc_restart
                      WRITE(*,*)
                      WRITE(*,*)' Check restart FILE at line 8'
                      WRITE(*,*)' File: ',restart_FILE
                      WRITE(*,*)
                      STOP
                   ENDIF

        READ(200,'(A)') string
        READ(string(30:),*,iostat=ierror) Nequilibration_restart
                   IF ( Nequilibration_restart <= 0 .OR. ierror /= 0 ) THEN
                      WRITE(*,*)
                      WRITE(*,*)' *** ERROR ***'
                      WRITE(*,*)' The total number of Equilibration'
                      WRITE(*,*)' must be positive integers, larger than zero'
                      WRITE(*,*) 'Equilibration: ',Nequilibration_restart
                      WRITE(*,*)
                      WRITE(*,*)' Check restart FILE at line 10'
                      WRITE(*,*)' File: ',restart_FILE
                      WRITE(*,*)
                      STOP
                   ENDIF

        READ(200,'(A)') string
        READ(string(30:),*,iostat=ierror) Ncheckpoints_eq_restart
                   IF ( Ncheckpoints_eq_restart <= 0 .OR. ierror /= 0 ) THEN
                      WRITE(*,*)
                      WRITE(*,*)' *** ERROR ***'
                      WRITE(*,*)' Checkpoints for Equilibration must be a positive integer'
                      WRITE(*,*) 'Checkpoint for Equilibration: ',Ncheckpoints_eq_restart
                      WRITE(*,*)
                      WRITE(*,*)' Check restart FILE at line 9'
                      WRITE(*,*)' File: ',restart_FILE
                      WRITE(*,*)
                      STOP
                   ENDIF

        READ(200,'(A)') string
        READ(string(30:),*,iostat=ierror) current_restart
                   IF (current_restart <= 0 .OR. ierror /= 0 ) THEN
                      WRITE(*,*)
                      WRITE(*,*)' *** ERROR ***'
                      WRITE(*,*)' Current point (MC step) must be positive integers, larger than zero'
                      WRITE(*,*) 'Current point: ',current_restart
                      WRITE(*,*)
                      WRITE(*,*)' Check restart FILE at line 11'
                      WRITE(*,*)' File: ',restart_FILE
                      WRITE(*,*)
                      STOP
                   ENDIF
        
        READ(200,'(A)') string
        READ(string(30:),*,iostat=ierror) eq_percent_restart 
                   IF ( eq_percent_restart < 0.0 .OR. ierror /= 0 ) THEN
                      WRITE(*,*)
                      WRITE(*,*)' *** ERROR ***'
                      WRITE(*,*)' The Evolution of the Equilibration'
                      WRITE(*,*)' must be positive float, larger than zero'
                      WRITE(*,*) 'Equilibration Evolution: ',eq_percent_restart
                      WRITE(*,*)
                      WRITE(*,*)' Check restart FILE at line 12'
                      WRITE(*,*)' File: ',restart_FILE
                      WRITE(*,*)
                      STOP
                   ENDIF

        READ(200,'(A)') string
        READ(string(30:),*,iostat=ierror) sampling_percent_restart
                   IF ( sampling_percent_restart < 0.0 .OR. ierror /= 0 ) THEN
                      WRITE(*,*)
                      WRITE(*,*)' *** ERROR ***'
                      WRITE(*,*)' The Evolution of the MC sampling'
                      WRITE(*,*)' must be positive float, larger than zero'
                      WRITE(*,*) 'Sampling Evolution ',sampling_percent_restart
                      WRITE(*,*)
                      WRITE(*,*)' Check restart FILE at line 11'
                      WRITE(*,*)' File: ',restart_FILE
                      WRITE(*,*)
                      STOP
                   ENDIF

        READ(200,'(A)') string
        READ(string(30:),*,iostat=ierror) Seed_restart
                   IF ( Seed_restart <= 0 .OR. ierror /= 0 ) THEN
                      WRITE(*,*)
                      WRITE(*,*)' *** ERROR ***'
                      WRITE(*,*)' The Seed for the Random number generator'
                      WRITE(*,*)' must be positive integers, larger than zero'
                      WRITE(*,*) 'Seed: ',Seed_restart
                      WRITE(*,*)
                      WRITE(*,*)' Check restart FILE at line 13'
                      WRITE(*,*)' File: ',restart_FILE
                      WRITE(*,*)
                      STOP
                   ELSEIF ( Seed_restart > 0  ) THEN
                      Seed = Seed_restart
                   ENDIF

        READ(200,*) !comment line
        READ(200,*) !comment line
        READ(200,*) !comment line
        READ(200,*) !comment line
        READ(200,*) !comment line
        READ(200,*) !comment line

        !- Initia States
        READ(200,'(A)') string
        READ(string(40:),*,iostat=ierror) Energy_Tot(0)
                                          Energy_Tot_fin(0)=0.0

                   IF ( ierror /= 0 ) THEN
                      WRITE(*,*)
                      WRITE(*,*)' *** ERROR ***'
                      WRITE(*,*)' Total Energy (initial state): ',Energy_Tot(0)
                      WRITE(*,*)
                      WRITE(*,*)' Check restart FILE at line 21'
                      WRITE(*,*)' File: ',restart_FILE
                      WRITE(*,*)
                      STOP
                   ENDIF

        READ(200,*) !comment line

        DO Molecule_a=1,NMolecules

           READ(200,'(A)') string
           READ(string(40:),*,iostat=ierror) Energy_2a_ini(0,Molecule_a)

           Energy_ab_ini(0,Molecule_a,Molecule_a)= 0.0    !E(a,a)=0

           DO Molecule_b=Molecule_a+1,NMolecules
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) Energy_ab_ini(0,Molecule_a,Molecule_b)
              Energy_ab_ini(0,Molecule_b,Molecule_a) = Energy_ab_ini(0,Molecule_a,Molecule_b)
           ENDDO
      
           READ(200,*) !comment line

        ENDDO

        !Was equilibration already done?
        !Equilibration evolution < 100%
        IF ( sampling_percent_restart == 0.0 .AND. eq_percent_restart > 0.0 ) THEN  !Equilibration

           DO t=1,NTrajectories
           
              READ(200,*) !comment line
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) Temperature(t)
              READ(200,*) !comment line
           
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) maxDisplacement(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) MaxRotation(t)   !radians
           
              READ(200,*) !comment line
           
              !Structural parameters 
              DO i=1,NMolecules
                 READ(200,*,iostat=ierror) ZNucleus(i,1),rf(t,i,1,1),rf(t,i,1,2),rf(t,i,1,3)
                 r(t,i,1,1) = rf(t,i,1,1)
                 r(t,i,1,2) = rf(t,i,1,2)
                 r(t,i,1,3) = rf(t,i,1,3)
!                 rex(t,i,1,1)=0.0 ; rex(t,i,1,2)=0.0 ; rex(t,i,1,3)=0.0
              ENDDO
           
              READ(200,*) !comment line
           
              !Acceptance criterion, adjusting step size
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) Mov_accepted_eq(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) adj_Mov_accepted_eq(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) Rot_accepted_eq(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) adj_Rot_accepted_eq(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) Exchange_tot(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) Exchange_acc(t)
           
              READ(200,*) !comment line
              READ(200,*) !comment line
              READ(200,*) !comment line
              READ(200,*) !comment line
           
              !Energies
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) Energy_Tot(t)
                                                Energy_Tot_fin(t) = Energy_Tot(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) MinEnergy(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) MaxEnergy(t)
           
              READ(200,*) !comment line
           
              DO Molecule_a=1,NMolecules
                 READ(200,'(A)') string
                 READ(string(40:),*,iostat=ierror) Energy_2a_ini(t,Molecule_a)
           
                 !-Initialization
                 Energy_2a_fin(t,Molecule_a) = Energy_2a_ini(t,Molecule_a)
                 Energy_ab_ini(t,Molecule_a,Molecule_a)= 0.0    !E(a,a)=0
                 Energy_ab_fin(t,Molecule_a,Molecule_a)= 0.0    !E(a,a)=0
           
                 DO Molecule_b=Molecule_a+1,NMolecules
                    READ(200,'(A)') string
                    READ(string(40:),*,iostat=ierror) Energy_ab_ini(t,Molecule_a,Molecule_b)
           
                    !-Initialization
                    Energy_ab_ini(t,Molecule_b,Molecule_a)= Energy_ab_ini(t,Molecule_a,Molecule_b)
                    Energy_ab_fin(t,Molecule_a,Molecule_b)= Energy_ab_ini(t,Molecule_a,Molecule_b)
                    Energy_ab_fin(t,Molecule_b,Molecule_a)= Energy_ab_ini(t,Molecule_a,Molecule_b)
                 ENDDO
             
                 READ(200,*) !comment line
           
              ENDDO

              !Equilibration Evolution <= 100%

              TotalEnergy(t)     =0.0                            !Sum of total energies per temperature
              TotalEnergy2(t)    =0.0                            !Sum of total energies**2 per temperature  
              Mov_accepted(t)    =0                              !After equilibration
              Rot_accepted(t)    =0                              !After equilibration
              adj_Mov_accepted(t)=0                              !Number of accepted configuration after linear MOVE, to update displacement size after 100 steps
              adj_Rot_accepted(t)=0                              !Number of accepted configuration after ROTATIONS, to update angle size after 100 steps

              DO hit=0,nBins
                 nhits(t,hit)=0
              ENDDO
              !-Radial Distribution Funtion
              DO indexRDF=0,BinsRDF
                 RDFHistogram(t,indexRDF)=0
                 RDFHistogramx(t,indexRDF)=0
                 RDFHistogramy(t,indexRDF)=0
                 RDFHistogramz(t,indexRDF)=0
              ENDDO

           ENDDO

           !Initialization Max & Min energy range for histogram analysis
           MinEnergyTotal = MinEnergy(1)
           MaxEnergyTotal = MaxEnergy(1)

           !Counting the energies OUTSIDE the range (Emax-Emin)
           nConfigELT=0
           nConfigEGT=0
          
           !Equilibration Evolution = 100%
           IF ( sampling_percent_restart == 0.0 .AND. eq_percent_restart == 100.0 ) THEN 

              READ(200,*) !comment line
  
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) MinEnergyTotal 
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) MaxEnergyTotal 
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) DeltaEHistogram

           ENDIF

        !during Samplig (after eq)
        ELSEIF ( sampling_percent_restart > 0.0 .AND. eq_percent_restart == 100.0 ) THEN  !Sampling Evolution > 0%

           DO t=1,NTrajectories

              READ(200,*) !comment line
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) Temperature(t)
              READ(200,*) !comment line
           
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) maxDisplacement(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) MaxRotation(t)   !radians
           
              READ(200,*) !comment line
           
              !Structural parameters 
              DO i=1,NMolecules
                 READ(200,*,iostat=ierror) ZNucleus(i,1),rf(t,i,1,1),rf(t,i,1,2),rf(t,i,1,3)
                 r(t,i,1,1) = rf(t,i,1,1)
                 r(t,i,1,2) = rf(t,i,1,2)
                 r(t,i,1,3) = rf(t,i,1,3)
              ENDDO
           
              READ(200,*) !comment line
           
              !Acceptance criterion, adjusting step size
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) Mov_accepted_eq(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) adj_Mov_accepted_eq(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) Rot_accepted_eq(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) adj_Rot_accepted_eq(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) Exchange_tot(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) Exchange_acc(t)
           
              READ(200,*) !comment line
              READ(200,*) !comment line
              READ(200,*) !comment line
              READ(200,*) !comment line
           
              !Energies
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) Energy_Tot(t)
                                                Energy_Tot_fin(t) = Energy_Tot(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) MinEnergy(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) MaxEnergy(t)
           
              READ(200,*) !comment line
           
              DO Molecule_a=1,NMolecules
                 READ(200,'(A)') string
                 READ(string(40:),*,iostat=ierror) Energy_2a_ini(t,Molecule_a)
           
                 !-Initialization
                 Energy_2a_fin(t,Molecule_a) = Energy_2a_ini(t,Molecule_a)
                 Energy_ab_ini(t,Molecule_a,Molecule_a)= 0.0    !E(a,a)=0
                 Energy_ab_fin(t,Molecule_a,Molecule_a)= 0.0    !E(a,a)=0
           
                 DO Molecule_b=Molecule_a+1,NMolecules
                    READ(200,'(A)') string
                    READ(string(40:),*,iostat=ierror) Energy_ab_ini(t,Molecule_a,Molecule_b)
           
                    !-Initialization
                    Energy_ab_ini(t,Molecule_b,Molecule_a)= Energy_ab_ini(t,Molecule_a,Molecule_b)
                    Energy_ab_fin(t,Molecule_a,Molecule_b)= Energy_ab_ini(t,Molecule_a,Molecule_b)
                    Energy_ab_fin(t,Molecule_b,Molecule_a)= Energy_ab_ini(t,Molecule_a,Molecule_b)
                 ENDDO
             
                 READ(200,*) !comment line
           
              ENDDO

              !during Samplig (after eq)
              READ(200,*) !comment line
              READ(200,*) !comment line
  
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) Mov_accepted(t) 
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) adj_Mov_accepted(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) Rot_accepted(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) adj_Rot_accepted(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) TotalEnergy(t)
              READ(200,'(A)') string
              READ(string(40:),*,iostat=ierror) TotalEnergy2(t)

              READ(200,*) !comment line
              READ(200,*) !comment line
              DO hit=0,nBins
                 READ(200,'(A)') string
                 READ(string,"(I20)",iostat=ierror) nhits(t,hit)
!                 READ(200,*,iostat=ierror) nhits(t,hit)
              ENDDO
      
              READ(200,*) !comment line
              READ(200,*) !comment line
              READ(200,*) !comment line
              DO indexRDF=0,BinsRDF
                 READ(200,"(I20)",iostat=ierror) RDFHistogramx(t,indexRDF)
              ENDDO
              READ(200,*) !comment line
              DO indexRDF=0,BinsRDF
                 READ(200,"(I20)",iostat=ierror) RDFHistogramy(t,indexRDF)
              ENDDO
              READ(200,*) !comment line
              DO indexRDF=0,BinsRDF
                 READ(200,"(I20)",iostat=ierror) RDFHistogramz(t,indexRDF)
              ENDDO
              READ(200,*) !comment line
              DO indexRDF=0,BinsRDF
                 READ(200,"(I20)",iostat=ierror) RDFHistogram(t,indexRDF)
              ENDDO

              READ(200,*) !comment line

           ENDDO
  
           READ(200,'(A)') string
           READ(string(40:),*,iostat=ierror) MinEnergyTotal 
           READ(200,'(A)') string
           READ(string(40:),*,iostat=ierror) MaxEnergyTotal 
           READ(200,'(A)') string
           READ(string(40:),*,iostat=ierror) DeltaEHistogram
           READ(200,'(A)') string
           READ(string(30:),*,iostat=ierror) nConfigEGT 
           READ(200,'(A)') string
           READ(string(30:),*,iostat=ierror) nConfigELT

        ENDIF

        READ(200,*) !comment line
        READ(200,*) !comment line
        READ(200,*) !comment line
        READ(200,*) !comment line

        CLOSE(UNIT=200) !restart FILE (previous one) to start from

     ELSE
        !Else FILE open failed. Tell user.
        WRITE (*,"(A)") 
        WRITE (*,"(A)")' *** ERROR ***'
        WRITE (*,"(A)")" RESTART FILE open failed. Restart flag activated from INPUT"
        WRITE (*,"(2A)")" ---  ",input
        WRITE (*,"(A)")" but no 'restart_"//output//"' FILE was found" 
        WRITE (*,"(A)") 
        STOP
     END IF
    

  !- no restarting
  ELSE !For any other case, we always start the MC sampling from the beginning
     !-Displacement of the CM
     DO t=1,NTrajectories
        maxDisplacement(t)=0.1*SQRT(MaxDisplacementRcm*Temperature(t))  !according to Flouront's code
     !   maxDisplacement(t)=MaxDisplacementRcm
        MaxRotation(t)=MaxAngle
     ENDDO
     
     !-Defining new atomic coordinates (r) to move. It is the same starting configuration for every temperature
     DO t=0,NTrajectories
        DO i=1,NMolecules
           DO j=1,NAtoms(i)
              DO k=1,3
                 r (t,i,j,k) =ro(i,j,k)   !coordinates for any configuration accepted or not
                 rf(t,i,j,k) =ro(i,j,k)   !coordinates for accepted configuration
                 rex(t,i,j,k)=0.0         !rex is a temporary variable
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     
     !Initialization. Total Energy per Temperature
     Energy_Tot(0)=0.0
     Energy_Tot_fin(0)=0.0
     
     DO Molecule_a=1,NMolecules
     
        !E(a,a)=0
        Energy_ab_ini(0,Molecule_a,Molecule_a)= 0.0
     
        DO Molecule_b=Molecule_a+1,NMolecules
           dist=SQRT( (r(0,Molecule_a,1,1) - r(0,Molecule_b,1,1))**2+&   !r_ab=SQRT[ (xa-xb)^2+
                      (r(0,Molecule_a,1,2) - r(0,Molecule_b,1,2))**2+&   !           (ya-yb)^2+
                      (r(0,Molecule_a,1,3) - r(0,Molecule_b,1,3))**2  )  !           (za-zb)^2 ]
     
           !From dot product, where |k|=1 and k . r12 = z1-z2, 0<=theta<=90
           theta=ACOS( ABS(r(0,Molecule_a,1,3) - r(0,Molecule_b,1,3))/dist )   !Cos(x)=(z1-z2)/|r12|
     
           indexAngle   =NINT( (theta-MinAngle)/SpacingAngles ) 
           indexDistance=NINT( (dist-MinDistance)/SpacingDistances ) 
     
           !computing the pair energy E(a,b)
           Energy_ab_ini(0,Molecule_a,Molecule_b)= Energy_LookUpTable(indexAngle,indexDistance)
     
           !E(a,b)=E(b,a)
           Energy_ab_ini(0,Molecule_b,Molecule_a)= Energy_ab_ini(0,Molecule_a,Molecule_b) 
     
           !Total energy
           Energy_Tot(0)= Energy_Tot(0) + Energy_ab_ini(0,Molecule_a,Molecule_b) 
        ENDDO
     
        !computing pair energy contribution for i-th
        Energy_2a_ini(0,Molecule_a)= 0.0                 !Initialization
        DO Molecule_b=1,NMolecules
           Energy_2a_ini(0,Molecule_a) = Energy_2a_ini(0,Molecule_a) + Energy_ab_ini(0,Molecule_a,Molecule_b)
        ENDDO
     
     ENDDO
     
     DO t=1,NTrajectories
        Energy_Tot(t)      =Energy_Tot(0)                  !The initial energy is the same for all temperatures
        Energy_Tot_fin(t)  =Energy_Tot(0)                  !The initial energy is the same for all temperatures
        TotalEnergy(t)     =0.0                            !Sum of total energies per temperature
        TotalEnergy2(t)    =0.0                            !Sum of total energies**2 per temperature  
        MinEnergy(t)       =Energy_Tot(0)                  !Initialization Max & Min energy range for histogram analysis
        MaxEnergy(t)       =Energy_Tot(0)
        Mov_accepted_eq(t) =0                              !Number of accepted configuration, to calculate the acceptance ratio
        Rot_accepted_eq(t) =0                              !Number of accepted configuration after ROTATIONS, to update step size maxDisplacement of CM
        adj_Mov_accepted_eq(t)=0 
        adj_Rot_accepted_eq(t)=0
        Mov_accepted(t)    =0                              !After equilibration
        Rot_accepted(t)    =0                              !After equilibration
        adj_Mov_accepted(t)=0                              !Number of accepted configuration after linear MOVE, to update displacement size after 100 steps
        adj_Rot_accepted(t)=0                              !Number of accepted configuration after ROTATIONS, to update angle size after 100 steps
        Exchange_tot(t)    =0                              !Number of exchanges 
        Exchange_acc(t)    =0                              !Number of exchanges accepted
        DO Molecule_a=1,NMolecules
           Energy_2a_ini(t,Molecule_a) = Energy_2a_ini(0,Molecule_a)
     
           DO Molecule_b=1,NMolecules
              Energy_ab_ini(t,Molecule_a,Molecule_b)= Energy_ab_ini(0,Molecule_a,Molecule_b) 
              !E(a,b)=E(b,a)
              Energy_ab_ini(t,Molecule_b,Molecule_a)= Energy_ab_ini(0,Molecule_a,Molecule_b) 
     
              Energy_ab_fin(t,Molecule_a,Molecule_b)= Energy_ab_ini(0,Molecule_a,Molecule_b)
              !E(a,b)=E(b,a)
              Energy_ab_fin(t,Molecule_b,Molecule_a)= Energy_ab_ini(0,Molecule_a,Molecule_b)
           ENDDO
        ENDDO
     ENDDO
     
     !Initialization Max & Min energy range for histogram analysis
     MinEnergyTotal =Energy_Tot(0)   !from J to kJ/mol
     MaxEnergyTotal =Energy_Tot(0)   !from J to kJ/mol
     !Counting the energies OUTSIDE the range (Emax-Emin)
     nConfigELT=0
     nConfigEGT=0
       
     DO t=1,NTrajectories
        DO hit=0,nBins
           nhits(t,hit)=0
        ENDDO
        !-Radial Distribution Funtion
        DO indexRDF=0,BinsRDF
           RDFHistogram(t,indexRDF)=0
           RDFHistogramx(t,indexRDF)=0
           RDFHistogramy(t,indexRDF)=0
           RDFHistogramz(t,indexRDF)=0
        ENDDO
     ENDDO
  ENDIF !end or restart option
     
!--------------------------------------------------------------------------------------------
  !-Write the Configuration FILE XYZ for the initial configuration
  OPEN(UNIT=300,FILE="initial_configuration.xyz",STATUS='REPLACE')

  WRITE(energy_tot_char,"(ES18.10E2)") Energy_Tot(0)*JtokJmol
  comment_line='Energy (kJ/mol): '//TRIM(ADJUSTL(energy_tot_char))
  FILEunit=300
  t=1  !for rf(t)

  CALL DrawMolecularConfiguration( ZNucleus,radiusSBC,NMolecules,NAtoms,  &
                                   Max_Atoms,MiddlePoint,NTrajectories,rf,&
                                   t,FILEunit,comment_line                 )

  CLOSE(300) !starting configuration

!--------------------------------------------------------------------------------------------
!-Writing the output FILE. An output FILE will be created for every checkpoint (unit=check+1)
  
  !-Removing blanks
  WRITE(MCcycles_char,"(I10)")   MCcycles  
  WRITE(MCstep_char,"(I10)")     MCstep    
  WRITE(Nequilibration_char,"(I10)") Nequilibration 
  WRITE(NMolecules_char,"(I10)") NMolecules
  WRITE(Max_Atoms_char,"(I10)")  Max_Atoms 
  WRITE(Ncheckpoints_mc_char,"(I10)") Ncheckpoints_mc 
  WRITE(Ncheckpoints_eq_char,"(I10)") Ncheckpoints_eq
  
  DO i=1,NMolecules
     WRITE(NAtoms_char(i),"(I10)") NAtoms(i)
  ENDDO

  OPEN(UNIT=10000,FILE=output//'.out',STATUS='REPLACE')
  
  WRITE(10000,'(A)')'****************************************************************'
  WRITE(10000,'(12A)')' Execution time: ',&
                                   today(1:4),"-",today(5:6),"-",today(7:8), ", ",&
                                   now(1:2),":",  now(3:4),":",  now(5:6)
  WRITE(10000,'(A)')'****************************************************************'
  WRITE(10000,'(A)')''
  WRITE(10000,'(A)')' Contents of the input FILE'
  WRITE(10000,'(A)')'----------------------------------------------------------------'
  WRITE(10000,"(A,T45,I20)")' Number of threads availables: ',aval_threads
  WRITE(10000,"(A,T45,I20)")' Number of threads used: ',threads
  WRITE(10000,"(A,T45,A20)")' Total MC Steps per Trajectory: ',TRIM(ADJUSTL(MCcycles_char))  
  WRITE(10000,"(A,T45,A20)")' Equilibration steps per Trajectory: ',TRIM(ADJUSTL(Nequilibration_char))  
  WRITE(10000,'(A,T45,F20.2)')' Minimum Temperature: ',To
  WRITE(10000,'(A,T45,F20.2)')' Maximum Temperature: ',Tf
  WRITE(10000,'(A,T45,I20)')' Number of Temperatures/Trajectories: ',NTrajectories
  WRITE(10000,"(A,T45,F20.2)")' Spherical BC, radius (Angst): ',radiusSBC
  WRITE(10000,"(A,T45,A20)")  ' The system has (Molecules/Atoms): ',TRIM(ADJUSTL(NMolecules_char))
  WRITE(10000,"(A,T45,F20.2)")' Max Displacement (Angst): ',MaxDisplacementRcm
  WRITE(10000,"(A,T45,F20.2)")' Max Rotation (degrees): ',MaxAngle*(180.0/PI) 
  WRITE(10000,'(A,T25,A40)')' Potential Energy: ',TRIM(ADJUSTL(EnergyName))
  WRITE(10000,"(A,T47,ES18.10E2)")' Initial Energy (kJ/mol): ',Energy_Tot(0)*JtokJmol
  WRITE(10000,"(A,T47,ES18.10E2)")' The Boltzmann constant (J/K): ',kB
  WRITE(10000,"(A,T45,I20)")' Seed for Random Number generator: ',Seed
  WRITE(10000,"(A,T45,I20)")' A configuration saved every (per traj): ',SaveConfig
  WRITE(10000,'(A)')''

  IF ( restart == 'yes' ) THEN
     WRITE(10000,"(A,T30,A35)") ' Simulation restarting from: ',TRIM(ADJUSTL(restart_FILE))
     WRITE(10000,'(A,T59,F6.2)') ' Equilibration evolution was (%): ',eq_percent_restart
     WRITE(10000,'(A,T45,I20)') '    of a total number of Eq steps: ',Nequilibration_restart
     WRITE(10000,'(A,T59,F6.2)') ' Sample evolution was (%): ',sampling_percent_restart
     WRITE(10000,'(A,T45,I20)') '    of a total number of MC steps: ',MCcycles_restart
     WRITE(10000,'(A)')
  ENDIF
  
  !coordinates
  WRITE(10000,'(A)')' The initial configuration is the same for every trajectory'
  WRITE(10000,'(A)')" Coordinates written at 'initial_configuration.xyz'"
  WRITE(10000,'(A)')''
  WRITE(10000,'(A)')'----------------------------------------------------------------'
  WRITE(10000,"(A)")' OUTPUT'
  WRITE(10000,'(A)')'----------------------------------------------------------------'
  
  !Flushing memory buffer
  CALL FLUSH(10000)  !writing the general output FILE 

  !----------------------------------------------------------------
  !STATUS, this FILE has the MC evolution
  OPEN(UNIT=1000,FILE='status_'//output,STATUS='REPLACE')                  !MC evolution
  !-writing the MC procedure evolution
  PartialStartRunTime=0.0
  PartialEndRunTime=0.0
  CALL DATE_AND_TIME(today,now)
  WRITE(1000,'(A)')
  WRITE(1000,"(12A)")'Date and time: ',&
                     today(1:4),"-",today(5:6),"-",today(7:8), ", ",&
                     now(1:2),":",  now(3:4),":",  now(5:6)
  WRITE(1000,'(A,T36,I10)')'Total Steps: ',MCcycles
  WRITE(1000,'(A,T36,I10)')'Equilibration:',Nequilibration
  WRITE(1000,'(A)')'Status evolution: '
  WRITE(1000,'(5X,A)') TRIM(ADJUSTL(Ncheckpoints_eq_char))//' checkpoints for equilibration and'
  WRITE(1000,'(5X,A)') TRIM(ADJUSTL(Ncheckpoints_mc_char))//' checkpoints for sampling'

  FLUSH(1000) !File "status"


!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
! -Monte Carlo Cycle: We define two independent loops, equilibration and sampling 
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------

  !- We have TWO restarting options: cold and warm 

  !- COLD restart: starting from the configurations obtained by a previous simulation run,
  !                but MC cycle (eq and sampling) run as a new/different simulation
  !-               For exemple, when you want to extend a simulation; equlibration, sampling,
  !-               or both

  !- WARM restart: basically a job recovery measure. If a job crashed, you can restart
  !-               it from where it was saved (any checkpoint), and obtain essentially the same result. 
  !-               Another alternative is to re-run a simulation increasing the number of checkpoints

!--------------------------------------------------------------------------------------------
!-Equilibration Loop
!--------------------------------------------------------------------------------------------

  check_restart_eq = 1 + NINT(eq_percent_restart*REAL(Ncheckpoints_eq)/100.0)

  !- COLD restart: extending the number of steps for equilibration or
  !- you can restart a NEW/different simulation from the last checkpoint
  !- saved for a previous simulation
!  IF ( sampling_percent_restart > 99.999 ) THEN
!     check_restart_eq = Ncheckpoints_eq 
!  ENDIF

  check_split_eq = Nequilibration/Ncheckpoints_eq

  DO check=check_restart_eq,Ncheckpoints_eq

     !skipping Eqilibration loop
     IF ( Nequilibration <= 1 .OR. sampling_percent_restart > 99.999 ) EXIT

     WRITE(check_char,"(I10)") check

     !----------------------------------------------------------------------------------------
     !-Start point of MC cycle
     check_MCstep = 1 + check_split_eq*(check - 1)
     check_MCcycles = check_split_eq*(check)

     WRITE(check_MCcycles_char,"(I10)") check_MCcycles

     !--------------------------------------------------------------------------------------------
     !-Writing the output FILE
   
     OPEN(UNIT=10000+check,FILE=output//"_"//TRIM(ADJUSTL(check_char))//"ckp_eq.out",STATUS='REPLACE')
   
     WRITE(10000+check,'(A)')'****************************************************************'
     WRITE(10000+check,'(12A)')' Execution time: ',&
                                      today(1:4),"-",today(5:6),"-",today(7:8), ", ",&
                                      now(1:2),":",  now(3:4),":",  now(5:6)
     WRITE(10000+check,'(A)')'****************************************************************'
     WRITE(10000+check,'(A)')''
     WRITE(10000+check,'(A)')' Contents of the input FILE'
     WRITE(10000+check,'(A)')'----------------------------------------------------------------'
     WRITE(10000+check,"(A,T45,I20)")' Number of threads availables: ',aval_threads
     WRITE(10000+check,"(A,T45,I20)")' Number of threads used: ',threads
     WRITE(10000+check,"(A,T45,A20)")' Total MC Steps per Trajectory: ',TRIM(ADJUSTL(MCcycles_char))  
     WRITE(10000+check,"(A,T45,A20)")' Equilibration steps per Trajectory: ',TRIM(ADJUSTL(Nequilibration_char))  
     WRITE(10000+check,'(A,T45,F20.2)')' Minimum Temperature: ',To
     WRITE(10000+check,'(A,T45,F20.2)')' Maximum Temperature: ',Tf
     WRITE(10000+check,'(A,T45,I20)')' Number of Temperatures/Trajectories: ',NTrajectories
     WRITE(10000+check,"(A,T45,F20.2)")' Spherical BC, radius (Angst): ',radiusSBC
     WRITE(10000+check,"(A,T45,A20)")  ' The system has (Molecules/Atoms): ',TRIM(ADJUSTL(NMolecules_char))
     WRITE(10000+check,"(A,T45,F20.2)")' Max Displacement (Angst): ',MaxDisplacementRcm
     WRITE(10000+check,"(A,T45,F20.2)")' Max Rotation (degrees): ',MaxAngle*(180.0/PI) 
     WRITE(10000+check,'(A,T25,A40)')' Potential Energy: ',TRIM(ADJUSTL(EnergyName))
     WRITE(10000+check,"(A,T47,ES18.10E2)")' Initial Energy (kJ/mol): ',Energy_Tot(0)*JtokJmol
     WRITE(10000+check,"(A,T47,ES18.10E2)")' The Boltzmann constant (J/K): ',kB
     WRITE(10000+check,"(A,T45,I20)")' Seed2 for Random Number generator: ',Seed2
     WRITE(10000+check,"(A,T45,I20)")' A configuration saved every (per traj): ',SaveConfig
     WRITE(10000+check,'(A)')''

     IF ( restart == 'yes' ) THEN
        WRITE(10000+check,"(A,T30,A35)") ' Simulation restarting from: ',TRIM(ADJUSTL(restart_FILE))
        WRITE(10000+check,'(A,T59,F6.2)') ' Equilibration evolution was (%): ',eq_percent_restart
        WRITE(10000+check,'(A,T45,I20)') '    of a total number of Eq steps: ',Nequilibration_restart
        WRITE(10000+check,'(A,T59,F6.2)') ' Sample evolution was (%): ',sampling_percent_restart
        WRITE(10000+check,'(A,T45,I20)') '    of a total number of MC steps: ',MCcycles_restart
        WRITE(10000+check,'(A)')
     ENDIF

     WRITE(10000+check,'(A)')          '----------------------------------------------------------------'
     WRITE(10000+check,'(A,T20,F5.1,A)')' Equilibration: ',(REAL(check)/REAL(Ncheckpoints_eq))*100.0,' %,  '//&
                                         TRIM(ADJUSTL(check_MCcycles_char))//' CYCLES DONE OF '//&
                                         TRIM(ADJUSTL(Nequilibration_char))
     WRITE(10000+check,'(A)')          '----------------------------------------------------------------'
   
     !Flushing memory buffer
     CALL FLUSH(10000+check)  !writing the general output FILE 

     !----------------------------------------------------------------------------------------

   
     DO MCstep=check_MCstep,check_MCcycles          !loop over MC steps for EQUILIBRATION
 
        
        !$OMP PARALLEL DO DEFAULT (SHARED) &
        !$OMP FIRSTPRIVATE (Seed) &
        !$OMP PRIVATE (t,SelectedMolecule,r_rot,dist,theta,indexAngle,indexDistance,Seed2 ) &
        !$OMP PRIVATE (Molecule_a,Molecule_b,DeltaEnergy,r_move,acceptance,k,i,n,Adjust)
        DO t=1,NTrajectories         !loop over temperatures/trajectories
            
           Seed2 = Seed + t

        DO n=1,NMolecules+1          !EXTRA loop over number of Molecules in the MC cluster. N+1 because the rotation of the whole cluster

           !Randomly selected the molecule which is translated and rotated
           SelectedMolecule=FLOOR( 1 + (NMolecules + Rotate)*RandomNumber(Seed2) )                  !if Ran()=0, then SM=1
   
           !if SM > N we rotate the whole cluster, otherwise a random atom is moved
           ! 1/(NMolecules + 1)  is the probability of attempting a rotation of the whole cluster
           IF ( SelectedMolecule > NMolecules ) THEN
   
              CALL RotateCluster( rf,NMolecules,MaxRotation(t),Seed2,NTrajectories,t,MiddlePoint,&   !IN
                                  r_rot                                                          )  !OUT
   
              DO k=1,3
                 DO i=1,NMolecules
                    r(t,i,1,k)=r_rot(i,k)
                 ENDDO
              ENDDO
   
              DO Molecule_a=1,NMolecules
                 DO Molecule_b=Molecule_a+1,NMolecules
                    dist=SQRT( (r(t,Molecule_a,1,1) - r(t,Molecule_b,1,1))**2+&   !r_ab=SQRT[ (xa-xb)^2+
                               (r(t,Molecule_a,1,2) - r(t,Molecule_b,1,2))**2+&   !           (ya-yb)^2+
                               (r(t,Molecule_a,1,3) - r(t,Molecule_b,1,3))**2  )  !           (za-zb)^2 ]
   
                    !From dot product, where |k|=1 and k . r12 = z1-z2, 0<=theta<=90
                    theta=ACOS( ABS(r(t,Molecule_a,1,3) - r(t,Molecule_b,1,3))/dist )   !Cos(x)=(z1-z2)/|r12|
   
                    indexAngle   =NINT( (theta-MinAngle)/SpacingAngles )
                    indexDistance=NINT( (dist-MinDistance)/SpacingDistances )
   
                    !computing the pair energy E(a,b)
                    Energy_ab_fin(t,Molecule_a,Molecule_b)= Energy_LookUpTable(indexAngle,indexDistance)
                 ENDDO
              ENDDO

              Energy_Tot_fin(t)=0.0                 !Initialization

              DO Molecule_a=1,NMolecules
                 !Initialization
                 Energy_2a_fin(t,Molecule_a)= 0.0
                 !E(a,a)=0
                 Energy_ab_fin(t,Molecule_a,Molecule_a)= 0.0

                 DO Molecule_b=Molecule_a+1,NMolecules
                    !E(a,b)=E(b,a)
                    Energy_ab_fin(t,Molecule_b,Molecule_a)= Energy_ab_fin(t,Molecule_a,Molecule_b)
   
                    !Total energy
                    Energy_Tot_fin(t)= Energy_Tot_fin(t) + Energy_ab_fin(t,Molecule_a,Molecule_b)
                 ENDDO

                 !computing pair energy contribution for i-th
                 DO Molecule_b=1,NMolecules
                    Energy_2a_fin(t,Molecule_a) = Energy_2a_fin(t,Molecule_a) + Energy_ab_fin(t,Molecule_a,Molecule_b)
                 ENDDO
              ENDDO
   
              DeltaEnergy=Energy_Tot_fin(t) - Energy_Tot(t)
   
           ELSE  !if the whole structure is not rotated, we move one atom randomly
   
              CALL DisplaceAtoms( rf,NMolecules,maxDisplacement(t),SelectedMolecule, &   !IN 
                                  Seed2,radiusSBC2,MiddlePoint,NTrajectories,t,       &   !IN
                                  r_move                                              )  !OUT
              
              
              DO j=1,NAtoms(SelectedMolecule)
                 DO k=1,3
                    r(t,SelectedMolecule,j,k)=r_move(k)
                 ENDDO
              ENDDO
              
              DO Molecule_a=1,NMolecules
                 IF ( Molecule_a == SelectedMolecule ) CYCLE             !Avoiding to compute self interaction
              
                 dist=SQRT( (r(t,Molecule_a,1,1) - r(t,SelectedMolecule,1,1))**2+&   !r_ab=SQRT[ (xa-xb)^2+
                            (r(t,Molecule_a,1,2) - r(t,SelectedMolecule,1,2))**2+&   !           (ya-yb)^2+
                            (r(t,Molecule_a,1,3) - r(t,SelectedMolecule,1,3))**2  )  !           (za-zb)^2 ]
              
                 !From dot product, where |k|=1 and k . r12 = z1-z2, 0<=theta<=90
                 theta=ACOS( ABS(r(t,Molecule_a,1,3) - r(t,SelectedMolecule,1,3))/dist )   !Cos(x)=(z1-z2)/|r12|
              
                 indexAngle   =NINT( (theta-MinAngle)/SpacingAngles )
                 indexDistance=NINT( (dist-MinDistance)/SpacingDistances )
              
                 !computing the pair energy E(a,b)
                 Energy_ab_fin(t,Molecule_a,SelectedMolecule)= Energy_LookUpTable(indexAngle,indexDistance)
              ENDDO
              
              Energy_2a_fin(t,SelectedMolecule)= 0.0                 !Initialization
              
              DO Molecule_a=1,NMolecules
                 IF ( Molecule_a == SelectedMolecule ) CYCLE             !Avoiding to compute self interaction

                 !E(a,a)=0
                 Energy_ab_fin(t,Molecule_a,Molecule_a)= 0.0

                 !E(a,b)=E(b,a)
                 Energy_ab_fin(t,SelectedMolecule,Molecule_a)= Energy_ab_fin(t,Molecule_a,SelectedMolecule)
              
                 !computing total 2-body energy contribution for i-th
                 Energy_2a_fin(t,SelectedMolecule) = Energy_2a_fin(t,SelectedMolecule) &
                                                     + Energy_ab_fin(t,Molecule_a,SelectedMolecule)
              ENDDO

              DeltaEnergy=Energy_2a_fin(t,SelectedMolecule) - Energy_2a_ini(t,SelectedMolecule)
   
           ENDIF
  
           accept(t)=.FALSE.

           !Metropolis acceptance criterion
           IF ( DeltaEnergy <= 0.0 ) THEN
              acceptance=1.0
           ELSE 
              acceptance=DEXP( -DeltaEnergy*beta(t) )
           ENDIF
   
           IF ( acceptance > RandomNumber(Seed2) ) THEN
    
              accept(t)=.TRUE.
   
              !saving the total energy.
              Energy_Tot(t)= Energy_Tot(t) + DeltaEnergy
   
              !Defining Bin width: ENERGY HISTOGRAM ANALYSIS 
              MaxEnergy(t)=Max( MaxEnergy(t),Energy_Tot(t) )
              MinEnergy(t)=Min( MinEnergy(t),Energy_Tot(t) )
   
              !if SM > N we rotate the whole cluster, otherwise a random atom is moved
              ! 1/(NMolecules + 1)  is the probability of attempting a rotation of the whole cluster
              IF ( SelectedMolecule > NMolecules ) THEN
   
                 DO k=1,3
                    DO i=1,NMolecules
                       rf(t,i,1,k)=r(t,i,1,k)
                    ENDDO
                 ENDDO
   
                 DO Molecule_b=1,NMolecules
                    !E(a,a)=0
                    Energy_ab_ini(t,Molecule_b,Molecule_b)= 0.0
   
                    !updating the pair energy contribution
                    Energy_2a_ini(t,Molecule_b) = Energy_2a_fin(t,Molecule_b)
   
                    DO Molecule_a=Molecule_b+1,NMolecules
                       !computing the pair energy E(a,b)
                       Energy_ab_ini(t,Molecule_a,Molecule_b)= Energy_ab_fin(t,Molecule_a,Molecule_b)
   
                       !E(a,b)=E(b,a)
                       Energy_ab_ini(t,Molecule_b,Molecule_a)= Energy_ab_ini(t,Molecule_a,Molecule_b)
                    ENDDO
                 ENDDO
   
                 !counting the number of accepted configuration after movement
                 Rot_accepted_eq(t)=Rot_accepted_eq(t)+1             !to calculate the acceptance ratio
                 adj_Rot_accepted_eq(t)=adj_Rot_accepted_eq(t)+1     !to update angle size, MaxRotation for the whole cluster
   
              ELSE  !if the whole structure is not rotated, we move one atom randomly
   
                 !Storing the new configuration 
                 DO k=1,3
                    rf(t,SelectedMolecule,1,k)=r(t,SelectedMolecule,1,k)
                 ENDDO
                 
                 !updating the pair energy contribution
                 Energy_2a_ini(t,SelectedMolecule) = Energy_2a_fin(t,SelectedMolecule)
                 
                 DO Molecule_a=1,NMolecules
                 
                    !updating the pair energy contribution for i-th
                    Energy_2a_ini(t,Molecule_a) = Energy_2a_ini(t,Molecule_a) + Energy_ab_fin(t,Molecule_a,SelectedMolecule) &
                                                                              - Energy_ab_ini(t,Molecule_a,SelectedMolecule)
                 
                    !computing the pair energy E(a,b)
                    Energy_ab_ini(t,Molecule_a,SelectedMolecule)= Energy_ab_fin(t,Molecule_a,SelectedMolecule)
                    
                    !E(a,b)=E(b,a)
                    Energy_ab_ini(t,SelectedMolecule,Molecule_a)= Energy_ab_ini(t,Molecule_a,SelectedMolecule)
                 
                 ENDDO
   
                 !counting the number of accepted configuration after movement
                 Mov_accepted_eq(t)=Mov_accepted_eq(t)+1             !to calculate the acceptance ratio
                 adj_Mov_accepted_eq(t)=adj_Mov_accepted_eq(t)+1             !to update step size maxDisplacement of CM
   
              ENDIF
   
              !-Write the Configuration FILE XYZ for the LOWEST energy config 
              !$OMP CRITICAL
              IF ( Energy_Tot(t) <= Energy_Tot(0) ) THEN                              !Energy_Tot(0) is the starting config Energy
   
                 Energy_Tot(0)=Energy_Tot(t)                                          !Saving the lowest energy value and coordinates
                 OPEN(UNIT=1200,FILE="lowest_Energy.xyz",STATUS='REPLACE')
   
                 WRITE(energy_tot_char,"(ES18.10E2)") Energy_Tot(0)*JtokJmol 
                 WRITE(t_char,"(F5.1)") Temperature(t)
   
                 comment_line='Energy (kJ/mol): '//TRIM(ADJUSTL(energy_tot_char))//' at '//TRIM(ADJUSTL(t_char))//' K'
                 FILEunit=1200
   
                 CALL DrawMolecularConfiguration( ZNucleus,radiusSBC,NMolecules,NAtoms,  &
                                                  Max_Atoms,MiddlePoint,NTrajectories,rf,&
                                                  t,FILEunit,comment_line                   )
   
                 CLOSE(1200)
   
             ENDIF
             !$OMP END CRITICAL   

           ELSE !if the NEW config is not accepted
   
              !if SM > N we rotate the whole cluster, otherwise a random atom is moved
              ! 1/(NMolecules + 1)  is the probability of attempting a rotation of the whole cluster
              IF ( SelectedMolecule > NMolecules ) THEN
   
                 !recovering the old configuration 
                 DO k=1,3
                    DO i=1,NMolecules
                       r(t,i,1,k)=rf(t,i,1,k)
                    ENDDO
                 ENDDO
   
                 DO Molecule_b=1,NMolecules
                    !E(a,a)=0
                    Energy_ab_fin(t,Molecule_b,Molecule_b)= 0.0
   
                    !recovering the pair energy contribution
                    Energy_2a_fin(t,Molecule_b) = Energy_2a_ini(t,Molecule_b)
   
                    DO Molecule_a=Molecule_b+1,NMolecules
                       !recovering the pair energy E(a,b)
                       Energy_ab_fin(t,Molecule_a,Molecule_b)= Energy_ab_ini(t,Molecule_a,Molecule_b)
                       !E(a,b)=E(b,a)
                       Energy_ab_fin(t,Molecule_b,Molecule_a)= Energy_ab_ini(t,Molecule_a,Molecule_b)
                    ENDDO
                 ENDDO
   
              ELSE  !if the whole structure is not rotated, we move one atom randomly
   
                 !recovering the old configuration 
                 DO k=1,3
                    r(t,SelectedMolecule,1,k)=rf(t,SelectedMolecule,1,k)
                 ENDDO
                 
                 !recovering the pair energy contribution
                 Energy_2a_fin(t,SelectedMolecule) = Energy_2a_ini(t,SelectedMolecule)
                 
                 DO Molecule_a=1,NMolecules
                    !recovering the pair energy E(a,b)
                    Energy_ab_fin(t,Molecule_a,SelectedMolecule)= Energy_ab_ini(t,Molecule_a,SelectedMolecule)
                    !E(a,b)=E(b,a)
                    Energy_ab_fin(t,SelectedMolecule,Molecule_a)= Energy_ab_ini(t,Molecule_a,SelectedMolecule)
                 ENDDO
   
              ENDIF
           ENDIF
   
        ENDDO !extra loop over the system size
   
        !-adjust step size 'maxDisplacement' after 100 steps
        IF ( MOD(MCstep,100) == 0 ) THEN                        !100 is the number of steps after which step size is updated
           Adjust=REAL(adj_Mov_accepted_eq(t))/REAL(100*NMolecules)       !Denominator is multiplied by NMolecules because the EXTRA loop over the number of molecules
           IF ( Adjust < 0.4 ) THEN                             !number of accepted steps/number of total steps should   
              maxDisplacement(t)=maxDisplacement(t)*0.9         !lie in [0.4,0.6]
           ELSEIF ( Adjust > 0.6 ) THEN
              maxDisplacement(t)=maxDisplacement(t)*1.1
           ENDIF
           !maxDisplacement must be less than Boundary Sphere diameter
           IF ( maxDisplacement(t) >= 2.0*radiusSBC ) THEN
              maxDisplacement(t) = 1.9*radiusSBC
           ENDIF
           adj_Mov_accepted_eq(t)=0                             !It restarts to count only the conf accepted every 100 steps
   
           !-adjust step size 'MaxRotation' after 100 steps
           Adjust=REAL(adj_Rot_accepted_eq(t))/100.0            !The probability of any rotation is 1/(n+1); i.e., NMolecule loop does not matter 
           IF ( Adjust < 0.4 ) THEN                             !number of accepted steps/number of total steps should   
              MaxRotation(t)=MaxRotation(t)*0.9                 !lie in [0.4,0.6]
           ELSEIF ( Adjust > 0.6 ) THEN
              MaxRotation(t)=MaxRotation(t)*1.1
           ENDIF
           !maxRotation must be less than 90 degrees 
           IF ( MaxRotation(t) >= 0.5*PI ) THEN
              MaxRotation(t) = 0.5*PI                  !Any rotation larger than 90 deg makes no sense in this contest
           ENDIF
           adj_Rot_accepted_eq(t)=0                    !It restarts to count only the conf accepted every 100 steps
        ENDIF

        ENDDO !loop over temperatures/trajectories
        !$OMP END PARALLEL DO
   
        !----------------------------------------------------------------------------------------
        !Replica exchange MC sampling (Parallel tempering). The exchange is between ExchangeTrajectory and ExchangeTrajectory+1
      
        IF ( 0.1 > RandomNumber(Seed) ) THEN                                                ! 0.1 is the probability of attempting an exchange
           ExchangeTrajectory=1+FLOOR( (NTrajectories-1)*RandomNumber(Seed) )               ! index, random trajectory to exchange in range [1:NTrajectories-1] 
   
           !-Total number of attempted exchanges
           Exchange_tot(ExchangeTrajectory)  =Exchange_tot(ExchangeTrajectory)  +1
           Exchange_tot(ExchangeTrajectory+1)=Exchange_tot(ExchangeTrajectory+1)+1
   
           DeltaBeta=beta(ExchangeTrajectory) - beta(ExchangeTrajectory+1)
           DeltaEnergy=Energy_Tot(ExchangeTrajectory) - Energy_Tot(ExchangeTrajectory+1)
           ExchangeAccept=MIN( 1.E0,DEXP(DeltaBeta*DeltaEnergy) )
   
           ! exchange accepted
           IF ( ExchangeAccept > RandomNumber(Seed) ) THEN
              !-Number of accepted exchanges
              Exchange_acc(ExchangeTrajectory)  =Exchange_acc(ExchangeTrajectory)  +1
              Exchange_acc(ExchangeTrajectory+1)=Exchange_acc(ExchangeTrajectory+1)+1
              ! swap configurations for an accepted configuration (rf)
              DO i=1,NMolecules
                 DO j=1,NAtoms(i)
                    DO k=1,3
                       !for r
                       rex(ExchangeTrajectory,i,j,k)=r(ExchangeTrajectory,i,j,k)          !rex is a temporary variable, so rex=rf(j)
                       r(ExchangeTrajectory,i,j,k)  =r(ExchangeTrajectory+1,i,j,k)        !r(j)  =r(j+1)
                       r(ExchangeTrajectory+1,i,j,k)=rex(ExchangeTrajectory,i,j,k)        !r(j+1)=rex
                       !for rf
                       rex(ExchangeTrajectory,i,j,k) =rf(ExchangeTrajectory,i,j,k)        !rex is a temporary variable, so rex=rf(j)
                       rf(ExchangeTrajectory,i,j,k)  =rf(ExchangeTrajectory+1,i,j,k)      !rf(j)  =rf(j+1)
                       rf(ExchangeTrajectory+1,i,j,k)=rex(ExchangeTrajectory,i,j,k)       !rf(j+1)=rex
                    ENDDO
                 ENDDO
              ENDDO
              ! swap Energies
              Energy_ex                       =Energy_Tot(ExchangeTrajectory)
              Energy_Tot(ExchangeTrajectory)  =Energy_Tot(ExchangeTrajectory+1)
              Energy_Tot(ExchangeTrajectory+1)=Energy_ex
   
              DO Molecule_a=1,NMolecules
               Energy_2a_ini_ex                               = Energy_2a_ini(ExchangeTrajectory,Molecule_a)
               Energy_2a_ini(ExchangeTrajectory,Molecule_a)   = Energy_2a_ini(ExchangeTrajectory+1,Molecule_a)
               Energy_2a_ini(ExchangeTrajectory+1,Molecule_a) = Energy_2a_ini_ex
               DO Molecule_b=Molecule_a+1,NMolecules
                Energy_ab_ini_ex                                         = Energy_ab_ini(ExchangeTrajectory,Molecule_a,Molecule_b)
                Energy_ab_ini(ExchangeTrajectory,Molecule_a,Molecule_b)  = Energy_ab_ini(ExchangeTrajectory+1,Molecule_a,Molecule_b)
                Energy_ab_ini(ExchangeTrajectory+1,Molecule_a,Molecule_b)= Energy_ab_ini_ex
   
                !E(a,b)=E(b,a) for t and t+1
                Energy_ab_ini(ExchangeTrajectory,Molecule_b,Molecule_a)  = Energy_ab_ini(ExchangeTrajectory,Molecule_a,Molecule_b)
                Energy_ab_ini(ExchangeTrajectory+1,Molecule_b,Molecule_a)= Energy_ab_ini(ExchangeTrajectory+1,Molecule_a,Molecule_b)
               ENDDO
              ENDDO
   
           ENDIF
        ENDIF
   
     ENDDO !loop over MC steps

     !----------------------------------------------------------------------------------------
   
     TotalAccepted_Mov =0
     TotalAccepted2_Mov=0     !after eq
     TotalAccepted_Rot =0
     TotalAccepted2_Rot=0     !after eq
     TotalExchanged=0
     DO t=1,NTrajectories        
        TotalAccepted_Mov =TotalAccepted_Mov +Mov_accepted_eq(t)
        TotalAccepted2_Mov=TotalAccepted2_Mov+Mov_accepted(t)    !after eq
        TotalAccepted_Rot =TotalAccepted_Rot +Rot_accepted_eq(t)
        TotalAccepted2_Rot=TotalAccepted2_Rot+Rot_accepted(t)    !after eq
        TotalExchanged=TotalExchanged+Exchange_acc(t)
     ENDDO

     DO t=1,NTrajectories
        MaxEnergyTotal=Max( MaxEnergyTotal,MaxEnergy(t) )
        MinEnergyTotal=Min( MinEnergyTotal,MinEnergy(t) )
     ENDDO


     WRITE(10000+check,'(A)')''
     WRITE(10000+check,'(A,/)')' The lowest and highest Energies are: '
     DO t=1,NTrajectories
        WRITE(10000+check,'(I3,3X,F5.1,T20,2(ES15.6E2,5X))')t,Temperature(t),&
                            MinEnergy(t)*JtokJmol,MaxEnergy(t)*JtokJmol
     ENDDO
     WRITE(10000+check,'(/,A,T20,2(ES15.6E2,5X))')' Global Min and Max Energy: ',&
                           MinEnergyTotal*JtokJmol,MaxEnergyTotal*JtokJmol
     WRITE(10000+check,'(A)')''

     WRITE(10000+check,'(A,T40,I20)')' Number of steps: ',MCcycles
     WRITE(10000+check,'(A,T40,I20)')' Equilibration: ',Nequilibration
     WRITE(10000+check,'(A,T40,I20)')' Total Movement / Trajectory: ',NMolecules*check_MCcycles
     WRITE(10000+check,'(A,T40,I20)')' Total Accepted: ',&
                               TotalAccepted_Mov+TotalAccepted_Rot+TotalAccepted2_Mov+TotalAccepted2_Rot
     WRITE(10000+check,'(A,T40,I20)')' Total Accepted (after eq): ',TotalAccepted2_Mov+TotalAccepted2_Rot
     WRITE(10000+check,'(A,T40,I20)')' Total Exchanged: ',TotalExchanged
     WRITE(10000+check,'(A)')''
     WRITE(10000+check,'(A)')" 1/(n+1)  is the probability of attempting a rotation"&
                            //" of the whole cluster, where 'n' is number of atoms"
     WRITE(10000+check,'(A,F5.1,A)')' The DISPLACEMENT inside a boundary sphere is within 0 and ',&
                                   2.0*radiusSBC,' Angtrom'
     WRITE(10000+check,'(A)')' The angle of ROTATION of the whole cluster is within 0 and 90 degree'
     WRITE(10000+check,'(A)')''

     WRITE(10000+check,'(A,/)') &
     ' #    T (K)  Max Disp. (Ang)   Accept      Ratio      '//&
                  'Max Rot. (deg)    Accept      Ratio    Tot accep rat   Exchang    Ex. ratio'
     DO t=1,NTrajectories
        WRITE(10000+check,"(T2,I2,T7,F5.1,T17,F9.5,T29,I9,T43,F7.4,T58,F10.4,T70,I9,T84,F7.4,T96,F7.4,T109,I8,T122,F7.4)")&
        t,Temperature(t),maxDisplacement(t),&
        Mov_accepted_eq(t),REAL(Mov_accepted_eq(t))/REAL(NMolecules*check_MCcycles),&
        MaxRotation(t)*(180.0/PI),&
        Rot_accepted_eq(t),REAL(Rot_accepted_eq(t))/REAL(check_MCcycles),&
        REAL(Mov_accepted_eq(t)+Rot_accepted_eq(t))/REAL(NMolecules*check_MCcycles),&
        Exchange_acc(t),REAL(Exchange_acc(t))/REAL(Exchange_tot(t))
     ENDDO


     WRITE(10000+check,'(A)')''
     WRITE(10000+check,'(A)')'****************************************************************'
     CALL DATE_AND_TIME(today,now)
     WRITE(10000+check,"(12A)")' Date and time: ',&
                     today(1:4),"-",today(5:6),"-",today(7:8), ", ",&
                     now(1:2),":",  now(3:4),":",  now(5:6)
     WRITE(10000+check,'(A)')''
     CALL CPU_TIME(EndTime)
     RunTime=EndTime-StartTime

     !$ wtime_end = omp_get_wtime()
     !$ RunTime = wtime_end - wtime_start

     days   =FLOOR( RunTime/86400.0 )
     hours  =FLOOR( (RunTime - REAL(days*86400)) /3600.0 )
     minutes=FLOOR( (RunTime - REAL(days*86400+hours*3600)) /60.0 )
     seconds=RunTime - REAL( days*86400+hours*3600+minutes*60 )
     WRITE(10000+check,'(A16,3(I4,X,A),X,F7.4,X,A)')' CPU time used: ',&
           days,'days',hours,'hours',minutes,'minutes',seconds,'seconds'
     WRITE(10000+check,'(A)')'****************************************************************'

     FLUSH(10000+check)  !OUTPUT FILE
     CLOSE(UNIT=10000+check)  !OUTPUT FILE

     !----------------------------------------------------------------------------------------
     !-STATUS, this FILE has the MC evolution
     WRITE(1000,'(A)')
     CALL CPU_TIME(PartialEndRunTime)
     RunTime=PartialEndRunTime-PartialStartRunTime   !Delta time
     PartialStartRunTime=PartialEndRunTime           !redefining starting time here to start the counting again

     !$ part_wtime_end = omp_get_wtime()
     !$ RunTime = part_wtime_end - part_wtime_start
     !$ part_wtime_start = part_wtime_end 

     CALL DATE_AND_TIME(today,now)
     WRITE(1000,"(12A)")'Date and time: ',&
                        today(1:4),"-",today(5:6),"-",today(7:8), ", ",&
                        now(1:2),":",  now(3:4),":",  now(5:6)

     WRITE(MCstep_char,"(I10)")     MCstep-1

     WRITE(1000,'(A,T25,I10,A6,I10)')'Steps for equilibration: '//TRIM(ADJUSTL(MCstep_char))//&
                                     ' of '//TRIM(ADJUSTL(Nequilibration_char))
     WRITE(1000,'(A,T25,F5.1,A)')'Equilibration: ',(REAL(check)/REAL(Ncheckpoints_eq))*100.0,' %'
     WRITE(1000,'(A,T25,ES12.5E2,A)') 'CPU time used: ',RunTime,' sec'
     WRITE(1000,'(A)')'Checkpoint FILEs:   checkpoint/equilibration/'//TRIM(ADJUSTL(check_char))//'/'

     IF ( check == Ncheckpoints_eq ) THEN
        CALL CPU_TIME(EndTime)
        RunTime=EndTime-StartTime
        
        !$ wtime_end = omp_get_wtime()
        !$ RunTime = wtime_end - wtime_start
        
        days   =FLOOR( RunTime/86400.0 )
        hours  =FLOOR( (RunTime - REAL(days*86400)) /3600.0 )
        minutes=FLOOR( (RunTime - REAL(days*86400+hours*3600)) /60.0 )
        seconds=RunTime - REAL( days*86400+hours*3600+minutes*60 )
        WRITE(1000,'(A)')
        WRITE(1000,'(A)')'*** EQUILIBRATION DONE ***'
        WRITE(1000,'(A)') '-----------------------------------------------------------------'
        WRITE(1000,'(A16,3(I4,X,A),X,F7.4,X,A)')' CPU time used: ',&
              days,'days',hours,'hours',minutes,'minutes',seconds,'seconds'
        WRITE(1000,'(A)') '-----------------------------------------------------------------'
     ENDIF
     FLUSH(1000)

     !----------------------------------------------------------------------------------------
     !- Saving FILEs to restart
     OPEN(UNIT=20000+check,FILE='restart_'//output//'_eq_'//TRIM(ADJUSTL(check_char)),STATUS='REPLACE')

     WRITE(20000+check,"(A)") '-----------------------------------------------------------------'
     WRITE(20000+check,"(12A)")'Date and time: ',&
                                today(1:4),"-",today(5:6),"-",today(7:8), ", ",&
                                now(1:2),":",  now(3:4),":",  now(5:6)
     WRITE(20000+check,"(A)") '-----------------------------------------------------------------'
     WRITE(20000+check,"(A)") 'RESTART SETTING'
     WRITE(20000+check,"(A)") '-----------------------------------------------------------------'
     WRITE(20000+check,"(A)")   
     WRITE(20000+check,"(A,T30,A20)") "MC cycles/sampling: ",TRIM(ADJUSTL(MCcycles_char))
     WRITE(20000+check,"(A,T30,I20)") "Checkpoints for sampling: ",Ncheckpoints_mc
     WRITE(20000+check,"(A,T30,A20)") "Equilibracion: ",TRIM(ADJUSTL(Nequilibration_char))
     WRITE(20000+check,"(A,T30,I20)") "Checkpoints for Equilibracion: ",Ncheckpoints_eq
     WRITE(20000+check,"(A,T30,I20)") "Current point: ",MCstep-1
     WRITE(20000+check,'(A,T30,F20.1)') "Equilibration Evolution (%): ",(REAL(check)/REAL(Ncheckpoints_eq))*100.0
     WRITE(20000+check,'(A,T30,A20)') "Sampling Evolution (%): ","0.0"

     WRITE(20000+check,"(A,T30,I20)") "Seed: ",Seed
     WRITE(20000+check,"(A)")   

     !- Initia States
     WRITE(20000+check,"(A)") '-----------------------------------------------------------------'
     WRITE(20000+check,"(A)") ' Initia States '
     WRITE(20000+check,"(A)") '-----------------------------------------------------------------'
     WRITE(20000+check,"(A)") '*** Energies in Joules '   
     WRITE(20000+check,"(A)")   
     WRITE(20000+check,"(A,T40,ES28.20E2)") 'Total Energy: ',Energy_Tot(0) 
     WRITE(20000+check,"(A)")   
     DO Molecule_a=1,NMolecules
        WRITE(20000+check,"(A32,I4,A4,T40,ES28.20E2)") 'Pair Energy for atom: ',&
                                                       Molecule_a,' is ',Energy_2a_ini(0,Molecule_a)
        DO Molecule_b=Molecule_a+1,NMolecules
           WRITE(20000+check,"(A25,2(I4,A5),T40,ES28.20E2)") '2b Energy for ',Molecule_a,' and ',&
                                                              Molecule_b,' is ',Energy_ab_ini(0,Molecule_a,Molecule_b)
        ENDDO
     WRITE(20000+check,"(A)")   
     ENDDO

     !- Trajectories
     DO t=1,NTrajectories
        WRITE(20000+check,"(A)") '-----------------------------------------------------------------'
        WRITE(20000+check,"(A,T20,I3,A4,T40,F10.5)") ' Temperature (K): ',t,' is ',Temperature(t)
        WRITE(20000+check,"(A)") '-----------------------------------------------------------------'
        WRITE(20000+check,"(A,T40,F30.25)") ' Max Disp. (Ang): ',maxDisplacement(t)
        WRITE(20000+check,"(A,T40,F30.25)") ' Max Rot. (rad): ',MaxRotation(t)!*(180.0/PI)
        WRITE(20000+check,"(A)")   
        !Structural parameters
        DO i=1,NMolecules
           WRITE(20000+check,FMT3) ZNucleus(i,1),rf(t,i,1,1),rf(t,i,1,2),rf(t,i,1,3)
        ENDDO           
        WRITE(20000+check,"(A)")   
        !Acceptance criterion, adjusting step size
        WRITE(20000+check,"(A,T40,I20)") 'Total Movements Accepted: ',Mov_accepted_eq(t)
        WRITE(20000+check,"(A,T40,I20)") 'Movements Accepted (per 100 Cycles): ',adj_Mov_accepted_eq(t)
        WRITE(20000+check,"(A,T40,I20)") 'Total Rotations Accepted: ',Rot_accepted_eq(t)
        WRITE(20000+check,"(A,T40,I20)") 'Rotations Accepted (per 100 Cycles): ',adj_Rot_accepted_eq(t)
        WRITE(20000+check,"(A,T40,I20)") 'Attempts for Exchange: ',Exchange_tot(t)
        WRITE(20000+check,"(A,T40,I20)") 'Number of Exchange Accepted: ',Exchange_acc(t)
        WRITE(20000+check,"(A)")   
        !Energy
        WRITE(20000+check,"(A)") '-----------------------'
        WRITE(20000+check,"(A)") '*** Energies in Joules '   
        WRITE(20000+check,"(A)")   
        WRITE(20000+check,"(A,T40,ES28.20E2)") 'Total Energy: ',Energy_Tot(t) 
        WRITE(20000+check,"(A,T40,ES28.20E2)") 'Min Energy: ',MinEnergy(t) 
        WRITE(20000+check,"(A,T40,ES28.20E2)") 'Max Energy: ',MaxEnergy(t) 
        WRITE(20000+check,"(A)")   
        DO Molecule_a=1,NMolecules
           WRITE(20000+check,"(A,I4,A4,T40,ES28.20E2)") 'Pair Energy for Atom ',Molecule_a,' is ',Energy_2a_ini(t,Molecule_a)
           DO Molecule_b=Molecule_a+1,NMolecules
              WRITE(20000+check,"(A23,I4,A5,I4,A4,T40,ES28.20E2)") '2b Energy for ',Molecule_a,' and ',&
                                                                    Molecule_b,' is ',Energy_ab_ini(t,Molecule_a,Molecule_b)
           ENDDO
           WRITE(20000+check,"(A)")   
        ENDDO

     ENDDO

     !----------------------------------------------------------------------------------------
     !Equilibration DONE
     IF ( check == Ncheckpoints_eq ) THEN
        !It is the 'Bin width'
        DeltaEHistogram=(MaxEnergyTotal - MinEnergyTotal)/(nBins)                          !nBins=500 by default for all systems

        WRITE(20000+check,"(A)")   
        WRITE(20000+check,"(A,T40,ES18.10E2)") 'Total Min Energy after Equilibration: ',MinEnergyTotal
        WRITE(20000+check,"(A,T40,ES18.10E2)") 'Total Max Energy after Equilibration: ',MaxEnergyTotal
        WRITE(20000+check,"(A,T40,ES18.10E2)") 'Bin Width for the Histogram: ',DeltaEHistogram
     ENDIF

     WRITE(20000+check,"(A)")   
     WRITE(20000+check,"(A)") '-----------------------------------------------------------------'
     WRITE(20000+check,"(A)") '*** END OF FILE ***'  
     WRITE(20000+check,"(A)") '-----------------------------------------------------------------'
     CLOSE(20000+check)

     !----------------------------------------------------------------------------------------
     !-Sorting FILEs
     IF ( check == 1 ) CALL SYSTEM("mkdir -p checkpoint/ ;"//&
                                   "mkdir -p checkpoint/equilibration/")     !this dir is only created in the beginning
     CALL SYSTEM("mkdir -p checkpoint/equilibration/"//TRIM(ADJUSTL(check_char)))
     CALL SYSTEM("cp -f "//input//" checkpoint/equilibration/"//TRIM(ADJUSTL(check_char)))
     CALL SYSTEM("mv -fu "//output//"_"//TRIM(ADJUSTL(check_char))//&
                 "ckp_eq.out checkpoint/equilibration/"//TRIM(ADJUSTL(check_char)))
!     CALL SYSTEM("mv -fu restart_"//output//"_eq_"//TRIM(ADJUSTL(check_char))//&
!                 " checkpoint/equilibration/"//TRIM(ADJUSTL(check_char)))

     CALL SYSTEM("mv -fu restart_"//output//"_eq_"//TRIM(ADJUSTL(check_char))//" checkpoint/")

     !- Deleting old checkpoint restart-files
     WRITE(check_char,"(I10)") check-1
     CALL SYSTEM("rm -f checkpoint/restart_"//output//"_eq_"//TRIM(ADJUSTL(check_char)) )

  ENDDO   !checkpoint loop

!--------------------------------------------------------------------------------------------

  !- We have TWO restarting options: cold and warm 

  !- COLD restart: starting from the configurations obtained by a previous simulation run,
  !                but MC cycle (eq and sampling) run as a NEW/different simulation
  !-               For exemple, when you want to extend a simulation; equlibration, sampling,
  !-               or both

  !- WARM restart: basically a job recovery measure. If a job crashed, you can restart
  !-               it from where it was saved (any checkpoint), and obtain essentially the same result. 
  !-               Another alternative is to re-run a simulation increasing the number of checkpoints

!---------------------------------------------------------------------------------------------------
!-Sampling Loop
!---------------------------------------------------------------------------------------------------

  check_restart_mc = 1 + NINT(sampling_percent_restart*REAL(Ncheckpoints_mc)/100.0)

  !- COLD restart: extending the number of steps for sampling or
  !- you can restart a NEW/different simulation from the last checkpoint
  !- saved for a previous simulation
  IF ( sampling_percent_restart > 99.999 ) THEN
     check_restart_mc = 1
  ENDIF

  check_split_mc = MCcycles/Ncheckpoints_mc

  DO check=check_restart_mc,Ncheckpoints_mc

     !skipping Sampling loop
     IF ( MCcycles <= 1 ) EXIT

     WRITE(check_char,"(I10)") check

     !--------------------------------------------------------------------------------------------
     !-Start point of MC cycle after equilibration
     check_MCstep = 1 + check_split_mc*(check - 1)
     check_MCcycles = check_split_mc*(check)

     WRITE(check_MCcycles_char,"(I10)") check_MCcycles

     !--------------------------------------------------------------------------------------------
     !-Writing the output FILE

     IF ( check < Ncheckpoints_mc ) THEN
        check2=check !doesn't change at all
        OPEN(UNIT=10000+check2,FILE=output//"_"//TRIM(ADJUSTL(check_char))//"ckp.out",STATUS='REPLACE')

        WRITE(10000+check2,'(A)')'****************************************************************'
        WRITE(10000+check2,'(12A)')' Execution time: ',&
                                         today(1:4),"-",today(5:6),"-",today(7:8), ", ",&
                                         now(1:2),":",  now(3:4),":",  now(5:6)
        WRITE(10000+check2,'(A)')'****************************************************************'
        WRITE(10000+check2,'(A)')''
        WRITE(10000+check2,'(A)')' Contents of the input FILE'
        WRITE(10000+check2,'(A)')'----------------------------------------------------------------'
        WRITE(10000+check2,"(A,T45,I20)")' Number of threads availables: ',aval_threads
        WRITE(10000+check2,"(A,T45,I20)")' Number of threads used: ',threads
        WRITE(10000+check2,"(A,T45,A20)")' Total MC Steps per Trajectory: ',TRIM(ADJUSTL(MCcycles_char))  
        WRITE(10000+check2,"(A,T45,A20)")' Equilibration steps per Trajectory: ',TRIM(ADJUSTL(Nequilibration_char))  
        WRITE(10000+check2,'(A,T45,F20.2)')' Minimum Temperature: ',To
        WRITE(10000+check2,'(A,T45,F20.2)')' Maximum Temperature: ',Tf
        WRITE(10000+check2,'(A,T45,I20)')' Number of Temperatures/Trajectories: ',NTrajectories
        WRITE(10000+check2,"(A,T45,F20.2)")' Spherical BC, radius (Angst): ',radiusSBC
        WRITE(10000+check2,"(A,T45,A20)")  ' The system has (Molecules/Atoms): ',TRIM(ADJUSTL(NMolecules_char))
        WRITE(10000+check2,"(A,T45,F20.2)")' Max Displacement (Angst): ',MaxDisplacementRcm
        WRITE(10000+check2,"(A,T45,F20.2)")' Max Rotation (degrees): ',MaxAngle*(180.0/PI) 
        WRITE(10000+check2,'(A,T25,A40)')' Potential Energy: ',TRIM(ADJUSTL(EnergyName))
        WRITE(10000+check2,"(A,T47,ES18.10E2)")' Initial Energy (kJ/mol): ',Energy_Tot(0)*JtokJmol
        WRITE(10000+check2,"(A,T47,ES18.10E2)")' The Boltzmann constant (J/K): ',kB
        WRITE(10000+check2,"(A,T45,I20)")' Seed for Random Number generator: ',Seed
        WRITE(10000+check2,"(A,T45,I20)")' A configuration saved every (per traj): ',SaveConfig
        WRITE(10000+check2,'(A)')''

        IF ( restart == 'yes' ) THEN
           WRITE(10000+check2,"(A,T30,A35)") ' Simulation restarting from: ',TRIM(ADJUSTL(restart_FILE))
           WRITE(10000+check2,'(A,T59,F6.2)') ' Equilibration evolution was (%): ',eq_percent_restart
           WRITE(10000+check2,'(A,T45,I20)') '    of a total number of Eq steps: ',Nequilibration_restart
           WRITE(10000+check2,'(A,T59,F6.2)') ' Sample evolution was (%): ',sampling_percent_restart
           WRITE(10000+check2,'(A,T45,I20)') '    of a total number of MC steps: ',MCcycles_restart
           WRITE(10000+check2,'(A)')
        ENDIF
        
        !coordinates
        WRITE(10000+check2,'(A)')' The initial configuration is the same for every trajectory'
        WRITE(10000+check2,'(A)')" Coordinates written at 'initial_configuration.xyz'"
        WRITE(10000+check2,'(A)')''
        WRITE(10000+check2,'(A)')'----------------------------------------------------------------'
        WRITE(10000+check2,'(A,F5.1,A)')' CHECKPOINT OUTPUT ',(REAL(check)/REAL(Ncheckpoints_mc))*100.0,' %,  '//&
                                       TRIM(ADJUSTL(check_MCcycles_char))//' CYCLES DONE OF '//&
                                       TRIM(ADJUSTL(MCcycles_char))
        WRITE(10000+check2,'(A)')'----------------------------------------------------------------'
        
        !Flushing memory buffer
        CALL FLUSH(10000+check2)  !writing the general output FILE 

     ELSE
        check2=0  !for the last checkpoint the output FILE is the 'FINAL' which has the results
     ENDIF

     !----------------------------------------------------------------------------------------

     DO MCstep=check_MCstep,check_MCcycles          !loop over MC steps


        !$OMP PARALLEL DO DEFAULT (SHARED) &
        !$OMP FIRSTPRIVATE (Seed) &
        !$OMP PRIVATE (t,SelectedMolecule,r_rot,dist,theta,indexAngle,indexDistance,Seed2 ) &
        !$OMP PRIVATE (Molecule_a,Molecule_b,DeltaEnergy,r_move,acceptance,k,i,n,Adjust)
        DO t=1,NTrajectories         !loop over temperatures/trajectories

           Seed2 = Seed + t
   
        DO n=1,NMolecules+1          !EXTRA loop over number of Molecules in the MC cluster. N+1 because the rotation of the whole cluster

           !Randomly selected the molecule which is translated and rotated
           SelectedMolecule=FLOOR( 1 + (NMolecules + Rotate)*RandomNumber(Seed2) )                  !if Ran()=0, then SM=1
   
           !if SM > N we rotate the whole cluster, otherwise a random atom is moved
           ! 1/(NMolecules + 1)  is the probability of attempting a rotation of the whole cluster
           IF ( SelectedMolecule > NMolecules ) THEN
   
              CALL RotateCluster( rf,NMolecules,MaxRotation(t),Seed2,NTrajectories,t,MiddlePoint,&   !IN
                                  r_rot                                                          )  !OUT
   
              DO k=1,3
                 DO i=1,NMolecules
                    r(t,i,1,k)=r_rot(i,k)
                 ENDDO
              ENDDO
   
              DO Molecule_a=1,NMolecules
                 DO Molecule_b=Molecule_a+1,NMolecules
                    dist=SQRT( (r(t,Molecule_a,1,1) - r(t,Molecule_b,1,1))**2+&   !r_ab=SQRT[ (xa-xb)^2+
                               (r(t,Molecule_a,1,2) - r(t,Molecule_b,1,2))**2+&   !           (ya-yb)^2+
                               (r(t,Molecule_a,1,3) - r(t,Molecule_b,1,3))**2  )  !           (za-zb)^2 ]

                    !From dot product, where |k|=1 and k . r12 = z1-z2, 0<=theta<=90
                    theta=ACOS( ABS(r(t,Molecule_a,1,3) - r(t,Molecule_b,1,3))/dist )   !Cos(x)=(z1-z2)/|r12|

                    indexAngle   =NINT( (theta-MinAngle)/SpacingAngles )
                    indexDistance=NINT( (dist-MinDistance)/SpacingDistances )

                    !computing the pair energy E(a,b)
                    Energy_ab_fin(t,Molecule_a,Molecule_b)= Energy_LookUpTable(indexAngle,indexDistance)
                 ENDDO
              ENDDO

              Energy_Tot_fin(t)=0.0                 !Initialization

              DO Molecule_a=1,NMolecules
                 !Initialization
                 Energy_2a_fin(t,Molecule_a)= 0.0
                 !E(a,a)=0
                 Energy_ab_fin(t,Molecule_a,Molecule_a)= 0.0

                 DO Molecule_b=Molecule_a+1,NMolecules
                    !E(a,b)=E(b,a)
                    Energy_ab_fin(t,Molecule_b,Molecule_a)= Energy_ab_fin(t,Molecule_a,Molecule_b)

                    !Total energy
                    Energy_Tot_fin(t)= Energy_Tot_fin(t) + Energy_ab_fin(t,Molecule_a,Molecule_b)
                 ENDDO

                 !computing pair energy contribution for i-th
                 DO Molecule_b=1,NMolecules
                    Energy_2a_fin(t,Molecule_a) = Energy_2a_fin(t,Molecule_a) + Energy_ab_fin(t,Molecule_a,Molecule_b)
                 ENDDO
              ENDDO

              DeltaEnergy=Energy_Tot_fin(t) - Energy_Tot(t)
   
           ELSE  !if the whole structure is not rotated, we move one atom randomly
   
              CALL DisplaceAtoms( rf,NMolecules,maxDisplacement(t),SelectedMolecule, &   !IN 
                                  Seed2,radiusSBC2,MiddlePoint,NTrajectories,t,       &   !IN
                                  r_move                                              )  !OUT
              
              
              DO j=1,NAtoms(SelectedMolecule)
                 DO k=1,3
                    r(t,SelectedMolecule,j,k)=r_move(k)
                 ENDDO
              ENDDO
              
              DO Molecule_a=1,NMolecules
                 IF ( Molecule_a == SelectedMolecule ) CYCLE             !Avoiding to compute self interaction

                 dist=SQRT( (r(t,Molecule_a,1,1) - r(t,SelectedMolecule,1,1))**2+&   !r_ab=SQRT[ (xa-xb)^2+
                            (r(t,Molecule_a,1,2) - r(t,SelectedMolecule,1,2))**2+&   !           (ya-yb)^2+
                            (r(t,Molecule_a,1,3) - r(t,SelectedMolecule,1,3))**2  )  !           (za-zb)^2 ]

                 !From dot product, where |k|=1 and k . r12 = z1-z2, 0<=theta<=90
                 theta=ACOS( ABS(r(t,Molecule_a,1,3) - r(t,SelectedMolecule,1,3))/dist )   !Cos(x)=(z1-z2)/|r12|

                 indexAngle   =NINT( (theta-MinAngle)/SpacingAngles )
                 indexDistance=NINT( (dist-MinDistance)/SpacingDistances )

                 !computing the pair energy E(a,b)
                 Energy_ab_fin(t,Molecule_a,SelectedMolecule)= Energy_LookUpTable(indexAngle,indexDistance)
              ENDDO

              Energy_2a_fin(t,SelectedMolecule)= 0.0                 !Initialization

              DO Molecule_a=1,NMolecules
                 IF ( Molecule_a == SelectedMolecule ) CYCLE             !Avoiding to compute self interaction

                 !E(a,a)=0
                 Energy_ab_fin(t,Molecule_a,Molecule_a)= 0.0

                 !E(a,b)=E(b,a)
                 Energy_ab_fin(t,SelectedMolecule,Molecule_a)= Energy_ab_fin(t,Molecule_a,SelectedMolecule)

                 !computing total 2-body energy contribution for i-th
                 Energy_2a_fin(t,SelectedMolecule) = Energy_2a_fin(t,SelectedMolecule) &
                                                     + Energy_ab_fin(t,Molecule_a,SelectedMolecule)
              ENDDO

              DeltaEnergy=Energy_2a_fin(t,SelectedMolecule) - Energy_2a_ini(t,SelectedMolecule)
   
           ENDIF
   
           accept(t)=.FALSE.
   
           !Metropolis acceptance criterion
           IF ( DeltaEnergy <= 0.0 ) THEN
              acceptance=1.0
           ELSE 
              acceptance=DEXP( -DeltaEnergy*beta(t) )
           ENDIF
   
           IF ( acceptance > RandomNumber(Seed2) ) THEN
    
              accept(t)=.TRUE.
   
              !saving the total energy.
              Energy_Tot(t)= Energy_Tot(t) + DeltaEnergy
   
              !if SM > N we rotate the whole cluster, otherwise a random atom is moved
              ! 1/(NMolecules + 1)  is the probability of attempting a rotation of the whole cluster
              IF ( SelectedMolecule > NMolecules ) THEN
   
                 DO k=1,3
                    DO i=1,NMolecules
                       rf(t,i,1,k)=r(t,i,1,k)
                    ENDDO
                 ENDDO
   
                 DO Molecule_b=1,NMolecules
                    !E(a,a)=0
                    Energy_ab_ini(t,Molecule_b,Molecule_b)= 0.0
   
                    !updating the pair energy contribution
                    Energy_2a_ini(t,Molecule_b) = Energy_2a_fin(t,Molecule_b)
   
                    DO Molecule_a=Molecule_b+1,NMolecules
                       !computing the pair energy E(a,b)
                       Energy_ab_ini(t,Molecule_a,Molecule_b)= Energy_ab_fin(t,Molecule_a,Molecule_b)
   
                       !E(a,b)=E(b,a)
                       Energy_ab_ini(t,Molecule_b,Molecule_a)= Energy_ab_ini(t,Molecule_a,Molecule_b)
                    ENDDO
                 ENDDO
   
                 !counting the number of accepted configuration after movement
                 Rot_accepted(t)=Rot_accepted(t)+1             !to calculate the acceptance ratio
                 adj_Rot_accepted(t)=adj_Rot_accepted(t)+1     !to update angle size, MaxRotation for the whole cluster
   
              ELSE  !if the whole structure is not rotated, we move one atom randomly
   
                 !Storing the new configuration 
                 DO k=1,3
                    rf(t,SelectedMolecule,1,k)=r(t,SelectedMolecule,1,k)
                 ENDDO
                 
                 !updating the pair energy contribution
                 Energy_2a_ini(t,SelectedMolecule) = Energy_2a_fin(t,SelectedMolecule)
                 
                 DO Molecule_a=1,NMolecules
                 
                    !updating the pair energy contribution for i-th
                    Energy_2a_ini(t,Molecule_a) = Energy_2a_ini(t,Molecule_a) + Energy_ab_fin(t,Molecule_a,SelectedMolecule) &
                                                                              - Energy_ab_ini(t,Molecule_a,SelectedMolecule)
                 
                    !computing the pair energy E(a,b)
                    Energy_ab_ini(t,Molecule_a,SelectedMolecule)= Energy_ab_fin(t,Molecule_a,SelectedMolecule)
                    
                    !E(a,b)=E(b,a)
                    Energy_ab_ini(t,SelectedMolecule,Molecule_a)= Energy_ab_ini(t,Molecule_a,SelectedMolecule)
                 
                 ENDDO
   
                 !counting the number of accepted configuration after movement
                 Mov_accepted(t)=Mov_accepted(t)+1             !to calculate the acceptance ratio
                 adj_Mov_accepted(t)=adj_Mov_accepted(t)+1     !to update step size maxDisplacement of CM
   
              ENDIF
   
              !-Write the Configuration FILE XYZ for the LOWEST energy config 
              !$OMP CRITICAL   
              IF ( Energy_Tot(t) <= Energy_Tot(0) ) THEN                              !Energy_Tot(0) is the starting config Energy
   
                 Energy_Tot(0)=Energy_Tot(t)                                          !Saving the lowest energy value and coordinates
                 OPEN(UNIT=1200,FILE="lowest_Energy.xyz",STATUS='REPLACE')
   
                 WRITE(energy_tot_char,"(ES18.10E2)") Energy_Tot(0)*JtokJmol 
                 WRITE(t_char,"(F5.1)") Temperature(t)
   
                 comment_line='Energy (kJ/mol): '//TRIM(ADJUSTL(energy_tot_char))//' at '//TRIM(ADJUSTL(t_char))//' K'
                 FILEunit=1200
   
                 CALL DrawMolecularConfiguration( ZNucleus,radiusSBC,NMolecules,NAtoms,  &
                                                  Max_Atoms,MiddlePoint,NTrajectories,rf,&
                                                  t,FILEunit,comment_line                   )
   
                 CLOSE(1200)
   
              ENDIF
              !$OMP END CRITICAL   
   
           ELSE !if the NEW config is not accepted
   
              !if SM > N we rotate the whole cluster, otherwise a random atom is moved
              ! 1/(NMolecules + 1)  is the probability of attempting a rotation of the whole cluster
              IF ( SelectedMolecule > NMolecules ) THEN
   
                 !recovering the old configuration 
                 DO k=1,3
                    DO i=1,NMolecules
                       r(t,i,1,k)=rf(t,i,1,k)
                    ENDDO
                 ENDDO
   
                 DO Molecule_b=1,NMolecules
                    !E(a,a)=0
                    Energy_ab_fin(t,Molecule_b,Molecule_b)= 0.0
   
                    !recovering the pair energy contribution
                    Energy_2a_fin(t,Molecule_b) = Energy_2a_ini(t,Molecule_b)
   
                    DO Molecule_a=Molecule_b+1,NMolecules
                       !recovering the pair energy E(a,b)
                       Energy_ab_fin(t,Molecule_a,Molecule_b)= Energy_ab_ini(t,Molecule_a,Molecule_b)
                       !E(a,b)=E(b,a)
                       Energy_ab_fin(t,Molecule_b,Molecule_a)= Energy_ab_ini(t,Molecule_a,Molecule_b)
                    ENDDO
                 ENDDO
   
              ELSE  !if the whole structure is not rotated, we move one atom randomly
   
                 !recovering the old configuration 
                 DO k=1,3
                    r(t,SelectedMolecule,1,k)=rf(t,SelectedMolecule,1,k)
                 ENDDO
                 
                 !recovering the pair energy contribution
                 Energy_2a_fin(t,SelectedMolecule) = Energy_2a_ini(t,SelectedMolecule)
                 
                 DO Molecule_a=1,NMolecules
                    !recovering the pair energy E(a,b)
                    Energy_ab_fin(t,Molecule_a,SelectedMolecule)= Energy_ab_ini(t,Molecule_a,SelectedMolecule)
                    !E(a,b)=E(b,a)
                    Energy_ab_fin(t,SelectedMolecule,Molecule_a)= Energy_ab_ini(t,Molecule_a,SelectedMolecule)
                 ENDDO
   
              ENDIF
           ENDIF
   
        ENDDO !extra loop over the system size
   
        !-adjust step size 'maxDisplacement' after 100 steps
        IF ( MOD(MCstep,100) == 0 ) THEN                        !100 is the number of steps after which step size is updated
           Adjust=REAL(adj_Mov_accepted(t))/REAL(100*NMolecules)          !Denominator is multiplied by NMolecules because the EXTRA loop over the number of molecules
           IF ( Adjust < 0.4 ) THEN                             !number of accepted steps/number of total steps should   
              maxDisplacement(t)=maxDisplacement(t)*0.9         !lie in [0.4,0.6]
           ELSEIF ( Adjust > 0.6 ) THEN
              maxDisplacement(t)=maxDisplacement(t)*1.1
           ENDIF
           !maxDisplacement must be less than Boundary Sphere diameter
           IF ( maxDisplacement(t) >= 2.0*radiusSBC ) THEN
              maxDisplacement(t) = 1.9*radiusSBC
           ENDIF
           adj_Mov_accepted(t)=0                                !It restarts to count only the conf accepted every 100 steps
   
           !-adjust step size 'MaxRotation' after 100 steps
           Adjust=REAL(adj_Rot_accepted(t))/100.0               !The probability of any rotation is 1/(n+1); i.e., NMolecule loop does not matter 
           IF ( Adjust < 0.4 ) THEN                             !number of accepted steps/number of total steps should   
              MaxRotation(t)=MaxRotation(t)*0.9                 !lie in [0.4,0.6]
           ELSEIF ( Adjust > 0.6 ) THEN
              MaxRotation(t)=MaxRotation(t)*1.1
           ENDIF
           !maxRotation must be less than 90 degrees 
           IF ( MaxRotation(t) >= 0.5*PI ) THEN
              MaxRotation(t) = 0.5*PI                  !Any rotation larger than 90 deg makes no sense in this contest
           ENDIF
           adj_Rot_accepted(t)=0                       !It restarts to count only the conf accepted every 100 steps
        ENDIF
   
        ENDDO !loop over temperatures/trajectories
        !$OMP END PARALLEL DO

   
        !----------------------------------------------------------------------------------------
        !Replica exchange MC sampling (Parallel tempering). The exchange is between ExchangeTrajectory and ExchangeTrajectory+1
      
        IF ( 0.1 > RandomNumber(Seed) ) THEN                                                ! 0.1 is the probability of attempting an exchange
           ExchangeTrajectory=1+FLOOR( (NTrajectories-1)*RandomNumber(Seed) )               ! index, random trajectory to exchange in range [1:NTrajectories-1] 
   
           !-Total number of attempted exchanges
           Exchange_tot(ExchangeTrajectory)  =Exchange_tot(ExchangeTrajectory)  +1
           Exchange_tot(ExchangeTrajectory+1)=Exchange_tot(ExchangeTrajectory+1)+1
   
           DeltaBeta=beta(ExchangeTrajectory) - beta(ExchangeTrajectory+1)
           DeltaEnergy=Energy_Tot(ExchangeTrajectory) - Energy_Tot(ExchangeTrajectory+1)
           ExchangeAccept=MIN( 1.E0,DEXP(DeltaBeta*DeltaEnergy) )
   
           ! exchange accepted
           IF ( ExchangeAccept > RandomNumber(Seed) ) THEN
              !-Number of accepted exchanges
              Exchange_acc(ExchangeTrajectory)  =Exchange_acc(ExchangeTrajectory)  +1
              Exchange_acc(ExchangeTrajectory+1)=Exchange_acc(ExchangeTrajectory+1)+1
              ! swap configurations for an accepted configuration (rf)
              DO i=1,NMolecules
                 DO j=1,NAtoms(i)
                    DO k=1,3
                       !for r
                       rex(ExchangeTrajectory,i,j,k)=r(ExchangeTrajectory,i,j,k)          !rex is a temporary variable, so rex=rf(j)
                       r(ExchangeTrajectory,i,j,k)  =r(ExchangeTrajectory+1,i,j,k)        !r(j)  =r(j+1)
                       r(ExchangeTrajectory+1,i,j,k)=rex(ExchangeTrajectory,i,j,k)        !r(j+1)=rex
                       !for rf
                       rex(ExchangeTrajectory,i,j,k) =rf(ExchangeTrajectory,i,j,k)        !rex is a temporary variable, so rex=rf(j)
                       rf(ExchangeTrajectory,i,j,k)  =rf(ExchangeTrajectory+1,i,j,k)      !rf(j)  =rf(j+1)
                       rf(ExchangeTrajectory+1,i,j,k)=rex(ExchangeTrajectory,i,j,k)       !rf(j+1)=rex
                    ENDDO
                 ENDDO
              ENDDO
              ! swap Energies
              Energy_ex                       =Energy_Tot(ExchangeTrajectory)
              Energy_Tot(ExchangeTrajectory)  =Energy_Tot(ExchangeTrajectory+1)
              Energy_Tot(ExchangeTrajectory+1)=Energy_ex
   
              DO Molecule_a=1,NMolecules
               Energy_2a_ini_ex                               = Energy_2a_ini(ExchangeTrajectory,Molecule_a)
               Energy_2a_ini(ExchangeTrajectory,Molecule_a)   = Energy_2a_ini(ExchangeTrajectory+1,Molecule_a)
               Energy_2a_ini(ExchangeTrajectory+1,Molecule_a) = Energy_2a_ini_ex
               DO Molecule_b=Molecule_a+1,NMolecules
                Energy_ab_ini_ex                                         = Energy_ab_ini(ExchangeTrajectory,Molecule_a,Molecule_b)
                Energy_ab_ini(ExchangeTrajectory,Molecule_a,Molecule_b)  = Energy_ab_ini(ExchangeTrajectory+1,Molecule_a,Molecule_b)
                Energy_ab_ini(ExchangeTrajectory+1,Molecule_a,Molecule_b)= Energy_ab_ini_ex
   
                !E(a,b)=E(b,a) for t and t+1
                Energy_ab_ini(ExchangeTrajectory,Molecule_b,Molecule_a)  = Energy_ab_ini(ExchangeTrajectory,Molecule_a,Molecule_b)
                Energy_ab_ini(ExchangeTrajectory+1,Molecule_b,Molecule_a)= Energy_ab_ini(ExchangeTrajectory+1,Molecule_a,Molecule_b)
               ENDDO
              ENDDO
   
           ENDIF
        ENDIF
   
        !----------------------------------------------------------------------------------------
        !This is the place to collect data (Energies) on the fly
        DO t=1,NTrajectories        
   
           TotalEnergy(t)=TotalEnergy(t)+Energy_Tot(t)
           TotalEnergy2(t)=TotalEnergy2(t)+Energy_Tot(t)**2
   
           !Drawing the conf accepted every SaveConfig MC cycles
           IF ( SaveConfig > 0 ) THEN                       !Do you wanna save the configs? Yes: SaveConfig>0 or No: SaveConfig<0
              IF ( accept(t) ) THEN                            !if you wanna save them, they'll only save when they are accepted for the Metropolis criterion
                 IF ( MOD(MCstep,SaveConfig) == 0 ) THEN    !No all of them will be saved, SaveConfig defines a cycle to save configs
   
                    WRITE(energy_tot_char,"(ES18.10E2)") Energy_Tot(t)*JtokJmol
                    comment_line='Energy (kJ/mol): '//TRIM(ADJUSTL(energy_tot_char))
                    FILEunit=500+t                                            !independent FILEs for each temperature
   
                    CALL DrawMolecularConfiguration( ZNucleus,radiusSBC,NMolecules,NAtoms,  &     
                                                     Max_Atoms,MiddlePoint,NTrajectories,rf,&
                                                     t,FILEunit,comment_line                 )   
                 ENDIF
              ENDIF
           ENDIF
   
           !-Histogram Analysis. Counting hits per Energy bin
           hit=NINT( (Energy_Tot(t) - MinEnergyTotal) / DeltaEHistogram )                !'hit' is the index for the frequency
           IF ( hit >= 0 .AND. hit <= nBins ) THEN
              nhits(t,hit)=nhits(t,hit)+1                                               !frequency counter for the energy in the bin
           ENDIF
           !Counting the energies OUTSIDE the range (Emax-Emin)
           IF ( hit < 0     ) nConfigELT=nConfigELT+1
           IF ( hit > nBins ) nConfigEGT=nConfigEGT+1
   
           !RDF histogram
           DO Molecule_a=1,NMolecules
              DO Molecule_b=Molecule_a+1,NMolecules
                 dist=SQRT(                                                    &   !r_ab=SQRT[ 
                            (rf(t,Molecule_a,1,1) - rf(t,Molecule_b,1,1))**2+&   !           (xa-xb)^2+
                            (rf(t,Molecule_a,1,2) - rf(t,Molecule_b,1,2))**2+&   !           (ya-yb)^2+
                            (rf(t,Molecule_a,1,3) - rf(t,Molecule_b,1,3))**2  )  !           (za-zb)^2 ]
                 indexRDF=NINT( (dist - MinR) / DeltaR )              
                 RDFHistogram(t,indexRDF)=RDFHistogram(t,indexRDF) + 1
   
                 !x-axis
                 distx= ABS(rf(t,Molecule_a,1,1) - rf(t,Molecule_b,1,1))
                 indexRDF=NINT( (distx - MinR) / DeltaR )               
                 RDFHistogramx(t,indexRDF)=RDFHistogramx(t,indexRDF) + 1
   
                 !y-axis
                 disty= ABS(rf(t,Molecule_a,1,2) - rf(t,Molecule_b,1,2))
                 indexRDF=NINT( (disty - MinR) / DeltaR )               
                 RDFHistogramy(t,indexRDF)=RDFHistogramy(t,indexRDF) + 1
   
                 !z-axis
                 distz= ABS(rf(t,Molecule_a,1,3) - rf(t,Molecule_b,1,3))
                 indexRDF=NINT( (distz - MinR) / DeltaR )               
                 RDFHistogramz(t,indexRDF)=RDFHistogramz(t,indexRDF) + 1
              ENDDO
           ENDDO
           
        ENDDO
   
     ENDDO !loop over MC steps
   
     !--------------------------------------------------------------------------------------------
     !OUTPUT
     
     TotalAccepted_Mov =0
     TotalAccepted2_Mov=0     !after eq
     TotalAccepted_Rot =0
     TotalAccepted2_Rot=0     !after eq
     TotalExchanged=0
     DO t=1,NTrajectories        
        TotalAccepted_Mov =TotalAccepted_Mov +Mov_accepted_eq(t)
        TotalAccepted2_Mov=TotalAccepted2_Mov+Mov_accepted(t)    !after eq
        TotalAccepted_Rot =TotalAccepted_Rot +Rot_accepted_eq(t)
        TotalAccepted2_Rot=TotalAccepted2_Rot+Rot_accepted(t)    !after eq
        TotalExchanged=TotalExchanged+Exchange_acc(t)
     ENDDO
     
     WRITE(10000+check2,'(A)')''
     WRITE(10000+check2,'(A,/)')' The lowest and highest Energies are: '
     DO t=1,NTrajectories
        WRITE(10000+check2,'(I3,3X,F5.1,T20,2(ES15.6E2,5X))')t,Temperature(t),&
                            MinEnergy(t)*JtokJmol,MaxEnergy(t)*JtokJmol
     ENDDO
     WRITE(10000+check2,'(/,A,T20,2(ES15.6E2,5X))')' Global Min and Max Energy: ',&
                           MinEnergyTotal*JtokJmol,MaxEnergyTotal*JtokJmol
     WRITE(10000+check2,'(A)')''
     
     WRITE(10000+check2,'(A,T40,I20)')' Number of steps: ',MCcycles
     WRITE(10000+check2,'(A,T40,I20)')' Equilibration: ',Nequilibration
     WRITE(10000+check2,'(A,T40,I20)')' Total Movement / Trajectory: ',NMolecules*check_MCcycles
     WRITE(10000+check2,'(A,T40,I20)')' Total Accepted: ',&
                               TotalAccepted_Mov+TotalAccepted_Rot+TotalAccepted2_Mov+TotalAccepted2_Rot
     WRITE(10000+check2,'(A,T40,I20)')' Total Accepted (after eq): ',TotalAccepted2_Mov+TotalAccepted2_Rot
     WRITE(10000+check2,'(A,T40,I20)')' Total Exchanged: ',TotalExchanged 
     WRITE(10000+check2,'(A)')''
     WRITE(10000+check2,'(A)')" 1/(n+1)  is the probability of attempting a rotation"&
                            //" of the whole cluster, where 'n' is number of atoms"
     WRITE(10000+check2,'(A,F5.1,A)')' The DISPLACEMENT inside a boundary sphere is within 0 and ',&
                                   2.0*radiusSBC,' Angtrom'
     WRITE(10000+check2,'(A)')' The angle of ROTATION of the whole cluster is within 0 and 90 degree'
     WRITE(10000+check2,'(A)')''
     
     
     WRITE(10000+check2,'(A,/)') &
     ' #    T (K)  Max Disp. (Ang)   Accept      Ratio      '//&
                  'Max Rot. (deg)    Accept      Ratio    Tot accep rat   Exchang    Ex. ratio'
     DO t=1,NTrajectories
        WRITE(10000+check2,"(T2,I2,T7,F5.1,T17,F9.5,T29,I9,T43,F7.4,T58,F10.4,T70,I9,T84,F7.4,T96,F7.4,T109,I8,T122,F7.4)")&
        t,Temperature(t),maxDisplacement(t),&
        Mov_accepted(t),REAL(Mov_accepted(t))/REAL(NMolecules*check_MCcycles),&
        MaxRotation(t)*(180.0/PI),&
        Rot_accepted(t),REAL(Rot_accepted(t))/REAL(check_MCcycles),&
        REAL(Mov_accepted(t)+Rot_accepted(t))/REAL(NMolecules*check_MCcycles),&
        Exchange_acc(t),REAL(Exchange_acc(t))/REAL(Exchange_tot(t))
     ENDDO
   
     !---------------------------------------------------------------------------------------------------
     !-Computing Average Energies Heat Capacities
   
     WRITE(10000+check2,'(A)')''
     WRITE(10000+check2,'(A,/)')'  #   T (K)   < E(T) > (kJ/mol)   Cv(T) (kJ/mol/K)'!   dCv/dT (kJ/mol/K/K)'
   
     DO t=1,NTrajectories
        AverageEnergy(t) =TotalEnergy(t) /check_MCcycles                                  !Average Energy per Temperature, < E(t) >   =E_tot/MCcycles
        AverageEnergy2(t)=TotalEnergy2(t)/check_MCcycles                                  !Square of the Average Energy,   < E(t)**2 >=E_tot**2/MCcycles
        Cv(t)=( AverageEnergy2(t) - AverageEnergy(t)**2 )/(kB*Temperature(t)**2)    !The Heat Capacity is Cv(t)=( <E(t)**2> - <E(t)>**2 )/(kB*T**2)
   
        !calculating dCv/dT, from DOI: 10.1021/ct200458m
   !     REWIND(9000+t)                                                  !Go back to initial point in the scratch FILE
   !     DO n=1,check_MCcycles
   !        READ(9000+t,'(ES30.20E3)') Energy
   !        Energy_Average(t) = ( Energy - AverageEnergy(t) )
   !        Energy_Average2(t)= ( Energy - AverageEnergy(t) )**2/check_MCcycles
   !        Energy_Average3(t)= ( Energy - AverageEnergy(t) )**3/check_MCcycles
   !        dCv(t)=  Energy_Average3(t)/( (kB**2)*(Temperature(t)**4) ) -&
   !               2*Energy_Average2(t)/( kB*Temperature(t)**3 )
   !     ENDDO
   
        WRITE(10000+check2,'(X,I2,3X,F5.1,T14,3(ES15.6E2,4X))')&
                  t,Temperature(t),AverageEnergy(t)*JtokJmol,Cv(t)*JtokJmol!,dCv(t)*JtokJmol
   
     ENDDO
   
     !---------------------------------------------------------------------------------------------------
     !-Results from the Multihistogram Analysis
   
     WRITE(10000+check2,'(A)')''
     WRITE(10000+check2,'(A,/,T5,A1,2(ES15.6E2,A))')' Energy Histrogram Analysis within',&
                                          '[',MinEnergyTotal*JtokJmol,',',MaxEnergyTotal*JtokJmol,' ] kJ/mol'
     WRITE(10000+check2,'(/,A)')' Number of Configuration with Energy'
     WRITE(10000+check2,'(T5,A,ES15.6E2,A9,I10,A)')'Below: ',MinEnergyTotal*JtokJmol,' kJ/mol, ',nConfigELT,' Config'
     WRITE(10000+check2,'(T5,A,ES15.6E2,A9,I10,A)')'Above: ',MaxEnergyTotal*JtokJmol,' kJ/mol, ',nConfigEGT,' Config'
     WRITE(10000+check2,'(/,A,I5)')' Energy Histogram Resolution: ',nBins
   
   
     OPEN(2000,FILE='histE.data',STATUS='REPLACE')   !summary results 
   
     WRITE(2000,*) 
     WRITE(2000,*) kB*JtokJmol
     WRITE(2000,'(I2)') NTrajectories
     WRITE(2000,*) (Temperature(t),t=1,NTrajectories)
     WRITE(2000,*) MinEnergyTotal*JtokJmol,MaxEnergyTotal*JtokJmol,nBins
   
     DO t=1,NTrajectories
   
        WRITE(2000,'(I2,2(5X,ES15.6E3),5X,F10.6)') t-1,MinEnergy(t)*JtokJmol,MaxEnergy(t)*JtokJmol,maxDisplacement(t)
   
        WRITE(t_char,"(I4)") t                                 !changing the variable for a character to concatenate it
        OPEN(2100+t,FILE='histE.'//TRIM(ADJUSTL(t_char)))         !detailed results for each temperature
        DO hit=0,nBins
           WRITE(2100+t,'(ES20.10E3,5X,I20)') MinEnergyTotal*JtokJmol+(hit)*DeltaEHistogram*JtokJmol,nhits(t,hit)
        ENDDO
        CLOSE(2100+t)
     ENDDO
     CLOSE(2000)
   
     !--------------------------------------------------------------------------------------------------
     !Radial Distribution Analysis
   
     DO t=1,NTrajectories
        WRITE(t_char,"(I4)") t
        OPEN(2300+t,FILE='RDFhist.'//TRIM(ADJUSTL(t_char)) )
        OPEN(2400+t,FILE='RDFhist_x.'//TRIM(ADJUSTL(t_char)) )
        OPEN(2600+t,FILE='RDFhist_y.'//TRIM(ADJUSTL(t_char)) )
        OPEN(2800+t,FILE='RDFhist_z.'//TRIM(ADJUSTL(t_char)) )
        DO indexRDF=0,BinsRDF
           WRITE(2300+t,'(ES20.10E3,5X,I20)') MinR+(indexRDF)*DeltaR,RDFHistogram(t,indexRDF)
           WRITE(2400+t,'(ES20.10E3,5X,I20)') MinR+(indexRDF)*DeltaR,RDFHistogramx(t,indexRDF)
           WRITE(2600+t,'(ES20.10E3,5X,I20)') MinR+(indexRDF)*DeltaR,RDFHistogramy(t,indexRDF)
           WRITE(2800+t,'(ES20.10E3,5X,I20)') MinR+(indexRDF)*DeltaR,RDFHistogramz(t,indexRDF)
        ENDDO
        CLOSE(2300+t)
        CLOSE(2400+t)
        CLOSE(2600+t)
        CLOSE(2800+t)
     ENDDO
   
     !---------------------------------------------------------------------------------------------------
     !-OUTPUT: Ending
     WRITE(10000+check2,'(A)')''
     WRITE(10000+check2,'(A)')''
     WRITE(10000+check2,'(A)')' NORMAL TERMINATION'
     WRITE(10000+check2,'(A)')'****************************************************************'
     CALL DATE_AND_TIME(today,now)
     WRITE(10000+check2,"(12A)")' Date and time: ',&
                     today(1:4),"-",today(5:6),"-",today(7:8), ", ",&
                     now(1:2),":",  now(3:4),":",  now(5:6)
     WRITE(10000+check2,'(A)')''
     CALL CPU_TIME(EndTime)
     RunTime=EndTime-StartTime
     
     !$ wtime_end = omp_get_wtime()
     !$ RunTime = wtime_end - wtime_start
     
     days   =FLOOR( RunTime/86400.0 )
     hours  =FLOOR( (RunTime - REAL(days*86400)) /3600.0 )
     minutes=FLOOR( (RunTime - REAL(days*86400+hours*3600)) /60.0 )
     seconds=RunTime - REAL( days*86400+hours*3600+minutes*60 )
     WRITE(10000+check2,'(A16,3(I4,X,A),X,F7.4,X,A)')' CPU time used: ',&
           days,'days',hours,'hours',minutes,'minutes',seconds,'seconds'
     WRITE(10000+check2,'(A)')'****************************************************************'
     
     FLUSH(10000+check2)  !OUTPUT FILE
     CLOSE(UNIT=10000+check2)  !OUTPUT FILE


     !----------------------------------------------------------------------------------------
     !- Saving FILEs to restart
     OPEN(UNIT=20000+check2,FILE='restart_'//output//'_'//TRIM(ADJUSTL(check_char)),STATUS='REPLACE')

     WRITE(20000+check2,"(A)") '-----------------------------------------------------------------'
     WRITE(20000+check2,"(12A)")'Date and time: ',&
                                today(1:4),"-",today(5:6),"-",today(7:8), ", ",&
                                now(1:2),":",  now(3:4),":",  now(5:6)
     WRITE(20000+check2,"(A)") '-----------------------------------------------------------------'
     WRITE(20000+check2,"(A)") 'RESTART SETTING'
     WRITE(20000+check2,"(A)") '-----------------------------------------------------------------'
     WRITE(20000+check2,"(A)")   
     WRITE(20000+check2,"(A,T30,A20)") "MC cycles/sampling: ",TRIM(ADJUSTL(MCcycles_char))
     WRITE(20000+check2,"(A,T30,I20)") "Checkpoints for sampling: ",Ncheckpoints_mc
     WRITE(20000+check2,"(A,T30,A20)") "Equilibracion: ",TRIM(ADJUSTL(Nequilibration_char))
     WRITE(20000+check2,"(A,T30,I20)") "Checkpoints for Equilibracion: ",Ncheckpoints_eq
     WRITE(20000+check2,"(A,T30,I20)") "Current point: ",MCstep-1
     WRITE(20000+check2,'(A,T30,A20)') "Equilibration Evolution (%): ","100.0"
     WRITE(20000+check2,'(A,T30,F20.1)') "Sampling Evolution (%): ",(REAL(check)/REAL(Ncheckpoints_mc))*100.0

     WRITE(20000+check2,"(A,T30,I20)") "Seed: ",Seed
     WRITE(20000+check2,"(A)")   

     !- Initia States
     WRITE(20000+check2,"(A)") '-----------------------------------------------------------------'
     WRITE(20000+check2,"(A)") ' Initia States '
     WRITE(20000+check2,"(A)") '-----------------------------------------------------------------'
     WRITE(20000+check2,"(A)") '*** Energies in Joules '   
     WRITE(20000+check2,"(A)")   
     WRITE(20000+check2,"(A,T40,ES28.20E2)") 'Total Energy: ',Energy_Tot(0) 
     WRITE(20000+check2,"(A)")   
     DO Molecule_a=1,NMolecules
        WRITE(20000+check2,"(A32,I4,A4,T40,ES28.20E2)") 'Pair Energy for atom: ',&
                                                       Molecule_a,' is ',Energy_2a_ini(0,Molecule_a)
        DO Molecule_b=Molecule_a+1,NMolecules
           WRITE(20000+check2,"(A25,2(I4,A5),T40,ES28.20E2)") '2b Energy for ',Molecule_a,' and ',&
                                                              Molecule_b,' is ',Energy_ab_ini(0,Molecule_a,Molecule_b)
        ENDDO
     WRITE(20000+check2,"(A)")   
     ENDDO

     !- Trajectories
     DO t=1,NTrajectories
        WRITE(20000+check2,"(A)") '-----------------------------------------------------------------'
        WRITE(20000+check2,"(A,T20,I3,A4,T40,F10.5)") ' Temperature (K): ',t,' is ',Temperature(t)
        WRITE(20000+check2,"(A)") '-----------------------------------------------------------------'
        WRITE(20000+check2,"(A,T40,F30.25)") ' Max Disp. (Ang): ',maxDisplacement(t)
        WRITE(20000+check2,"(A,T40,F30.25)") ' Max Rot. (rad): ',MaxRotation(t)!*(180.0/PI)
        WRITE(20000+check2,"(A)")   
        !Structural parameters
        DO i=1,NMolecules
           WRITE(20000+check2,FMT3) ZNucleus(i,1),rf(t,i,1,1),rf(t,i,1,2),rf(t,i,1,3)
        ENDDO           
        WRITE(20000+check2,"(A)")   
        !Acceptance criterion, adjusting step size
        WRITE(20000+check2,"(A,T40,I20)") 'Total Movements Accepted: ',Mov_accepted_eq(t)
        WRITE(20000+check2,"(A,T40,I20)") 'Movements Accepted (per 100 Cycles): ',adj_Mov_accepted_eq(t)
        WRITE(20000+check2,"(A,T40,I20)") 'Total Rotations Accepted: ',Rot_accepted_eq(t)
        WRITE(20000+check2,"(A,T40,I20)") 'Rotations Accepted (per 100 Cycles): ',adj_Rot_accepted_eq(t)
        WRITE(20000+check2,"(A,T40,I20)") 'Attempts for Exchange: ',Exchange_tot(t)
        WRITE(20000+check2,"(A,T40,I20)") 'Number of Exchange Accepted: ',Exchange_acc(t)
        WRITE(20000+check2,"(A)")   
        !Energy
        WRITE(20000+check2,"(A)") '-----------------------'
        WRITE(20000+check2,"(A)") '*** Energies in Joules '   
        WRITE(20000+check2,"(A)")   
        WRITE(20000+check2,"(A,T40,ES28.20E2)") 'Total Energy: ',Energy_Tot(t) 
        WRITE(20000+check2,"(A,T40,ES28.20E2)") 'Min Energy: ',MinEnergy(t) 
        WRITE(20000+check2,"(A,T40,ES28.20E2)") 'Max Energy: ',MaxEnergy(t) 
        WRITE(20000+check2,"(A)")   
        DO Molecule_a=1,NMolecules
           WRITE(20000+check2,"(A,I4,A4,T40,ES28.20E2)") 'Pair Energy for Atom ',Molecule_a,' is ',Energy_2a_ini(t,Molecule_a)
           DO Molecule_b=Molecule_a+1,NMolecules
              WRITE(20000+check2,"(A23,I4,A5,I4,A4,T40,ES28.20E2)") '2b Energy for ',Molecule_a,' and ',&
                                                                    Molecule_b,' is ',Energy_ab_ini(t,Molecule_a,Molecule_b)
           ENDDO
           WRITE(20000+check2,"(A)")   
        ENDDO
        !during Samplig (after eq)
        WRITE(20000+check2,"(A)") '-----------------------'
        WRITE(20000+check2,"(A)") 'After Equilibration: '
        WRITE(20000+check2,"(A,T40,I20)") 'Total Movements Accepted: ',Mov_accepted(t)
        WRITE(20000+check2,"(A,T40,I20)") 'Movements Accepted (per 100 Cycles): ',adj_Mov_accepted(t)
        WRITE(20000+check2,"(A,T40,I20)") 'Total Rotations Accepted: ',Rot_accepted(t)
        WRITE(20000+check2,"(A,T40,I20)") 'Rotations Accepted (per 100 Cycles): ',adj_Rot_accepted(t)
        WRITE(20000+check2,"(A,T40,ES28.20E2)") 'Sum of E_tot per temperature: ',TotalEnergy(t)
        WRITE(20000+check2,"(A,T40,ES28.20E2)") 'Sum of E_tot**2 per temperature: ',TotalEnergy2(t)

        WRITE(20000+check2,"(A)")   
        WRITE(20000+check2,"(A)") 'Multihistogram Analysis (hits)'
        WRITE(20000+check2,"(I20)") ( nhits(t,hit), hit=0,nBins )
        WRITE(20000+check2,"(A)")   
        WRITE(20000+check2,"(A)") 'Radial Distribution Analysis:'
        WRITE(20000+check2,"(A)") 'x-axis: '
        WRITE(20000+check2,"(I20)") (RDFHistogramx(t,indexRDF), indexRDF=0,BinsRDF )
        WRITE(20000+check2,"(A)") 'y-axis: '
        WRITE(20000+check2,"(I20)") (RDFHistogramy(t,indexRDF), indexRDF=0,BinsRDF )
        WRITE(20000+check2,"(A)") 'z-axis: '
        WRITE(20000+check2,"(I20)") (RDFHistogramz(t,indexRDF), indexRDF=0,BinsRDF )
        WRITE(20000+check2,"(A)") 'Total: '
        WRITE(20000+check2,"(I20)") (RDFHistogram(t,indexRDF), indexRDF=0,BinsRDF )

        WRITE(20000+check2,"(A)")   
     ENDDO


     WRITE(20000+check2,"(A,T40,ES18.10E2)") 'Total Min Energy after Equilibration: ',MinEnergyTotal
     WRITE(20000+check2,"(A,T40,ES18.10E2)") 'Total Max Energy after Equilibration: ',MaxEnergyTotal
     WRITE(20000+check2,"(A,T40,ES18.10E2)") 'Bin Width for the Histogram: ',DeltaEHistogram
     WRITE(20000+check2,"(A,T30,I20)") 'Configs out (upper): ',nConfigEGT
     WRITE(20000+check2,"(A,T30,I20)") 'Configs out (lower): ',nConfigELT


     WRITE(20000+check2,"(A)")   
     WRITE(20000+check2,"(A)") '-----------------------------------------------------------------'
     WRITE(20000+check2,"(A)") '*** END OF FILE ***'  
     WRITE(20000+check2,"(A)") '-----------------------------------------------------------------'
     CLOSE(20000+check2)



     !----------------------------------------------------------------------------------------
     !-Sorting out FILEs
     CALL SYSTEM("mkdir -p histograms ; mv -fu histE* histograms/")
     CALL SYSTEM("mkdir -p radial_distribution ; mv -fu RDFhist* radial_distribution/")
     CALL SYSTEM("mkdir -p configurations")
     IF ( check == 1 ) CALL SYSTEM("mkdir -p checkpoint/")     !this dir is only created in the beginning

     !----------------------------------------------------------------------------------------
     !-STATUS, this FILE has the MC evolution
     WRITE(1000,'(A)')
     CALL CPU_TIME(PartialEndRunTime)
     RunTime=PartialEndRunTime-PartialStartRunTime   !Delta time
     PartialStartRunTime=PartialEndRunTime           !redefining starting time here to start the counting again

     !$ part_wtime_end = omp_get_wtime()
     !$ RunTime = part_wtime_end - part_wtime_start
     !$ part_wtime_start = part_wtime_end 

     CALL DATE_AND_TIME(today,now)
     WRITE(1000,"(12A)")'Date and time: ',&
                        today(1:4),"-",today(5:6),"-",today(7:8), ", ",&
                        now(1:2),":",  now(3:4),":",  now(5:6)

     WRITE(MCstep_char,"(I10)")     MCstep-1

     WRITE(1000,'(A,T19,A)')'Sample steps: ',TRIM(ADJUSTL(MCstep_char))//' of '//TRIM(ADJUSTL(MCcycles_char))
     WRITE(1000,'(A,T18,F5.1,A)')'MC Evolution: ',(REAL(check)/REAL(Ncheckpoints_mc))*100.0,' %'
     WRITE(1000,'(A,T18,ES12.5E2,A)') 'CPU time used: ',RunTime,' sec'

     IF ( check < Ncheckpoints_mc ) THEN
        CALL SYSTEM("mkdir -p checkpoint/"//TRIM(ADJUSTL(check_char)))

        INQUIRE (FILE="initial_configuration.xyz",EXIST=filexist)
        IF ( FILExist ) CALL SYSTEM("cp -f initial_configuration.xyz configurations/")

        INQUIRE (FILE="lowest_Energy.xyz",EXIST=filexist)
        IF ( FILExist ) CALL SYSTEM("cp -f lowest_Energy.xyz configurations/")

        IF ( SaveConfig > 0 ) THEN          !Save configs? Yes: SaveConfig>0 or No: SaveConfig<=0, SaveConfig defines a cycle to save the configurations         
           DO t=1,NTrajectories
              WRITE(t_char,"(I10)") t
              INQUIRE (FILE="configuration_"//TRIM(ADJUSTL(t_char))//".xyz",EXIST=filexist)
              IF ( FILExist ) CALL SYSTEM("cp -f configuration_"//TRIM(ADJUSTL(t_char))//".xyz configurations/")
           ENDDO
        ENDIF

        CALL SYSTEM("cp -f "//input//" checkpoint/"//TRIM(ADJUSTL(check_char)))
        CALL SYSTEM("mv -fu "//output//"_"//TRIM(ADJUSTL(check_char))//"ckp.out checkpoint/"//TRIM(ADJUSTL(check_char)))
        CALL SYSTEM("mv -fu histograms/ radial_distribution/ configurations/ checkpoint/"//TRIM(ADJUSTL(check_char)))

        WRITE(1000,'(A)')'Checkpoint FILEs:   checkpoint/'//TRIM(ADJUSTL(check_char))//'/'

!        CALL SYSTEM("mv -fu restart_"//output//"_"//TRIM(ADJUSTL(check_char))//&
!                    " checkpoint/"//TRIM(ADJUSTL(check_char)) )

     ELSE !in the end

        CALL SYSTEM("mkdir -p checkpoint/")  !when we have only 1 checkpoint, it is necessary

        INQUIRE (FILE="lowest_Energy.xyz",EXIST=filexist)
        IF ( FILExist ) CALL SYSTEM("mv -fu lowest_Energy.xyz configurations/")

        INQUIRE (FILE="initial_configuration.xyz",EXIST=filexist)
        IF ( FILExist ) CALL SYSTEM("mv -fu initial_configuration.xyz configurations/")

        IF ( SaveConfig > 0 ) THEN          !Save configs? Yes: SaveConfig>0 or No: SaveConfig<=0, SaveConfig defines a cycle to save the configurations         
           DO t=1,NTrajectories
              WRITE(t_char,"(I10)") t
              INQUIRE (FILE="configuration_"//TRIM(ADJUSTL(t_char))//".xyz",EXIST=filexist)
              IF ( FILExist ) CALL SYSTEM("mv -fu configuration_"//TRIM(ADJUSTL(t_char))//".xyz configurations/")
           ENDDO
        ENDIF

        WRITE(1000,'(A)')
        WRITE(1000,'(A)')'*** END OF MONTE CARLO CYCLE ***'

!        CALL SYSTEM("mv -fu restart_"//output//"_"//TRIM(ADJUSTL(check_char))//" checkpoint/")

     ENDIF

     CALL SYSTEM("mv -fu restart_"//output//"_"//TRIM(ADJUSTL(check_char))//" checkpoint/")

     !- Deleting old checkpoint restart-files
     WRITE(check_char,"(I10)") check-1
     CALL SYSTEM("rm -f checkpoint/restart_"//output//"_"//TRIM(ADJUSTL(check_char)) )
  
     WRITE(1000,'(A)')
     FLUSH(1000)

  ENDDO   !checkpoint loop
!---------------------------------------------------------------------------------------------------

  CLOSE(UNIT=1000) !Status FILE

  INQUIRE (FILE="./checkpoint/",EXIST=direxist)
  IF ( direxist ) THEN
     CALL SYSTEM("tar -czf _saveData.tar.gz checkpoint/")
     CALL SYSTEM("rm -rf checkpoint/")
  ENDIF

  !-Closing all the FILEs
  IF ( SaveConfig > 0 ) THEN          !Save configs? Yes: SaveConfig>0 or No: SaveConfig<=0, SaveConfig defines a cycle to save the configurations  
     DO t=1,NTrajectories
        CLOSE(UNIT=500+t) !Configuration XYZ FILE
     ENDDO
  ENDIF

!STOP  !comment this to avoid warnings such as: "IEEE_UNDERFLOW_FLAG"  and  "IEEE_DENORMAL"
END PROGRAM atomic_PTMC
