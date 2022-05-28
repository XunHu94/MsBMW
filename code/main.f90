    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    Module ran_mod
      Implicit None
    ! ran return a uniform random number between 0-1  
    ! norma return a normal distribution  
    contains
      function ran()   !returns random number between 0 - 1  
        implicit none
        integer , save :: flag = 0
        real :: ran        
        if(flag==0) then
          call random_seed()
          flag = 1
        endif
        call random_number(ran)
      end function ran

      function normal(mean,sigma)
        implicit none
        integer :: flag 
        real, parameter :: pi = 3.141592653589793239
        real :: u1, u2, y1, y2, normal, mean, sigma
        save flag 
        data flag /0/
        u1 = ran(); u2 = ran()
        if (flag.eq.0) then
          y1 = sqrt(-2.0d0*log(u1))*cos(2.0d0*pi*u2)
          normal = mean + sigma*y1
          flag = 1
        else
          y2 = sqrt(-2.0d0*log(u1))*sin(2.0d0*pi*u2)
          normal = mean + sigma*y2
          flag = 0
        endif 
      end function normal 
    End Module ran_mod

    !***********************************************************************************
    module varmod
    implicit none
    
    ! Declare constant variables

    ! Input/output units used
    INTEGER, PARAMETER :: lin=1,lind=2,lins=3,ParametersFile=15,ObFunpar=1001
    ! Program version
    REAL, PARAMETER :: VERSION=1.0

    ! Maximum number of categories (classes of values)
    INTEGER, PARAMETER :: MAXCUT=7

    ! Maximum dimensions of the simulation grid
    INTEGER, PARAMETER :: MAXX=150, MAXY=1, MAXZ=150

    ! Maximum dimensions of the training image
    INTEGER, PARAMETER :: MAXXTR=300, MAXYTR=1, MAXZTR=300

    ! Maximum dimensions of the data search neighborhood (odd numbers)
    INTEGER, PARAMETER :: MAXCTX =150, MAXCTY=1, MAXCTZ=150

    INTEGER, PARAMETER :: MAXXY=MAXX*MAXY, MAXXYTR=MAXXTR*MAXYTR
    INTEGER, PARAMETER :: MAXXYZ=MAXXY*MAXZ
    INTEGER, PARAMETER :: MAXXYZTR=MAXXYTR*MAXZTR
    INTEGER, PARAMETER :: MAXCTXY=MAXCTX*MAXCTY
    INTEGER, PARAMETER :: MAXCTXYZ=MAXCTXY*MAXCTZ

    ! Maximum number of original sample data
    INTEGER, PARAMETER :: MAXDAT=1000000

    ! Maximum number of conditioning nodes
    INTEGER, PARAMETER :: MAXNOD=5000

    ! Maximum number of multiple grids
    INTEGER, PARAMETER :: MAXMULT=5

    ! Minimum correction factor in the servosystem (if used)
    REAL, PARAMETER :: MINCOR=1.0

    REAL, PARAMETER :: EPSILON=1.0e-20, DEG2RAD=3.141592654/180.0
    INTEGER, PARAMETER :: UNEST=-99


    ! Variable declarations

    ! Data locations
    REAL, DIMENSION(MAXDAT) :: x, y, z

    ! Data values
    INTEGER, DIMENSION(MAXDAT) :: vr
    real,dimension(MAXDAT)::vr_d,vr_s,vr_ss

    ! Number of original sample data
    INTEGER :: nd

    ! Number of categories
    INTEGER :: ncut

    ! Category thresholds
    INTEGER, DIMENSION(MAXCUT) :: thres

    ! Dimension specifications of the simulation grid
    INTEGER :: nx, ny, nz, nxy, nxyz
    REAL :: xmn, xsiz, ymn, ysiz, zmn, zsiz

    ! Training images
    INTEGER, DIMENSION (MAXXTR,MAXYTR,MAXZTR,MAXMULT) :: trainim

    ! Dimensions of the training images
    INTEGER, DIMENSION (MAXMULT) :: nxtr, nytr, nztr, nxytr, nxyztr

    ! Integer debugging level (0=none,1=normal,3=serious)
    INTEGER :: idbg

    ! Number of realizations to generate
    INTEGER :: nsim

    ! Realization and number of conditioning data retained
    INTEGER, DIMENSION (MAXX,MAXY,MAXZ) :: simim, numcd,FaciesSimTmp,FaciesSimTmp1
    real, DIMENSION (MAXX,MAXY) :: ESnySim, EActualSeismic, CovSny, VarSnySim, VarActualSeismic, Rsny,ESnySimBMcMCIte, CovSnyBMcMCIte, VarSnySimBMcMCIte, RsnyBMcMCIte

    ! Parameters defining the search ellipse
    REAL, DIMENSION (MAXMULT) :: radius, radius1, radius2
    REAL, DIMENSION (MAXMULT) :: sanis1, sanis2
    REAL, DIMENSION (MAXMULT) :: sang1, sang2, sang3
    REAL, DIMENSION (3,3) :: rotmat

    ! Number of grid node locations in data search
    ! neighborhood (no search tree)
    INTEGER :: nlsearch

    ! Relative grid node coordinates in data search
    ! neighborhood
    INTEGER, DIMENSION (MAXCTXYZ) :: ixnode,iynode,iznode

    ! Number of conditioning data
    INTEGER :: nodmax

    ! Total number of multiple grids, number of multiple grids
    ! simulated using a search tree
    INTEGER :: nmult, streemult

    ! Current multiple grid number
    INTEGER :: ncoarse
    integer::ncoarset

    ! Dimensions of the current multiple grid
    INTEGER :: nxcoarse, nycoarse, nzcoarse
    INTEGER :: nxycoarse, nxyzcoarse

    ! Number of nodes in the data template
    INTEGER :: nltemplate

    ! Relative node coordinates in the data template (with search tree)
    real, DIMENSION (MAXNOD) :: ixtemplate,iytemplate,iztemplate

    ! Relative coordinates and values of conditioning data
    INTEGER, DIMENSION (MAXNOD) :: cnodex,cnodey,cnodez,cnodev,cnodev1

    ! Number of conditioning data retained
    INTEGER :: ncnode

    ! Index of the nodes successively visited in the current multiple grid
    real, DIMENSION(MAXXY) :: order 
    
    !The number of chains
    integer::NumChains

    !ObjectiveFuction£º1--RMSE, and 2--Correlation Coefficient
    integer::ObjectiveFuction
    
    !assignway£º1--single grid, and 2--patch
    integer::assignway

    !Seismic parameters
    integer::wcy
    INTEGER::NumSampledWalve,HalfNumSampledWalve
    REAL,DIMENSION(0:500)::wavelet
    real, DIMENSION (MAXX,MAXY,MAXZ) :: snyor
    real, DIMENSION (MAXX,MAXY,MAXZ) :: ActualSeismic
    real,dimension(MAXX,MAXY,MAXZ)::DenSim,VpSim,VpSimtmp,DenSimtmp
    real,dimension(MAXX,MAXY,MAXZ)::SnySim

    real,dimension(MAXCUT)::tesro,tedro,tessro
    character(len=2) :: cFilename
    integer::sandncode,mudncode,unimud
    Real  time
    integer :: itera,IterationPrevious
    real :: Temperature
    real(kind=4):: SNVar
    real,DIMENSION  (MAXCUT) ::EDen,VDen,EVp,VVp,slopeVp,interceptVp,ResMeanVp,RMSEVp
    real :: Allcon,altcon,alttoAllcon
    integer :: Maxiterations,MaxMultiGridIterations
    real :: TIScanning,DMathcingRate,CoolingRate,CoolingCriteria,GridCriteria,InitialTemperature

    end module
    !***********************************************************************************
    
    
    
    !***********************************************************************************
    program main
    !c----------------------------------------------------------------------
    !c
    !c
    !c
    !c----------------------------------------------------------------------
    use varmod
    use ran_mod
    implicit none
    integer::imult
    integer::ix,iy,iz
    integer::i,j,k
    character(100)::namefl,ParametersFileName,iFacieOutPath,iDensityOutPath,iVelocityOutPath,iSyntheticRecordOutPath
    character(len=100)::FaciesOutfl,SeismicOutfl,DenOutfl,VpOutfl,AIoutfl,syoutfl,ChainFile
    real::Tsnycc,AverageCC
    character(len=4)::ctemp
    integer::iChain
    
    call SYSTEM('md '//trim(adjustl('Results')))
    call SYSTEM('del/q '//'Results\')
    
    call readparm

    WRITE(*,*) ' ObjectiveFuction =',ObjectiveFuction
    WRITE(*,*) ' Assignway =',Assignway
    
    EActualSeismic(1:nx,1:ny) = 0.0
    VarActualSeismic(1:nx,1:ny) = 0.0
    
    do ix = 1,nx
        do iy = 1,ny
            do iz = nz, 1,-1
                EActualSeismic(ix,iy)=EActualSeismic(ix,iy)+ActualSeismic(ix,iy,iz)
            end do                   
            EActualSeismic(ix,iy)=real(EActualSeismic(ix,iy)/nz)
            do iz=1,nz
                VarActualSeismic(ix,iy)=VarActualSeismic(ix,iy)+(ActualSeismic(ix,iy,iz)-EActualSeismic(ix,iy))**2
            end do
            VarActualSeismic(ix,iy)=real(VarActualSeismic(ix,iy)/nz)
        end do
    end do
    
    do iChain=1,NumChains

        write(ctemp,'(i4)')iChain
        ChainFile="Results\chain"// trim(adjustl(ctemp))
        call SYSTEM('md '//trim(adjustl(ChainFile)))
        
        numcd(1:nx,1:ny,1:nz)=UNEST
        DenSim(1:nx,1:ny,1:nz)=0.1
        VpSim(1:nx,1:ny,1:nz)=0.1
        SnySim(1:nx,1:ny,1:nz)=UNEST
        Tsnycc=0.0
    
        !open(ParametersFile,file="Results\Parameters.txt",status='UNKNOWN')
        ParametersFileName = "Results\Parameters"//"Chain"// trim(adjustl(ctemp))//".txt"
        open(ParametersFile,file=ParametersFileName,status='UNKNOWN')
    
        !!!
        if (ObjectiveFuction==1) then
        !OPEN(ObFunpar,file='Results\02RMSE.txt',status='UNKNOWN')
        WRITE(ParametersFile,*) "MP_ite  ","BMcMC_ite  ","Grid_Level      ","RMSE      "," AR   ","Temperatures      ","Time  "
        !WRITE(ObFunpar,*) "BMcMC_ite      ","RMSE      "
        else
        !OPEN(ObFunpar,file='Results\02CCAverage.txt',status='UNKNOWN')
        WRITE(ParametersFile,*)"MP_ite  ","BMcMC_ite  ","Grid_Level      ","CC       "," AR   ","Temperatures      ","Time "
        !WRITE(ObFunpar,*) "BMcMC_ite     ","CCAverage      "
        end if
    
        WCY=1
        itera=0
        imult=nmult
        alttoAllcon=1.0
        Temperature = InitialTemperature
        FaciesSimTmp(1:nx,1:ny,1:nz) = UNEST
    
        !call SYSTEM('md '//trim(adjustl('Results\FaciesOut')))
        !call SYSTEM('md '//trim(adjustl('Results\DensityOut')))
        !call SYSTEM('md '//trim(adjustl('Results\VelocityOut')))
        !call SYSTEM('md '//trim(adjustl('Results\SyntheticRecordOut')))
        
        iFacieOutPath = trim(adjustl(ChainFile))//"\FaciesOut"
        iDensityOutPath = trim(adjustl(ChainFile))//"\DensityOut"
        iVelocityOutPath = trim(adjustl(ChainFile))//"\VelocityOut"
        iSyntheticRecordOutPath = trim(adjustl(ChainFile))//"\SyntheticRecordOut"
        call SYSTEM('md '//trim(adjustl(iFacieOutPath)))
        call SYSTEM('md '//trim(adjustl(iDensityOutPath)))
        call SYSTEM('md '//trim(adjustl(iVelocityOutPath)))
        call SYSTEM('md '//trim(adjustl(iSyntheticRecordOutPath)))
    
        DO WCY = 1,1000
            IF(itera.GT.Maxiterations) goto 50
            write(ctemp,'(i4)')WCY
        
            !FaciesOutfl="Results\FaciesOut\Facies"// trim(adjustl(cTemp)) 
            !DenOutfl="Results\DensityOut\Density"// trim(adjustl(cTemp))
            !VpOutfl="Results\VelocityOut\Velocity"// trim(adjustl(cTemp))
            !SeismicOutfl="Results\SyntheticRecordOut\SyntheticRecord"// trim(adjustl(cTemp))
            !syoutfl="Results\CCOut"// trim(adjustl(cTemp))
            !AIoutfl="Results\AIOut"// trim(adjustl(cTemp))
            
            FaciesOutfl=trim(adjustl(iFacieOutPath))//"\Facies"// trim(adjustl(cTemp)) 
            DenOutfl=trim(adjustl(iDensityOutPath))//"\DensityOut"// trim(adjustl(cTemp))
            VpOutfl=trim(adjustl(iVelocityOutPath))//"\VelocityOut"// trim(adjustl(cTemp))
            SeismicOutfl=trim(adjustl(iSyntheticRecordOutPath))//"\SyntheticRecordOut"// trim(adjustl(cTemp))
            !syoutfl="Results\CCOut"// trim(adjustl(cTemp))
            !AIoutfl="Results\AIOut"// trim(adjustl(cTemp))
        
            OPEN(1,file=FaciesOutfl,status='UNKNOWN')
            WRITE(1,511)
511         format('PETREL:Property',/,'1')
            namefl='Facies_'// trim(adjustl(cTemp))
            WRITE(1,*)namefl

            OPEN(2,file=SeismicOutfl,status='UNKNOWN')
            WRITE(2,512)
512         format('PETREL:Property',/,'1')
            namefl='Seismic_SNY_'// trim(adjustl(cTemp))
            WRITE(2,*)namefl
!           OPEN(3,file=syoutfl,status='UNKNOWN')
!           WRITE(3,5112)
!5112       format('PETREL:Property',/,'1')
!           namefl='Seismic_r_'// trim(adjustl(cTemp))
!           WRITE(3,*)namefl
            OPEN(5,file=DenOutfl,status='UNKNOWN')
            WRITE(5,515)
515         format('PETREL:Property',/,'1')
            namefl='Density_'// trim(adjustl(cTemp))
            WRITE(5,*)namefl
            OPEN(6,file=VpOutfl,status='UNKNOWN')
            WRITE(6,516)
516         format('PETREL:Property',/,'1')
            namefl='P-velocity '// trim(adjustl(cTemp))
            WRITE(6,*)namefl
!           OPEN(7,file=AIoutfl,status='UNKNOWN')
!           WRITE(7,517)
!517        format('PETREL:Property',/,'1')
!           namefl='Seismic_PIM_'// trim(adjustl(cTemp))
!           WRITE(7,*)namefl
             
            simim(1:nx,1:ny,1:nz) = UNEST
            ESnySim(1:nx,1:ny) = 0.0
            CovSny(1:nx,1:ny) = 0.0
            VarSnySim(1:nx,1:ny) = 0.0
            Rsny(1:nx,1:ny) = 0.0
            AverageCC = 0.0
            nsim = 0
            Allcon = 0.0
            altcon = 0.0
        
            if(alttoAllcon<GridCriteria.or.((itera-IterationPrevious)>MaxMultiGridIterations).and.imult>1) then
                imult = imult-1
                IterationPrevious = itera
            end if
            ncoarse = 2**(imult-1)
            nxcoarse = MAX(1,(nx-1)/ncoarse+1)
            nycoarse = MAX(1,(ny-1)/ncoarse+1)
            nzcoarse = MAX(1,(nz-1)/ncoarse+1)
            nxycoarse = nxcoarse*nycoarse
            nxyzcoarse = nxycoarse*nzcoarse
            ncoarset = ncoarse
            
            WRITE(*,*) 'WCY=',WCY,'     Grid_Level:', imult

            ! Set up the spiral search:
            CALL SortTemplate(imult)

            ! Rescale data template used to construct the seach tree
            ixnode(1:nltemplate) = ncoarse*ixtemplate(1:nltemplate)
            iynode(1:nltemplate) = ncoarse*iytemplate(1:nltemplate)
            iznode(1:nltemplate) = ncoarse*iztemplate(1:nltemplate)

            CALL AssignData()

            ! Work out a random path for this realization:
            CALL RandomPath()

            CALL Simulation(imult)

            alttoAllcon=altcon/Allcon
              
            do iz = nz, 1,-1
                do iy = 1, ny
                    do ix = 1, nx                  
                        if(simim(ix,iy,iz).ge.1)then
                            write(1, '(2i6)') thres(simim(ix,iy,iz))
                        else
                            Call upscale(ix,iy,iz)
                            if(simim(ix,iy,iz)<1) then 
                                write(*,*) 'Error:simim(ix,iy,iz)<1'
                                stop
                            end if
                            write(1, '(2i6)') thres(simim(ix,iy,iz))
                        end if
                        FaciesSimTmp(ix,iy,iz) = simim(ix,iy,iz)
                        DenSimtmp(ix,iy,iz) = DenSim(ix,iy,iz)
                        VpSimtmp(ix,iy,iz) = VpSim(ix,iy,iz)
                    
                        if(simim(ix,iy,iz)==trainim(ix,iy,iz,1)) nsim=nsim+1
                        write(5, *) DenSim(ix,iy,iz)
                        write(6, *) VpSim(ix,iy,iz)
                        !write(7, *) DenSim(ix,iy,iz)*VpSim(ix,iy,iz)
                    end do
                end do
            end do
        
            do ix = 1,nx
                do iy = 1,ny
                    call convolution(nz,ix,iy,1,DenSim(ix,iy,1:nz),VpSim(ix,iy,1:nz),SnySim(ix,iy,1:nz))
                end do
            end do
        
            do iz = nz, 1,-1
                do iy = 1, ny
                    do ix = 1, nx                  
                        write(2, *) SnySim(ix,iy,iz)
                        Tsnycc = Tsnycc+ABS(SnySim(ix,iy,iz)-ActualSeismic(ix,iy,iz))
                    end do
                end do
            end do        

            do ix = 1,nx
                do iy = 1,ny
                    do iz = nz, 1,-1
                        ESnySim(ix,iy)=ESnySim(ix,iy)+SnySim(ix,iy,iz)
                        CovSny(ix,iy)=CovSny(ix,iy)+SnySim(ix,iy,iz)*ActualSeismic(ix,iy,iz)
                    end do                   
                    ESnySim(ix,iy)=real(ESnySim(ix,iy)/nz)
                    CovSny(ix,iy)=real(CovSny(ix,iy)/nz)
                    do iz=1,nz
                        VarSnySim(ix,iy)=VarSnySim(ix,iy)+(SnySim(ix,iy,iz)-ESnySim(ix,iy))**2
                    end do
                    VarSnySim(ix,iy)=real(VarSnySim(ix,iy)/nz)
                    Rsny(ix,iy)=(CovSny(ix,iy)-ESnySim(ix,iy)*EActualSeismic(ix,iy))/((VarSnySim(ix,iy)*VarActualSeismic(ix,iy))**0.5)
                    AverageCC = AverageCC+Rsny(ix,iy)
                    !write(3, *) Rsny(ix,iy)
                end do
            end do
            AverageCC = AverageCC/nx/ny
        
            close(1)
            close(2)
            !close(3)
            close(5)
            close(6)
            !close(7)
        
            call CPU_TIME(time)
        
            if (ObjectiveFuction==1) then
                WRITE(ParametersFile,1002) WCY,int(itera+Allcon),imult,Tsnycc,real(altcon/Allcon),Temperature,time
1002            format(I8,I16,I12,F20.1,F10.3,E14.4,F13.1)
                WRITE(ParametersFile,1004) WCY,int(itera+Allcon),imult,Tsnycc,real(altcon/Allcon),Temperature,time
1004            format(I8,I16,I12,F20.1,F10.3,E14.4,F13.1)
            else
                WRITE(ParametersFile,1003) WCY,int(itera+Allcon),imult,AverageCC,real(altcon/Allcon),Temperature,time
1003            format(I8,I16,I12,F16.3,F9.3,E14.4,F13.1)
            end if
        
            !!!Adaptive cooling schame
            if(real(altcon/Allcon)>CoolingCriteria) then
                Temperature=Temperature*CoolingRate
            end if

        end do

50      continue
        !close(ObFunpar)
        close(ParametersFile)
    end do
    
    
    end
    !***********************************************************************************
    