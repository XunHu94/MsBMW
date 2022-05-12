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
        call random_number(ran)     ! built in fortran 90 random number function  
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
    INTEGER, PARAMETER :: lin=1,lind=2,lins=3,dzsnycc=15,multigrid=16,timpar=7,RMSEpar=1001,wellpar=1002,recordpar=1005
    ! Program version
    REAL, PARAMETER :: VERSION=1.0

    ! Maximum number of categories (classes of values)
    INTEGER, PARAMETER :: MAXCUT=7

    ! Maximum dimensions of the simulation grid
    INTEGER, PARAMETER :: MAXX=160, MAXY=1, MAXZ=80

    ! Maximum dimensions of the training image
    INTEGER, PARAMETER :: MAXXTR=300, MAXYTR=1, MAXZTR=160

    ! Maximum dimensions of the data search neighborhood (odd numbers)
    INTEGER, PARAMETER :: MAXCTX =160, MAXCTY=1, MAXCTZ=80

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
    real, DIMENSION (MAXX,MAXY) :: ESnySim, EAngleSeismic, CovSny, VarSnySim, VarAngleSeismic, Rsny

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

    !assignway：1-single grid, and 2--patch
    integer::assignway

    !Seismic parameters
    integer::ncy=1,wcy=6
    INTEGER::nnd
    REAL,DIMENSION(0:500)::wavelet
    real, DIMENSION (MAXX,MAXY,MAXZ) :: snyor
    real, DIMENSION (MAXX,MAXY,MAXZ) :: AngleSeismic
    real,dimension(MAXX,MAXY,MAXZ)::DenSim,VpSim,VpSimtmp,DenSimtmp
    real,dimension(MAXX,MAXY,MAXZ)::SnySim

    real,dimension(MAXCUT)::tesro,tedro,tessro
    character(len=2) :: cFilename
    integer::sandncode,mudncode,unimud
    real,DIMENSION  (MAXXYZ) ::unif,uni
    Real time_begin , time_end,time_end1,time1,time2,time3
    integer :: itera
    real :: TK
    real(kind=4):: SNVar,seisRMS
    real,DIMENSION  (MAXCUT) ::EDen,VDen,EVp,VVp,slopeVp,interceptVp,ResMeanVp,RMSEVp
    real :: Allcon,altcon,alttoAllcon
    integer :: Maxiterations
    real :: TIScanning,DMathcingRate,CoolingRate,CoolingCriteria,GridCriteria,T0

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
    integer::ind
    integer::ix,iy,iz
    integer::i,j,k
    character(20)::namefl
    character(len=30),dimension(250)::FaciesOutfl,SeismicOutfl,DenOutfl,VpOutfl,AIoutfl,Seisrofl,syoutfl
    real,dimension(0:250)::Tsnycc,Tsnyccj,TVp,TDen
    character(len=4)::ctemp
    
    call SYSTEM('del/q '//'SIMOUTFL\') 
    
    open(timpar,file="SIMOUTFL\time.txt",status="UNKNOWN")
    OPEN(RMSEpar,file='SIMOUTFL\SeismicRMSE',status='UNKNOWN')
    OPEN(wellpar,file='SIMOUTFL\wells1',status='UNKNOWN')
    OPEN(recordpar,file='SIMOUTFL\recordpar1',status='UNKNOWN')
    call CPU_TIME(time_begin)
    write(timpar,*) 0,"          ",time_begin
     
    call readparm

    WRITE(*,*) ' Assignway =',Assignway

    numcd(1:nx,1:ny,1:nz)=UNEST
    DenSim(1:nx,1:ny,1:nz)=0.1
    VpSim(1:nx,1:ny,1:nz)=0.1
    SnySim(1:nx,1:ny,1:nz)=UNEST
    Tsnycc(0:250)=0.0
    TDen(0:250)=0.0
    TVp(0:250)=0.0

    open(dzsnycc,file="SIMOUTFL\zsnycc.txt",status='UNKNOWN')
    open(multigrid,file="SIMOUTFL\multigrid.txt",status='UNKNOWN')
    !
    WRITE(dzsnycc,*)"MP ite.    ","RMSE    ","BMcMC ite.    ","suggested models in a mp ite.    ","the number of the acceptance of suggested models in a mp ite.    "," Beta parameter   ","run time    ","temperature    "
    WRITE(multigrid,*) "grid level ","MP iterations ","BMcMC iterations "
    WCY=1
    Tsnyccj(0)=-10000000
    Tsnyccj(1)=-10000000
    itera=0
    imult=nmult
    alttoAllcon=1.0
    TK = T0
    FaciesSimTmp(1:nx,1:ny,1:nz) = UNEST
    
    DO WCY = 1,250
        call CPU_TIME(time1)
        IF(WCY.GT.250) EXIT
        if(wcy.ge.250)then
            write(*,*)" wcy>250"
            stop
        else
            write(ctemp,'(i4)')WCY
        end if
        write(*,*) WCY
        
        FaciesOutfl(WCY)="SIMOUTFL\FaciesOut"// trim(adjustl(cTemp))
        SeismicOutfl(WCY)="SIMOUTFL\SeismicOut"// trim(adjustl(cTemp))
        syoutfl(WCY)="SIMOUTFL\rsnyOut"// trim(adjustl(cTemp))
        DenOutfl(WCY)="SIMOUTFL\DenOut"// trim(adjustl(cTemp))
        VpOutfl(WCY)="SIMOUTFL\VpOut"// trim(adjustl(cTemp))
        AIoutfl(WCY)="SIMOUTFL\AIOut"// trim(adjustl(cTemp))
        Seisrofl(WCY)="SIMOUTFL\Seisro"// trim(adjustl(cTemp))
        
        OPEN(1,file=FaciesOutfl(WCY),status='UNKNOWN')
        WRITE(1,511)
511     format('PETREL:Property',/,'1')
        namefl='Facies_'// trim(adjustl(cTemp))
        WRITE(1,*)namefl

        OPEN(2,file=SeismicOutfl(WCY),status='UNKNOWN')
        WRITE(2,512)
512     format('PETREL:Property',/,'1')
        namefl='Seismic_SNY_'// trim(adjustl(cTemp))
        WRITE(2,*)namefl
        OPEN(3,file=syoutfl(WCY),status='UNKNOWN')
        WRITE(3,5112)
5112     format('PETREL:Property',/,'1')
        namefl='Seismic_r_'// trim(adjustl(cTemp))
        WRITE(3,*)namefl
        OPEN(5,file=DenOutfl(WCY),status='UNKNOWN')
        WRITE(5,515)
515     format('PETREL:Property',/,'1')
        namefl='Density_'// trim(adjustl(cTemp))
        WRITE(5,*)namefl
        OPEN(6,file=VpOutfl(WCY),status='UNKNOWN')
        WRITE(6,516)
516     format('PETREL:Property',/,'1')
        namefl='P-velocity '// trim(adjustl(cTemp))
        WRITE(6,*)namefl
        OPEN(7,file=AIoutfl(WCY),status='UNKNOWN')
        WRITE(7,517)
517     format('PETREL:Property',/,'1')
        namefl='Seismic_PIM_'// trim(adjustl(cTemp))
        WRITE(7,*)namefl
        OPEN(8,file=Seisrofl(WCY),status='UNKNOWN')
        WRITE(8,518)
518     format('PETREL:Property',/,'1')
        namefl='Seismic_RO'// trim(adjustl(cTemp))
        WRITE(8,*)namefl
             
        simim(1:nx,1:ny,1:nz) = UNEST
        ESnySim(1:nx,1:ny) = 0.0
        EAngleSeismic(1:nx,1:ny) = 0.0
        CovSny(1:nx,1:ny) = 0.0
        VarSnySim(1:nx,1:ny) = 0.0
        VarAngleSeismic(1:nx,1:ny) = 0.0
        Rsny(1:nx,1:ny) = 0.0
        nsim = 0
        Allcon = 0.0
        altcon = 0.0
        
        WRITE(*,*) 'WCY ',WCY
        
        !if(alttoAllcon<GridCriteria.or.(wcy>(nmult-imult+1)*100.or.itera>(nmult-imult+1)*10000)) then!多点迭代和反演迭代任选其一约束
        if(alttoAllcon<GridCriteria) then
            imult = imult-1
            WRITE(multigrid,*) imult,wcy,itera
        end if    
        if(imult<1) imult = 1
        ncoarse = 2**(imult-1)
        nxcoarse = MAX(1,(nx-1)/ncoarse+1)
        nycoarse = MAX(1,(ny-1)/ncoarse+1)
        nzcoarse = MAX(1,(nz-1)/ncoarse+1)
        nxycoarse = nxcoarse*nycoarse
        nxyzcoarse = nxycoarse*nzcoarse
        ncoarset = ncoarse
        WRITE(*,*) 'working on grid: ', imult

        ! Set up the spiral search:
        CALL SortTemplate(imult)

        ! Rescale data template used to construct the seach tree
        ixnode(1:nltemplate) = ncoarse*ixtemplate(1:nltemplate)
        iynode(1:nltemplate) = ncoarse*iytemplate(1:nltemplate)
        iznode(1:nltemplate) = ncoarse*iztemplate(1:nltemplate)

        CALL AssignData1()

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
                            write(*,*) 'Error:simim(ix,iy,iz)<1--line≈472'
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
                    write(7, *) DenSim(ix,iy,iz)*VpSim(ix,iy,iz)
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
                    write(8, *) SnySim(ix,iy,iz)-AngleSeismic(ix,iy,iz)
                    Tsnycc(WCY) = Tsnycc(WCY)+ABS(SnySim(ix,iy,iz)-AngleSeismic(ix,iy,iz))
                end do
            end do
        end do        

        do ix = 1,nx
            do iy = 1,ny
                do iz = nz, 1,-1
                    ESnySim(ix,iy)=ESnySim(ix,iy)+SnySim(ix,iy,iz)
                    EAngleSeismic(ix,iy)=EAngleSeismic(ix,iy)+AngleSeismic(ix,iy,iz)
                    CovSny(ix,iy)=CovSny(ix,iy)+SnySim(ix,iy,iz)*AngleSeismic(ix,iy,iz)
                end do                   
                ESnySim(ix,iy)=real(ESnySim(ix,iy)/nz)
                EAngleSeismic(ix,iy)=real(EAngleSeismic(ix,iy)/nz)
                CovSny(ix,iy)=real(CovSny(ix,iy)/nz)
                do iz=1,nz
                    VarSnySim(ix,iy)=VarSnySim(ix,iy)+(SnySim(ix,iy,iz)-ESnySim(ix,iy))**2
                    VarAngleSeismic(ix,iy)=VarAngleSeismic(ix,iy)+(AngleSeismic(ix,iy,iz)-EAngleSeismic(ix,iy))**2
                end do
                VarSnySim(ix,iy)=real(VarSnySim(ix,iy)/nz)
                VarAngleSeismic(ix,iy)=real(VarAngleSeismic(ix,iy)/nz)
                Rsny(ix,iy)=(CovSny(ix,iy)-ESnySim(ix,iy)*EAngleSeismic(ix,iy))/((VarSnySim(ix,iy)*VarAngleSeismic(ix,iy))**0.5)
                write(3, *) Rsny(ix,iy)
            end do
        end do 
        
        close(1)
        close(2)
        close(3)
        close(5)
        close(6)
        close(7)
        close(8)
        
        call CPU_TIME(time2)
        
        WRITE(dzsnycc,*) WCY,Tsnycc(WCY),itera+Allcon,Allcon,altcon,real(altcon/Allcon),time2-time1,TK
        !annealing temperature
        if(real(altcon/Allcon)>CoolingCriteria) then
            TK=TK*CoolingRate
        end if

        call cpu_time(time_end1)    
        write(timpar,*) wcy,"          ",time_end1 - time_begin
        time_begin=time_end1
        if(itera==Maxiterations) go to 50
    end do

50  continue
    call cpu_time(time_end)
    write(timpar,*) '   total time',"     ",time_end
    close(timpar)
    close(RMSEpar)
    close(wellpar)
    close(recordpar)
    close(dzsnycc)
    close(multigrid)
    
    end
    !***********************************************************************************
     

    !***********************************************************************************
    subroutine readparm
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    implicit none
    real,dimension(10)::var,var_d,var_s,var_ss
    INTEGER, DIMENSION(2) :: seed
    character::ConFaciesFile*40,ConDenFile*40,ConVpFile*40,str*40,TemplateFile*40,datafl1*40,AcSeismicFile*40
    CHARACTER (LEN=20), DIMENSION(MAXMULT) :: trainfl
    INTEGER, DIMENSION (MAXMULT) :: ivrltr
    logical::testfl
    integer::ioerror
    integer::ixl,iyl,izl,ivrl,icut
    integer::ix,iy,iz,imult,ntr
    integer::i,j,k,nvari
    integer::ixv,Facie,TatolFacie
    
    WRITE(*,999) VERSION
999 format(/' MsBMW Version: ',f8.3/)

    WRITE(*,*) 'Which parameter file do you want to use?'
    INQUIRE(file="MsBMW.par",exist=testfl)
    OPEN(lin,file="MsBMW.par",status='OLD')

12  READ(lin,'(a4)',end=198) str(1:4)
    IF(str(1:4)/='STAR') GO TO 12

    READ(lin,'(a20)',err=198) ConFaciesFile
    CALL chknam(ConFaciesFile,40)
    WRITE(*,*) ' data file = ',ConFaciesFile

    READ(lin,'(a20)',err=198) ConDenFile
    CALL chknam(ConDenFile,40)
    WRITE(*,*) ' data file = ',ConDenFile

    READ(lin,'(a20)',err=198) ConVpFile
    CALL chknam(ConVpFile,40)
    WRITE(*,*) ' data file = ',ConVpFile

    READ(lin,'(a20)',err=198) datafl1
    CALL chknam(datafl1,40)
    WRITE(*,*) ' data file = ',datafl1

    READ(lin,'(a20)',err=198) AcSeismicFile
    CALL chknam(AcSeismicFile,40)
    WRITE(*,*) ' data file = ',AcSeismicFile

    READ(lin,*,err=198) assignway
    WRITE(*,*) ' assignway = ',assignway

    READ(lin,*,err=198) ixl,iyl,izl,ivrl
    WRITE(*,*) ' input columns = ',ixl,iyl,izl,ivrl

    READ(lin,*,err=198) ncut
    WRITE(*,*) ' number of categories = ',ncut
    IF(ncut>MAXCUT) THEN
        WRITE(*,*) 'ERROR: maximum number of categories: ',MAXCUT
        WRITE(*,*) '       you have asked for : ',ncut
    END IF

    READ(lin,*,err=198) (thres(icut), icut=1,ncut)
    WRITE(*,*) ' categories = ', (thres(icut), icut=1,ncut)

    READ(lin,*,err=198) nsim
    WRITE(*,*) ' number of realizations = ',nsim

    READ(lin,*,err=198) nx,xmn,xsiz
    WRITE(*,*) ' X grid specification = ',nx,xmn,xsiz

    READ(lin,*,err=198) ny,ymn,ysiz
    WRITE(*,*) ' Y grid specification = ',ny,ymn,ysiz

    READ(lin,*,err=198) nz,zmn,zsiz
    WRITE(*,*) ' Z grid specification = ',nz,zmn,zsiz

    IF(nx>MAXX.or.ny>MAXY.or.nz>MAXZ) THEN
        WRITE(*,*) 'ERROR: available grid size: ',MAXX,MAXY,MAXZ
        WRITE(*,*) '       you have asked for : ',nx,ny,nz
        STOP
    END IF
    nxy  = nx*ny
    nxyz=nxy*nz

    READ(lin,*,err=198) ixv
    WRITE(*,*) ' random number seed = ',ixv

    ! Initialize the random seed of the simulation:
    seed(1)=ixv
    seed(2)=ixv+1
    CALL random_seed(PUT=seed(1:2))

    READ(lin,'(a20)',err=198) TemplateFile
    CALL chknam(TemplateFile,20)
    WRITE(*,*) ' data template file = ',TemplateFile

    READ(lin,*,err=198) nodmax
    WRITE(*,*) ' maximum conditioning data = ',nodmax
    IF(nodmax.gt.MAXNOD) THEN
        WRITE(*,*) 'ERROR: maximum available cond. data: ',MAXNOD
        WRITE(*,*) '       you have asked for :  ',nodmax
        STOP
    END IF
    
    READ(lin,*,err=198) Maxiterations
    WRITE(*,*) ' Maxiterations = ',Maxiterations
    
    READ(lin,*,err=198) GridCriteria
    WRITE(*,*) ' The grid degradation criteria = ',GridCriteria
    
    READ(lin,*,err=198) TIScanning
    WRITE(*,*) ' the search scope of TI = ',TIScanning
    
    READ(lin,*,err=198) DMathcingRate
    WRITE(*,*) ' Data events match rate = ',DMathcingRate
    
    READ(lin,*,err=198) T0
    WRITE(*,*) ' Initial temperature = ',T0
    
    READ(lin,*,err=198) CoolingRate
    WRITE(*,*) ' Cooling rate = ',CoolingRate
    
    READ(lin,*,err=198) CoolingCriteria
    WRITE(*,*) ' The cooling criteria = ',CoolingCriteria

    READ(lin,*,err=198) nmult, streemult
    WRITE(*,*) ' multiple grid simulation = ',nmult,streemult
    IF(nmult>MAXMULT) THEN
        WRITE(*,*) 'ERROR: maximum number of mult. grids: ',MAXMULT
        WRITE(*,*) '       you have asked for :  ',nmult
        STOP
    END IF

    IF(streemult>nmult) THEN
        WRITE(*,*) 'ERROR: the number of grids using the search tree'
        WRITE(*,*) 'must be less than the total number of mult. grids.'
        STOP
    END IF

    ! Now read the information related to each multiple grid.
    ! If only one training image is provided in the parameter file,
    ! this single image will be used for all multiple grids.

    DO imult=nmult,1,-1
        ncoarse=2**(imult-1)
        WRITE(*,*) 'Multiple grid ', imult
        READ(lin,'(a20)',IOSTAT=ioerror) trainfl(imult)
        IF(ioerror<0.AND.imult==nmult) THEN
            STOP 'ERROR in parameter file!'
        ELSE IF(ioerror<0) THEN
            trainfl(imult)=trainfl(nmult)
            WRITE(*,*) ' training image file = ',trainfl(imult)
            nxtr(imult)=nxtr(nmult)
            nytr(imult)=nytr(nmult)
            nztr(imult)=nztr(nmult)
            WRITE(*,*) ' training grid dimensions = ',nxtr(imult),nytr(imult),nztr(imult)
            ivrltr(imult)=ivrltr(nmult)
            WRITE(*,*) ' column for variable = ',ivrltr(imult)
            radius(imult)=radius(nmult)
            radius1(imult)=radius1(nmult)
            radius2(imult)=radius2(nmult)
            WRITE(*,*) ' data search neighborhood radii = ',radius(imult),radius1(imult),radius2(imult)
            sang1(imult)=sang1(nmult)
            sang2(imult)=sang2(nmult)
            sang3(imult)=sang3(nmult)
            WRITE(*,*) ' search anisotropy angles = ',sang1(imult),sang2(imult),sang3(imult)
        ELSE
            CALL chknam(trainfl(imult),20)
            WRITE(*,*) ' training image file = ',trainfl(imult)

            READ(lin,*,err=198) nxtr(imult), nytr(imult), nztr(imult)
            WRITE(*,*) ' training grid dimensions = ',nxtr(imult),nytr(imult),nztr(imult)

            READ(lin,*,err=198) ivrltr(imult)
            WRITE(*,*) ' column for variable = ',ivrltr(imult)

            READ(lin,*,err=198) radius(imult),radius1(imult),radius2(imult)
            WRITE(*,*) ' data search neighborhood radii = ',radius(imult),radius1(imult),radius2(imult)

            READ(lin,*,err=198) sang1(imult),sang2(imult),sang3(imult)
            WRITE(*,*) ' search anisotropy angles = ',sang1(imult),sang2(imult),sang3(imult)
        END IF

        IF(nxtr(imult)>MAXXTR.or.nytr(imult)>MAXYTR) THEN
            WRITE(*,*) 'ERROR: available train. grid size: ',MAXXTR,MAXYTR, MAXZTR
            WRITE(*,*) '       you have asked for : ',nxtr(imult),nytr(imult),nztr(imult)
            STOP
        END IF
        nxytr(imult) = nxtr(imult)*nytr(imult)
        nxyztr(imult)=nxytr(imult)*nztr(imult)

        IF(radius(imult)<EPSILON.OR.radius1(imult)<EPSILON.OR.radius2(imult)<EPSILON) STOP 'radius must be greater than zero'
        sanis1(imult)=radius1(imult)/radius(imult)
        sanis2(imult)=radius2(imult)/radius(imult)

    END DO

    !
    ! Now, read the data if the file exists
    !
    INQUIRE(file=ConFaciesFile,exist=testfl)
    IF(.NOT.testfl) THEN
        WRITE(*,*) 'WARNING data file ',ConFaciesFile,' does not exist!'
        WRITE(*,*) '   - Hope your intention was to create an ','unconditional simulation'
        nd=0
    ELSE
        WRITE(*,*) 'Reading input data'
        OPEN(lin,file=ConFaciesFile,status='OLD')
        READ(lin,*,err=199)
        READ(lin,*,err=199) nvari
        !*******************************************************************
        OPEN(lind,file=ConDenFile,status='OLD')
        READ(lind,*,err=199)
        READ(lind,*,err=199)
        OPEN(lins,file=ConVpFile,status='OLD')
        READ(lins,*,err=199)
        READ(lins,*,err=199)
        !***************************************************************
        DO i=1,nvari
            READ(lin,*,err=199)
            READ(lind,*,err=199)
            READ(lins,*,err=199)
        END DO
        IF(ixl>nvari.OR.iyl>nvari.OR.izl>nvari.OR.ivrl>nvari) THEN
            WRITE(*,*) 'ERROR: you have asked for a column number'
            WRITE(*,*) '       greater than available in file'
            STOP
        END IF
        !
        ! Read all the data until the end of the file:
        ! nd: number of data read in the data file.
        !
        nd=0
        DO
            READ(lin,*,IOSTAT=ioerror) (var(j),j=1,nvari)
            IF(ioerror<0) EXIT
            READ(lind,*) (var_d(j),j=1,nvari)
            READ(lins,*) (var_s(j),j=1,nvari)
            nd=nd+1
            IF(nd>MAXDAT) STOP' ERROR exceeded MAXDAT - check source file'
            !
            ! Acceptable data, assign the value, X, Y coordinates:
            ! x,y: vectors of data coordinates
            ! vr: vector of data values
            !
            DO icut=1,ncut
                IF(nint(var(ivrl))==thres(icut)) THEN
                    vr(nd)=icut
                    EXIT
                END IF
            END DO
            vr_d(nd)=var_d(ivrl)
            vr_s(nd)=var_s(ivrl)
            vr_ss(nd)=var_ss(ivrl)
            IF(ixl>0) x(nd)=var(ixl)
            IF(iyl>0) y(nd)=var(iyl)
            IF(izl>0) z(nd)=var(izl)
        END DO
        CLOSE(lin)
    END IF

    !
    ! Now, read the training images if the files exist:
    !
    DO imult=nmult,1,-1
        WRITE(*,*) trainfl(imult)
        INQUIRE(file=trainfl(imult),exist=testfl)
        IF(.NOT.testfl) THEN
            STOP 'Error in training image file'
        ELSE
            WRITE(*,*) 'Reading training image ',imult
            OPEN(lin,file=trainfl(imult),status='OLD')
            READ(lin,*,err=197)
            READ(lin,*,err=197) nvari
            DO i=1,nvari
                READ(lin,*,err=197)
            END DO
            IF(ivrltr(imult)>nvari) THEN
                WRITE(*,*) 'ERROR: you have asked for a column number'
                WRITE(*,*) '       greater than available in file'
                STOP
            END IF
            !
            ! Read all the data until the end of the training file:
            ! ntr: number of data read in the training file.
            !
            ntr=0

            DO
                READ(lin,*, IOSTAT=ioerror) (var(j),j=1,nvari)
                IF(ioerror<0) EXIT
                ntr = ntr+1
                IF(ntr>nxyztr(imult)) STOP ' ERROR exceeded nxyztr -check inc file'
                iz = 1+(ntr-1)/nxytr(imult)
                iy = 1+(ntr-(iz-1)*nxytr(imult)-1)/nxtr(imult)
                ix = ntr-(iz-1)*nxytr(imult)-(iy-1)*nxtr(imult)
                DO icut=1,ncut
                    IF(nint(var(ivrltr(imult)))==thres(icut)) THEN
                        trainim(ix,iy,nztr(imult)-iz+1,imult) = icut
                        EXIT
                    END IF
                END DO
            END DO
            CLOSE(lin)
        END IF
    END DO

    !
    ! Now, read the file defining the data template if it exists:
    !
    INQUIRE(file=TemplateFile,exist=testfl)
    IF(.NOT.testfl) THEN
        STOP 'Error in data template file'
    ELSE
        WRITE(*,*) 'Reading data template file'
        OPEN(lin,file=TemplateFile,status='OLD')
        READ(lin,*,err=196)
        READ(lin,*,err=196) nvari
        DO i=1,nvari
            READ(lin,*,err=196)
        END DO
        IF(nvari/=3) STOP 'ERROR: the number of columns should be 3'
        !
        ! Read all the data locations until the end of the training file:
        ! nltemplate: number of data locations in the template.
        !
        nltemplate=0
        DO
            READ(lin,*, IOSTAT=ioerror) (var(j),j=1,nvari)
            IF(ioerror<0) EXIT
            nltemplate=nltemplate+1
            IF(nltemplate>MAXNOD) THEN
                STOP 'ERROR exceeded MAXNOD - check source code'
            END IF
            ixtemplate(nltemplate)=nint(var(1))
            iytemplate(nltemplate)=nint(var(2))
            iztemplate(nltemplate)=nint(var(3))
        END DO
        CLOSE(lin)
    END IF
    
    INQUIRE(file=datafl1,exist=testfl)
    OPEN(lin,file=datafl1,status='OLD')
    IF(.NOT.testfl) THEN
        WRITE(*,*) 'WARNING data file ',datafl1,' does not exist!'
    ELSE
        WRITE(*,*) 'Reading input data'
        nnd=0
        DO
            READ(lin,*,IOSTAT=ioerror) wavelet(nnd)
            IF(ioerror<0) EXIT
            nnd=nnd+1
        end do
        CLOSE(lin)
        nnd=nnd-1
    END IF

    INQUIRE(file=AcSeismicFile,exist=testfl)
    INQUIRE(file='SeismicData',exist=testfl)
    IF(.NOT.testfl) THEN
        STOP 'Error in Seismic data file'
    ELSE
        OPEN(lin,file='SeismicData',status='OLD')
        READ(lin,*,err=197)
        READ(lin,*,err=197) nvari
        DO i=1,nvari
            READ(lin,*,err=197)
        END DO
        ntr=0
        DO
            READ(lin,*, IOSTAT=ioerror) (var(j),j=1,nvari)
            IF(ioerror<0) EXIT
            ntr=ntr+1
            IF(ntr>nx*ny*nz) STOP ' ERROR exceeded nxyztr -check inc file'
            iz=1+(ntr-1)/(nx*ny)
            iy=1+(ntr-(iz-1)*(nx*ny)-1)/nx
            ix=ntr-(iz-1)*(nx*ny)-(iy-1)*nx
            AngleSeismic(ix,iy,nz-iz+1)=var(nvari)
        END DO
        CLOSE(lin)
    end if
    
    INQUIRE(file="pdf",exist=testfl)
    OPEN(1,file="pdf",status='OLD')
    
    read(1,*) TatolFacie
    do i=1,TatolFacie
        read(1,*) Facie
        READ(1,*,err=198) EDen(Facie+1),VDen(Facie+1)
        WRITE(*,*) ' EDen2,VDen2 = ',EDen(Facie+1),VDen(Facie+1) 
        READ(1,*,err=198) EVp(Facie+1),VVp(Facie+1)
        WRITE(*,*) ' EVp2,VVp2 = ',EVp(Facie+1),VVp(Facie+1)
        READ(1,*,err=198) slopeVp(Facie+1),interceptVp(Facie+1),ResMeanVp(Facie+1),RMSEVp(Facie+1)
    end do
    
    READ(1,*,err=198) SNVar
    WRITE(*,*) ' input columns = ',SNVar
    close(1)

    RETURN
    !
    ! Error in an Input File Somewhere:
    !
192 STOP 'ERROR in vertical proportions file!'
196 STOP 'ERROR in data template file!'
197 STOP 'ERROR in training image file!'
198 STOP 'ERROR in parameter file!'
199 STOP 'ERROR in data file!'
    end subroutine readparm
    !***********************************************************************************
    
    
    !***********************************************************************************
    SUBROUTINE RandomPath()
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    implicit none
    ! Declare local variables
    REAL , DIMENSION  (MAXXYZ) :: sim
    real , DIMENSION  (MAXXYZ) :: c,d,e,f,g,h
    integer::jx,jy,jz,ix,iy,iz,iix,iiy,iiz
    INTEGER :: ingd,ind

    CALL Random_Number(sim(:nxyzcoarse))

    DO ingd=1, nxycoarse
        jz=1
        jy=1+int((ingd-(jz-1)*nxycoarse-1)/nxcoarse)
        jx=ingd-(jz-1)*nxycoarse-(jy-1)*nxcoarse
        ix=(jx-1)*ncoarse+1
        iy=(jy-1)*ncoarse+1
        iz=(jz-1)*ncoarse+1

        ncnode=0
        do ind=1,nltemplate
            iix=ix+ixnode(ind)
            iiy=iy+iynode(ind)
            iiz=iz+iznode(ind)
            IF(iix>=1.and.iix<=nx.and.iiy>=1.and.iiy<=ny.and.iiz>=1.and.iiz<=nz) THEN
                IF(simim(iix,iiy,iiz)>UNEST) ncnode=ncnode+1
            end if
        end do
        order(ingd)=ingd
        sim(ingd)=sim(ingd)-ncnode
    end do
    CALL sortem(1,nxycoarse,sim,1,order,c,d,e,f,g,h)
    

    END SUBROUTINE RandomPath
    !***********************************************************************************


    !***********************************************************************************
    SUBROUTINE Simulation (imult)
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    use ran_mod
    implicit none
    ! Declare dummy arguments
    INTEGER, INTENT(IN) :: imult

    ! Declare local variables
    INTEGER :: ingd, ic, ix, iy, iz, irepo,ixr,iyr,izr
    INTEGER :: jx, jy, jz
    real::pptmp
    integer::pptmpi
    real,dimension(nz)::rangz,rangtm
    real , DIMENSION  (MAXXYZ) :: c,d,e,f,g,h
    integer::i,NcoarseNz
    integer :: kz
    real,dimension(MAXX,MAXY,MAXZ)::pden,pspeed,pspeeds
    
    irepo=max(1,min((nxycoarse/100),10000))
    !
    ! Main loop over all simulation grid nodes:
    !

    DO ingd=1, nxycoarse
        IF(mod(ingd,irepo)==0) WRITE(*,104) ingd
104     format('   currently on node ',i9)
        !
        ! Figure out the location of this point and make sure it has
        ! not been assigned a value already:
        !
        jz=1
        jy=1+int((order(ingd)-(jz-1)*nxycoarse-1)/nxcoarse)
        jx=order(ingd)-(jz-1)*nxycoarse-(jy-1)*nxcoarse
        ix=(jx-1)*ncoarse+1
        iy=(jy-1)*ncoarse+1
        iz=(jz-1)*ncoarse+1

        rangtm(1:nz)=2
        rangz(1:nz)=UNEST
        do i=1,(nz-1)/ncoarse+1
            rangz((i-1)*ncoarse+1)=(i-1)*ncoarse+1
            call RANDOM_NUMBER(pptmp)
            rangtm((i-1)*ncoarse+1)=pptmp
        end do
        CALL sortem(1,nz,rangtm,1,rangz,c,d,e,f,g,h)
        do i=1,(nz-1)/ncoarse+1
            iz=rangz(i)
            
            IF(simim(ix, iy, iz)<0) THEN
                
                ! First, get close conditioning data:
                CALL SearchClosestNodes(ix,iy,iz)

                ! If no conditioning data were found, use the marginal pdf as
                ! local cpdf, otherwise infer the local cpdf from the training image:
                IF(ncnode==0) THEN
                    IF(WCY.GT.1)THEN
                        simim(ix,iy,iz)=FaciesSimTmp(ix,iy,iz)
                        DenSim(ix,iy,iz)=DenSimtmp(ix,iy,iz)         
                        VpSim(ix,iy,iz)=VpSimtmp(ix,iy,iz)   
                    ELSE
                        numcd(ix,iy,iz)=0
                        call RANDOM_NUMBER(pptmp)
                        do pptmpi=1,ncut
                            if(pptmp.le.real(pptmpi)/real(ncut)) exit
                        end do
                        simim(ix,iy,iz)=pptmpi
                        DenSim(ix,iy,iz)=normal(EDen(simim(ix,iy,iz)),VDen(simim(ix,iy,iz)))      
                        VpSim(ix,iy,iz)=slopeVp(simim(ix,iy,iz))*DenSim(ix,iy,iz)+interceptVp(simim(ix,iy,iz))+normal(ResMeanVp(simim(ix,iy,iz)),RMSEVp(simim(ix,iy,iz)))
                         
                    end if
                    if(simim(ix,iy,iz)<1) write(*,*) '1'
                ELSE
                    itera=itera+1
                    CALL InferCpdfdz(ix,iy,iz,imult)
                    do ixr=1,nx
                        do iyr=1,ny
                            call convolution(nz,ixr,iyr,1,DenSim(ixr,iyr,1:nz),VpSim(ixr,iyr,1:nz),SnySim(ixr,iyr,1:nz))           
                        end do
                    end do
                
                    seisRMS=0.0
                    do izr=nz, 1,-1
                        do iyr=1, ny
                            do ixr=1, nx
                                seisRMS=seisRMS+ABS(SnySim(ixr,iyr,izr)-AngleSeismic(ixr,iyr,izr))
                                if(ixr == 20) then
                                    if(simim(ixr,iyr,izr)<0) then
                                        write(wellpar,*) FaciesSimTmp(ixr,iyr,izr)-1,DenSim(ixr,iyr,izr),VpSim(ixr,iyr,izr)
                                    else
                                        write(wellpar,*) simim(ixr,iyr,izr)-1,DenSim(ixr,iyr,izr),VpSim(ixr,iyr,izr)
                                    end if 
                                end if
                            end do   
                        end do
                    end do
                    
                    call CPU_TIME(time3)
                    write(RMSEpar,*) itera,seisRMS,time3,TK
                    write(recordpar,*) WCY,ix,iy,iz,simim(ix,iy,iz)
                END IF 
            END IF
        end do
    END DO

    END SUBROUTINE Simulation
    !***********************************************************************************
    

    !***********************************************************************************
    SUBROUTINE SearchClosestNodes(ix,iy,iz)
    !c----------------------------------------------------------------------
    !c
    !c
    !c
    !c----------------------------------------------------------------------
    use varmod
    implicit none
    ! Declare dummy arguments
    INTEGER, INTENT(IN) :: ix, iy,iz

    ! Declare local variables
    INTEGER :: i, j, k, ind, icut

    !
    ! First, spiral away from the node being simulated and node all
    ! the nearby nodes that have been simulated
    !
    ncnode=0

    ! Consider all the nearby nodes until enough have been found:

    cnodev(1:MAXNOD)=UNEST
    DO ind=1,nltemplate
        IF(ncnode==nodmax) EXIT
        i=ix+ixnode(ind)
        j=iy+iynode(ind)
        k=iz+iznode(ind)
     
        IF(i>=1.and.i<=nx.and.j>=1.and.j<=ny.and.k>=1.and.k<=nz) THEN
            cnodev1(ind)=simim(i,j,k)
            IF(cnodev1(ind)>UNEST) THEN
                ncnode=ncnode+1
                cnodex(ncnode) = ixnode(ind)
                cnodey(ncnode) = iynode(ind)
                cnodez(ncnode) = iznode(ind)
                cnodev(ncnode)= simim(i,j,k)
            END IF
            
        END IF
    END DO

    END SUBROUTINE SearchClosestNodes
    !***********************************************************************************


    !***********************************************************************************
    SUBROUTINE InferCpdfdz(ix,iy,iz,imult)
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    use ran_mod
    implicit none
    INTEGER, INTENT(IN) :: ix, iy, iz, imult
    INTEGER :: icd, ixtr, iytr, iztr, i, j, k,icd1
    INTEGER :: icut
    integer,dimension(MAXNOD,10000)::tmpx,tmpy,tmpz
    INTEGER, DIMENSION (MAXCTX,MAXCTY,MAXZ)::SimTmp1,simtmp2
    REAL,dimension(MAXCTX,MAXCTY,MAXZ)::DenTmp1,VpTmp1,DenTmp2,VpTmp2
    REAL,dimension(MAXZ)::snytmp1,snytmp2
    integer::nre,incy
    integer::iix,iiy,iiz,nctx,ncty,nctz,nctx1,nctx2,ncty1,ncty2,nctz1,nctz2
    integer::iiix,iiiy,iiiz
    real::pptmp
    real(8)::lktmp1,lktmp2
    integer::zzt,zzb,kk
    integer:: zzk,ixTem,iyTem,izTem,upsnn,countEventData,CountSearchTI
    real , DIMENSION  (MAXXYZTR) :: c,d,e,f,g,h
    real,dimension(MAXXYZTR)::simsj,ordertmp
    integer::nxyzti,nxyti,nxti,inct
    real(8) :: kz,NcnodeRatio1,NcnodeRatio,SearchTiRatio
    real,dimension(MAXCUT) :: DenTmpCenter

    logical::better
    better=.FALSE.
    tmpx(1:MAXNOD,1:10000)=0
    tmpy(1:MAXNOD,1:10000)=0
    tmpz(1:MAXNOD,1:10000)=0
    uni(1:MAXXYZ) = 0
    unif(1:MAXXYZ) = 0
    SimTmp1(MAXCTX,MAXCTY,MAXZ) = UNEST

    nxyzti=nxtr(imult)*nytr(imult)*nztr(imult)
    nxyti=nxtr(imult)*nytr(imult)
    nxti=nxtr(imult)
    call RANDOM_NUMBER(simsj(:nxyzti))
    do inct=1,nxyzti
        ordertmp(inct)=inct  
    end do
    CALL sortem(1,nxyzti,simsj,1,ordertmp,c,d,e,f,g,h)
    nctx1=-int(min(ix-1,abs(ixnode(nltemplate))))
    nctx2=int(min(nx-ix,abs(ixnode(nltemplate))))
    ncty1=-int(min(iy-1,abs(iynode(nltemplate))))
    ncty2=int(min(ny-iy,abs(iynode(nltemplate))))
    nctz1=-int(min(iz-1,abs(iznode(nltemplate))))
    nctz2=int(min(nz-iz,abs(iznode(nltemplate))+ncoarset-1))  
    
    do inct=1,nxyzti
        iztr=1+int((ordertmp(inct)-1)/nxyti)
        iytr=1+int((ordertmp(inct)-(iztr-1)*nxyti-1)/nxti)
        ixtr=ordertmp(inct)-(iztr-1)*nxyti-(iytr-1)*nxti
        if((ixtr+nctx1)<1.or.(ixtr+nctx2)>nxtr(imult).or.(iytr+ncty1)<1.or.(iytr+ncty2)>nytr(imult).or.(iztr+nctz1)<1.or.(iztr+nctz2)>nztr(imult)) go to 45
        countEventData=0
        CountSearchTI=0
        SearchTiRatio=0.0
        NcnodeRatio=0.0
        DO icd1=1,ncnode
            i=ixtr+cnodex(icd1)
            j=iytr+cnodey(icd1)
            k=iztr+cnodez(icd1)
            IF(trainim(i,j,k,imult)==cnodev(icd1)) THEN
                countEventData=countEventData+1
                CountSearchTI=CountSearchTI+1
            end if
            SearchTiRatio=real(CountSearchTI)/real(nxyzti)
            NcnodeRatio1=real(countEventData)/real(ncnode)
        end do
        if(NcnodeRatio1>NcnodeRatio) then
            NcnodeRatio=NcnodeRatio1
            icd=icd1
            tmpx(icd,1)=ixtr
            tmpy(icd,1)=iytr
            tmpz(icd,1)=iztr
            if(SearchTiRatio>=TIScanning.or.NcnodeRatio1>=DMathcingRate) then
                continue
                exit
            end if
        end if
45      continue
    end do
    
    IF(WCY.GT.1)THEN
        lktmp1=0.0
        do ixTem=nctx1,nctx2,ncoarset
            iiix=ix+ixTem
            i=ixTem-nctx1+1
            do iyTem=ncty1,ncty2,ncoarset
                iiiy=iy+iyTem
                j=iyTem-ncty1+1
                call convolution(nz,iiix,iiiy,1,DenSim(iiix,iiiy,1:nz),VpSim(iiix,iiiy,1:nz),snytmp1(1:nz))
                do iiiz = nz, 1,-1
                    if(SNVar==0.0) then
                        lktmp1=lktmp1+abs(snytmp1(iiiz)-AngleSeismic(iiix,iiiy,iiiz))
                    else
                        lktmp1=lktmp1+(snytmp1(iiiz)-AngleSeismic(iiix,iiiy,iiiz))**2/(2*SNVar*SNVar)
                    end if
                end do
            end do
        end do
    ELSE
        lktmp1=1000000000.0
    END IF
                                                   
    do nre=1,1 
        do incy=1,ncy
            lktmp2=0.0
            Rsny(1:nx,1:ny) = 0.0
            simtmp2(1:nx,1:ny,1:nz) = UNEST
            do i=1,ncut
                DenTmpCenter(i) = normal(EDen(i),VDen(i))
            end do
            
            do ixTem=nctx1,nctx2,ncoarset
                iix=tmpx(icd,nre)+ixTem
                iiix=ix+ixTem
                i=ixTem-nctx1+1
                do iyTem=ncty1,ncty2,ncoarset
                    iiy=tmpy(icd,nre)+iyTem
                    iiiy=iy+iyTem
                    j=iyTem-ncty1+1
                    do kz=1,nz
                        DenTmp2(i,j,kz)=DenSim(iiix,iiiy,kz)
                        VpTmp2(i,j,kz)=VpSim(iiix,iiiy,kz)
                    end do
                    do izTem=nctz1,nctz2,ncoarset
                        iiz=tmpz(icd,nre)+izTem
                        iiiz=iz+izTem
                        simtmp2(i,j,iiiz)=trainim(iix,iiy,iiz,imult)
                            
                        call RANDOM_NUMBER(pptmp)
                        DenTmp2(i,j,iiiz) = DenTmpCenter(simtmp2(i,j,iiiz))+(pptmp-0.5)*VDen(simtmp2(i,j,iiiz))
                        VpTmp2(i,j,iiiz)=slopeVp(simtmp2(i,j,iiiz))*DenTmp2(i,j,iiiz)+interceptVp(simtmp2(i,j,iiiz))+normal(ResMeanVp(simtmp2(i,j,iiiz)),RMSEVp(simtmp2(i,j,iiiz)))
                            
                        if(simim(iiix,iiiy,iiiz).gt.0) then
                            simtmp2(i,j,iiiz)=simim(iiix,iiiy,iiiz)
                                
                            DenTmp2(i,j,iiiz)=DenSim(iiix,iiiy,iiiz)
                            VpTmp2(i,j,iiiz)=VpSim(iiix,iiiy,iiiz)
                        end if
                        do upsnn=1,ncoarset-1
                            if(iiiz+upsnn<=iz+nctz2) then
                                simtmp2(i,j,iiiz+upsnn)=simtmp2(i,j,iiiz)
                                    
                                DenTmp2(i,j,iiiz+upsnn)=DenTmp2(i,j,iiiz)
                                VpTmp2(i,j,iiiz+upsnn)=VpTmp2(i,j,iiiz)
                            end if
                        end do
                    end do
                        
                    call convolution(nz,i,j,1,DenTmp2(i,j,1:nz),VpTmp2(i,j,1:nz),snytmp2(1:nz))
                    do iiiz = nz, 1,-1
                        if(SNVar==0.0) then
                            lktmp2=lktmp2+abs(snytmp2(iiiz)-AngleSeismic(iiix,iiiy,iiiz))
                        else
                            lktmp2=lktmp2+(snytmp2(iiiz)-AngleSeismic(iiix,iiiy,iiiz))**2/(2*SNVar*SNVar)
                        end if
                    end do
                end do
            end do
            
            call RANDOM_NUMBER(pptmp)
            Allcon=Allcon+1.0
            if(pptmp<exp((lktmp1-lktmp2)/TK)) then
                altcon=altcon+1.0
                BETTER=.TRUE.
                lktmp1=lktmp2
                do ixTem=nctx1,nctx2,ncoarset
                    iiix=ix+ixTem
                    i=ixTem-nctx1+1
                    do iyTem=ncty1,ncty2,ncoarset
                        iiiy=iy+iyTem
                        j=iyTem-ncty1+1
                        do izTem=nctz1,nctz2
                            iiiz=iz+izTem
                            SimTmp1(i,j,iiiz)=simtmp2(i,j,iiiz)
                            DenTmp1(i,j,iiiz)=DenTmp2(i,j,iiiz)
                            VpTmp1(i,j,iiiz)=VpTmp2(i,j,iiiz)
                        end do
                    end do
                end do
            end if
        end do
    end do
   
    if(.NOT.better)THEN
        do ixTem=nctx1,nctx2,ncoarset
            iiix=ix+ixTem
            i=ixTem-nctx1+1
            do iyTem=ncty1,ncty2,ncoarset
                iiiy=iy+iyTem
                j=iyTem-ncty1+1
                do izTem=nctz1,nctz2
                    iiiz=iz+izTem
                    if(simim(iiix,iiiy,iiiz)>0) then
                        SimTmp1(i,j,iiiz)=simim(iiix,iiiy,iiiz)
                    else
                        SimTmp1(i,j,iiiz)=FaciesSimTmp(iiix,iiiy,iiiz)
                    end if
                    DenTmp1(i,j,iiiz)=DenSim(iiix,iiiy,iiiz)
                    VpTmp1(i,j,iiiz)=VpSim(iiix,iiiy,iiiz)
                end do
            end do
        end do
    END IF
    if(assignway.eq.2)then
        do ixTem=nctx1,nctx2,ncoarset
            iiix=ix+ixTem
            i=ixTem-nctx1+1
            do iyTem=ncty1,ncty2,ncoarset
                iiiy=iy+iyTem
                j=iyTem-ncty1+1
                do izTem=nctz1,nctz2
                    iiiz=iz+izTem
                    simim(iiix,iiiy,iiiz)=SimTmp1(i,j,iiiz)
                    DenSim(iiix,iiiy,iiiz)=DenTmp1(i,j,iiiz)
                    VpSim(iiix,iiiy,iiiz)=VpTmp1(i,j,iiiz)
                end do
            end do
        end do
    else if(assignway.eq.1)then
        i=nctx+1
        j=ncty+1
        simim(ix,iy,iz)=SimTmp1(i,j,iz)
        DenSim(ix,iy,iz)=DenTmp1(i,j,iz)
        VpSim(ix,iy,iz)=VpTmp1(i,j,iz)
    end if
    
    !if(wcy>1) then
    !simim(1,1,80)=1
    !open(1,file="SIMOUTFL\faciestest.txt",status="UNKNOWN")
    !write(1,*) 'PETREL:Property'
    !write(1,*) '1'
    !write(1,*) 'Facies_1'
    !do iiz=nz,1,-1
    !    do iiy=1,ny
    !        do iix=1,nx
    !            write(1,*) simim(iix,iiy,iiz)
    !        end do
    !    end do
    !end do
    !close(1)
    !end if
    
    END SUBROUTINE InferCpdfdz
    !***********************************************************************************



    !***********************************************************************************
    subroutine Convolution(num1,xtmp,ytmp,ztmp,den,speed,rarry) 
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    implicit none
    integer::num1
    integer::i,j,k,xtmp,ytmp,ztmp
    real::dt,f0,sinang
    real,dimension(0:nz-1)::reflect
    real,dimension(MAXZ)::den,speed,speeds,EI,kk
    real,dimension(0:nz+nnd-1)::sny
    real,dimension(num1)::rarry

    sny(0:nz+nnd-1)=0.0
    reflect(0)=0.0
    EI(1)=den(1)*speed(1)
    do k=1,nz-1
        EI(k+1)=den(k+1)*speed(k+1)
        reflect(k)=(EI(k+1)-EI(k))/(EI(k+1)+EI(k))
    end do
    
    if(nnd.gt.nz)then
        do i=0,nz+nnd-1
            if(i.lt.nz)then
                sny(i)=0.0
                do j=0,i
                    sny(i)=sny(i)+wavelet(j)*reflect(i-j)
                end do
            else if(i.ge.nz.and.i.lt.nnd)then
                sny(i)=0.0
                do j=i-nz+1,i
                    sny(i)=sny(i)+wavelet(j)*reflect(i-j)
                end do
            else if(i.ge.nnd)then
                sny(i)=0.0
                do j=i-nz+1,nnd
                    sny(i)=sny(i)+wavelet(j)*reflect(i-j)
                end do
            end if
        end do
    else
        do i=0,nz+nnd-1
            if(i.lt.nnd)then
                sny(i)=0.0
                do j=0,i
                    sny(i)=sny(i)+wavelet(j)*reflect(i-j)
                end do
            else if(i.ge.nnd.and.i.lt.nz-1)then
                sny(i)=0.0
                do j=0,nnd
                    sny(i)=sny(i)+wavelet(j)*reflect(i-j)
                end do
            else if(i.ge.nz-1)then
                sny(i)=0.0
                do j=i-nz+1,nnd 
                    sny(i)=sny(i)+wavelet(j)*reflect(i-j)
                end do
            end if
        end do
    end if

    do i=1,num1
        rarry(i)=sny(i+33+ztmp)
    end do
    sny(0:num1+nnd-1)=0.0
    end subroutine Convolution

    
    
    SUBROUTINE RelocateData()
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    IMPLICIT NONE

    ! Declare local variables
    INTEGER :: ix, iy, iz, id, id2
    REAL :: xx, yy, zz, test, test2
    LOGICAL :: testind
    real::pptmp
    integer::i

    !
    ! Loop over all sample data:
    !
    DO id=1,nd
        !
        ! Calculate the coordinates of the closest simulation grid node:
        !
        CALL getindx(nx,xmn,xsiz,x(id),ix,testind)
        CALL getindx(ny,ymn,ysiz,y(id),iy,testind)
        CALL getindx(nz,zmn,zsiz,z(id),iz,testind)
        xx=xmn+real(ix-1)*xsiz
        yy=ymn+real(iy-1)*ysiz
        zz=zmn+real(iz-1)*zsiz
        test=abs(xx-x(id))+abs(yy-y(id))+abs(zz-z(id))
        !
        ! Assign this data to the node (unless there is a closer data):
        !
        IF(simim(ix,iy,iz)>0) THEN
            id2 = simim(ix,iy,iz)
            test2 = abs(xx-x(id2))+abs(yy-y(id2))+abs(zz-z(id2))
            IF(test<test2) simim(ix,iy,iz)=id
        ELSE
            simim(ix,iy,iz)=id
        END IF
    END DO

102 format('Warning data values ',2i5,' are both assigned to ',/,'the same node - taking the closest')

    !
    ! Now, enter data values into the simulation grid:
    !

    DO iz=1,nz
        DO iy=1,ny
            DO ix=1,nx
                id=simim(ix,iy,iz)
                IF(id>0) THEN
                    simim(ix,iy,iz) = vr(id)
                    DenSim(ix,iy,iz) = vr_d(id)
                    VpSim(ix,iy,iz)=vr_s(id)
                END IF
            END DO
        END DO
    END DO
    
    stop
    
    END SUBROUTINE RelocateData

    !***********************************************************************************

    !***********************************************************************************
    SUBROUTINE SortTemplate (imult)
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    implicit none

    ! Declare dummy arguments
    INTEGER, INTENT(IN) :: imult

    ! Declare local variables
    INTEGER :: ind,n
    REAL :: xx,yy,zz,hsqd,cont
    REAL, DIMENSION(MAXNOD) :: tmp
    real, DIMENSION(MAXNOD) :: e,f,g,h

    CALL SetRotMat(imult)

    DO ind=1,nltemplate
        xx=ixtemplate(ind)
        yy=iytemplate(ind)
        zz=iztemplate(ind)
        hsqd = 0.0
        DO n=1,3
            cont=rotmat(n,1)*xx+rotmat(n,2)*yy+rotmat(n,3)*zz
            hsqd = hsqd + cont*cont
        END DO
        tmp(ind) = hsqd
    END DO
    
    call sortem(1,nltemplate,tmp,3,ixtemplate,iytemplate,iztemplate,e,f,g,h)
    
    end subroutine SortTemplate
    !***********************************************************************************


    !***********************************************************************************
    SUBROUTINE SetRotMat(imult)
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    implicit none
    ! Declare dummy arguments
    INTEGER, INTENT(IN) :: imult

    ! Declare local variables
    REAL :: afac1,afac2,sina,sinb,sint,cosa,cosb,cost
    REAL :: alpha,beta,theta
    !
    ! Converts the input angles to three angles which make more
    ! mathematical sense:
    !
    !         alpha   angle between the major axis of anisotropy and the
    !                 E-W axis. Note: Counter clockwise is positive.
    !         beta    angle between major axis and the horizontal plane.
    !                 (The dip of the ellipsoid measured positive down)
    !         theta   Angle of rotation of minor axis about the major axis
    !                 of the ellipsoid.
    !
    IF(sang1(imult)>=0.0.and.sang1(imult)<270.0) then
        alpha=(90.0-sang1(imult))*DEG2RAD
    ELSE
        alpha=(450.0-sang1(imult))*DEG2RAD
    ENDIF
    beta=-1.0*sang2(imult)*DEG2RAD
    theta=sang3(imult)*DEG2RAD
    !
    ! Get the required sines and cosines:
    !
    sina  = sin(alpha)
    sinb  = sin(beta)
    sint  = sin(theta)
    cosa  = cos(alpha)
    cosb  = cos(beta)
    cost  = cos(theta)
    !
    ! Construct the rotation matrix in the required memory:
    !
    afac1 = 1.0 / max(sanis1(imult),EPSILON)
    afac2 = 1.0 / max(sanis2(imult),EPSILON)
    rotmat(1,1) = cosb * cosa
    rotmat(1,2) = cosb * sina
    rotmat(1,3) = -sinb
    rotmat(2,1) = afac1*(-cost*sina + sint*sinb*cosa)
    rotmat(2,2) = afac1*(cost*cosa + sint*sinb*sina)
    rotmat(2,3) = afac1*( sint * cosb)
    rotmat(3,1) = afac2*(sint*sina + cost*sinb*cosa)
    rotmat(3,2) = afac2*(-sint*cosa + cost*sinb*sina)
    rotmat(3,3) = afac2*(cost * cosb)

    END SUBROUTINE SetRotMat
    !***********************************************************************************




    !***********************************************************************************
    SUBROUTINE AssignData()
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    implicit none
    ! Declare local variables
    INTEGER :: ix, iy, iz, id, id2
    INTEGER :: jx, jy, jz
    REAL :: xx, yy, zz, test, test2
    REAL :: xmncoarse, ymncoarse, xsizcoarse, ysizcoarse
    REAL :: zmncoarse, zsizcoarse
    LOGICAL :: testind
    INTEGER, DIMENSION (MAXX,MAXY,MAXZ) :: simtemp

    !
    ! Define specifications of the current multiple grid.
    !
    xmncoarse=xmn*ncoarse
    ymncoarse=ymn*ncoarse
    zmncoarse=zmn*ncoarse
    xsizcoarse=xsiz*ncoarse
    ysizcoarse=ysiz*ncoarse
    zsizcoarse=zsiz*ncoarse
    simtemp(1:nxcoarse,1:nycoarse,1:nzcoarse)=UNEST

    !
    ! Loop over all the original sample data
    !
    DO id=1,nd
        !
        ! Calculate the coordinates of the closest simulation grid node:
        !
        x(id) = x(id)-xmncoarse
        y(id) = y(id)-ymncoarse
        z(id) = z(id)-zmncoarse
        CALL getindx(nxcoarse,xmncoarse,xsizcoarse,x(id),ix,testind)
        CALL getindx(nycoarse,ymncoarse,ysizcoarse,y(id),iy,testind)
        CALL getindx(nzcoarse,zmncoarse,zsizcoarse,z(id),iz,testind)
        xx  = xmncoarse + real(ix-1)*xsizcoarse
        yy  = ymncoarse + real(iy-1)*ysizcoarse
        zz  = zmncoarse + real(iz-1)*zsizcoarse
        test = abs(xx-x(id)) + abs(yy-y(id)) + abs(zz-z(id))
        !
        ! Assign this data to the node (unless there is a closer data):
        !
        IF(simtemp(ix,iy,iz)>0) THEN
            id2 = simtemp(ix,iy,iz)
            test2 = abs(xx-x(id2)) + abs(yy-y(id2)) + abs(zz-z(id2))
            IF(test<test2) simtemp(ix,iy,iz)=id
        ELSE
            simtemp(ix,iy,iz)=id
        END IF
    END DO
    !
    ! Now, enter data values into the simulated grid:
    !
    DO iz=1,nzcoarse
        DO iy=1,nycoarse
            DO ix=1,nxcoarse
                id=simtemp(ix,iy,iz)
                IF(id>0) THEN
                    jz=(iz-1)*ncoarse+1
                    jy=(iy-1)*ncoarse+1
                    jx=(ix-1)*ncoarse+1

                    ! Check if there is already a simulated value; if yes, replace it.            
                    simim(jx,jy,jz) = vr(id)     
                    DenSim(jx,jy,jz) = vr_d(id)
                    VpSim(jx,jy,jz)=vr_s(id)
                    
                    ! Indicates with a special value assigned to numcd that a sample data
                    ! has been assigned to the node.
                    numcd(jx,jy,jz)=10*UNEST
                END IF
            END DO
        END DO
    END DO
    
    write(*,*) 'check out!'

    END SUBROUTINE AssignData
    !***********************************************************************************



    !***********************************************************************************
    SUBROUTINE UnassignData()
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    implicit none
    ! Declare local variables
    INTEGER :: ix,iy,iz,jx,jy,jz

    ! Loop over all the nodes of the current simulation grid:

    DO jz=1, nzcoarse
        iz=(jz-1)*ncoarse+1
        DO jy=1, nycoarse
            iy=(jy-1)*ncoarse+1
            DO jx=1, nxcoarse
                ix=(jx-1)*ncoarse+1

                ! Check if an original sample data has been assigned to the current node.
                ! If yes, remove it.

                IF(numcd(ix,iy,iz)==UNEST*10) THEN
                    simim(ix,iy,iz)=UNEST
                    numcd(ix,iy,iz)=UNEST
                END IF
            END DO
        END DO
    END DO

    END SUBROUTINE UnassignData
    !***********************************************************************************
    
    !***********************************************************************************
    SUBROUTINE upscale(ix,iy,iz)
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    implicit none
    integer :: upix,upiy,upiz,Fineix,Fineiy,Fineiz
    INTEGER, INTENT(IN) :: ix,iy,iz

    upix = (ix-1)/ncoarse+1
    upiy = (iy-1)/ncoarse+1
    upiz = (iz-1)/ncoarse+1
    
    do Fineix = max((upix-1)*ncoarse+1,1),min(upix*ncoarse,nx)
        do Fineiy = max((upiy-1)*ncoarse+1,1),min(upiy*ncoarse,ny)
            do Fineiz = max((upiz-1)*ncoarse+1,1),min(upiz*ncoarse,nz)
                
                if(simim(Fineix,Fineiy,Fineiz)>0) then
                    simim(ix,iy,iz) = simim(Fineix,Fineiy,Fineiz)
                    DenSim(ix,iy,iz) = DenSim(Fineix,Fineiy,Fineiz)
                    VpSim(ix,iy,iz) = VpSim(Fineix,Fineiy,Fineiz)                     
                    SnySim(ix,iy,iz) = SnySim(Fineix,Fineiy,Fineiz)
                    go to 56
                end if
            end do
        end do
    end do 
    
56  continue
    
    END SUBROUTINE upscale
    
    
    !***********************************************************************************
    SUBROUTINE AssignData1()
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    implicit none
    ! Declare local variables
    INTEGER :: ix, iy, iz, id, id2
    INTEGER :: jx, jy, jz
    REAL :: xx, yy, zz, test, test2
    REAL :: xmncoarse, ymncoarse, xsizcoarse, ysizcoarse
    REAL :: zmncoarse, zsizcoarse
    LOGICAL :: testind
    INTEGER, DIMENSION (MAXX,MAXY,MAXZ) :: simtemp

    !
    ! Define specifications of the current multiple grid.
    !
    xsizcoarse=xsiz*ncoarse
    ysizcoarse=ysiz*ncoarse
    zsizcoarse=zsiz*ncoarse
    simtemp(1:nxcoarse,1:nycoarse,1:nzcoarse)=UNEST

    DO id=1,nd
        jx=int((x(id)-xmn)/xsizcoarse)+1
        jx=jx*xsizcoarse
        if(jx>nx) jx=nx
        jy=int((y(id)-ymn)/ysizcoarse)+1
        jy=jy*ysizcoarse
        if(jy>ny) jy=ny
        jz=int((z(id)-zmn)/zsizcoarse)+1
        jz=jz*zsizcoarse
        if(jz>nz) jz=nz
        simim(jx,jy,jz) = vr(id)
        DenSim(jx,jy,jz) = vr_d(id)
        VpSim(jx,jy,jz)=vr_s(id)
        numcd(jx,jy,jz)=10*UNEST
    END DO
    
    do ix=1,nx
        do iy=1,ny
            do iz=1,nz
                Call upscale(ix,iy,iz)
            end do
        end do
    end do

    END SUBROUTINE AssignData1
    !***********************************************************************************