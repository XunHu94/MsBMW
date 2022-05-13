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