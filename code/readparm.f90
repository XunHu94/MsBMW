!***********************************************************************************
    subroutine readparm
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    implicit none
    real,dimension(10)::var,var_d,var_s,var_ss
    INTEGER, DIMENSION(2) :: seed
    character::ConFaciesFile*40,ConDenFile*40,ConVpFile*40,str*40,TemplateFile*40,SeismicWaveletFile*40,SeismicFile*40,InputFoderName*40,MsBMWPar*40,RockPhysicsModelFile*40
    CHARACTER (LEN=20), DIMENSION(MAXMULT) :: trainfl
    INTEGER, DIMENSION (MAXMULT) :: ivrltr
    logical::testfl
    integer::ioerror
    integer::ixl,iyl,izl,ivrl,icut
    integer::ix,iy,iz,imult,ntr
    integer::i,j,k,nvari
    integer::ixv,Facie,TatolFacie
    
    InputFoderName='RealData/'
    MsBMWPar='MsBMW.par'
    
    WRITE(*,999) VERSION
999 format(/' MsBMW Version: ',f8.3/)

    MsBMWPar=trim(adjustl(InputFoderName))// trim(adjustl(MsBMWPar))
    WRITE(*,*) 'Which parameter file do you want to use?'
    INQUIRE(file=MsBMWPar,exist=testfl)
    OPEN(lin,file=MsBMWPar,status='OLD')

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

    READ(lin,'(a20)',err=198) SeismicWaveletFile
    CALL chknam(SeismicWaveletFile,40)
    WRITE(*,*) ' data file = ',SeismicWaveletFile

    READ(lin,'(a20)',err=198) SeismicFile
    CALL chknam(SeismicFile,40)
    WRITE(*,*) ' data file = ',SeismicFile
    
    READ(lin,'(a20)',err=198) RockPhysicsModelFile
    CALL chknam(RockPhysicsModelFile,40)
    WRITE(*,*) ' data file = ',RockPhysicsModelFile
    
    READ(lin,*,err=198) NumChains
    WRITE(*,*) ' assignway = ',NumChains

    READ(lin,*,err=198) ObjectiveFuction
    WRITE(*,*) ' assignway = ',ObjectiveFuction
    
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
    
    READ(lin,*,err=198) MaxMultiGridIterations
    WRITE(*,*) ' Maxiterations = ',MaxMultiGridIterations
    
    READ(lin,*,err=198) GridCriteria
    WRITE(*,*) ' The grid degradation criteria = ',GridCriteria
    
    READ(lin,*,err=198) TIScanning
    WRITE(*,*) ' the search scope of TI = ',TIScanning
    
    READ(lin,*,err=198) DMathcingRate
    WRITE(*,*) ' Data events match rate = ',DMathcingRate
    
    READ(lin,*,err=198) InitialTemperature
    WRITE(*,*) ' Initial temperature = ',InitialTemperature
    
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

    ConFaciesFile = trim(adjustl(InputFoderName))// trim(adjustl(ConFaciesFile))
    ConDenFile = trim(adjustl(InputFoderName))// trim(adjustl(ConDenFile))
    ConVpFile = trim(adjustl(InputFoderName))// trim(adjustl(ConVpFile))
    SeismicWaveletFile = trim(adjustl(InputFoderName))// trim(adjustl(SeismicWaveletFile))
    SeismicFile = trim(adjustl(InputFoderName))// trim(adjustl(SeismicFile))
    TemplateFile = trim(adjustl(InputFoderName))// trim(adjustl(TemplateFile))
    RockPhysicsModelFile = trim(adjustl(InputFoderName))// trim(adjustl(RockPhysicsModelFile))
    
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
        trainfl(imult) = trim(adjustl(InputFoderName))// trim(adjustl(trainfl(imult)))
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
    
    INQUIRE(file=SeismicWaveletFile,exist=testfl)
    OPEN(lin,file=SeismicWaveletFile,status='OLD')
    IF(.NOT.testfl) THEN
        WRITE(*,*) 'WARNING data file ',SeismicWaveletFile,' does not exist!'
    ELSE
        WRITE(*,*) 'Reading input data'
        NumSampledWalve=0
        DO
            READ(lin,*,IOSTAT=ioerror) wavelet(NumSampledWalve)
            IF(ioerror<0) EXIT
            NumSampledWalve=NumSampledWalve+1
        end do
        CLOSE(lin)
        NumSampledWalve=NumSampledWalve-1
    END IF
    HalfNumSampledWalve=(NumSampledWalve+1)/2

    INQUIRE(file=SeismicFile,exist=testfl)
    IF(.NOT.testfl) THEN
        STOP 'Error in Seismic data file'
    ELSE
        OPEN(lin,file=SeismicFile,status='OLD')
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
            ActualSeismic(ix,iy,nz-iz+1)=var(nvari)
        END DO
        CLOSE(lin)
    end if
    
    INQUIRE(file=RockPhysicsModelFile,exist=testfl)
    OPEN(1,file=RockPhysicsModelFile,status='OLD')
    
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
    WRITE(*,*) ' Noises = ',SNVar
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