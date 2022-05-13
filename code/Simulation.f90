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