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
    INTEGER :: ingd,ic,ix,iy,iz,ixr,iyr,izr
    INTEGER :: jx, jy, jz
    real::pptmp
    integer::pptmpi
    real,dimension(nz)::rangz,rangtm
    real , DIMENSION  (MAXXYZ) :: c,d,e,f,g,h
    integer::i,NcoarseNz
    integer :: kz
    real :: TsnyccBMcMCIte,AverageCCBMcMCIte
    
    ! Main loop over all simulation grid nodes:
    DO ingd=1, nxycoarse
        ! Figure out the location of this point and make sure it has
        ! not been assigned a value already:
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
                ELSE
                    if (ObjectiveFuction==1) then
                        CALL InferCpdfRMSE(ix,iy,iz,imult)
                        
!                        do ixr=1,nx
!                            do iyr=1,ny
!                                call convolution(nz,ixr,iyr,1,DenSim(ixr,iyr,1:nz),VpSim(ixr,iyr,1:nz),SnySim(ixr,iyr,1:nz))
!                            end do
!                        end do
!                        seisRMS=0.0
!                        do izr=nz, 1,-1
!                            do iyr=1, ny
!                                do ixr=1, nx
!                                    seisRMS=seisRMS+ABS(SnySim(ixr,iyr,izr)-ActualSeismic(ixr,iyr,izr))
!                                end do
!                            end do
!                        end do
!                        write(ObFunpar,1004) itera,seisRMS
!1004                    format(I14,F16.3) 
                    else
                        itera=itera+1
                        CALL InferCpdfCC(ix,iy,iz,imult)
                        
!                        ESnySimBMcMCIte(1:nx,1:ny) = 0.0
!                        CovSnyBMcMCIte(1:nx,1:ny) = 0.0
!                        VarSnySimBMcMCIte(1:nx,1:ny) = 0.0
!                        RsnyBMcMCIte(1:nx,1:ny) = 0.0
!                        AverageCCBMcMCIte=0.0
!                        do ixr=1,nx
!                            do iyr=1,ny
!                                call convolution(nz,ixr,iyr,1,DenSim(ixr,iyr,1:nz),VpSim(ixr,iyr,1:nz),SnySim(ixr,iyr,1:nz))
!                            end do
!                        end do
!                        do izr = nz, 1,-1
!                            do iyr = 1, ny
!                                do ixr = 1, nx                  
!                                    write(2, *) SnySim(ixr,iyr,izr)
!                                    TsnyccBMcMCIte = TsnyccBMcMCIte+ABS(SnySim(ixr,iyr,izr)-ActualSeismic(ixr,iyr,izr))
!                                end do
!                            end do
!                        end do
!                        do ixr = 1,nx
!                            do iyr = 1,ny
!                                do izr = nz, 1,-1
!                                    ESnySimBMcMCIte(ixr,iyr)=ESnySimBMcMCIte(ixr,iyr)+SnySim(ixr,iyr,izr)
!                                    CovSnyBMcMCIte(ixr,iyr)=CovSnyBMcMCIte(ixr,iyr)+SnySim(ixr,iyr,izr)*ActualSeismic(ixr,iyr,izr)
!                                end do                   
!                                ESnySimBMcMCIte(ixr,iyr)=real(ESnySimBMcMCIte(ixr,iyr)/nz)
!                                CovSnyBMcMCIte(ixr,iyr)=real(CovSnyBMcMCIte(ixr,iyr)/nz)
!                                do izr=1,nz
!                                    VarSnySimBMcMCIte(ixr,iyr)=VarSnySimBMcMCIte(ixr,iyr)+(SnySim(ixr,iyr,izr)-ESnySimBMcMCIte(ixr,iyr))**2
!                                end do
!                                VarSnySimBMcMCIte(ixr,iyr)=real(VarSnySimBMcMCIte(ixr,iyr)/nz)
!                                if(VarSnySimBMcMCIte(ixr,iyr)/=0.0) RsnyBMcMCIte(ixr,iyr)=(CovSnyBMcMCIte(ixr,iyr)-ESnySimBMcMCIte(ixr,iyr)*EActualSeismic(ixr,iyr))/((VarSnySimBMcMCIte(ixr,iyr)*VarActualSeismic(ixr,iyr))**0.5)
!                                AverageCCBMcMCIte = AverageCCBMcMCIte+RsnyBMcMCIte(ixr,iyr)
!                            end do
!                        end do
!                        AverageCCBMcMCIte = AverageCCBMcMCIte/nx/ny
!                        write(ObFunpar,1005) itera,AverageCCBMcMCIte
!1005                    format(I14,F16.3) 
                    end if
                END IF 
            END IF
        end do
    END DO

    END SUBROUTINE Simulation
    !***********************************************************************************