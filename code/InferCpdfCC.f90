  !***********************************************************************************
    SUBROUTINE InferCpdfCC(ix,iy,iz,imult)
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    use ran_mod
    implicit none
    INTEGER, INTENT(IN) :: ix, iy, iz, imult
    INTEGER :: ixtr, iytr, iztr, i, j, k,icd
    INTEGER :: icut
    integer::tmpx,tmpy,tmpz
    INTEGER, DIMENSION (MAXCTX,MAXCTY,MAXZ)::SimTmp1,simtmp2
    REAL,dimension(MAXCTX,MAXCTY,MAXZ)::DenTmp1,VpTmp1,DenTmp2,VpTmp2
    REAL,dimension(MAXZ)::snytmp1,snytmp2
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
    real, DIMENSION (MAXX,MAXY) :: ESnySim1,CovSny1,VarSnySim1,ESnySim2,CovSny2, VarSnySim2

    logical::better
    better=.FALSE.
    tmpx=0
    tmpy=0
    tmpz=0
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
        DO icd=1,ncnode
            i=ixtr+cnodex(icd)
            j=iytr+cnodey(icd)
            k=iztr+cnodez(icd)
            IF(trainim(i,j,k,imult)==cnodev(icd)) THEN
                countEventData=countEventData+1
                CountSearchTI=CountSearchTI+1
            end if
            SearchTiRatio=real(CountSearchTI)/real(nxyzti)
            NcnodeRatio1=real(countEventData)/real(ncnode)
        end do
        if(NcnodeRatio1>NcnodeRatio) then
            NcnodeRatio=NcnodeRatio1
            tmpx=ixtr
            tmpy=iytr
            tmpz=iztr
            if(SearchTiRatio>=TIScanning.or.NcnodeRatio1>=DMathcingRate) then
                continue
                exit
            end if
        end if
45      continue
    end do
    
    IF(WCY.GT.1)THEN
        lktmp1=0.0
        ESnySim1(1:nx,1:ny) = 0.0
        CovSny1(1:nx,1:ny) = 0.0
        VarSnySim1(1:nx,1:ny) = 0.0
        do ixTem=nctx1,nctx2,ncoarset
            iiix=ix+ixTem
            do iyTem=ncty1,ncty2,ncoarset
                iiiy=iy+iyTem
                
                ! correlation coefficient
                call convolution(nz,iiix,iiiy,1,DenSim(iiix,iiiy,1:nz),VpSim(iiix,iiiy,1:nz),snytmp1(1:nz))
                do iiiz = nz, 1,-1
                    ESnySim1(iiix,iiiy)=ESnySim1(iiix,iiiy)+snytmp1(iiiz)
                    CovSny1(iiix,iiiy)=CovSny1(iiix,iiiy)+snytmp1(iiiz)*ActualSeismic(iiix,iiiy,iiiz)
                    
                end do                   
                ESnySim1(iiix,iiiy)=real(ESnySim1(iiix,iiiy)/nz)
                CovSny1(iiix,iiiy)=real(CovSny1(iiix,iiiy)/nz)
                do iiiz=1,nz
                    VarSnySim1(iiix,iiiy)=VarSnySim1(iiix,iiiy)+(snytmp1(iiiz)-ESnySim1(iiix,iiiy))**2
                end do
                
                VarSnySim1(iiix,iiiy)=real(VarSnySim1(iiix,iiiy)/nz)
                lktmp1=lktmp1+(CovSny1(iiix,iiiy)-ESnySim1(iiix,iiiy)*EActualSeismic(iiix,iiiy))/((VarSnySim1(iiix,iiiy)*VarActualSeismic(iiix,iiiy))**0.5)          
            end do
        end do
    ELSE
        lktmp1=-1000000000.0
    END IF
                                                   
    lktmp2=0.0
    ESnySim2(1:nx,1:ny) = 0.0
    CovSny2(1:nx,1:ny) = 0.0
    VarSnySim2(1:nx,1:ny) = 0.0
    Rsny(1:nx,1:ny) = 0.0
    simtmp2(1:nx,1:ny,1:nz) = UNEST
    do i=1,ncut
        DenTmpCenter(i) = normal(EDen(i),VDen(i))
    end do
            
    do ixTem=nctx1,nctx2,ncoarset
        iix=tmpx+ixTem
        iiix=ix+ixTem
        i=ixTem-nctx1+1
        do iyTem=ncty1,ncty2,ncoarset
            iiy=tmpy+iyTem
            iiiy=iy+iyTem
            j=iyTem-ncty1+1
                    
            do kz=1,nz
                DenTmp2(i,j,kz)=DenSim(iiix,iiiy,kz)
                VpTmp2(i,j,kz)=VpSim(iiix,iiiy,kz)
            end do
                    
            do izTem=nctz1,nctz2,ncoarset
                iiz=tmpz+izTem
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
                    
            ! correlation coefficient
            call convolution(nz,i,j,1,DenTmp2(i,j,1:nz),VpTmp2(i,j,1:nz),snytmp2(1:nz))
            do iiiz = nz, 1,-1
                ESnySim2(i,j)=ESnySim2(i,j)+snytmp2(iiiz)
                CovSny2(i,j)=CovSny2(i,j)+snytmp2(iiiz)*ActualSeismic(iiix,iiiy,iiiz)
            end do                   
                    
            ESnySim2(i,j)=real(ESnySim2(i,j)/nz)
            CovSny2(i,j)=real(CovSny2(i,j)/nz)
            do iiiz=1,nz
                VarSnySim2(i,j)=VarSnySim2(i,j)+(snytmp2(iiiz)-ESnySim2(i,j))**2
            end do
            VarSnySim2(i,j)=real(VarSnySim2(i,j)/nz)
            lktmp2=lktmp2+(CovSny2(i,j)-ESnySim2(i,j)*EActualSeismic(iiix,iiiy))/((VarSnySim2(i,j)*VarActualSeismic(iiix,iiiy))**0.5)
        end do
    end do
            
    call RANDOM_NUMBER(pptmp)
    Allcon=Allcon+1.0
    if(pptmp<exp((lktmp2-lktmp1)/Temperature)) then
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
    
    END SUBROUTINE InferCpdfCC
    !***********************************************************************************