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