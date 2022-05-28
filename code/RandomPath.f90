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