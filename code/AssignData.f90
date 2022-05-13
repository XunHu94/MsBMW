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