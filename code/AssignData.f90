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

    END SUBROUTINE AssignData
    !***********************************************************************************