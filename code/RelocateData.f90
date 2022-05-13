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