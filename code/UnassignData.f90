 !***********************************************************************************
    SUBROUTINE UnassignData()
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    implicit none
    ! Declare local variables
    INTEGER :: ix,iy,iz,jx,jy,jz

    ! Loop over all the nodes of the current simulation grid:

    DO jz=1, nzcoarse
        iz=(jz-1)*ncoarse+1
        DO jy=1, nycoarse
            iy=(jy-1)*ncoarse+1
            DO jx=1, nxcoarse
                ix=(jx-1)*ncoarse+1

                ! Check if an original sample data has been assigned to the current node.
                ! If yes, remove it.

                IF(numcd(ix,iy,iz)==UNEST*10) THEN
                    simim(ix,iy,iz)=UNEST
                    numcd(ix,iy,iz)=UNEST
                END IF
            END DO
        END DO
    END DO

    END SUBROUTINE UnassignData
    !***********************************************************************************