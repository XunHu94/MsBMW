!***********************************************************************************
    SUBROUTINE SearchClosestNodes(ix,iy,iz)
    !c----------------------------------------------------------------------
    !c
    !c
    !c
    !c----------------------------------------------------------------------
    use varmod
    implicit none
    ! Declare dummy arguments
    INTEGER, INTENT(IN) :: ix, iy,iz

    ! Declare local variables
    INTEGER :: i, j, k, ind, icut

    !
    ! First, spiral away from the node being simulated and node all
    ! the nearby nodes that have been simulated
    !
    ncnode=0

    ! Consider all the nearby nodes until enough have been found:

    cnodev(1:MAXNOD)=UNEST
    DO ind=1,nltemplate
        IF(ncnode==nodmax) EXIT
        i=ix+ixnode(ind)
        j=iy+iynode(ind)
        k=iz+iznode(ind)
     
        IF(i>=1.and.i<=nx.and.j>=1.and.j<=ny.and.k>=1.and.k<=nz) THEN
            cnodev1(ind)=simim(i,j,k)
            IF(cnodev1(ind)>UNEST) THEN
                ncnode=ncnode+1
                cnodex(ncnode) = ixnode(ind)
                cnodey(ncnode) = iynode(ind)
                cnodez(ncnode) = iznode(ind)
                cnodev(ncnode)= simim(i,j,k)
            END IF
            
        END IF
    END DO

    END SUBROUTINE SearchClosestNodes
    !***********************************************************************************