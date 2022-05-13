 !***********************************************************************************
    SUBROUTINE SortTemplate (imult)
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    implicit none

    ! Declare dummy arguments
    INTEGER, INTENT(IN) :: imult

    ! Declare local variables
    INTEGER :: ind,n
    REAL :: xx,yy,zz,hsqd,cont
    REAL, DIMENSION(MAXNOD) :: tmp
    real, DIMENSION(MAXNOD) :: e,f,g,h

    CALL SetRotMat(imult)

    DO ind=1,nltemplate
        xx=ixtemplate(ind)
        yy=iytemplate(ind)
        zz=iztemplate(ind)
        hsqd = 0.0
        DO n=1,3
            cont=rotmat(n,1)*xx+rotmat(n,2)*yy+rotmat(n,3)*zz
            hsqd = hsqd + cont*cont
        END DO
        tmp(ind) = hsqd
    END DO
    
    call sortem(1,nltemplate,tmp,3,ixtemplate,iytemplate,iztemplate,e,f,g,h)
    
    end subroutine SortTemplate
    !***********************************************************************************