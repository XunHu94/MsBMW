  !***********************************************************************************
    SUBROUTINE SetRotMat(imult)
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    implicit none
    ! Declare dummy arguments
    INTEGER, INTENT(IN) :: imult

    ! Declare local variables
    REAL :: afac1,afac2,sina,sinb,sint,cosa,cosb,cost
    REAL :: alpha,beta,theta
    !
    ! Converts the input angles to three angles which make more
    ! mathematical sense:
    !
    !         alpha   angle between the major axis of anisotropy and the
    !                 E-W axis. Note: Counter clockwise is positive.
    !         beta    angle between major axis and the horizontal plane.
    !                 (The dip of the ellipsoid measured positive down)
    !         theta   Angle of rotation of minor axis about the major axis
    !                 of the ellipsoid.
    !
    IF(sang1(imult)>=0.0.and.sang1(imult)<270.0) then
        alpha=(90.0-sang1(imult))*DEG2RAD
    ELSE
        alpha=(450.0-sang1(imult))*DEG2RAD
    ENDIF
    beta=-1.0*sang2(imult)*DEG2RAD
    theta=sang3(imult)*DEG2RAD
    !
    ! Get the required sines and cosines:
    !
    sina  = sin(alpha)
    sinb  = sin(beta)
    sint  = sin(theta)
    cosa  = cos(alpha)
    cosb  = cos(beta)
    cost  = cos(theta)
    !
    ! Construct the rotation matrix in the required memory:
    !
    afac1 = 1.0 / max(sanis1(imult),EPSILON)
    afac2 = 1.0 / max(sanis2(imult),EPSILON)
    rotmat(1,1) = cosb * cosa
    rotmat(1,2) = cosb * sina
    rotmat(1,3) = -sinb
    rotmat(2,1) = afac1*(-cost*sina + sint*sinb*cosa)
    rotmat(2,2) = afac1*(cost*cosa + sint*sinb*sina)
    rotmat(2,3) = afac1*( sint * cosb)
    rotmat(3,1) = afac2*(sint*sina + cost*sinb*cosa)
    rotmat(3,2) = afac2*(-sint*cosa + cost*sinb*sina)
    rotmat(3,3) = afac2*(cost * cosb)

    END SUBROUTINE SetRotMat
    !***********************************************************************************