PROGRAM wave2D

!*****************************************!
!* 2d shallow-Water model                *!
!*                                       *!
!* including:                            *!
!* - horizontal pressure-gradient force  *!
!* - Shapiro filter                      *!
!* - flooding algorithm                  *!
!*                                       *!
!* Author: J. Kaempf, 2008               *!
!*****************************************!

    USE param
    USE sub
    USE io

    IMPLICIT NONE

! local parameters
    REAL :: time
    INTEGER :: n, ntot, nout
    CHARACTER(len=50) :: outfile

!**********
    CALL INIT  ! initialisation
!**********

! set local parameters

! set epsilon for Shapiro filter
    eps = 0.05

! runtime parameters
    ntot = 1000

! output parameter
    nout = 10

    ! create output file
    outfile = "output.nc"
    CALL create_nc(outfile)

    ! write out initial conditions
    CALL write_nc(outfile, 0.)
!---------------------------
! simulation loop
!---------------------------
    DO n = 1, ntot

        time = REAL(n)*dt

! forcing
        DO j = 0, ny + 1
            eta(j, 1) = 0.2*sin(time*2.0*PI/10.)
        END DO

! call predictor
        CALL dyn

! updating including Shapiro filter

        CALL shapiro

        DO j = 0, ny + 1
        DO k = 0, nx + 1
            h(j, k) = hzero(j, k) + eta(j, k)
            wet(j, k) = 1
            IF (h(j, k) < hmin) wet(j, k) = 0
            u(j, k) = un(j, k)
            v(j, k) = vn(j, k)
        END DO
        END DO

! data output
        IF (MOD(n, nout) == 0) THEN
            CALL write_nc(outfile, 0.)
            WRITE (6, *) "Data output at time = ", time
        END IF

    END DO ! end of iteration loop

END PROGRAM wave2D
