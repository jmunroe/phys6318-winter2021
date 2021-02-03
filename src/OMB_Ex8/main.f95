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

! local parameters
    REAL :: hmax, time, dtmax
    REAL :: c, lambda
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
    nout = 5

! determine maximum water depth
    hmax = 0.

    DO j = 1, ny
    DO k = 1, nx
        hmax = MAX(hmax, h(j, k))
    END DO
    END DO

! maximum phase speed
    c = SQRT(2*g*hmax)

! determine stability parameter
    lambda = dt*SQRT(g*hmax)/MIN(dx, dy)

    IF (lambda > 1) THEN
        WRITE (6, *) "This will not work. Do you know why?"
        STOP
    END IF

    ! create output file
    outfile = 'output.nc'
    CALL create_nc(outfile)
    
    DO j = 26, 26
    DO k = 26, 26
        eta(j, k) = 1.0
    END DO
    END DO

    ! write out initial conditions
    CALL write_nc(outfile, 0.)

!---------------------------
! simulation loop
!---------------------------
    DO n = 1, ntot

        time = REAL(n)*dt

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
            CALL write_nc(outfile, time)

            WRITE (6, *) "Data output at time = ", time
        END IF

    END DO ! end of iteration loop

END PROGRAM wave2D
