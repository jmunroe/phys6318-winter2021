PROGRAM multi

!*****************************************!
!* 1d multi-layer shallow-Water model    *!
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
    REAL :: time
    REAL :: amplitude, period
    REAL :: hmax
    INTEGER :: n, ntot, nout
    CHARACTER(LEN=50) :: outfile

    period = 8.0 ! forcing period in seconds
    amplitude = 1.0 ! forcing amplitude

!---- initialisation ----
    CALL INIT
!------------------------

    hmax = 100.0 ! total water depth
    dt = 0.05

    wl = period*sqrt(g*hmax)
    write (6, *) "shallow-water wavelength (m) is ", wl
    ps = wl/period
    write (6, *) "shallow-water phase speed (m/s) is ", ps


! set epsilon for Shapiro filter
    eps = 0.05


! runtime parameters
    ntot = INT(100./dt)

! output parameter
    nout = INT(1./dt)

    ! create output file
    outfile = 'output.nc'
    CALL create_nc(outfile)

    ! write out initial conditions
    CALL write_nc(outfile, 0.)

!---- simulation loop ----
    DO n = 1, ntot
!-------------------------
        time = REAL(n) * dt

        DO i = 1, nz
            eta(i, 1) = amplitude*SIN(2.*pi*time/period)
        END DO

!---- prognostic equations ----
        CALL dyn
!------------------------------
! data output
        IF (MOD(n, nout) == 0) THEN

            CALL write_nc(outfile, time)

            WRITE (6, *) "Data output at time = ", time
        END IF

!---- end of iteration loop ----
    END DO
!-------------------------------

END PROGRAM multi
