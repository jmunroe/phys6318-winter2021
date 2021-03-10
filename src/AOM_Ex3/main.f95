!=========================================
! Exercise 3: Short Surface Gravity Waves
!=========================================
! Author: J. Kaempf, August 2009

PROGRAM slice

    USE param
    USE sub
    USE io

! local parameters
    INTEGER :: n, ntot, nout
    REAL :: wl, ps
    CHARACTER(LEN=50) :: outfile

    period = 8.0 ! forcing period in seconds
    amplitude = 1.0 ! forcing amplitude

    wl = G*period*period/(2.*PI)
    write (6, *) "deep-water wavelength (m) is ", wl
    ps = wl/period
    write (6, *) "deep-water phase speed (m/s) is ", ps

!**********
    CALL INIT  ! initialisation
!**********

! runtime parameters
    ntot = INT(100./dt)

! output parameter
    nout = INT(1./dt)

    ! create output file
    outfile = 'output.nc'
    CALL create_nc(outfile)

    ! write out initial conditions
    CALL write_nc(outfile, 0.)
!---------------------------
! simulation loop
!---------------------------
    DO n = 1, ntot

        time = REAL(n) * dt

! variation of forcing
        ad = amplitude*SIN(2.*PI*time/period)

! prognostic equations
        CALL dyn

! data output
        IF (MOD(n, nout) == 0) THEN
            CALL write_nc(outfile, time)
            
            WRITE (6, *) "Data output at time = ", time

        END IF

    END DO ! end of iteration loop

END PROGRAM slice
