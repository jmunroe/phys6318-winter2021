!=============================
! Exercise 8: Free convection
!=============================
! Author: J. Kaempf, August 2009

PROGRAM slice

    USE param
    USE sub
    USE io

! local parameters
    INTEGER :: ntot, nout
    CHARACTER(LEN=50) :: outfile

!**********
    CALL INIT  ! initialisation
!**********

! runtime
    ntot = INT(6.*3600./dt)
    time = 0.0

! output frequency
    nout = INT(5.*60./dt)

    ! create output file
    outfile = 'output.nc'
    CALL create_nc(outfile)

    ! write out initial conditions
    CALL write_nc(outfile, 0.)

!---------------------------
! simulation loop
!---------------------------
    DO n = 1, ntot

        time = REAL(n)*dt

! prognostic equations
        CALL dyn

! data output
        IF (MOD(n, nout) == 0) THEN
            CALL write_nc(outfile, time)
            WRITE (6, *) "Data output at time (hours) = ", time/(3600.)
        END IF
!ENDIF

    END DO ! end of iteration loop

END PROGRAM slice
