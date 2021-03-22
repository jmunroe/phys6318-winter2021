!===========================================
! Exercise 13: Strfatified flows on a slope
!===========================================
! Author: J. Kaempf, August 2009

PROGRAM slice

    USE param
    USE sub
    USE random
    USE io

! local parameters
    INTEGER :: ntot, nout
    CHARACTER(LEN=50) :: outfile

!**********
    CALL INIT  ! initialisation
!**********

! runtime parameters
    ntot = INT(120.*60./dt)
    time = 0.0

! output frequency
    nout = INT(60./dt)

    ! create output file
    outfile = 'output.nc'
    CALL create_nc(outfile)

    ! write out initial conditions
    CALL write_nc(outfile, 0.)

!initialisation of random function
    ist = -1
    randm = ran3(ist)

    DO k = 0, nx + 1
    DO i = 41, nz + 1
        rho(i, k) = 0.2 + ran3(ist)*0.00001
    END DO
    END DO

! bottom inclination (5 degrees converted to radians)
    alpha = 5.0*PI/180.

!---------------------------
! simulation loop
!---------------------------
    DO n = 1, ntot

        time = REAL(n)*dt
        !write (6, *) "time (hours)", time/(3600.)

! prognostic equations
        CALL dyn

! data output
        IF (MOD(n, nout) == 0) THEN

            CALL write_nc(outfile, time)

           WRITE (6, *) "Data output at time (mins) = ", time/(60.)
        END IF
!ENDIF

    END DO ! end of iteration loop

END PROGRAM slice
