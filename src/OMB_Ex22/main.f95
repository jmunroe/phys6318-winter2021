PROGRAM plume

!*****************************************!
!* 2d reduced-gravity plume model        *!
!*                                       *!
!* including:                            *!
!* - semi-implicit bottom friction       *!
!* - nonlinear terms                     *!
!* - horizontal pressure-gradient force  *!
!* - Coriolis force (f plane)                 *!
!* - flooding algorithm                  *!
!* - TDV advection scheme                  *!
!* - Eulerian tracer prediction          *!
!*                                       *!
!* Author: J. Kaempf, 2008               *!
!*****************************************!

    USE param
    USE sub
    USE io

! local parameters
    INTEGER :: n, ntot, nout
    CHARACTER(len=50) :: outfile

!**********
    CALL INIT  ! initialisation
!**********

! runtime parameters
    ntot = INT(5.*24.*3600./dt)

! output parameter
    nout = INT(1.*3600./dt)

! initial tracer distributions (disabled here)
    DO i = 1, nz
    DO j = 0, ny + 1
    DO k = 0, nx + 1
        T(i, j, k) = 0.
        TN(i, j, k) = 0.
    END DO
    END DO
    END DO

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
        write (6, *) "time (hours)", time/(3600.)

        ad = 0.0 ! adjustment not used in this exercise

! call prognostic equations
        CALL dyn

! data output
        IF (MOD(n, nout) == 0) THEN
            CALL write_nc(outfile, time)

            WRITE (6, *) "Data output at time = ", time/(24.*3600.)
        END IF

    END DO ! end of iteration loop

END PROGRAM plume
