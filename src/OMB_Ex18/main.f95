PROGRAM wind

!*****************************************!
!* 2d shallow-Water model                *!
!*                                       *!
!* including:                            *!
!* - wind-stress forcing                 *!
!* - semi-implicit bottom friction       *!
!* - nonlinear terms                     *!
!* - horizontal pressure-gradient force  *!
!* - Coriolis force (beta plane)         *!
!* - flooding algorithm                  *!
!* - TDV advection scheme                  *!
!* - Eulerian tracer prediction          *!
!* - Lagrangian float prediction         *!
!*                                       *!
!* Author: J. Kaempf, 2008               *!
!*****************************************!

    USE param
    USE sub
    USE io

! local parameters
    REAL :: time, day
    INTEGER :: n, ntot, nout, nnout
    CHARACTER(len=50) :: outfile

!**********
    CALL INIT  ! initialisation
!**********

! initialisation of Eulerian tracer

    DO j = 0, ny + 1
    DO k = 0, nx + 1
        T(j, k) = 0.
        TN(j, k) = 0.
    END DO
    END DO

! set epsilon for Shapiro filter
    eps = 0.00

! runtime parameters
    ntot = INT(150.*24.*3600./dt)
    time = 0.0

    nnout = ntot/150

    ! output parameter
    nout = INT(60.*3600./dt)

    ! create output file
    outfile = "output.nc"
    CALL create_nc(outfile)

    ! write out initial conditions
    CALL write_nc(outfile, 0.)

! open files for output (float locations)
    OPEN (101, file='TRx.dat', form='formatted', recl=10000000)
    OPEN (102, file='TRy.dat', form='formatted', recl=10000000)

!---------------------------
! simulation loop
!---------------------------
    DO n = 1, ntot

        time = REAL(n)*dt

        day = time/(24.*3600)

        write (6, *) "time (days)", time/(24.*3600.)

        ad = MIN(time/(50.*24.*3600.), 1.0)

! call predictor
        CALL dyn

! Lagrangian float prediction

        DO i = 1, ntra

! locate grid cell of tracer
            jpos = INT((tra(i, 1) + 0.5*dy)/dy) + 1
            kpos = INT((tra(i, 2) + 0.5*dx)/dx) + 1
            uu = 0.5*(u(jpos, kpos) + u(jpos, kpos - 1))
            vv = 0.5*(v(jpos, kpos) + v(jpos - 1, kpos))

! change of location
            tra(i, 1) = tra(i, 1) + dt*vv
            tra(i, 2) = tra(i, 2) + dt*uu

        END DO

! data output (float locations)
        IF (INT(day) > 49) THEN
            IF (MOD(n, nnout) == 0) THEN
                WRITE (101, '(5000F12.6)') (tra(i, 2)/1000., i=1, ntra)
                WRITE (102, '(5000F12.6)') (tra(i, 1)/1000., i=1, ntra)
                WRITE (6, *) "Float output at time = (days)", day
            END IF
        END IF

! data output (snapshot distributions)
        IF (MOD(n, nout) == 0) THEN
            CALL write_nc(outfile, time)

            WRITE (6, *) "Data output at time = ", time/(24.*3600.)
        END IF

    END DO ! end of iteration loop

END PROGRAM wind
