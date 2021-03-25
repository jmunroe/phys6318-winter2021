!==================================
! Exercise 14: Positive estuaries
!==================================
! Author: J. Kaempf, August 2009

PROGRAM slice

    USE param
    USE sub
    USE io

! local parameters
    INTEGER :: ntot, nout, nout2
    REAL :: ttt
    CHARACTER(LEN=50) :: outfile

!**********
    CALL INIT  ! initialisation
!**********

! runtime parameters
    ntot = INT(20.*24.*60.*60./dt)
    time = 0.0

! output frequency
    nout = INT(2.*3600./dt) ! snapshots
    nout2 = INT(12.*3600.0/dt) ! tidally averaged

! open output files
    OPEN (10, file='q.dat', form='formatted')
    OPEN (20, file='u.dat', form='formatted')
    OPEN (30, file='w.dat', form='formatted')
    OPEN (40, file='eta.dat', form='formatted')
    OPEN (50, file='rho.dat', form='formatted')
    OPEN (60, file='rhoM.dat', form='formatted')
    OPEN (61, file='etaR.dat', form='formatted')
    OPEN (62, file='uM.dat', form='formatted')
    OPEN (63, file='age.dat', form='formatted')
    OPEN (64, file='wM.dat', form='formatted')

! output of initial distributions
    DO i = 1, nz
        WRITE (10, '(151F12.6)') (q(i, k)/(RHOREF*G), k=1, nx)
        WRITE (20, '(151F12.6)') (u(i, k), k=1, nx)
        WRITE (30, '(151F12.6)') (w(i, k), k=1, nx)
        WRITE (50, '(151F12.6)') (rho(i, k), k=1, nx)
    END DO

    WRITE (40, '(151F12.6)') (q(0, k)/(RHOREF*G), k=1, nx)

    ! create output file
    outfile = 'output.nc'
    CALL create_nc(outfile)

    ! write out initial conditions
    CALL write_nc(outfile, 0.)

! tidal period
    period = 12.*3600.

!---------------------------
! simulation loop
!---------------------------
    DO n = 1, ntot

        time = time + dt

        ttt = time/(6.*3600.0)
        ad = min(1., ttt) ! adjustment parameter

        write (6, *) "time (days)", time/(24.*3600.), ad

! prognostic equations
        CALL dyn

! snapshot data output
        IF (MOD(n, nout) == 0) THEN
            CALL write_nc(outfile, time)

            DO i = 1, nz
                WRITE (10, '(151F12.6)') (q(i, k)/(RHOREF*G), k=1, nx)
                WRITE (20, '(151F12.6)') (u(i, k), k=1, nx)
                WRITE (30, '(151F12.6)') (w(i, k), k=1, nx)
                WRITE (50, '(151F12.6)') (rho(i, k), k=1, nx)
            END DO
            WRITE (40, '(151F12.6)') (q(0, k)/(RHOREF*G), k=1, nx)
            WRITE (6, *) "Data output at time = ", time/(24.*3600.)
        END IF

! output of tidally averaged fields
        IF (MOD(n, nout2) == 0) THEN
            DO i = 1, nz
                WRITE (60, '(151F12.6)') (rhoM(i, k)/(12.*3600.), k=1, nx)
                WRITE (62, '(151F12.6)') (uM(i, k)/(12.*3600.), k=1, nx)
                WRITE (64, '(151F12.6)') (wM(i, k)/(12.*3600.), k=1, nx)
                WRITE (63, '(151F12.6)') (c1(i, k)/(24.*3600), k=1, nx)
! reset fields
                DO k = 0, nx + 1
                    rhoM(i, k) = 0.0
                    uM(i, k) = 0.0
                    wM(i, k) = 0.0
                END DO
            END DO
! output of tidal range
            WRITE (61, '(151F12.6)') (etaMAX(k) - etaMIN(k), k=1, nx)
            DO k = 0, nx + 1
                etaMAX(k) = 0.0
                etaMIN(k) = 0.0
            END DO
            WRITE (6, *) "Data output at time = ", time/(24.*3600.)
        END IF

    END DO ! end of iteration loop

END PROGRAM slice
