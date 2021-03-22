!==============================================
! Exercise 10: Slope Convection near the Shore
!==============================================
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

! runtime parameters
    ntot = INT(12.*3600./dt)
    time = 0.0

! output frequency
    nout = INT(6.*60./dt)

! open files for data output
    OPEN (70, file='TRx.dat', form='formatted', recl=10000000)
    OPEN (80, file='TRz.dat', form='formatted', recl=10000000)

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
        !write (6, *) "time (hours)", time/(3600.)

! prognostic equations
        CALL dyn

! float prediction scheme
        DO ii = 1, ntra

! locate grid cell of tracer
            ipos = INT(tra(ii, 1)/dz) + 1
            kpos = INT(tra(ii, 2)/dx + 0.5) + 1
            uu = 0.5*(u(ipos, kpos) + u(ipos, kpos - 1))
            ww = 0.5*(w(ipos, kpos) + w(ipos + 1, kpos))

! avoid stranding
            IF (.not. wet(ipos, kpos) .AND. uu > 0) uu = 0.0
            IF (.not. wet(ipos, kpos - 1) .AND. uu < 0) uu = 0.0
            IF (.not. wet(ipos + 1, kpos) .AND. ww < 0) ww = 0.0
            IF (ipos == 1 .AND. ww > 0) ww = 0.0
            IF (ipos == nz .AND. ww < 0) ww = 0.0

! change of location
            tra(ii, 1) = tra(ii, 1) - dt*ww
            tra(ii, 2) = tra(ii, 2) + dt*uu

            tra(ii, 1) = MAX(tra(ii, 1), 0.0)

        END DO

! data output
        IF (MOD(n, nout) == 0) THEN

            CALL write_nc(outfile, time)

            WRITE (70, '(3000F12.6)') (tra(ii, 2), ii=1, ntra)
            WRITE (80, '(3000F12.6)') (tra(ii, 1), ii=1, ntra)

            WRITE (6, *) "Data output at time (hours) = ", time/(3600.)
        END IF

    END DO ! end of iteration loop

END PROGRAM slice
