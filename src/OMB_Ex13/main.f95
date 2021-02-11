PROGRAM wind

!*****************************************!
!* 2d shallow-Water model                *!
!*                                       *!
!* including:                            *!
!* - wind-stress forcing                 *!
!* - semi-implicit bottom friction       *!
!* - nonlinear terms                     *!
!* - horizontal pressure-gradient force  *!
!* - Shapiro filter                      *!
!* - flooding algorithm                  *!
!* - TDV advection scheme                *!
!*                                       *!
!* Author: J. Kaempf, 2008               *!
!*****************************************!

!******************************************************
! This program requires the bathymetry data "topo.dat"
! from Exercise 10 as an input file!
!******************************************************

    USE param
    USE sub
    USE io

! local parameters
    REAL :: time
    INTEGER :: n, ntot, nout
    CHARACTER(len=50) :: outfile

! Flux limiters available
!======================================
! mode = 1 => upstream scheme
! mode = 2 => Lax-Wendroff scheme
! mode = 3 => Superbee scheme
! mode = 4 => Super-C scheme
!=======================================
! Selections are made in the "dyn" subroutine

!**********
    CALL INIT  ! initialisation
!**********

    ! set local parameters

    ! set epsilon for Shapiro filter
    eps = 0.05

    ! runtime parameters
    ntot = INT(10.*24.*3600./dt)

    ! output parameter
    nout = INT(2.*3600./dt)

    ! create output file
    outfile = OUTPUT_FILE
    !outfile = "output_2.nc"
    CALL create_nc(outfile)

    ! write out initial conditions
    CALL write_nc(outfile, 0.)

!---------------------------
! simulation loop
!---------------------------
    DO n = 1, ntot

        time = REAL(n)*dt

        !write (6, *) "time (days)", time/(24.*3600.)

! gradual adjustment of wind forcing

        taux = 0.0
        tauy = 0.2*MIN(time/(1.*24.*3600.), 1.0)

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
            WRITE (6, *) "Data output at time = ", time/(24.*3600.)
        END IF

    END DO ! end of iteration loop

END PROGRAM wind
