PROGRAM kelvin

!*****************************************!
!* 2d shallow-Water model                *!
!*                                       *!
!* including:                            *!
!* - wind-stress forcing                 *!
!* - semi-implicit bottom friction       *!
!* - nonlinear terms                     *!
!* - horizontal pressure-gradient force  *!
!* - Coriolis force (semi-implicit)      *!
!* - Shapiro filter                      *!
!* - flooding algorithm                  *!
!* - TDV advection scheme                  *!
!* - Eulerian tracer prediction          *!
!*                                       *!
!* Author: J. Kaempf, 2008               *!
!*****************************************!
!*
!* Lateral momentum diffusion/friction
!* is not included.
!*
!*****************************************

    USE param
    USE sub
    USE io

! local parameters
    REAL :: time
    REAL :: ad
    INTEGER :: n, ntot, nout
    CHARACTER(len=50) :: outfile

!======================================
! mode = 1 => upstream scheme
! mode = 2 => Lax-Wendroff scheme
! mode = 3 => Superbee scheme
! mode = 4 => Super-C scheme
!=======================================

    mode = 3

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

! set local parameters

! set epsilon for Shapiro filter
    eps = 0.05

    ! runtime parameters
    ntot = INT(12.*3600./dt)

    ! output parameter
    nout = INT(10.*60./dt)

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

        write (6, *) "time (days)", time/(24.*3600.)

        taux = 0.0
        tauy = 0.0

! sea-level forcing
        ad = 1.0*SIN(time/(2.*3600.)*2.*PI)
        eta(6, 2) = ad

! call predictor
        CALL dyn

! data output
        IF (MOD(n, nout) == 0) THEN
            CALL write_nc(outfile, time)

            WRITE (6, *) "Data output at time = ", time/(24.*3600.)

        END IF

    END DO ! end of iteration loop

END PROGRAM kelvin

