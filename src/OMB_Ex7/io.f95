MODULE io
    USE param

CONTAINS

!=======================
    SUBROUTINE create_nc(outfile)
        USE netcdf
        IMPLICIT NONE
        CHARACTER(LEN=50), INTENT(IN) :: outfile

        INTEGER :: ncid, x_dimid, z_dimid, time_dimid
        INTEGER :: time_id
        INTEGER :: h_id, hzero_id
        INTEGER :: u_id, eta_id
        INTEGER, DIMENSION(3) :: dimids

        !create netCDF dataset: enter define mode
        CALL check(nf90_create(outfile, NF90_CLOBBER, ncid))

        ! define dimensions: from name and length
        CALL check(nf90_def_dim(ncid, "time", NF90_UNLIMITED, time_dimid))
        CALL check(nf90_def_dim(ncid, "x", nx, x_dimid))
        CALL check(nf90_def_dim(ncid, "z", nz, z_dimid))

        !Define coordinate variables
        CALL check(nf90_def_var(ncid, "time", NF90_FLOAT, time_dimid, time_id))
        CALL check(nf90_put_att(ncid, time_id, "units", "seconds since 2000-01-01"))
        CALL check(nf90_put_att(ncid, time_id, "calendar", "gregorian"))

        !CALL check(nf90_def_var(ncid, "x", NF90_FLOAT, x_dimid, x_varid))
        !CALL check(nf90_put_att(ncid, x_varid, "units", "m"))
        !CALL check(nf90_def_var(ncid, "z", NF90_FLOAT, z_dimid, z_varid))
        !CALL check(nf90_put_att(ncid, z_varid, "units", "m"))

        dimids = (/z_dimid, x_dimid, time_dimid/)

        !define variables: from name, type, dims
        CALL check(nf90_def_var(ncid, "hzero", NF90_FLOAT, (/z_dimid, x_dimid/), hzero_id))
        CALL check(nf90_def_var(ncid, "h", NF90_FLOAT, dimids, h_id))
        CALL check(nf90_def_var(ncid, "u", NF90_FLOAT, dimids, u_id))
        CALL check(nf90_def_var(ncid, "eta", NF90_FLOAT, dimids, eta_id))

        !nf90_put_att  ! assign attribute values
        CALL check(nf90_put_att(ncid, h_id, "units", "m"))
        CALL check(nf90_put_att(ncid, h_id, "description", "interval thickness"))
        CALL check(nf90_put_att(ncid, u_id, "units", "m s-1"))
        CALL check(nf90_put_att(ncid, u_id, "description", "velocity"))
        CALL check(nf90_put_att(ncid, eta_id, "units", "m"))
        CALL check(nf90_put_att(ncid, eta_id, "description", "displacement"))

        !end definitions: leave define mode
        CALL check(nf90_enddef(ncid))

        !Write Data
        CALL check(nf90_put_var(ncid, hzero_id, hzero(:, 1:nx)))

        ! close: save new netcdf dataset
        CALL check(nf90_close(ncid))

    END SUBROUTINE create_nc

!:=========================================================================
    SUBROUTINE write_nc(outfile, time)
        USE netcdf
        IMPLICIT NONE
        CHARACTER(LEN=50), INTENT(IN) :: outfile
        REAL, INTENT(IN) :: time
        INTEGER :: ncid, time_dimid
        INTEGER :: u_id, eta_id, h_id, time_id
        INTEGER :: n

        !Open the netCDF file for writing
        CALL check(nf90_open(outfile, NF90_WRITE, ncid))

        ! Get ID of unlimited (time) dimension
        CALL check(nf90_inquire(ncid, unlimiteddimid=time_dimid))

        ! How many time records are there so far?
        CALL check(nf90_inquire_dimension(ncid, time_dimid, len=n))

        !Look up the variable id
        CALL check(nf90_inq_varid(ncid, "h", h_id))

        ! write the data to the file
        CALL check(nf90_put_var(ncid, h_id, h(1:nz, 1:nx), start=(/1, 1, n + 1/), count=(/nz, nx, 1/)))

        CALL check(nf90_inq_varid(ncid, "u", u_id))
        CALL check(nf90_put_var(ncid, u_id, u(1:nz, 1:nx), start=(/1, 1, n + 1/), count=(/nz, nx, 1/)))

        CALL check(nf90_inq_varid(ncid, "eta", eta_id))
        CALL check(nf90_put_var(ncid, eta_id, eta(1:nz, 1:nx), start=(/1, 1, n + 1/), count=(/nz, nx, 1/)))

        CALL check(nf90_inq_varid(ncid, "time", time_id))
        CALL check(nf90_put_var(ncid, time_id, time, start=(/ n+1 /)))

        ! close the netCDF file
        CALL check(nf90_close(ncid))
    END SUBROUTINE write_nc
!:=========================================================================t

!Check (ever so slightly modified from www.unidata.ucar.edu)
!:======================================================================
    SUBROUTINE check(istatus)
        USE netcdf
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: istatus
        IF (istatus /= nf90_noerr) THEN
            write (*, *) TRIM(ADJUSTL(nf90_strerror(istatus)))
        END IF
    END SUBROUTINE check
!:======================================================================

END MODULE io
