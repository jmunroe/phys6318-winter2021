MODULE io
    USE param

CONTAINS

!=======================
    SUBROUTINE create_nc(outfile)
        USE netcdf

        IMPLICIT NONE

        CHARACTER(LEN=50), INTENT(IN) :: outfile
        INTEGER :: ncid, time_dimid
        INTEGER :: xh_dimid, yh_dimid
        INTEGER :: xu_dimid, yu_dimid
        INTEGER :: xv_dimid, yv_dimid
        INTEGER :: z_dimid
        INTEGER :: time_id
        INTEGER :: h_id, eta_id, T_id
        INTEGER :: u_id, v_id
        INTEGER :: xh_id, yh_id
        INTEGER :: xu_id, yu_id
        INTEGER :: xv_id, yv_id
        INTEGER :: z_id

        !create netCDF dataset: enter define mode
        CALL check(nf90_create(outfile, NF90_CLOBBER, ncid))

        ! define dimensions: from name and length
        CALL check(nf90_def_dim(ncid, "time", NF90_UNLIMITED, time_dimid))
        CALL check(nf90_def_dim(ncid, "xh", nx + 2, xh_dimid))
        CALL check(nf90_def_dim(ncid, "yh", ny + 2, yh_dimid))
        CALL check(nf90_def_dim(ncid, "xu", nx + 2, xu_dimid))
        CALL check(nf90_def_dim(ncid, "yu", ny + 2, yu_dimid))
        CALL check(nf90_def_dim(ncid, "xv", nx + 2, xv_dimid))
        CALL check(nf90_def_dim(ncid, "yv", ny + 2, yv_dimid))
        CALL check(nf90_def_dim(ncid, "z", nz, z_dimid))

        !Define coordinate variables
        CALL check(nf90_def_var(ncid, "time", NF90_FLOAT, time_dimid, time_id))
        CALL check(nf90_put_att(ncid, time_id, "units", "seconds since 2000-01-01"))
        CALL check(nf90_put_att(ncid, time_id, "calendar", "gregorian"))

        CALL check(nf90_def_var(ncid, "xh", NF90_FLOAT, xh_dimid, xh_id))
        CALL check(nf90_put_att(ncid, xh_id, "units", "m"))
        CALL check(nf90_def_var(ncid, "yh", NF90_FLOAT, yh_dimid, yh_id))
        CALL check(nf90_put_att(ncid, yh_id, "units", "m"))

        CALL check(nf90_def_var(ncid, "xu", NF90_FLOAT, xu_dimid, xu_id))
        CALL check(nf90_put_att(ncid, xu_id, "units", "m"))
        CALL check(nf90_def_var(ncid, "yu", NF90_FLOAT, yu_dimid, yu_id))
        CALL check(nf90_put_att(ncid, yu_id, "units", "m"))

        CALL check(nf90_def_var(ncid, "xv", NF90_FLOAT, xv_dimid, xv_id))
        CALL check(nf90_put_att(ncid, xv_id, "units", "m"))
        CALL check(nf90_def_var(ncid, "yv", NF90_FLOAT, yv_dimid, yv_id))
        CALL check(nf90_put_att(ncid, yv_id, "units", "m"))

        CALL check(nf90_def_var(ncid, "z", NF90_FLOAT, z_dimid, z_id))
        CALL check(nf90_put_att(ncid, z_id, "units", "m"))

        !define variables: from name, type, dims
        CALL check(nf90_def_var(ncid, "h", NF90_FLOAT, (/z_dimid, yh_dimid, xh_dimid, time_dimid/), h_id))
        CALL check(nf90_def_var(ncid, "eta", NF90_FLOAT, (/z_dimid, yh_dimid, xh_dimid, time_dimid/), eta_id))
        CALL check(nf90_def_var(ncid, "u", NF90_FLOAT, (/z_dimid, yu_dimid, xu_dimid, time_dimid/), u_id))
        CALL check(nf90_def_var(ncid, "v", NF90_FLOAT, (/z_dimid, yv_dimid, xv_dimid, time_dimid/), v_id))
        CALL check(nf90_def_var(ncid, "T", NF90_FLOAT, (/z_dimid, yh_dimid, xh_dimid, time_dimid/), T_id))

        !nf90_put_att  ! assign attribute values
        CALL check(nf90_put_att(ncid, h_id, "units", "m"))
        CALL check(nf90_put_att(ncid, h_id, "description", "interval thickness"))
        CALL check(nf90_put_att(ncid, u_id, "units", "m s-1"))
        CALL check(nf90_put_att(ncid, u_id, "description", "velocity"))
        CALL check(nf90_put_att(ncid, v_id, "units", "m s-1"))
        CALL check(nf90_put_att(ncid, v_id, "description", "velocity"))
        CALL check(nf90_put_att(ncid, eta_id, "units", "m"))
        CALL check(nf90_put_att(ncid, eta_id, "description", "displacement"))
        CALL check(nf90_put_att(ncid, T_id, "units", ""))
        CALL check(nf90_put_att(ncid, T_id, "description", "Eulerian tracer"))

        !end definitions: leave define mode
        CALL check(nf90_enddef(ncid))

        !Write Data
        CALL check(nf90_put_var(ncid, xh_id, xh))
        CALL check(nf90_put_var(ncid, yh_id, yh))
        CALL check(nf90_put_var(ncid, xu_id, xu))
        CALL check(nf90_put_var(ncid, yu_id, yu))
        CALL check(nf90_put_var(ncid, xv_id, xv))
        CALL check(nf90_put_var(ncid, yv_id, yv))

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
        INTEGER :: u_id, v_id, eta_id, h_id, T_id, time_id
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
        CALL check(nf90_put_var(ncid, h_id, h, start=(/1, 1, 1, n + 1/), count=(/nz, ny + 2, nx + 2, 1/)))

        CALL check(nf90_inq_varid(ncid, "u", u_id))
        CALL check(nf90_put_var(ncid, u_id, u, start=(/1, 1, 1, n + 1/), count=(/nz, ny + 2, nx + 2, 1/)))

        CALL check(nf90_inq_varid(ncid, "v", v_id))
        CALL check(nf90_put_var(ncid, v_id, v, start=(/1, 1, 1, n + 1/), count=(/nz, ny + 2, nx + 2, 1/)))

        CALL check(nf90_inq_varid(ncid, "eta", eta_id))
        CALL check(nf90_put_var(ncid, eta_id, eta, start=(/1, 1, 1, n + 1/), count=(/nz, ny + 2, nx + 2, 1/)))

        CALL check(nf90_inq_varid(ncid, "T", T_id))
        CALL check(nf90_put_var(ncid, T_id, T, start=(/1, 1, 1, n + 1/), count=(/nz, ny + 2, nx + 2, 1/)))

        CALL check(nf90_inq_varid(ncid, "time", time_id))
        CALL check(nf90_put_var(ncid, time_id, time, start=(/n + 1/)))

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
