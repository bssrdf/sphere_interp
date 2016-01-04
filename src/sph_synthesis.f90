MODULE SPH_SYNTHESIS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A module which contains the function
! to perform spherical synthesis and
! (sph coefficients to function values)
! its helper functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE VAND_SPH
IMPLICIT NONE

CONTAINS

FUNCTION synthesis(nlat, nlon, coef) result(vals)
    ! Calculates function values arising
    ! along an equally spaced grid (in lat/lon)
    ! from a spherical harmonic list of
    ! coefficients
    !  ! inputs:
    !   nlat -- integer --
    !       number of latitudes to use
    !   nlon -- integer --
    !       number of longitudes to use
    !   coef -- real, dimension(lmax**2 + 1) --
    !       coefficients of spherical harmonics
    !
    ! result:
    !   2-d matrix of values -- real, dimension(0:nlat-1, 0:nlon-1) --
    !       Let the output be V, then
    !       V(i,j) = value of function at
    !           phi = i*pi/(nlat-1)
    !           theta = j*2*pi/nlon
    !       NOTE: Since i=0 and i=nlat correspond to the poles,
    !           V is constant for all these values

    ! inputs
    integer, intent(in) :: nlat
    integer, intent(in) :: nlon
    real, dimension(:) :: coef

    ! output
    real, dimension(0:nlat-1, 0:nlon-1) :: vals

    ! work variables
    real, allocatable, dimension(:,:) :: M_T
    real, dimension(2,0:nlon*(nlat - 2) + 1) :: locs
    real, dimension(0:nlon*(nlat - 2) + 1) :: lin_vals
    integer :: nlocs, i, j

    ! form grid 
    nlocs = nlon*(nlat - 2) + 2
    DO j=0,nlon-1
        DO i=1,nlat-2
            !the i index is kinda funny because we only want 1 solution at poles
            locs(1,(i-1)*nlon + j) = j*2*pi/nlon
            locs(2,(i-1)*nlon + j) = i*pi/(nlat-1)
        END DO
    END DO
    !south pole
    locs(1,nlocs-2) = 0.0
    locs(2,nlocs-2) = pi
    !north pole
    locs(1,nlocs-1) = 0.0
    locs(2,nlocs-1) = 0.0

    ! form vand_sph for regular grid
    M_T = form_vand_sph_mat(int(sqrt(size(coef) - 1.0)), locs, nlocs)

    ! calculate values
    lin_vals = MATMUL(coef(2:), M_T) + coef(1)

    ! reshape into output
    ! North pole
    vals(0,:) = lin_vals(nlocs-1)
    DO i=1,nlat-2
        vals(i,:) = lin_vals((i-1)*nlon:i*nlon-1)
    END DO
    ! South Pole
    vals(nlat-1,:) = lin_vals(nlocs-2)

    return

END FUNCTION synthesis

END MODULE SPH_SYNTHESIS
