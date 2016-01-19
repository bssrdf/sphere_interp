MODULE SPH_INTERPOLATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A Module which provides subroutines
! to perform interpolation using
! spherical harmonics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE SPH_ANALYSIS
USE SPH_SYNTHESIS
USE VAND_SPH
IMPLICIT NONE

CONTAINS

SUBROUTINE interpolation(nlocs, in_data, nlat, nlon, lmax_in, lsq, interpolated_values, lambda, sv)
    ! Interpolates in_data onto a sphere using spherical harmonics
    ! see documentation for algorithm and its derivation
    !
    ! inputs:
    !   nlocs -- integer --
    !     the number of locations. This is used to enforce that the first dimension
    !     of in_data is 3
    !   in_data -- real, dimension(3,nlocs) -- 
    !     in_data(1:2,i) should be the location (theta, phi) 
    !       with 0 < theta < 2*pi, 0 < phi < pi
    !     in_data(3) should be the corresponding known value at that location
    !   nlat -- integer --
    !       number of latitudes to use
    !   nlon -- integer --
    !       number of longitudes to use
    !   lmax_in -- integer -- 
    !     maximum degree of spherical harmonics.
    !     if lmax < 1, then defaults to CEILING(sqrt(2*nlocs))
    !     Default value garuntees at least twice as many basis functions
    !     as data points are used - giving plenty of high frequency
    !     representations which the interpolation algorithm will smooth out
    !   lsq -- logical -- 
    !     if .true. uses a least squares fit rather than the algorithm
    !     to perform the analysis
    !
    ! outputs:
    !   interpolated_values -- real, dimension(0:nlat-1, 0:nlon-1) --
    !     Let the output be V, then
    !     V(i,j) = value of function at
    !         phi = i*pi/(nlat-1)
    !         theta = j*2*pi/nlon
    !     NOTE: Since i=0 and i=nlat correspond to the poles,
    !         V is constant for all these values
    !   lambda -- real --
    !     If least_squares is .false., then
    !       value of lambda for the result. See documentation.
    !       Recall the following meaning for lambda:
    !          0  -- harm. sph. representation result is constant (mean)
    !          1  -- harm. sph. representation balances lsq and constant fits
    !         inf -- harm. sph. representation is least squares fit
    !     If least_squares is .true., then this is the rank of M
    !   sv -- real, dimension(nlocs) --
    !     If least_squares is .true., then
    !       singular values of vandermonde-like matrix, extra values are set to 0
    !     otherwise an unallocated array
    !
    !   optional inputs:

    ! inputs
    integer, intent(in) :: nlocs
    real, dimension(3,nlocs), intent(in) :: in_data
    integer, intent(in) :: nlat
    integer, intent(in) :: nlon
    integer, intent(in) :: lmax_in
    logical, intent(in) :: lsq

    ! outputs
    real, dimension(0:nlat-1, 0:nlon-1), intent(out) :: interpolated_values
    real, intent(out) :: lambda
    real, dimension(nlocs), intent(out) :: sv

    ! work variables
    real, allocatable, dimension(:, :) :: M_T
    real, allocatable, dimension(:) :: coef
    integer :: info
    integer :: rank
    integer :: lmax
    real, allocatable, dimension(:) :: sv_tmp

    ! define default lmax if needed
    ! This will allow a lot of high frequency noise which should
    ! get smoothed out by algorithm
    if (lmax_in .LE. 0) then
        ! TODO: figure out a good default value for lmax
        lmax = CEILING(SQRT(2.0*nlocs))
    else
        lmax = lmax_in
    endif
    
    
    ! perform interpolation
    sv = 0
    M_T = form_vand_sph_mat(lmax, in_data(1:2,:), nlocs)
    if (lsq) then
        call analysis_lsq(M_T, in_data(3,:), coef, sv_tmp, rank, info)
        if (info .NE. 0) then
            write(*,*) 'svd failed with info = ', info
        endif
        lambda = rank
        sv(1:size(sv_tmp)) = sv(:)
    else
        call analysis(M_T, in_data(3,:), coef, lambda)
    endif

    interpolated_values = synthesis(nlat, nlon, coef)

END SUBROUTINE interpolation

SUBROUTINE geo_interpolation(nlocs, in_data, nlat, nlon, lmax_in, lsq, interpolated_values, lambda, sv)
    ! Interpolates in_data onto a sphere using spherical harmonics
    ! see documentation for algorithm and its derivation
    !
    ! This is a wrapper function to interpolation which translates lat/lons
    ! to the theta/phi used in throughout the algorithm
    !
    ! inputs:
    !   nlocs -- integer --
    !     the number of locations. This is used to enforce that the first dimension
    !     of locs is 2
    !   in_data -- real, dimension(3,nlocs) -- 
    !     in_data(1:2,i) should be the location (lat, lon) 
    !       with -90 < lat < 90, 0 < lon < 360
    !     in_data(3) should be the corresponding known value at that location
    !   nlat -- integer --
    !       number of latitudes to use
    !   nlon -- integer --
    !       number of longitudes to use
    !   lmax_in -- integer -- 
    !     maximum degree of spherical harmonics.
    !     if lmax < 1, then defaults to CEILING(sqrt(2*nlocs))
    !     Default value garuntees at least twice as many basis functions
    !     as data points are used - giving plenty of high frequency
    !     representations which the interpolation algorithm will smooth out
    !   lsq -- logical -- 
    !     if .true. uses a least squares fit rather than the algorithm
    !     to perform the analysis
    !
    ! outputs:
    !   interpolated_values -- real, dimension(0:nlat-1, 0:nlon-1) --
    !     Let the output be V, then
    !     V(i,j) = value of function at
    !         phi = 90 - i*180/(nlat-1)
    !         theta = j*360/nlon
    !     NOTE: Since i=0 and i=nlat correspond to the poles,
    !         V is constant for all these values
    !   lambda -- real --
    !     If least_squares is .false., then
    !       value of lambda for the result. See documentation.
    !       Recall the following meaning for lambda:
    !          0  -- harm. sph. representation result is constant (mean)
    !          1  -- harm. sph. representation balances lsq and constant fits
    !         inf -- harm. sph. representation is least squares fit
    !     If least_squares is .true., then this is the rank of M
    !   sv -- real, allocatable, dimension(rank(M)) --
    !     If least_squares is .true., then
    !       singular values of vandermonde-like matrix, extra values are set to 0
    !     otherwise an unallocated array

    ! inputs
    integer, intent(in) :: nlocs
    real, dimension(3,nlocs), intent(in) :: in_data
    integer, intent(in) :: nlat
    integer, intent(in) :: nlon
    integer, intent(in) :: lmax_in
    logical, intent(in) :: lsq

    ! outputs
    real, intent(out) :: lambda
    real, dimension(0:nlat-1, 0:nlon-1), intent(out) :: interpolated_values
    real, dimension(nlocs), intent(out) :: sv

    ! work variables
    real, dimension(3,nlocs) :: in_data_adj

    ! convert lat/lons into theta/phi
    in_data_adj(3,:) = in_data(3,:)
    in_data_adj(1,:) = in_data(2,:) * 2 * pi / 360.0
    in_data_adj(2,:) = (90.0 - in_data(1,:)) * pi / 180.0

    ! perform interpolation
    call interpolation(nlocs, in_data_adj, nlat, nlon, lmax_in, lsq, interpolated_values, lambda, sv)

END SUBROUTINE geo_interpolation

END MODULE SPH_INTERPOLATION
