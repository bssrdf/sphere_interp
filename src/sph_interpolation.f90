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

SUBROUTINE interpolation(nlocs, in_data, nlat, nlon, interpolated_values, lambda, lmax_in)
    ! Interpolates in_data onto a sphere using spherical harmonics
    ! see documentation for algorithm and its derivation
    !
    ! inputs:
    !   nlocs -- integer --
    !     the number of locations. This is used to enforce that the first dimension
    !     of locs is 2
    !   in_data -- real, dimension(3,nlocs) -- 
    !     in_data(1:2,i) should be the location (theta, phi) 
    !       with 0 < theta < 2*pi, 0 < phi < pi
    !     in_data(3) should be the corresponding known value at that location
    !   nlat -- integer --
    !       number of latitudes to use
    !   nlon -- integer --
    !       number of longitudes to use
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
    !     value of lambda for the result. See documentation.
    !     Recall the following meaning for lambda:
    !        0  -- harm. sph. representation result is constant (mean)
    !        1  -- harm. sph. representation balances lsq and constant fits
    !       inf -- harm. sph. representation is least squares fit
    !
    !   optional inputs:
    !   lmax -- integer -- 2*nlocs --
    !     maximum degree of spherical harmonics
    !     Default value should capture lots of high frequencies
    !     which the interpolation algorithm will smooth out

    ! inputs
    integer, intent(in) :: nlocs
    real, dimension(3,nlocs), intent(in) :: in_data
    integer, intent(in) :: nlat
    integer, intent(in) :: nlon
    integer, optional, intent(in) :: lmax_in

    ! outputs
    real, intent(out) :: lambda
    real, dimension(0:nlat-1, 0:nlon-1), intent(out) :: interpolated_values

    ! work variables
    real, allocatable, dimension(:, :) :: M_T
    real, allocatable, dimension(:) :: coef
    integer :: lmax

    ! define default lmax if needed
    ! TODO: figure out a good default value for lmax
    lmax = 2*nlocs
    ! This will allow a lot of high frequency noise which should
    ! get smoothed out by algorithm
    if (present(lmax_in)) lmax = lmax_in
    

    ! perform interpolation
    M_T = form_vand_sph_mat(lmax, in_data(1:2,:), nlocs)
    call analysis(M_T, in_data(3,:), coef, lambda)
    interpolated_values = synthesis(nlat, nlon, coef)

END SUBROUTINE interpolation

SUBROUTINE geo_interpolation(nlocs, in_data, nlat, nlon, interpolated_values, lambda, lmax_in)
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
    !     value of lambda for the result. See documentation.
    !     Recall the following meaning for lambda:
    !        0  -- harm. sph. representation result is constant (mean)
    !        1  -- harm. sph. representation balances lsq and constant fits
    !       inf -- harm. sph. representation is least squares fit
    !
    !   optional inputs:
    !     lmax -- integer -- 2*nlocs --
    !       maximum degree of spherical harmonics
    !       The default value should capture lots of high frequencies
    !       which the interpolation algorithm will smooth out

    ! inputs
    integer, intent(in) :: nlocs
    real, dimension(3,nlocs), intent(in) :: in_data
    integer, intent(in) :: nlat
    integer, intent(in) :: nlon
    integer, optional, intent(in) :: lmax_in

    ! outputs
    real, intent(out) :: lambda
    real, dimension(0:nlat-1, 0:nlon-1), intent(out) :: interpolated_values

    ! work variables
    real, dimension(3,nlocs) :: in_data_adj

    ! convert lat/lons into theta/phi
    in_data_adj(3,:) = in_data(3,:)
    in_data_adj(1,:) = in_data(2,:) * 2 * pi / 360.0
    in_data_adj(2,:) = (90.0 - in_data(1,:)) * pi / 180.0

    ! perform interpolation
    call interpolation(nlocs, in_data_adj, nlat, nlon, interpolated_values, lambda, lmax_in)

END SUBROUTINE geo_interpolation

END MODULE SPH_INTERPOLATION
