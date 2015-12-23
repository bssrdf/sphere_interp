MODULE VAND_SPH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A module containing all the
! functions for forming a
! vandermonde-like matrix for
! the spherical harmonics
!
! OpenMP parallelizes across locations
! so each processor takes care of a subset
! of the locations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE legendre, ONLY: legendre_nextl_cos
IMPLICIT NONE

real, parameter :: pi = 4*atan(1.0)
real, dimension(:), allocatable :: facts

CONTAINS

! procedures for forming vand_sph
FUNCTION form_vand_sph_mat(lmax,locs,nlocs) result(M_T)
    ! Forms the vandermonde-like matrix
    ! This is a convenience function that calls into
    ! initialize_vand_sph and set_basis_vand_sph
    !
    ! inputs:
    !   lmax -- integer -- 
    !     maximum degree of spherical harmonics
    !   locs -- real, dimension(2,nlocs) -- 
    !     list of locations (theta, phi) with 0 < theta < 2*pi, 0 < phi < pi
    !   nlocs -- integer --
    !     the number of locations. This is used to enforce that the first dimension
    !     of locs is 2
    !
    ! result:
    !   vandermonde-like matrix -- dimension(lmax**2, nlocs) --
    !     entry (i,j) corresponds to basis function i, location j

    ! inputs
    integer, intent(in) :: lmax
    real, intent(in), dimension(2,*) :: locs
    integer, intent(in) :: nlocs

    ! output
    real, allocatable, dimension(:, :) :: M_T

    !temporary variable
    integer ierror

    ! initialize M_T
    M_T = initialize_vand_sph(locs, nlocs)

    ! build it out for full basis
    call set_basis_vand_sph(M_T, lmax, ierror)
    IF (ierror /= 0) THEN
        write(*,*) 'Issue building M_T'
    END IF
END FUNCTION form_vand_sph_mat

FUNCTION initialize_vand_sph(locs, nlocs) result(M_T)
    ! Initializes the vandermonde-like matrix (using only l=1)
    !
    ! inputs:
    !   locs -- real, dimension(2,nlocs) -- 
    !     list of locations (theta, phi) with 0 < theta < 2*pi, 0 < phi < pi
    !   nlocs -- integer --
    !     the number of locations. This is used to enforce that the first dimension
    !     of locs is 2
    !
    ! result:
    !   vandermonde-like matrix -- dimension(3, nlocs) --
    !     entry (i,j) corresponds to basis function i, location j

    ! inputs
    real, intent(in), dimension(2,*) :: locs
    integer, intent(in) :: nlocs

    ! output
    real, allocatable, dimension(:,:) :: M_T

    ! temporary variables
    integer :: p
    real, allocatable, dimension(:) :: leg

    allocate(M_T(3,nlocs))
   
!$OMP PARALLEL DO PRIVATE(leg)
    DO p=1,nlocs
        allocate(leg(1))
        leg = 1
        call eval_sph(cos(locs(1,p)), locs(2,p), leg, M_T(:,p))
        deallocate(leg)
    END DO
!$OMP END PARLLEL DO
END FUNCTION initialize_vand_sph

SUBROUTINE set_basis_vand_sph(M_T, lmax, ierror)
    ! Appends basis functions as appropriate to create
    ! vandermonde-like matrix with given lmax
    !
    ! inputs:
    !   M_T -- real, dimension(lmax_prev*(lmax_prev+2), nlocs) --
    !     previously computed M_T, must be allocatable
    !  lmax -- integer --
    !     desired maximum l of basis functions
    !
    ! outputs
    !   M_T -- real, dimension(lmax*(lmax+2), nlocs) --
    !     M_T of desired dimensions
    !   ierror -- integer --
    !     error indicator. Possible values are:
    !       0 - calculation finished
    !       1 - M_T has not been initialized
    !       2 - deallocation error
    !       3 - allocation error

    ! input
    integer, intent(in) :: lmax
    real, intent(inout), allocatable, dimension(:,:) :: M_T

    ! output
    integer :: ierror
    
    ! temporary variables
    real :: cos_theta, phi
    real :: t_r, c1, c2
    integer :: m, l, loc, ALLOC_ERR
    integer :: prev_lmax, nloc
    real, dimension(size(M_T,1), size(M_T,2)) :: M_T_copy
    real, allocatable, dimension(:) :: leg_val
    real, allocatable, dimension(:) :: leg
   
    IF (.NOT. allocated(M_T)) THEN
        ierror = 1     
        return
    END IF

    ! extract needed parameters
    prev_lmax = int(sqrt(1.0 + size(M_T,1))) - 1
    nloc = size(M_T,2)
    M_T_copy = M_T

    ! if prev_lmax == lmax, do nothing
    IF (prev_lmax == lmax) THEN
        ierror = 0
        return
    END IF

    ! resize M_T
    deallocate(M_T, STAT = ALLOC_ERR)
    IF ( ALLOC_ERR .NE. 0 ) THEN
        ierror = 2
        return
    END IF
    allocate(M_T(lmax*(lmax+2),nloc), STAT = ALLOC_ERR)
    IF ( ALLOC_ERR .NE. 0 ) THEN
        ierror = 3
        return
    END IF

    ! the trivial case where prev_lmax > lmax, so truncating
    IF (prev_lmax > lmax) THEN
       M_T(:,:) = M_T_copy(1:lmax*(lmax+2),:)
       ierror = 0
       return
    END IF

    ! if lmax > prev_lmax, need to calculate further
    M_T(1:prev_lmax*(prev_lmax+2),:) = M_T_copy(:,:)
!$OMP PARLLEL DO PRIVATE(leg_val, phi, cos_theta)
    DO loc=1,nloc
        ! instead of calculating all the legendre polynomials from
        ! scratch, extract them from M_T and calcuate the new ones
        
        call extract_basis_info(M_T(:,loc), prev_lmax, leg_val, phi, cos_theta)
        ! cannot perform nested parallelization because leg_val depends on its previous value
        DO l=prev_lmax+1,lmax
            call eval_sph(cos_theta, phi, leg_val, M_T( ml2sIdx(-l,l):ml2sIdx(l,l) ,loc ) )
        END DO
    END DO
!$OMP END PARALLEL DO
END SUBROUTINE set_basis_vand_sph

! helper routines that do most of the computation 
SUBROUTINE eval_sph(cos_theta, phi, leg, sph_vals)
    ! calculates values for spherical harmonics
    ! of order l at theta_phi
    !
    ! inputs:
    !   cos_theta -- real -- 
    !     cosine of theta, 0 < theta < 2*pi
    !   phi -- real --
    !     angle from North (0 < phi < pi)
    !   leg -- real, allocatable, dimension(0:[l or (l-1)]) --
    !     value of legendre polyniomials for [l or (l-1)]
    !     this is determined by the size of sph_vals
    !   sph_vals -- real, dimension(2*l + 1) --
    !     the size of this indicates which l is to
    !     be calculated
    !
    !  outputs:
    !    leg -- real, dimension(0:l+1) --
    !      value of legendre polynimials for l+1
    !    sph_vals -- real, dimension(-l:l) -- 
    !      sph_vals(m) -> Y_l^m(theta, phi)

    ! inputs
    real, intent(in) :: cos_theta
    real, intent(in) :: phi
    
    ! outputs
    real, intent(inout), allocatable, dimension(:) :: leg
    real, intent(inout), dimension(0:) :: sph_vals

    ! temporary variables
    integer :: m,l
    real :: c1, c2

    ! figure out which l we are trying to fill
    l = int(size(sph_vals) - 1)/2

    ! if leg is of dimension(0:l-1), update leg
    IF(l == size(leg)) THEN
        call legendre_nextl_cos(cos_theta, leg) 
    ELSE IF (l+1 /= size(leg)) THEN
        ! THIS IS BAD - error should be thrown, but I'm lazy
        ! and don't know how to throw errors in fortran
        sph_vals(:) = 0
        return
    END IF

    ! calculate values
    c1 = sqrt((2*l + 1)/(4*pi))
    ! Y_l^0
    sph_vals(0 + l) = c1 * leg(0)
    DO m=1,l
        ! normalization
        c2 = factorial(l - m)
        c2 = c2 / factorial(l + m)
        c2 = sqrt(2*c2)
        ! legendre polynomial value
        c2 = c1*c2*leg(m)
        ! Y_l^-m
        sph_vals(-m+l) = c2*sin(m*phi)
        !Y_l^m
        sph_vals(m+l) = c2*cos(m*phi)
    END DO
END SUBROUTINE eval_sph

SUBROUTINE extract_basis_info(M_T_col, prev_lmax, leg, phi, cos_theta)
    !Given a column of M_T, will extract the
    !values of the legendre polynomial in this row
    
    ! inputs
    real, intent(in), dimension(:) :: M_T_col
    integer, intent(in) :: prev_lmax

    ! outputs
    real, intent(out), allocatable, dimension(:) :: leg
    real, intent(out) :: phi
    real, intent(out) :: cos_theta

    ! temporary variables
    real :: t_r, c1, c2
    integer :: m

    allocate( leg(0:prev_lmax) )

    ! get angles
    ! calculate angle phi
    t_r = 3/(4*pi)
    t_r = t_r - M_T_col(ml2sIdx(0,1))**2
    t_r = sqrt(t_r)
    t_r = - M_T_col(ml2sIdx(1,1))/t_r
    phi = acos( t_r )

    ! cos of theta
    t_r = sqrt(4*pi/3)
    cos_theta = t_r*M_T_col(ml2sIdx(1,0))

    ! extraction of leg_l from M_T
    c1 = sqrt((2*prev_lmax + 1)/(4*pi))
    leg(0) = M_T_col(ml2sIdx(0,prev_lmax)) / c1
    DO m=1,prev_lmax
        c2 = factorial(prev_lmax - abs(m))
        c2 = c2 / factorial(prev_lmax + abs(m))
        c2 = sqrt(2*c2)
        ! legendre polynomial value
        c2 = c1*c2
        leg(m) = M_T_col(ml2sIdx(m,prev_lmax))/(cos(m*phi)*c2)
    END DO
END SUBROUTINE extract_basis_info

FUNCTION factorial(n) result(val)
    !Calculates the factorial
    !
    !This is done by filling a storage array,
    !so that it is quick to retrieve values if called repeadedly
    integer, intent(in) :: n
    integer :: n_old, i
    real :: val
    real, dimension(0:n) :: tmp

    !retrieve value from memory if available
    if (allocated(facts) .and. size(facts) >= n+1) then
        val = facts(n)
        return
    end if

    !copy over old values
    IF (allocated(facts)) THEN
        n_old = size(facts)-1
        tmp(0:n_old) = facts(0:n_old) 
        deallocate(facts)
    ELSE
        n_old = 0
        tmp(0) = 1
    END IF

    !allocate space for new values
    allocate(facts(0:n))
    facts(0:n_old) = tmp(0:n_old)

    DO i=n_old+1,n
        facts(i) = i*facts(i-1)
    END DO
    
    val = facts(n)
    return
END FUNCTION factorial

! helper function for indexing
FUNCTION ml2sIdx(m,l) result(idx)
    ! Computes array index for spherical
    ! harmonics given m,l
    !
    ! Note that it is assumed (0,0) is
    ! assumed to not be in array
    integer m,l,idx

    idx = l*(l+1) + m
END FUNCTION ml2sIdx

END MODULE VAND_SPH
