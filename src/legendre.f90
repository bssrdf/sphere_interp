MODULE LEGENDRE
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A module containing all the
! functions needed for dealing
! with legendre polynomials
!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

! functions to compute complete set of legendre polynomial values
FUNCTION compute_legendre(theta, maxl) result(vals)
    IMPLICIT NONE

    !inputs
    real, intent(in) :: theta
    integer, intent(in) :: maxl

    !outputs
    real, allocatable, dimension(:) :: vals
    
    ! temporary variables
    integer :: l
    real, dimension(1) :: val_0

    val_0 = 1

    allocate(vals(0:(maxl+1)*(maxl+2)/2))

    vals(0:1) = legendre_nextl(theta, val_0)
    DO l=2,maxl
        vals(ml2lIdx(0,l):ml2lIdx(l,l)) = &
            legendre_nextl(theta, vals(ml2lIdx(0,l-1):ml2lIdx(l-1,l-1)))
    END DO
END FUNCTION compute_legendre

FUNCTION legendre_nextl(theta, leg_l) result(leg_nextl)
    IMPLICIT NONE
    ! Computes values of Legendre polynomial P_l+1^m (cos(theta))
    ! for m = 0,l+1
    !
    ! inputs:
    !   theta -- angle at which to evaluate legendre polynomial
    !   leg_l -- list of values such that leg_l(m) = P_l^m (cos(theta)).
    !       Thus should be of dimension(0:l)
    !       Note that if l=1, leg_l = \( 1 \)
    !
    ! return:
    !   real dimension(0:l+1) such that return(l+1,m) = P_l+1^m (cos(theta))

    ! inputs
    real, intent(in) :: theta
    real, intent(in), dimension(0:) :: leg_l

    ! outputs
    real, dimension(0:size(leg_l)+1) :: leg_nextl

    ! define convenience variables
    real :: c,s
    integer :: l,m

    l = size(leg_l)

    c = cos(theta)
    DO m=0,l
        leg_nextl(m) = (2*l+1)*c*leg_l(m) - (l+m)*leg_l(m)
        leg_nextl(m) = leg_nextl(m) / (l-m+1)
    END DO
    leg_nextl(l+1) = -(2*l+1)*sqrt(1-c**2)*leg_l(l)
END FUNCTION legendre_nextl

FUNCTION ml2lIdx(m,l) result(idx)
    IMPLICIT NONE
    ! Computes array index for legendre
    ! polynomial given m,l
    !
    ! Note that it is assumed (0,0) is
    ! assumed to not be in array
    integer m,l,idx

    idx = l*(l+1)/2 + abs(m) - 1
END FUNCTION ml2lIdx

! subroutine doing one set of legendre polynomials
SUBROUTINE legendre_nextl_cos(cos_theta, leg) 
    IMPLICIT NONE
    ! Computes values of Legendre polynomial P_l+1^m (cos(theta))
    ! for m = 0,l+1
    !
    ! inputs:
    !   cos_theta -- real --
    !     cosine of angle at which to evaluate legendre polynomial
    !   leg -- real, allocatable, dimension(0:l)
    !       list of values such that leg_l(m) = P_l^m (cos(theta)).
    !       Note that if l=1, leg = \( 1 \)
    !
    ! output:
    !   leg -- real, allocatable, dimension(0:l+1)
    !       list of values such that leg_l(m) = P_l^m (cos(theta)).

    ! inputs
    real, intent(in) :: cos_theta
    real, allocatable, intent(inout), dimension(:) :: leg

    ! define convenience variables
    real :: c
    integer :: l,m
    real, dimension(0:size(leg)-1) :: leg_copy

    l = size(leg)
    leg_copy = leg

    deallocate(leg)
    allocate(leg(0:l+1))

    c = cos_theta
    DO m=0,l
        leg(m) = (2*l+1)*c*leg_copy(m) - (l+m)*leg_copy(m)
        leg(m) = leg(m) / (l-m+1)
    END DO
    leg(l+1) = -(2*l+1)*sqrt(1-c**2)*leg_copy(l)
END SUBROUTINE legendre_nextl_cos

END MODULE LEGENDRE
