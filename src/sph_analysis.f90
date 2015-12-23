MODULE SPH_ANALYSIS
    USE VAND_SPH
    IMPLICIT NONE

CONTAINS

SUBROUTINE toab(coef, mean, a, b)
    ! This subroutine takes in a coefficient array
    ! and a function mean and translates it into
    ! coefficient matrices a,b which can be used in spherepack
    !
    ! inputs:
    !   coef -- real, dimension(0:(l+1)**2) --
    !       coefficient array for spherical harmonics
    !
    ! outputs:
    !   a -- allocatable array, can be fed into spherepack
    !   b -- allocatable array, can be fed into spherepack

    ! inputs
    real, dimension(0:),intent(in) :: coef
    real, intent(in) :: mean

    ! outputs
    real, allocatable, dimension(:,:), intent(out) :: a
    real, allocatable, dimension(:,:), intent(out) :: b

    ! helper variables
    integer :: maxl, l, m
    
    ! Since coef is a linear array ranging from 0:(maxl+1)**2
    maxl = int(sqrt(size(coef)-1.0) - 1)

    allocate( a(0:maxl,0:maxl) )
    allocate( b(0:maxl,0:maxl) )

    do l=0,maxl
    do m=0,l
        a(m,l) = coef(ml2sIdx(m,l))
        if (m > 0) then
            b(m,l) = coef(ml2sIdx(-1*m,l))
        else
            b(m,l) = 0
        endif
    end do
    end do
END SUBROUTINE toab

SUBROUTINE scalef(f, mean, scalefact)
    ! Scales f to values of (-1,1)
    !
    ! inputs:
    !   f -- real, dimension(1:nloc) --
    !     input array
    !
    ! outputs:
    !   f -- real, dimension(1:nloc) --
    !     rescaled input array
    !   mean -- real --
    !     mean of input f
    !   scalefact -- real --
    !     scale factor for f
    !
    !  input f can be formed by taking
    !  scalefact*(f + mean)
    
    !inputs
    real, dimension(:), intent(inout) :: f

    !outputs
    real, intent(out) :: mean
    real, intent(out) :: scalefact

    !temporary variables
    integer :: nloc

    ! rescale and save values
    nloc = size(f)
    mean = sum(f)
    mean = mean/float(nloc)
    scalefact = maxval(abs(f))
    f = f/scalefact

END SUBROUTINE scalef

END MODULE SPH_ANALYSIS
