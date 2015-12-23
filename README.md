# sphere_interp
A Fortran package that performs interpolation on a sphere using spherical harmonics

The theory and main algorithm can be found in the doc folder
as a latex document.

Implementation is done in the src folder.
This includes modules to build a vandermonde like matrix
given data locations (vand_sph.f90), calculate legendre
polynomials (legendre.f90), and perform spherical analysis
(sph_analysis).

In the future this will call into spherepack subroutines in
order to perform the synthesis step.
Furthermore, there will be hooks for python to call into this
package.

Note on language:
Following the convention of spherepack,
spherical analysis calculates coefficients of spherical harmonics
given data points.
Whereas spherical synthesis uses those coefficients to calculate
interpolated values.
