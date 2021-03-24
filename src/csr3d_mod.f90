
module csr3d_mod

use, intrinsic :: iso_fortran_env
use elliptic_integral_mod 
use fft_mod 

implicit none

integer, parameter, private :: dp = REAL64

contains




subroutine csr3d_calc()
real(dp) :: phi, m, mc, elb, eld
real(dp) :: ellipkinc, ellipeinc

print *, 'Done!'



phi = 0.1
m = 0.5

call ellipinc(phi, m, ellipkinc, ellipeinc)
print *, 'ellipinc(phi, m, ellipkinc, ellipeinc)', phi, m, ellipkinc, ellipeinc

print *, ' ------ negative m test ---------' 
phi = 0.1761732772710074
m = -2001999.9999999998

call ellipinc(phi, m, ellipkinc, ellipeinc)
print *, 'ellipinc(phi, m, ellipkinc, ellipeinc)', phi, m, ellipkinc, ellipeinc


!     Inputs: phi = amplitude, mc = complementary parameter, 0 <= m < 1
!
!     Outputs: elb = B(phi|m), eld = D(phi|m)





end subroutine









end module


