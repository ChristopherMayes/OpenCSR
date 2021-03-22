
module csr3d_mod

use, intrinsic :: iso_fortran_env
use external_stuff_mod ! For FFTW, GSL special functions

implicit none

integer, parameter, private :: dp = REAL64

contains

subroutine csr3d_calc()
print *, 'Done!'
end subroutine

end module


