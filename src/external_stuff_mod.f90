module external_stuff_mod

implicit none

contains

subroutine ccfft3d(a,b,idir,n1,n2,n3,iskiptrans)
use, intrinsic :: iso_c_binding
use omp_lib
implicit none
include 'fftw3.f03'

integer :: idir(3)
type(C_PTR) :: plan
complex(C_DOUBLE_COMPLEX), dimension(:,:,:) ::a, b
integer :: dir, fdir, n1, n2, n3, iskiptrans, n_threads

if (idir(1) == 1) then
  fdir = FFTW_BACKWARD
else
  fdir = FFTW_FORWARD
endif

print *, 'fftw_execute_dft...'
!$ n_threads = omp_get_max_threads()
!$ if (n_threads > 1) then
!$   print *, 'n_threads: ', n_threads
!$   call fftw_plan_with_nthreads(n_threads)
!$ endif

plan = fftw_plan_dft_3d(n3,n2,n1, a,b, fdir,FFTW_ESTIMATE)
call fftw_execute_dft(plan,a, b)
call fftw_destroy_plan(plan)

!$ if (n_threads > 1) call fftw_cleanup_threads()

print *, '...done'

end subroutine

end module