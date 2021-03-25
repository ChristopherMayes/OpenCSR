program test_csr3d

use, intrinsic :: iso_fortran_env
use csr3d_mod

implicit none

type(mesh3d_struct) :: mesh3d
integer, parameter :: dp = REAL64



real(dp) :: gamma, rho, sigmas(3), gauss_cutoff
integer ::  sizes(3)


! File handling
character(200) :: in_file
integer :: open_status, namelist_file
namelist / csr3d_test_params / &
  sizes, gamma, rho, &
  sigmas, gauss_cutoff

!Namelist defaults
gamma = 500
rho = 1
sizes = [64, 32, 128]
sigmas = [10.0e-3_dp, 10.0e-3_dp, 10.0e-3_dp]
gauss_cutoff = 5.0

!Read namelist
in_file = 'test_opensc.in'
if (command_argument_count() > 0) call get_command_argument(1, in_file)
open(newunit=namelist_file, file = in_file, status = 'old', iostat=open_status)
if (open_status /= 0) then
  print *, 'Input file missing: ', in_file
  print *, 'Using defaults'
else 
  read (namelist_file, nml = csr3d_test_params)
  close (namelist_file)
endif

print *, '------------------------'
write(*, csr3d_test_params)
print *, '------------------------'



mesh3d%size = sizes
call gaussian_mesh(mesh3d, sigmas, gauss_cutoff)


call csr3d_calc()
call ellipinc_test()




contains






! Fills mesh with Gaussian data
subroutine gaussian_mesh(mesh3d, sigmas, cutoff)
type(mesh3d_struct) :: mesh3d
real(dp) :: x, y, z, sigmas(3), cutoff
real(dp) :: dx, dy, dz, xmin, ymin, zmin, norm
integer :: i, j, k, isize, jsize, ksize

if (.not. allocated(mesh3d%density)) then
  allocate(mesh3d%density(mesh3d%size(1), mesh3d%size(2), mesh3d%size(3)))
  allocate(mesh3d%density_prime(mesh3d%size(1), mesh3d%size(2), mesh3d%size(3)))
  allocate(mesh3d%wake(mesh3d%size(1), mesh3d%size(2), mesh3d%size(3), 3))   
endif


! Fill mesh
mesh3d%min   = -cutoff*sigmas
mesh3d%max   =  cutoff*sigmas
mesh3d%delta = (mesh3d%max - mesh3d%min)/(mesh3d%size-1)

! Convenient local variables. 
isize = mesh3d%size(1)
jsize = mesh3d%size(2)
ksize = mesh3d%size(3)

dx = mesh3d%delta(1)
dy = mesh3d%delta(2)
dz = mesh3d%delta(3)
xmin = mesh3d%min(1)
ymin = mesh3d%min(2)
zmin = mesh3d%min(3)

do k = 0, ksize-1
  z = k*dz + zmin

  do j=0, jsize-1
    y=j*dy + ymin

   do i=0, isize-1
     x = i*dx + xmin
     
     mesh3d%density(i,j,k)= gauss3_prime(x, y, z, sigmas(1), sigmas(2), sigmas(3))
     
    enddo
  enddo
enddo

! Normalize
norm = sum(mesh3d%density) * dx * dy * dz
mesh3d%density  = mesh3d%density / norm


end subroutine 

! Gaussian
elemental real(dp) function gauss(x, sigma)
real(dp), intent(in):: x, sigma
real(dp), parameter :: pi = 4.0_dp*atan2(1.0_dp,1.0_dp)
gauss =  1/(sigma * sqrt(2 * pi)) * exp( - (x)**2 / (2 * sigma**2) )
end function

! 3D Gaussian
elemental real(dp) function gauss3(x, y, z, sigma_x, sigma_y, sigma_z)
real(dp), intent(in):: x, y, z, sigma_x, sigma_y, sigma_z
real(dp), parameter :: pi = 4.0_dp*atan2(1.0_dp,1.0_dp)
gauss3 =  gauss(x, sigma_x)*gauss(y, sigma_y)*gauss(z, sigma_z)
end function

! 3D d/dz Gaussian
elemental real(dp) function gauss3_prime(x, y, z, sigma_x, sigma_y, sigma_z)
real(dp), intent(in):: x, y, z, sigma_x, sigma_y, sigma_z
real(dp), parameter :: pi = 4.0_dp*atan2(1.0_dp,1.0_dp)
gauss3_prime =  gauss(x, sigma_x)*gauss(y, sigma_y)*gauss(z, sigma_z)*(-z/sigma_z**2)
end function



end program