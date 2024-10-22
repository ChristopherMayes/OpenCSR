program test_csr3d

use, intrinsic :: iso_fortran_env
use csr3d_mod

implicit none
integer, parameter :: dp = REAL64

!type mesh3d_struct
!  integer :: size(3) = [64, 64, 64]     ! Grid shape
!  real(dp) :: min(3)                    ! Minimim in each dimension
!  real(dp) :: max(3)                    ! Maximum in each dimension
!  real(dp) :: delta(3)                  ! Grid spacing
!  real(dp) :: gamma                     ! Relativistic gamma
!  real(dp) :: charge                    ! Total charge on mesh
!  real(dp), allocatable, dimension(:,:,:) :: density        ! Charge density grid
!  real(dp), allocatable, dimension(:,:,:,:) :: wake        ! wake grid
!end type

! Compatibility with OpenSpaceCharge
type mesh3d_struct
  integer :: nlo(3) = [ 1,  1,  1]       ! Lowest  grid index in x, y, z (m) of rho and the quantity being computed (phi or E)
  integer :: nhi(3) = [64, 64, 64]       ! Highest grid index in x, y, z (m) of rho and the quantity being computed (phi or E)
  integer :: npad(3) = [ 1,  1,  1]      ! Array padding for cyclic convolution
  real(dp) :: min(3)                    ! Minimim in each dimension
  real(dp) :: max(3)                    ! Maximum in each dimension
  real(dp) :: delta(3)                  ! Grid spacing
  real(dp) :: gamma                     ! Relativistic gamma
  real(dp) :: charge                    ! Total charge on mesh
  real(dp), allocatable, dimension(:,:,:) :: rho        ! Charge density grid
  real(dp), allocatable, dimension(:,:,:) :: phi        ! electric potential grid
  real(dp), allocatable, dimension(:,:,:,:) :: efield   ! electric field grid
  real(dp), allocatable, dimension(:,:,:,:) :: bfield   ! magnetic field grid
end type



type(mesh3d_struct) :: mesh3d




real(dp) :: gamma, rho, sigmas(3), gauss_cutoff
real(dp) :: center(3), x0, y0, z0, dummy
integer ::  sizes(3)
logical :: normalize

! File handling
character(200) :: in_file
integer :: open_status, namelist_file
namelist / csr3d_test_params / &
  sizes, gamma, rho, &
  sigmas, gauss_cutoff, center, &
  normalize

!--------------
!Namelist defaults
gamma = 500
rho = 1
sizes = [4,4,4]
sigmas = [10.0e-6_dp, 10.0e-6_dp, 10.0e-6_dp]
gauss_cutoff = 5.0
center = [0.0_dp, 0.0_dp, 0.0_dp]
normalize = .true.  ! Use normalized units 1/m^2, otherwise V/m

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

!
!x0 = -1e-6_dp
!y0 = -1e-6_dp
!z0 = -1e-6_dp
!gamma = 500.0_dp
!
!
!print *, x0, y0, z0, gamma
!
!!dummy = alpha(x0, y0, z0, gamma)
!!print *, 'Alpha: ', dummy
!
!dummy = psi(x0, y0, z0, gamma, 2)
!print *, 'psi_y: ', dummy
!
!stop


mesh3d%nhi = sizes
print *, 'Gaussian mesh'
call gaussian_mesh(mesh3d, sigmas, gauss_cutoff)

!call calc_density_derivative(mesh3d%rho, mesh3d%rho_prime, mesh3d%delta(3))

!print *, '------------------------'
!print *, 'Writing 2d density'
!call write_2d('density0', mesh3d%rho(:,sizes(2)/2,:))
!print *, 'sizes(2)/2', sizes(2)/2
!call write_2d('density0_prime', mesh3d%rho_prime(:,sizes(2)/2,:))
!


mesh3d%gamma = gamma
call print_mesh3d(mesh3d)


print *, '------------------------'
print *, 'csr3d_steady_state'
call csr3d_steady_state(mesh3d, rho)


print *, '------------------------'
print *, 'Writing lines, planes'
!call write_lines(mesh3d, center(1), center(2), center(3))

call write_line('x_line.dat', mesh3d, center, 'x')
call write_line('y_line.dat', mesh3d, center, 'y')
call write_line('z_line.dat', mesh3d, center, 'z')

center = 0
center(1) = -sigmas(1)
call write_line('z_line_-sigma_x.dat', mesh3d, center, 'z')
center(1) = sigmas(1)
call write_line('z_line_+sigma_x.dat', mesh3d, center, 'z')

center = 0
center(2) = -sigmas(2)
call write_line('z_line_-sigma_y.dat', mesh3d, center, 'z')
center(2) = sigmas(2)
call write_line('z_line_+sigma_y.dat', mesh3d, center, 'z')

!call write_plane(mesh3d)

center = 0
call write_plane('xz_plane.dat', mesh3d, center, 'xz')
center(2) = sigmas(2)
call write_plane('xz_plane_+sigma_y.dat', mesh3d, center, 'xz')
center(2) = -sigmas(2)
call write_plane('xz_plane_-sigma_y.dat', mesh3d, center, 'xz')


!call csr3d_calc()
!call ellipinc_test()






contains






! Fills mesh with Gaussian data
subroutine gaussian_mesh(mesh3d, sigmas, cutoff)
type(mesh3d_struct) :: mesh3d
real(dp) :: x, y, z, sigmas(3), cutoff
real(dp) :: dx, dy, dz, xmin, ymin, zmin, norm
integer :: i, j, k, isize, jsize, ksize

if (.not. allocated(mesh3d%rho)) then
  allocate(mesh3d%rho(mesh3d%nhi(1), mesh3d%nhi(2), mesh3d%nhi(3)))
  allocate(mesh3d%efield(mesh3d%nhi(1), mesh3d%nhi(2), mesh3d%nhi(3), 3))   
endif


! Fill mesh
mesh3d%min   = -cutoff*sigmas
mesh3d%max   =  cutoff*sigmas
mesh3d%delta = (mesh3d%max - mesh3d%min)/(mesh3d%nhi-1)

! Convenient local variables. 
isize = mesh3d%nhi(1)
jsize = mesh3d%nhi(2)
ksize = mesh3d%nhi(3)

dx = mesh3d%delta(1)
dy = mesh3d%delta(2)
dz = mesh3d%delta(3)
xmin = mesh3d%min(1)
ymin = mesh3d%min(2)
zmin = mesh3d%min(3)

do k = 1, ksize
  z = (k-1)*dz + zmin

  do j=1, jsize
    y=(j-1)*dy + ymin

   do i=1, isize
     x = (i-1)*dx + xmin
     
     mesh3d%rho(i,j,k)= gauss3(x, y, z, sigmas(1), sigmas(2), sigmas(3))
    ! mesh3d%rho_prime(i,j,k)= gauss3_prime(x, y, z, sigmas(1), sigmas(2), sigmas(3))
     
    enddo
  enddo
enddo

! Normalize
norm = sum(mesh3d%rho) * dx * dy * dz
mesh3d%rho  = mesh3d%rho / norm


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


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine interpolate_field(x, y, z, mesh3d, wake)
!
! Interpolate field on mesh
!
! Input:
!   x, y, z   -- REAL64: coordinates to interpolate
!   mesh3d    -- mesh3d_struct:  contains efield, bfield
!
! Output:
!   wake(3)      -- REAL64 : interpolated wake field at x, y, z   
!
!-

subroutine interpolate_field(x, y, z, mesh3d, wake, density)
type(mesh3d_struct) ::  mesh3d
real(dp) :: x, y, z
real(dp), optional :: wake(3), density
real(dp) :: hxi,hyi,hzi,ab,de,gh
integer :: ip,jp,kp,ip1,jp1,kp1
integer :: nflag
nflag=0
hxi=1.d0/mesh3d%delta(1); hyi=1.d0/mesh3d%delta(2); hzi=1.d0/mesh3d%delta(3)      
ip=floor((x-mesh3d%min(1))*hxi+1)
jp=floor((y-mesh3d%min(2))*hyi+1)
kp=floor((z-mesh3d%min(3))*hzi+1)

if(ip<1 .or. ip>mesh3d%nhi(1)-1)then
  nflag=1
  write(6,*)'ierror: ip=', ip, (x-mesh3d%min(1)), (x-mesh3d%min(1))/mesh3d%delta(1)
  if(ip<1)then
    ip=1
  else
    ip=mesh3d%nhi(1)-1
  endif
endif
if(jp<1 .or. jp>mesh3d%nhi(2)-1)then
  nflag=1
  write(6,*)'jerror: jp=', jp, (y-mesh3d%min(2)), (y-mesh3d%min(2))/mesh3d%delta(2)
  write(6,*)ab,de,gh
  if(jp<1)then
    jp=1
  else
    jp=mesh3d%nhi(2)-1
  endif
endif
  
if(kp<1 .or. kp>mesh3d%nhi(3)-1)then
    nflag=1
    write(6,*)'kerror:  kp=',kp
    write(6,*)'z=',z
    write(6,*)'mesh3d%min(3)=',mesh3d%min(3)
!!!!!!!!!!write(6,*)ab,de,gh
  if(kp<1)then
    kp=1
  else
    kp=mesh3d%nhi(3)-1
  endif
endif
ab=((mesh3d%min(1)-x)+ip*mesh3d%delta(1))*hxi
de=((mesh3d%min(2)-y)+jp*mesh3d%delta(2))*hyi
gh=((mesh3d%min(3)-z)+kp*mesh3d%delta(3))*hzi
if(nflag.eq.1)then
  write(6,*)ab,de,gh
  nflag=0
endif

ip1=ip+1
jp1=jp+1
kp1=kp+1

if (present(wake)) then
    wake=mesh3d%efield(ip, jp,  kp,  :)*ab*de*gh                 &
     +mesh3d%efield(ip, jp1, kp,  :)*ab*(1.-de)*gh            &
     +mesh3d%efield(ip, jp1, kp1, :)*ab*(1.-de)*(1.-gh)       &
     +mesh3d%efield(ip, jp,  kp1, :)*ab*de*(1.-gh)            &
     +mesh3d%efield(ip1,jp,  kp1, :)*(1.-ab)*de*(1.-gh)       &
     +mesh3d%efield(ip1,jp1, kp1, :)*(1.-ab)*(1.-de)*(1.-gh)  &
     +mesh3d%efield(ip1,jp1, kp,  :)*(1.-ab)*(1.-de)*gh       &
     +mesh3d%efield(ip1,jp,  kp,  :)*(1.-ab)*de*gh
endif

if (present(density)) then
    density=mesh3d%rho(ip, jp,  kp  )*ab*de*gh                 &
          +mesh3d%rho(ip, jp1, kp  )*ab*(1.-de)*gh            &
          +mesh3d%rho(ip, jp1, kp1 )*ab*(1.-de)*(1.-gh)       &
          +mesh3d%rho(ip, jp,  kp1 )*ab*de*(1.-gh)            &
          +mesh3d%rho(ip1,jp,  kp1 )*(1.-ab)*de*(1.-gh)       &
          +mesh3d%rho(ip1,jp1, kp1 )*(1.-ab)*(1.-de)*(1.-gh)  &
          +mesh3d%rho(ip1,jp1, kp  )*(1.-ab)*(1.-de)*gh       &
          +mesh3d%rho(ip1,jp,  kp  )*(1.-ab)*de*gh
endif


end subroutine



! -----------------
! -----------------
!! Write line utility
!!
subroutine write_line(fname, mesh3d, center, axis)
type(mesh3d_struct) :: mesh3d
real(dp) :: center(3), x, y, z, wake(3)
integer :: i, outfile
character(*) :: fname
character(1) :: axis

open(newunit=outfile, file = trim(fname))

select case(axis)
case('x')
    z = center(3)
    y = center(2)
    do i = 1, mesh3d%nhi(1) -1 ! skip last point
      x = (i-1)*mesh3d%delta(1) + mesh3d%min(1) 
      call interpolate_field(x, y, z, mesh3d, wake=wake)
      write(outfile, '(6(1pe14.7,1x))') x,y,z, wake(1:3)
    enddo  

case('y')
x = center(1)
    z = center(3)
    do i = 1, mesh3d%nhi(2) -1 ! skip last point
      y = (i-1)*mesh3d%delta(2) + mesh3d%min(2) 
      call interpolate_field(x, y, z, mesh3d, wake=wake)
      write(outfile, '(6(1pe14.7,1x))') x,y,z, wake(1:3)
    enddo  

case('z')
    x = center(1)
    y = center(2)
    do i = 1, mesh3d%nhi(3) -1 ! skip last point
      z = (i-1)*mesh3d%delta(3) + mesh3d%min(3) 
      call interpolate_field(x, y, z, mesh3d, wake=wake)
      write(outfile, '(6(1pe14.7,1x))') x,y,z,wake(1:3)
    enddo
end select

close(outfile) 
end subroutine


! -----------------
! -----------------
!! Write plane utility
!!
subroutine write_plane(fname, mesh3d, center, axes)
type(mesh3d_struct) :: mesh3d
real(dp) :: center(3), x, y, z, wake(3), density
integer :: i, k, outfile
character(*) :: fname
character(2) :: axes

open(newunit=outfile,  file = trim(fname))

write(outfile, '(6(a, 1x))') 'x', 'z', 'Wx', 'Wy', 'Ws', 'density'

select case(axes)
case('xz')
y = center(2)
do k = 1, mesh3d%nhi(3) -1 ! skip last point
  z = (k-1)*mesh3d%delta(3) + mesh3d%min(3) 
  do i = 1, mesh3d%nhi(1) -1 ! skip last point
    x = (i-1)*mesh3d%delta(1) + mesh3d%min(1) 
    call interpolate_field(x, y, z, mesh3d, wake=wake, density=density)
    write(outfile, *) x, z, wake(1), wake(2), wake(3), density
  enddo
enddo  


end select

close(outfile)  
end subroutine


! -----------------
! -----------------
!! Write grid utility
!!
!subroutine write_2d(fname, grid)
!real(dp) :: grid(:,:)
!integer :: i, j, outfile
!character(*) :: fname
!
!open(newunit=outfile, file = trim(fname))
!
!write(outfile, *)  size(grid, 1)
!write(outfile, *)  size(grid, 2)
!do j = 1, size(grid, 2)
!  do i = 1, size(grid, 1)
!    write(outfile, *) grid(i, j)
!  enddo
!enddo  
!close(outfile)  
!end subroutine
!! -----------------
! -----------------




!------------------------------------------------------------------------
!+
subroutine print_mesh3d(mesh3d)
type(mesh3d_struct) :: mesh3d
print *, '------------------------'
print *, 'Mesh: '
print *, 'size: ', mesh3d%nhi
print *, 'min: ', mesh3d%min
print *, 'max: ', mesh3d%max
print *, 'delta: ', mesh3d%delta
print *, 'gamma: ', mesh3d%gamma
print *, 'charge: ', mesh3d%charge
if (allocated(mesh3d%rho)) print *, 'density allocated'
if (allocated(mesh3d%efield)) print *, 'wake allocated'
print *, '------------------------'
end subroutine


!------------------------------------------------------------------------
!+
! Subroutine csr3d_steady_state(mesh3d, offset)
!
! Performs the space charge calculation using the integrated Green function method
! and FFT-based convolutions. 
!
! Input:
!   mesh3d        -- mesh3d_struct: populated with %rho
!
!   offset        -- real(3), optional: Offset coordinates x0, y0, z0 to evaluate the field,
!                    relative to rho. 
!                    Default: (0,0,0)
!                        For example, an offset of (0,0,10) can be used to compute 
!                        the field at z=+10 m relative to rho. 
!
!
!
! Output:
!   mesh3d        -- mesh3d_struct: populated with %wake                            
!
!
!
! Notes: 
!
!     
!
!-
subroutine csr3d_steady_state(mesh3d, rho, offset, normalize)
type(mesh3d_struct), intent(inout):: mesh3d
real(dp), intent(in) :: rho
logical, optional :: normalize
real(dp), allocatable, dimension(:,:,:,:) :: image_efield   ! electric field grid
real(dp), optional :: offset(3)
real(dp) :: offset0(3)
real(dp), parameter :: c_light = 299792458.0


print *, 'csr3d_steady_state'

! Free space field
call csr3d_steady_state_solver(mesh3d%rho, mesh3d%gamma, rho, mesh3d%delta, &
   mesh3d%efield, offset=offset0, normalize=normalize)


end subroutine csr3d_steady_state




end program