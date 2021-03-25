
module csr3d_mod

use, intrinsic :: iso_fortran_env
use elliptic_integral_mod, only : ellipinc 
use fft_mod, only : ccfft3d


implicit none

! Fortran 2008
!integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
!integer, parameter, private :: qp = REAL128

type mesh3d_struct
  integer :: size(3) = [64, 64, 64]     ! Grid shape
  real(dp) :: min(3)                    ! Minimim in each dimension
  real(dp) :: max(3)                    ! Maximum in each dimension
  real(dp) :: delta(3)                  ! Grid spacing
  real(dp) :: gamma                     ! Relativistic gamma
  real(dp) :: rho                       ! bending radius (positive)
  real(dp) :: charge                    ! Total charge on mesh
  real(dp), allocatable, dimension(:,:,:) :: density        ! Charge density grid
  real(dp), allocatable, dimension(:,:,:) :: density_prime        ! d/dz density
  real(dp), allocatable, dimension(:,:,:,:) :: wake        ! wake grid

end type





!------------------------------------------------------------------------
!------------------------------------------------------------------------
contains




!------------------------------------------------------------------------
!+
subroutine print_mesh3d(mesh3d)
type(mesh3d_struct) :: mesh3d
print *, '------------------------'
print *, 'Mesh: '
print *, 'size: ', mesh3d%size
print *, 'min: ', mesh3d%min
print *, 'max: ', mesh3d%max
print *, 'delta: ', mesh3d%delta
print *, 'gamma: ', mesh3d%gamma
print *, 'rho: ', mesh3d%rho
print *, 'charge: ', mesh3d%charge
if (allocated(mesh3d%density)) print *, 'density allocated'
if (allocated(mesh3d%density_prime)) print *, 'density_prime allocated'
if (allocated(mesh3d%wake)) print *, 'wake allocated'
print *, '------------------------'
end subroutine




subroutine csr3d_calc()


print *, 'Done!'
end subroutine










subroutine ellipinc_test()
real(dp) :: phi, m, mc, elb, eld
real(dp) :: ellipkinc, ellipeinc
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
!   at_cathode    -- logical, optional: Maintain constant voltage at the cathode 
!                       using image charges. Default is False. 
!
!   calc_bfield    -- logical, optional: Calculate the magnetic field mesh3d%bfield
!
!                     Default: False
!
!
! Output:
!   mesh3d        -- mesh3d_struct: populated with %efield, and optionally %bfield                             
!
!
!
! Notes: 
!   The magnetic field components can be calculated by:
!     Bx = -(beta/c) * Ey
!     By =  (beta/c) * Ex
!     Bz = 0
!   The image charges move in the opposite direction, so the signs are flipped. 
!     
!
!-
subroutine csr3d_steady_state(mesh3d, offset)
type(mesh3d_struct) :: mesh3d
real(dp), allocatable, dimension(:,:,:,:) :: image_efield   ! electric field grid
real(dp), optional :: offset(3)
real(dp) :: offset0(3), beta
real(dp), parameter :: c_light = 299792458.0

if (.not. present(offset)) then
  offset0 = 0
else
  offset0 = offset
endif

! Free space field
call csr3d_steady_state_solver(mesh3d%density, mesh3d%gamma, mesh3d%rho, mesh3d%delta, mesh3d%wake, offset=offset0)


end subroutine csr3d_steady_state




!------------------------------------------------------------------------
!+
! Subroutine calc_density_derivative(density, density_prime, dz)
!
! Input:
!
!
! Output:
!
!-
subroutine calc_density_derivative(density, density_prime, dz)
real(dp), intent(in), dimension(:,:,:) :: density
real(dp), intent(in) :: dz
complex(dp), allocatable, dimension(:,:,:), intent(out) :: density_prime
integer :: nz

! Sizes
nz = size(density, 3)

! Allocate double-sized arrays. This array will be placed in the first octant
allocate( density_prime( 2*size(density, 1), 2*size(density, 2), 2*size(density, 3) ))
density_prime = 0 ! Zero out everything

! Endpoints: Simple backward and forward differences
density_prime(:,:,1)  = (density(:,:,2)  - density(:,:,1))/dz
density_prime(:,:,nz) = (density(:,:,nz) - density(:,:,nz-1))/dz

! Central differences
density_prime(:,:,2:nz-1) =  (density(:,:,3:nz) - density(:,:,1:nz-2))/(2*dz)
!
end subroutine


!------------------------------------------------------------------------
!+
! Subroutine csr3d_steady_state_solver(rho, gamma, delta, efield, phi, offset)
!
! Deposits particle arrays onto mesh
!
! Input:
!   rho          -- REAL64(:,:,:): charge density array in x, y, z
!   delta        -- REAL64(3): vector of grid spacings dx, dy, dz
!   gamma        -- REAL64: relativistic gamma
!   icomp        -- integer: Field component requested:
!                        0: phi (scalar potential)
!                       
!
!   efield        -- REAL64(:,:,:,3), optional: allocated electric field array to populate.
!                      
!                                     The final index corresponds to components
!                                     1: Ex
!                                     2: Ey
!                                     3: Ez                                   
!                                     If present, all components will be computed.    
!
!   phi           -- REAL64(:,:,:), optional: allocated potential array to populate
!
!   offset        -- real(3), optional: Offset coordinates x0, y0, z0 to evaluate the field,
!                    relative to rho. 
!                    Default: (0,0,0)
!                        For example, an offset of (0,0,10) can be used to compute 
!                        the field at z=+10 m relative to rho. 
!
! Output:
!   efield        -- REAL64(:,:,:,:) : electric field                                 
!   phi           -- REAL64(:,:,:)   : potential
!
!
! Notes: 
!
!-
subroutine csr3d_steady_state_solver(density, gamma, rho, delta, wake, offset)


real(dp), intent(in), dimension(:,:,:) :: density
real(dp), intent(in) :: gamma, rho, delta(3)
real(dp), optional, intent(out), dimension(:,:,:,:) :: wake
real(dp), intent(in), optional :: offset(3)
! internal arrays
complex(dp), allocatable, dimension(:,:,:) :: complex_density_prime, cgrn
real(dp) :: factr, offset0=0
real(dp), parameter :: clight=299792458.0
real(dp), parameter :: fpei=299792458.0**2*1.00000000055d-7  ! this is 1/(4 pi eps0) after the 2019 SI changes

integer :: nx, ny, nz, nx2, ny2, nz2
integer :: icomp, ishift, jshift, kshift

! Sizes
nx = size(density, 1); ny = size(density, 2); nz = size(density, 3)
nx2 = 2*nx; ny2 = 2*ny; nz2 = 2*nz; 

! Allocate complex scratch arrays
allocate(cgrn(nx2, ny2, nz2))

! density -> d/dz density -> FFT 
! This will calculate d/dz density and place in the complex array
call calc_density_derivative(density, complex_density_prime, delta(3))
call ccfft3d(complex_density_prime, complex_density_prime, [1,1,1], nx2, ny2, nz2, 0) 

! Loop over Gx, Gy, Gs
do icomp=1, 3

  call get_cgrn_csr3d(cgrn, delta, gamma, rho, icomp, offset=offset)
  
  !  cgrn -> FFT(cgrn)
  call ccfft3d(cgrn, cgrn, [1,1,1], nx2, ny2, nz2, 0)  
  
  ! Multiply FFT'd arrays, re-use cgrn
  cgrn=complex_density_prime*cgrn

  ! Inverse FFT
  !call ccfft3d(cgrn, cgrn, [-1,-1,-1], nx2, ny2, nz2, 0)  
  
  ! This is where the output is shifted to
  ishift = nx-1
  jshift = ny-1
  kshift = nz-1
  
  ! Extract field
  ! ???? Factor
  factr = fpei/(nx2*ny2*nz2)
  
  wake(:,:,:,icomp) = factr * real(cgrn(1+ishift:nx+ishift, 1+jshift:ny+jshift, 1+kshift:nz+kshift), dp)

    
enddo

end subroutine csr3d_steady_state_solver


!------------------------------------------------------------------------
!+
! Subroutine get_cgrn_csr3d(cgrn, delta, gamma, rho, icomp, offset)
!
! Computes the free space Green function on a mesh with given spacings in the lab frame.
! The computation is performed in the rest fram by boosting the coordinates by gamma.
!
!
! Input:
!   cgrn         -- COMPLEX128(:,:,:): pre-allocated array 
!   delta        -- REAL64(3): vector of grid spacings dx, dy, dz
!   gamma        -- REAL64: relativistic gamma
!   rho          -- REAL64: bending radius (must be positive)
!   icomp        -- integer: Wake component requested:
!
!                        1: Wx
!                        2: Wy
!                        3: Ws
!   offset        -- real(3), optional: Offset coordinates for the center of the grid in [m]. 
!                    Default: (0,0,0)
!                        For example, an offset of (0,0,10) can be used to compute 
!                        the field at z=+10 m relative to the rho_mesh center. 
!                              
! Output:
!   cgrn         -- COMPLEX128(:,:,:): Green function array
!   
!                
! Notes:
!   Internally, mesh spacings in scaled units
!   x - > x/rho
!   y ->  y/rho
!   z ->  z/(2*rho)   
!
!
!-
subroutine get_cgrn_csr3d(cgrn, delta, gamma, rho, icomp, offset)

complex(dp), intent(out), dimension(0:,0:,0:) :: cgrn ! Convenient indexing 
real(dp), intent(in), dimension(3) :: delta
integer, intent(in) :: icomp
real(dp), intent(in) :: gamma, rho
real(dp), intent(in), optional :: offset(3)
! Local
real(dp) :: dx,dy,dz
real(dp) :: u,v,w, umin, vmin, wmin
real(dp) :: gval, factor
integer :: imin, imax, jmin, jmax, kmin, kmax
integer :: i,j,k, isize, jsize, ksize

! Mesh spacings in scaled units
! x - > x/rho
! y ->  y/rho
! z ->  z/(2*rho)

dx=delta(1)/rho; dy=delta(2)/rho; dz=delta(3)/(2*rho)

! ????
factor = 1.0 /(dx*dy*dz)

! Grid min
isize = size(cgrn,1); jsize=size(cgrn,2); ksize=size(cgrn,3)
umin = (1-isize/2) *dx
vmin = (1-jsize/2) *dy
wmin = (1-ksize/2) *dz

! Add optional offset
if (present(offset)) then
  umin = umin + offset(1)/rho
  vmin = vmin + offset(2)/rho
  wmin = wmin + offset(3)/(2*rho)
endif

! !$ print *, 'OpenMP Green function calc osc_get_cgrn_freespace'
!$OMP PARALLEL DO &
!$OMP DEFAULT(FIRSTPRIVATE), &
!$OMP SHARED(cgrn)
do k = 0, ksize-1
  w = k*dz + wmin

  do j=0, jsize-1
    v=j*dy + vmin

   do i=0, isize-1
     u = i*dx + umin
     
     !gval=psi(u, w, v, gamma, icomp)*factor
     !cgrn(i,j,k)= cmplx(gval, 0, dp)
     
     cgrn(i,j,k)=  psi(u, w, v, gamma, icomp)*factor
     
    enddo
  enddo
enddo
!$OMP END PARALLEL DO



end subroutine get_cgrn_csr3d





!------------------------------------------------------------------------
!+
! elemental real(dp) function psi_s(x, y, z, gamma)
!
!
! Eq. 24 from Ref[X] without the prefactor e beta^2 / (2 rho^2)
!-
elemental real(dp) function psi_s(x, y, z, gamma)
real(dp), intent(in) :: x, y, z, gamma
real(dp) :: beta, beta2, kap, alp
    
beta2 = 1-1/gamma*2
beta = sqrt(beta2)    
    
! Check for the origin
if ((x == 0) .and. (y == 0) .and. (z == 0)) then
    psi_s = 0
else
    beta2 = beta**2
    
    alp = alpha(x, y, z, beta2)
    kap = 2*(alp - z)/beta ! Simpler form of kappa
    !kap = sqrt(x**2 + y**2 + 4*(1+x) * sin(alp)**2) 
    
    psi_s =  (cos(2*alp) - 1/(1+x)) / (kap - beta*(1+x)*sin(2*alp) )
endif
end function psi_s


!------------------------------------------------------------------------
!+
! elemental real(dp) function psi_x(x, y, z, gamma)
!
!
! Eq. 24 from Ref[X] without the prefactor e beta^2 / (2 rho^2)
! This is actually psi_x_hat
!
!
!
! icomp : integer
!          1: x
!          2: y
!          3: s
!-
elemental real(dp) function psi(x, y, z, gamma, icomp)
real(dp), intent(in) :: x, y, z, gamma
integer, intent(in) :: icomp
real(dp) :: beta, beta2, kap, alp
real(dp) :: sin2a, cos2a, kap2, sin2a2, x2, y2, y4, xp, xp2, xy, xy2, f1, f2, arg2, f, e

beta2 = 1-1/gamma*2
beta = sqrt(beta2)

alp = alpha(x, y, z, gamma)
kap = 2*(alp - z)/beta ! Simpler form of kappa
!kap = sqrt(x**2 + y**2 + 4*(1+x) * sin(alp)**2) 

if (icomp == 3) then
    psi = (cos(2*alp) - 1/(1+x)) / (kap - beta*(1+x)*sin(2*alp) )
    return
endif

! Common patterns
sin2a = sin(2*alp)
cos2a = cos(2*alp)

kap2 = kap**2
sin2a2 = sin2a**2

x2 = x**2 
y2 = y**2
y4 = y2**2
xp = x + 1
xp2 = xp**2
xy2 = x2 + y2
xy = sqrt(xy2)

! More complicated pattens
f1 = 2 + 2*x +x2
f2 = (2+x)**2
arg2 = -4 * xp / xy2 

! Elliptic integrals of the first and second kind F(phi|m), E(phi|m)
call ellipinc(alp, arg2, F, E)

if (icomp == 1) then
! psi_x    
! psi_x is actually psi_x_hat that includes the psi_phi term
! There is an extra ] in the numerator of the second term. All terms should multiply E. 

    psi = f1*F / (xp*xy) - (x2*f2 + y2*f1)*E / (xp*(y2+f2)*xy)  &
        + ( kap2 - 2*beta2*xp2 + beta2*xp*f1*cos2a  ) / (beta *xp*(kap2 - beta2*xp2*sin2a2)) &
        + kap*( y4 - x2*f2 - 2*beta2*y2*xp2 )*sin2a / ( xy2*(y2 + f2)*(kap2-beta2*xp2*sin2a2)  ) &
        + kap*beta2*xp*( x2*f2 + y2*f1 )*sin2a*cos2a / ( xy2*(y2+f2)*(kap2-beta2*xp2*sin2a2)  ) &
        - (2/beta2)* F/xy ! Include the phi term      

elseif (icomp == 2) then
! psi_y
    psi = y * ( &
        F/xy - (x*(2+x)+y2)*E / ((y2+f2)*xy) &
        - beta*(1-xp*cos2a) / (kap2-beta2*xp2*sin2a2) &
        + kap*xp*( -(2+beta2)*y2 + (-2+beta2)*x*(2+x) ) * sin2a / ( (y4 + x2*f2 + 2*y2*f1)*( kap2-beta2*xp2*sin2a2 ) ) &
        + kap*beta2*xp2*(y2 + x*(2+x))*sin2a*cos2a / ( ( y4 + x2*f2 + 2*y2*f1)*(kap2 -beta2*xp2*sin2a2)  ) &
        )        
endif


end function psi





!------------------------------------------------------------------------
!+
! elemental real(dp) function alpha(x, y, z, gamma)
!
!
!
!-
elemental real(dp) function alpha(x, y, z, gamma)
real(dp), intent(in) :: x, y, z, gamma
real(dp) :: beta2, b, c, eta, nu, zeta,  omega3, m
real(dp) :: arg1, arg2, arg3, temp

beta2 = 1-1/gamma*2

if (z == 0.0) then
    ! Quadratic solution 
    
    b = 3 * (1 - beta2 - beta2*x) / beta2 / (1+x)    
    c = -3*(x**2 + y**2)/(4*(1+x))

    alpha = sqrt(-b + sqrt(b**2 - 4*c))/2
    
else    
    !Quartic solution 
        
    ! Terms of the depressed quartic equation
    eta = -6 * z / (beta2 * (1+x))
    nu = 3 * (1/beta2 - 1 - x) / (1+x)
    zeta = (3/4) * (4* z**2 /beta2 - x**2 - y**2) / (1+x)
    
    ! Omega calc and cube root
    temp = (eta**2/16 - zeta * nu/6 + nu**3/216)  
    omega3 =  (temp + sqrt(temp**2 - (zeta/3 + nu**2/36)**3))**(1/3)
    
    ! Eq. (A2) from Ref[1]
    m = -nu/3 + (zeta/3 + nu**2/36) /omega3 + omega3
     
    arg1 = sqrt(2 * abs(m))
    arg2 = -2 * (m + nu)
    arg3 = 2 * eta / arg1
    
    if (z < 0) then
        alpha =  (arg1 + sqrt(abs(arg2 - arg3)))/2
    else
        alpha = (-arg1 + sqrt(abs(arg2 + arg3)))/2
    endif

endif
end function alpha


end module


