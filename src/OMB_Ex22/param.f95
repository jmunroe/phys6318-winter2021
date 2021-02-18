MODULE param

INTEGER, PARAMETER :: nx = 101	
INTEGER, PARAMETER :: ny = 51	
INTEGER, PARAMETER :: nz = 2
REAL, PARAMETER :: PI = 3.1416

REAL :: htotal(0:ny+1,0:nx+1) 
REAL :: hzero(nz,0:ny+1,0:nx+1), h(nz,0:ny+1,0:nx+1)
REAL :: dp(0:nz,0:ny+1,0:nx+1) 
REAL :: eta(nz+1,0:ny+1,0:nx+1), etan(nz,0:ny+1,0:nx+1)
REAL :: eta0(nz+1,0:ny+1,0:nx+1)  
REAL :: dhdt(nz,0:ny+1,0:nx+1)
REAL :: u(nz,0:ny+1,0:nx+1), un(nz,0:ny+1,0:nx+1)
REAL :: v(nz,0:ny+1,0:nx+1), vn(nz,0:ny+1,0:nx+1)
REAL :: T(nz,0:ny+1,0:nx+1), TN(nz,0:ny+1,0:nx+1)
REAL :: dt, dx, dy, g, r
REAL :: xh(0:nx + 1), yh(0:ny + 1)
REAL :: xu(0:nx + 1), yu(0:ny + 1)
REAL :: xv(0:nx + 1), yv(0:ny + 1)
REAL :: f(0:ny+1), ad
REAL :: rho(0:nz)
REAL :: hmin, time

INTEGER :: i,j,k
INTEGER :: wet(nz,0:ny+1,0:nx+1)

REAL :: CuP(0:ny+1,0:nx+1), CuN(0:ny+1,0:nx+1)
REAL :: CvP(0:ny+1,0:nx+1), CvN(0:ny+1,0:nx+1)
REAL :: Cu(0:ny+1,0:nx+1), Cv(0:ny+1,0:nx+1)
REAL :: B(0:ny+1,0:nx+1), BN(0:ny+1,0:nx+1)

INTEGER :: mode

END MODULE param
