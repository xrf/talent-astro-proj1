subroutine init
  use global
  use zone
  implicit none

  integer :: i, j, k
  real    :: ridt, xvel, yvel, zvel, width, widthz, widthy
  real    :: xmin, xmax, ymin, ymax, zmin, zmax
  real    :: r
  real :: sun_radius=.2
  real :: sun_origin_x=.5
  real :: sun_origin_y=.5

  ! --------------------------------------------------------------------------
  ! Set up geometry and boundary conditions of grid
  !
  ! Boundary condition flags: nleft, nright
  !   0: reflecting boundary condition
  !   1: inflow/outflow boundary condition (zero gradients)
  !   2: fixed inflow boundary condition (values set by dinflo, pinflo, etc.)
  !   3: periodic (nmax + 1 = nmin; nmin - 1 = nmax)
  !
  ! Geometry flag : ngeom        |  Cartesian:
  !   0: planar                  |    gx = 0, gy = 0, gz = 0   (x, y, z)
  !   1: cylindrical radial      |  Cylindrical:
  !   2: spherical   radial  3D= {     gx = 1, gy = 3, gz = 0   (s, phi, z)
  !   3: cylindrical angle       |
  !   4: polar angle (theta)     |  Spherical:
  !   5: azimuth (phi)           |    gx = 2, gy = 4, gz = 5   (r, theta, phi)

  ! Define the computational grid...

  ngeomx = 0
  ngeomy = 0
  ngeomz = 0

  nleftx = 1
  nrightx= 1
  nlefty = 1
  nrighty= 1
  nleftz = 1
  nrightz= 1

  xmin   = 0.0
  xmax   = 1.0
  ymin   = 0.0
  ymax   = 1.0
  zmin   = 0.0
  zmax   = 1.0

  ! if any dimension is angular, multiply coordinates by pi
  if(ngeomy >= 3) then
     ymin = ymin * pi
     ymax = ymax * pi
  endif
  if(ngeomz >= 3) then
     zmin = zmin * pi
     zmax = zmax * pi
  endif

  ! --------------------------------------------------------------------------
  ! set up parameters from the problem

  gam    = 5. / 3.

  gamm = gam - 1.0

  ! --------------------------------------------------------------------------
  ! set time and cycle counters

  time   = 0.0
  timep  = 0.0
  timem  = 0.0
  ncycle = 0
  ncycp  = 0
  ncycd  = 0
  ncycm  = 0
  nfile  = 0

  ! Set up grid coordinates

  call grid(imax,xmin,xmax,zxa,zxc,zdx)
  call grid(jmax,ymin,ymax,zya,zyc,zdy)
  call grid(kmax,zmin,zmax,zza,zzc,zdz)

  if (ndim <= 2) zzc(1) = 0.0
  if (ndim == 1) zyc(1) = 0.0

  ! --------------------------------------------------------------------------
  ! Log parameters of problem in history file

  write (8, "('Oblique Sod shock tube in ',i1,' dimensions.')") ndim
  if (ndim == 1) then
     write (8, "('Grid dimensions: ',i4)") imax
  else if (ndim == 2) then
     write (8, "('Grid dimensions: ',i4,' x ',i4)") imax, jmax
  else
     write (8, "('Grid dimensions: ',i4,' x ',i4,' x ',i4)") imax, jmax, kmax
  endif
  write (8, *)
  write (8, "('Ratio of specific heats = ',f5.3)") gam
  write (8, *)

  ! initialize grid:

  do k = 1, kmax
     do j = 1, jmax
        do i = 1, imax
           r = (zxc(i) - sun_origin_x) ** 2 &
             + (zyc(j) - sun_origin_y) ** 2
!            + zzc(k) ** 2
           if (r >= sun_radius ** 2) then ! outside
              zro(i,j,k) = 1
              zpr(i,j,k) = 1
           else
              zro(i,j,k) = 10
              zpr(i,j,k) = 1
           endif
           zux(i,j,k) = 0.
           zuy(i,j,k) = 0.
           zuz(i,j,k) = 0.
           zfl(i,j,k) = 0.
        enddo
     enddo
  enddo

  do i = 1, imax
     zpr(i, 1, 1) = 1e7
  enddo

  ! --------------------------------------------------------------------------
  ! Compute Courant-limited timestep

  ridt = 0.
  if(ndim == 1) then
     do i = 1, imax
        svel = sqrt(gam*zpr(i,1,1)/zro(i,1,1))/zdx(i)
        xvel = abs(zux(i,1,1)) / zdx(i)
        ridt = max(xvel,svel,ridt)
     enddo
  else if(ndim == 2) then
     do j = 1, jmax
        do i = 1, imax
           widthy = zdy(j)
           if(ngeomy > 2) widthy = widthy*zxc(i)
           width  = min(zdx(i),widthy)
           svel = sqrt(gam*zpr(i,j,1)/zro(i,j,1))/width
           xvel = abs(zux(i,j,1)) / zdx(i)
           yvel = abs(zuy(i,j,1)) / widthy
           ridt = max(xvel,yvel,svel,ridt)
        enddo
     enddo
  else if(ndim == 3) then
     do k = 1, kmax
        do j = 1, jmax
           do i = 1, imax
              widthy = zdy(j)
              widthz = zdz(k)
              if(ngeomy > 2) widthy = widthy*zxc(i)
              if(ngeomz > 2) widthz = widthz*zxc(i)
              if(ngeomz == 5) widthz = widthz*sin(zyc(j))
              width  = min(zdx(i),widthy,widthz)
              svel = sqrt(gam*zpr(i,j,k)/zro(i,j,k))/width
              xvel = abs(zux(i,j,k)) / zdx(i)
              yvel = abs(zuy(i,j,k)) / widthy
              zvel = abs(zuz(i,j,k)) / widthz
              ridt = max(xvel,yvel,zvel,svel,ridt)
           enddo
        enddo
     enddo
  endif

  dt = courant / ridt

  return
end subroutine init

! ----------------------------------------------------------------------------

subroutine grid(nzones, xmin, xmax, xa, xc, dx)
  !
  ! Create grid to cover physical size from xmin to xmax
  ! number of physical grid zones is nzones
  !
  ! xa(1) is left boundary location - at xmin
  ! xa(nzones+1) is right boundary location - at xmax
  ! --------------------------------------------------------------------------
  implicit none

  ! LOCALS
  integer :: nzones, n
  real, dimension(nzones) :: xa, dx, xc
  real :: dxfac, xmin, xmax

  ! --------------------------------------------------------------------------

  dxfac = (xmax - xmin) / float(nzones)
  do n = 1, nzones
     xa(n) = xmin + (n-1) * dxfac
     dx(n) = dxfac
     xc(n) = xa(n) + 0.5 * dx(n)
  enddo

  return
end subroutine grid
