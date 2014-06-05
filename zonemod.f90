module zone

  integer :: ngeomx,  ngeomy,  ngeomz       ! xyz geometry flag
  integer :: nleftx,  nlefty,  nleftz       ! xyz lower boundary condition
  integer :: nrightx, nrighty, nrightz      ! xyz upper boundary condition

  integer, parameter :: imax = 100
  integer, parameter :: jmax = 200
  integer, parameter :: kmax = 1

  real, dimension(imax, jmax, kmax) :: zro, zpr, zux, zuy, zuz, zfl

  real, dimension(imax) :: zxa, zdx, zxc
  real, dimension(jmax) :: zya, zdy, zyc
  real, dimension(kmax) :: zza, zdz, zzc

  ! Distance: AU
  ! Mass: 10^12 kg
  ! Time: year

  real, parameter :: box_xmin = -10
  real, parameter :: box_xmax = 10

  real, parameter :: box_ymin = -20
  real, parameter :: box_ymax = 20

  real, parameter :: ambient_density = 6
  real, parameter :: ambient_pressure = 32.8

  real, parameter :: injection_radius = .5

  real, parameter :: wind_density = 94
  real, parameter :: wind_pressure = 3e5
!  real, parameter :: wind_speed = 1000.

  real, parameter :: shock_density = 20
  real, parameter :: shock_pressure = 1e5
  real, parameter :: shock_velocity = 100
  real, parameter :: shock_start_time = 2

  real, parameter :: sun_origin_x = 0
  real, parameter :: sun_origin_y = 0
  real, parameter :: sun_origin_z = 0

  integer, parameter :: shock_ny = 10

end module zone
