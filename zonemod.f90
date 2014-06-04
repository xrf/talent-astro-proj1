module zone

  integer :: ngeomx,  ngeomy,  ngeomz       ! xyz geometry flag
  integer :: nleftx,  nlefty,  nleftz       ! xyz lower boundary condition
  integer :: nrightx, nrighty, nrightz      ! xyz upper boundary condition

  integer, parameter :: imax = 100
  integer, parameter :: jmax = 100
  integer, parameter :: kmax = 1

  real, dimension(imax, jmax, kmax) :: zro, zpr, zux, zuy, zuz, zfl

  real, dimension(imax) :: zxa, zdx, zxc
  real, dimension(jmax) :: zya, zdy, zyc
  real, dimension(kmax) :: zza, zdz, zzc

  ! Distance: AU
  ! Mass: 10^12 kg
  ! Time: year

  real, parameter :: box_xmin = -4
  real, parameter :: box_xmax = 4
  real, parameter :: box_ymin = -4
  real, parameter :: box_ymax = 4

  real, parameter :: ambient_density = 6
  real, parameter :: ambient_pressure = 32.8

  real, parameter :: injection_radius = .5

  real, parameter :: wind_density = 94
  real, parameter :: wind_pressure = 3e5
!  real, parameter :: wind_speed = 1000.

  real, parameter :: shock_density = 6
  real, parameter :: shock_velocity = 250
  real, parameter :: shock_start_time = .3

  real, parameter :: sun_origin_x = 0
  real, parameter :: sun_origin_y = 0
  real, parameter :: sun_origin_z = 0

end module zone
