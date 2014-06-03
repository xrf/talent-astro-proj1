module zone

  integer, parameter :: imax=100, jmax=100, kmax=1   ! memory dimensions

  real, dimension(imax, jmax, kmax) :: zro, zpr, zux, zuy, zuz, zfl

  real, dimension(imax) :: zxa, zdx, zxc
  real, dimension(jmax) :: zya, zdy, zyc
  real, dimension(kmax) :: zza, zdz, zzc

  integer :: ngeomx,  ngeomy,  ngeomz       ! xyz geometry flag
  integer :: nleftx,  nlefty,  nleftz       ! xyz lower boundary condition
  integer :: nrightx, nrighty, nrightz      ! xyz upper boundary condition

  real :: ambient_density = 1.
  real :: ambient_pressure = 1.


  real :: injection_radius = .05
  real :: solar_density = 3.
  real :: solar_pressure = 8e6
  real :: solar_speed = 1000.

  real :: shock_density = 1
  real :: shock_pressure = 5e7

  real :: sun_origin_x = .5
  real :: sun_origin_y = .5
  real :: sun_origin_z = 0

  real :: shock_start_time = 0.0015

end module zone
