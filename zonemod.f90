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

  real, parameter :: box_xmin = -1.5e14
  real, parameter :: box_xmax = 1.5e14
  real, parameter :: box_ymin = -1.5e14
  real, parameter :: box_ymax = 1.5e14

  real, parameter :: ambient_density = 1.67e-24
  real, parameter :: ambient_pressure = 2.2e-12

  real, parameter :: injection_radius = 7.48e9 ! .5 AU

  real, parameter :: wind_density = 28e-24 ! ?
  real, parameter :: wind_pressure = 8e-8
!  real, parameter :: wind_speed = 1000.

  real, parameter :: shock_density = 1.67e-14 ! ?
  real, parameter :: shock_energy = 0
  real            :: shock_pressure
  real, parameter :: shock_start_time = 0.0015

  real, parameter :: sun_origin_x = 0
  real, parameter :: sun_origin_y = 0
  real, parameter :: sun_origin_z = 0

end module zone
