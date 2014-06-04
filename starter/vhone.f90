program VHone
  use global
  use sweeps
  implicit none
  character(len=8)  :: todayis
  character(len=2)  :: rstrt
  character(len=4)  :: tmp1
  character(len=50) :: hstfile
  character(len=50) :: filename
  integer :: ncycend, nprin, nmovie, ndump, imax, n
  real :: endtime, tprin, tmovie, dxfac, xmin, xmax, ridt, xvel, dt3, dtx
  real :: box_size        = 4.62e19
  real :: ambient_density = 1e-24
  real :: blast_pressure
  real :: blast_energy    = 1e51
  namelist /hinput/ rstrt, prefix, ncycend, ndump, nprin, nmovie, endtime, tprin, tmovie

  ! begin by reading input file (namelist) for job control values
  open (unit=15, file='indat', status='old', form='formatted')
  read (15, nml=hinput)
  close (15)

  ! name and open a file for recording the simulation metadata
  hstfile = 'output/' // trim(prefix) // '.hst'
  open (unit=8, file=hstfile, form='formatted')
  call date_and_time(todayis)
  write (8,*) 'History File for VH-1 simulation run on ', todayis(5:6), ' / ', &
              todayis(7:8), ' / ', todayis(1:4)
  write (8,*)

  ! set some flags for geometry and boundary conditions
  sweep  = 'x'   ! default is x-sweep (only option for 1D)
  ngeom  = 2     ! 2 = spherical geometry
  nleft  = 0     ! reflecting at xmin
  nright = 0     ! reflecting at xmax

  ! create a grid of imax zones, making room for 6 'ghost zones' on each end
  imax = 100           ! total number of real zones on grid
  nmin = 7             ! first real zone
  nmax = imax + 6      ! last real zone
  xmin = 0.            ! x value at left edge of grid
  xmax = box_size      ! x value at right edge of grid

  ! initialize grid coordinates: xa0(n) is coordinate of left edge of zone n
  dxfac = (xmax - xmin) / float(imax)  ! width of each zone
  do n = nmin, nmax
     xa0(n) = xmin + (n - nmin) * dxfac
     dx0(n) = dxfac
  enddo

  gam  = 1.4   ! We always need a ratio of specific heats, gamma
  gamm = gam - 1.0

  ! initial conditions for Sedov blast (supernova)
  blast_pressure = blast_energy / (4. * pi * xa0(nmin + 1) ** 3 / 3.) * gamm
  do n = nmin, nmax
     r(n) = ambient_density
     p(n) = blast_pressure / 1e10
     u(n) = 0.0
     v(n) = 0.0            ! note that we have to carry around the transverse
     w(n) = 0.0            ! velocities even though this is a 1D code
     f(n) = 0.0            ! set initial flattening to zero
  enddo
  p(nmin) = blast_pressure

  ! Write out initial data to a file
  nfile = 0
  write (tmp1,"(i4)") nfile + 1000 ; nfile = nfile + 1
  filename = 'output/' // trim(prefix) // tmp1 // '.dat'
  open (unit=3,file=filename,form='formatted')
  do n = nmin, nmax
     write (3, 1003) xa0(n), r(n), p(n), u(n)
  enddo
  close (3)

  ! log parameters of problem in history file
  write (8, "('Grid dimensions: ',i4)") imax
  write (8, *)
  write (8, "('Ratio of specific heats = ',f5.3)") gam
  write (8, *)

  ! compute initial time step based on Courant condition to ensure stability
  ridt = 0.
  do n = nmin, nmax
     svel = sqrt(gam*p(n)/r(n))/dx0(n)
     xvel = abs(u(n)) / dx0(n)
     ridt = max(xvel,svel,ridt)
  enddo
  dt = courant / ridt

  ! set time and cycle counters
  time   = 0.0
  timep  = 0.
  ncycle = 0
  ncycp  = 0

  ! main computational loop
  do while (ncycle < ncycend)
     if (time + dt > endtime) then
        dt = endtime - time       ! set dt so time lands on endtime
        ncycend = ncycle - 1      ! set ncycend so simulation stops after this step
     endif

     ! construct total energy; reset Lagrangian coordinates to Eulerian grid
     !    Eulerian coordinates, xa0(n), will stay fixed;
     !    Lagrangian coordinates, xa(n), will move with the flow each time step
     do n = nmin, nmax
        e (n) = p(n) / (r(n) * gamm) + .5 * (u(n) ** 2 + v(n) ** 2 + w(n) ** 2)
        xa(n) = xa0(n)
        dx(n) = dx0(n)
     enddo

     ! perform 1D hydro update with PPMLR algorithm
     call ppmlr

     ! update time and cycle counters
     time   = time  + dt
     timep  = timep + dt
     ncycle = ncycle + 1
     ncycp  = ncycp  + 1

     ! time step control
     ridt = 0.         ! searching for maximum, so start out with zero
     do n = nmin, nmax
        svel = sqrt(gam * p(n) / r(n)) / dx0(n)
        xvel = abs(u(n)) / dx0(n)
        ridt = max(xvel, svel, ridt)
     enddo
     dtx = courant / ridt     ! global time constraint for given courant parameter
     dt3 = 1.1 * dt           ! limiting constraint on rate of increase of dt
     dt  = min(dt3, dtx)      ! use smallest required timestep

     ! data output
     if (ncycp >= nprin .or. timep >= tprin) then    ! print out a data file
        write (tmp1,"(i4)") nfile + 1000; nfile = nfile + 1
        filename = 'output/' // trim(prefix) // tmp1 // '.dat'
        write (8,6001) trim(prefix) // tmp1, time, ncycle
        open (unit=3, file=filename, form='formatted')
        do n = nmin, nmax
           write (3, 1003) xa0(n), r(n), p(n), u(n)
        enddo
        close (3)
        timep = 0.
        ncycp = 0
     endif
  enddo
  close (8)

1003 format(1pe13.5,4e13.5)
6001 format('Wrote ',a18,' to disk at time =',1pe12.5,' (ncycle =', i7,')')

end program VHone
