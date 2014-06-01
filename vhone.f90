program VHone
  use global
  use zone
  use sweepsize
  implicit none

  character(len=8)  :: todayis
  character(len=2)  :: rstrt
  character(len=3)  :: dumpfile
  character(len=50) :: hstfile

  integer :: ncycend, nprin, ndump, nmovie, start_cycl
  real    :: endtime, tprin, tmovie, olddt
  real    :: start_time, end_time, run_time, zones

  namelist /hinput/ rstrt, prefix, ncycend, ndump, nprin, nmovie, endtime, &
                    tprin, tmovie

  ! --------------------------------------------------------------------------
  ! read input file

  open (unit=15, file='indat', status='old', form='formatted')
  read (15, nml=hinput)
  close (15)

  ! create some file names for metadata (hstfile) and restart file (dumpfile)
  hstfile  = 'output/' // trim(prefix) // '.hst'
  dumpfile = 'daa'

  ! check that arrays are large enough for desired number of physical zones
  if (max(imax, jmax, kmax) + 12 > maxsweep) then
     write (*, *) 'maxsweep too small'
     stop
  endif

  ! set the number of dimensions based on array sizes
  if (jmax * kmax == 1) then
     ndim = 1
  else if (kmax == 1) then
     ndim = 2
  else
     ndim = 3
  endif

  if ((rstrt == 'NO') .or. (rstrt == 'no')) then
     ! initialize variables for new problem
     open (unit=8, file=hstfile, form='formatted')
     call date_and_time(todayis)
     write (8,*) 'History File for VH-1 simulation run on ', &
                 todayis(5:6), ' / ', todayis(7:8), ' / ', todayis(1:4)
     write (8,*)
     call init
     call prin
  else
     ! restart from old dumpfile
     dumpfile = 'd' // rstrt
     open (unit=8, file=hstfile, form='formatted', position='append')
     call undump(dumpfile)
     if(timep >= tprin)  timep = tprin  - 2. * dt
     ! adjust so that dt is not 0 on 1st step
     if(timem >= tmovie) timem = tmovie - 2. * dt
  endif

  call cpu_time(start_time)
  start_cycl = ncycle

  ! --------------------------------------------------------------------------
  ! main loop
  do while (ncycle < ncycend)

     !  write(*, *) 'STEP = ', ncycle

     ncycle = ncycle + 2
     ncycp  = ncycp  + 2
     ncycd  = ncycd  + 2
     ncycm  = ncycm  + 2
     olddt  = dt
     svel   = 0.

     if (time + 2. * dt > endtime) then ! set dt to land on endtime
        write (8, *) 'cutting to the end...', ncycle, ncycend
        dt = .5 * (endtime - time)
        ncycend = ncycle - 1
        ncycp   = nprin
        ncycd   = ndump
     else if (timep + 2. * dt > tprin) then ! set dt to land on tprin
        dt = .5 * (tprin - timep)
        ncycp = nprin
     else if (timem + 2. * dt > tmovie) then ! set dt to land on tmovie
        dt = .5 * (tmovie - timem)
        ncycm = nmovie
     endif

     ! alternate sweeps to approximate 2nd order operator splitting

     call sweepx
     if(ndim >  1) call sweepy
     if(ndim == 3) call sweepz

     time  = time  + dt
     timep = timep + dt
     timem = timem + dt

     if(ndim == 3) call sweepz
     if(ndim >  1) call sweepy
     call sweepx

     time  = time  + dt
     timep = timep + dt
     timem = timem + dt
     dt    = olddt

     ! check constraints on the timestep
     call dtcon

     ! output movie images/datasets/dumpfiles based on simulation time or
     ! cycle number

     if (ncycm >= nmovie .or. timem >= tmovie) then
        call images
        timem = 0.
        ncycm = 0
     endif

     if (ncycp >= nprin .or. timep >= tprin) then
        call prin
        timep = 0.
        ncycp = 0
     endif

     if (ncycd >= ndump) then
        ncycd = 0
        call dump(dumpfile)
     endif

  enddo
  ! --------------------------------------------------------------------------

  call cpu_time(end_time)
  run_time = end_time - start_time
  zones    = (imax * jmax * kmax) * (ncycle - start_cycl)
  write (8, *) 'successful stop at cycle number', ncycle
  if (run_time > 4000.) then
     write (8, *) 'elapsed time = ', run_time / 3600., ' hours'
  else
     write (8, *) 'elapsed time = ', run_time / 60., ' minutes'
  endif
  write (8, "('speed = ',f5.1,' kz/s')") 1.0e-3 * zones / run_time

  close (8)

end program VHone
