!**************************************************************************
!
! Program SLFWD_LTT
! -----------------
!
! This is a GPS slant delay forward model developed at FMI, in co-operation
! with ECMWF, in 2006. This slant delay forward model is based on ray
! tracing for determination of the GPS signal path through the NWP model
! grid.
!
! authors Reima Eresmaa and Sean Healy
! email   Reima.Eresmaa@fmi.fi
! date    24.11.2006
!
! OpenIFS implementation: Lauri Tuppi
! Implementation of dn/dpsi term, geometrical correction and condensate fields: Maksym Vasiuta
! date 16.12.2021
!
! Example run call:
! ./slfwd_v8.exe oifs2016120100+000 p000 $inputdir $outputdir 1 66 t1279
!
!--------------------------------------------------------------------------
!
!**************************************************************************

!--------------------------------------------------------------------------
!
! ------------ MAIN PROGRAM ------------
!
!--------------------------------------------------------------------------

! Variables:
!
! commandline:
! year_it: year of the analysis
! month_it: month of the analysis
! day_it: day of the analysis
! hour_it: hour of the analysis
! fclen: forecast length
! year_vt: year at the slant delay calculation
! month_v: month tat the slant delay calculation
! day_vt: day at the slant delay calculation
! hour_vt: hour at the slant delay calculation
! gribdir: path to the directory where the input grib file is
! resdir: output directory
! gfname: name of the input grib file
! ccyc: 10 character identifier of the current analysis cycle
! fcid: identifier constructed from forecast length
! nx: number of longitudes in grib file
! ny: number of latitudes in grib file
! nlev: number of vertical levels in grib file
! dlon: spacing between longitudes in degrees
! dlat: spacing between latitudes in degrees
! afull: hybrid A coefficient, vertical sigma coordinate
! bfull: hybrid B coefficient, vertical sigma coordinate
! fis_bg: 2D geopotential of the Earth's surface
! lnps_bg: 2D ln(pressure at the Earth's surface)
! t_bg: 3D atmospheric temperature
! q_bg: 3D specific humidity
! dvecrad: latitude increment in radians, constant, used in one of the subroutines
! maxobs: limit for number of rays
! nhor: maximum number of columns used in the 2D horizontal projection
! nobs_tot: number of rays
! stid: vector of 4 letter station identifier
! lat: vector of latitude for one station
! lon: vector of longitude for one station
! hei: vector of height for one station
! aa_geom: vector of azimuth angles
! za_geom: vector of zenith angles
! lood: identifier if the column in question has been used in the projection to the 2D plane (not sure)
! sdelay_nwp: vector of slant delays
! ierr: error code if reading the grib file encounters problems
! dpm: number of days in each month
! dpmleap: number of days in each month in the leap years
! t: time indicator used to time calculation of one slant delay field
! sec1: start time
! sec2: end time
! expid: experiment identifier
! rceq: radius of the Earth
! stationfile: text file containing a list of stations
! ycoordsfile: text file containing latitudes in Gaussian regular grid
! stations: 4 letter identifiers of all stations
! longitudes: longitudes of all stations
! latitudes: latitudes of all stations
! heights: altitudes of all stations
! m_lats: vector of latitudes in Gaussian regular grid
! date: month and day of the analysis
! nobs: number of rays in one sky view
! i,j,k: loop indices
! nm_lats: number of latitudes in Gaussian regular grid
! hour: hour of the analysis
! dvec_lon: distance between two longitudes			!!!
! dvec_lat: distance between two latitudes
! zen: number of zenith angles
! azi: number of azimuth angles

PROGRAM slfwd_ltt_v8
use eccodes
!--------------------------------------------------------------------------
!
! 1. DECLARATION OF VARIABLES FOR THE MAIN PROGRAM
!
!--------------------------------------------------------------------------

  IMPLICIT NONE

! Definitions for the time (to be read from the command line)
  CHARACTER(LEN=90) :: commandline
  INTEGER :: year_it, month_it, day_it, hour_it, fclen
  INTEGER :: year_vt, month_vt, day_vt, hour_vt

! Definitions for the directories and file names
  CHARACTER(LEN=90) :: gribdir, resdir
  CHARACTER(LEN=90) :: gfname
  CHARACTER(LEN=10) :: ccyc
  CHARACTER(LEN=2) :: fcid

! Definitions for the NWP domain and the needed fields
  INTEGER :: nx, ny, nlev
  REAL :: dlon, dlat
  REAL, ALLOCATABLE :: afull(:), bfull(:)
  REAL, ALLOCATABLE :: fis_bg(:,:), lnps_bg(:,:), t_bg(:,:,:), q_bg(:,:,:)
  ! REAL, ALLOCATABLE :: clwc_bg(:,:,:), ciwc_bg(:,:,:), crwc_bg(:,:,:), cswc_bg(:,:,:)

! Horizontal resolution of the 2d projected array
  REAL :: dvecrad

! Definitions for the observational data
  INTEGER, PARAMETER  :: maxobs=37801, nhor=120
  CHARACTER(LEN=4), DIMENSION(maxobs) :: stid
  REAL, DIMENSION(maxobs) :: lat, lon, hei, aa_geom, za_geom

! Definitions for the modelled slant delays
  INTEGER, DIMENSION(maxobs) :: lood
  REAL, DIMENSION(maxobs) :: sdelay_nwp=0.0

! Local work variables
  INTEGER :: ierr
  INTEGER, DIMENSION(12), PARAMETER :: &
       dpm=    (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /), &
       dpmleap=(/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  INTEGER, DIMENSION(8) :: t
  REAL :: sec1, sec2
  CHARACTER(LEN=3) :: expid
  REAL, PARAMETER :: rceq=6378137.0

  character(len=90) :: stationfile, ycoordsfile
  character(len=4),dimension(:),allocatable :: stations
  real,dimension(:),allocatable :: longitudes, latitudes, height, m_lats
  character(len=25) :: outfile1, outfile2
  integer :: date, nobs, i,j, nm_lats
  character(len=3) :: hour
  real :: dvec_lon
  real, allocatable :: dvec_lat(:)
  integer :: zen, azi, start, sstop

  real,dimension(1) :: zen_sdelay_nwp
  real,dimension(:),allocatable :: z_sdelay_nwp

  real,dimension(:),allocatable :: zen_angles, azi_angles
  character(len=4) :: pert
  integer :: nancount=0

  character(len=4) :: start_station, end_station
  character(len=5) :: resolution
  integer :: start_station_i, end_station_i
!--------------------------------------------------------------------------
!
! 2. THE MAIN PROGRAM CODE
!
!--------------------------------------------------------------------------
  print*,'###### Welcome to using slfwd slant delay operator. ######'

! 2.1 Define the time of the observations and forecast field
! The name of the input file is cut into pieces, and the values assigned to the variables.
  CALL get_command_argument(1,commandline)
  read(commandline(5:8),'(i4)') year_it
  read(commandline(9:10),'(i2)') month_it
  read(commandline(11:12),'(i2)') day_it
  read(commandline(13:14),'(i2)') hour_it
  read(commandline(16:18),'(i3)') fclen
!--------------------------------------------------------------------------
!
! 2.2 Read the paths for the input grib file and output directory.
! Then read the list of GNSS stations and the list of latitudes in
! regular Gaussian grid of the input file.
  print*,'Reading the configurations etc ... '
  call get_command_argument(2,pert)
  call get_command_argument(3,gribdir)
  call get_command_argument(4,resdir)
  call get_command_argument(5,start_station)
  call get_command_argument(6,end_station)
  call get_command_argument(7,resolution)
  read(start_station,*) start_station_i
  read(end_station,*) end_station_i
  print*,' ... done'
  ierr=0
  !stationfile='LTT_conf_files/stationCoordinates_extd.txt'
  stationfile='LTT_conf_files/stationCoordinates.txt'
  if (resolution=='t639') then
    ycoordsfile='LTT_conf_files/lats_t639.txt'
  elseif (resolution=='t1279') then
    ycoordsfile='LTT_conf_files/lats_t1279.txt'
  else
    print*,'Damn, that resolution is undefined!'
  endif
  call readStations(stationfile, ierr, stations, longitudes, latitudes, height)
    if (ierr/=0) then
      write (*,'(A,I6)') 'Cannot read the list of stations'
      stop
    endif
  call read_y_coords(ycoordsfile,nm_lats,m_lats)
!--------------------------------------------------------------------------
!
! 2.3 Derive the name for the GRIB-file and retrieve the needed NWP fields
!
  fcid=ch2(fclen)
  CALL cycleid( year_it, month_it, day_it, hour_it, ccyc)
  gfname=commandline
  WRITE (*,'(A)',ADVANCE='NO') &
       'Reading the GRIB file '//TRIM(gribdir)//TRIM(gfname)//' ... '
  CALL getDomain2(gribdir, gfname, nx, ny, nlev, dlon, dlat, ierr)
  ALLOCATE (afull(nlev+1), bfull(nlev+1))
  ALLOCATE (fis_bg(nx,ny), lnps_bg(nx,ny))
  ALLOCATE (t_bg(nx,ny,nlev), q_bg(nx,ny,nlev))
  ! ALLOCATE (clwc_bg(nx,ny,nlev), ciwc_bg(nx,ny,nlev))
  ! ALLOCATE (cswc_bg(nx,ny,nlev), crwc_bg(nx,ny,nlev))
  allocate (dvec_lat(nx-1))
  CALL getFields2(&
       gribdir, gfname, resolution, nx, ny, nlev, &
       afull, bfull, &
       fis_bg, lnps_bg, t_bg, q_bg, ierr)
  IF (ierr/=0) THEN
     WRITE (*,'(A)') 'Cannot read GRIB-file: ' // gfname
     STOP
  END IF
  WRITE (*,'(A)') ' ... done'

!--------------------------------------------------------------------------
!
! 2.4 Calculate the valid time of the forecast, derive the name for the
!     observation file and read the observations
!
  year_vt  = year_it
  month_vt = month_it
  day_vt   = day_it
  hour_vt  = hour_it+fclen
  DO WHILE (hour_vt>=24)
     hour_vt=hour_vt-24
     day_vt=day_vt+1
  END DO
  IF ( (MOD(year_it,4)==0 .AND. MOD(year_it,100)/=0) &
       .OR. MOD(year_it,400)==0) THEN
     IF (day_vt>dpm(month_it)) THEN
        day_vt=day_vt-dpmleap(month_vt)
        month_vt=month_vt+1
     END IF
  ELSE
     IF (day_vt>dpm(month_it)) THEN
        day_vt=day_vt-dpm(month_vt)
        month_vt=month_vt+1
     END IF
  END IF

!--------------------------------------------------------------------------
!
! 2.5 Set the horizontal resolution of the projected grid.
! Set the numbers of zenith and azimuth angles.
! Begin loop over all GNSS stations and initialize the slant delay sky view.
!

  dvec_lon = rceq * degtor(dlon)
  do j=1,nm_lats-1
    dvec_lat(j) = rceq * degtor(m_lats(j)-m_lats(j+1))
  enddo
  dvecrad = degtor(dlat) ! this should be replaced sometimes

  zen=85 !105!75!89
  azi=360
  allocate(zen_angles(zen))
  allocate(azi_angles(azi))
  allocate(z_sdelay_nwp(azi))
  nobs=zen*azi+1
  do i=start_station_i,end_station_i !1,size(stations)
    print*,'####### Initializing sky view for station ',i,stations(i),' #######'
    call createSlantGPSfield(stations(i), longitudes(i), latitudes(i), height(i), nobs, stid, lat, lon, hei, &
        zen, azi, aa_geom, za_geom, zen_angles, azi_angles)
    print*,'  ... Done!'

!--------------------------------------------------------------------------
!
! 2.6 Call the actual forward modelling subroutine.
! Also, take time how long it takes to calculate the sky view.
!
    print*,'####### Calculating actual slant delays for station ',stations(i),' #######'

    CALL DATE_AND_TIME(values=t)
    sec1 = t(5)*3600.0+t(6)*60.0+t(7)*1.0+t(8)*0.001

    CALL slantGPS_ltt( &
         1, 1, stid(1), lat(1), lon(1), hei(1), &
         aa_geom(1), za_geom(1), &
         nx, ny, nlev, &
         dlon, &
         afull, bfull, &
         lnps_bg, fis_bg, t_bg, q_bg, &
        !  clwc_bg, ciwc_bg, cswc_bg, crwc_bg, &
         nhor, dvecrad, lood, zen_sdelay_nwp, nm_lats, dvec_lon, dvec_lat, m_lats )

    !print*, zen_sdelay_nwp
    print*,'Zenith delay = ', zen_sdelay_nwp
    sdelay_nwp(1) = zen_sdelay_nwp(1)

    do j=1,zen
      start=azi*(j-1)+2
      sstop=start+azi-1
!     do j=1,3060
!      start=(j-1)*1000+1
!      sstop=j*1000
      !print*,start,sstop,size(sdelay_nwp(start:sstop))
      z_sdelay_nwp = sdelay_nwp(start:sstop)

      CALL slantGPS_ltt( &
           azi, azi, stid(start:sstop), lat(start:sstop), lon(start:sstop), hei(start:sstop), &
           aa_geom(start:sstop), za_geom(start:sstop), &
           nx, ny, nlev, &
           dlon, &
           afull, bfull, &
           lnps_bg, fis_bg, t_bg, q_bg, &
          !  clwc_bg, ciwc_bg, cswc_bg, crwc_bg, &
           nhor, dvecrad, lood, z_sdelay_nwp, nm_lats, dvec_lon, dvec_lat, m_lats )

      sdelay_nwp(start:sstop) = z_sdelay_nwp
    enddo

!#########################################################################################

    CALL DATE_AND_TIME(values=t)
    sec2 = t(5)*3600.0+t(6)*60.0+t(7)*1.0+t(8)*0.001
    WRITE (*,*) '  ... Done! Time:', sec2-sec1
!--------------------------------------------------------------------------
!
! 2.7 Write the modelled slant delays in an ascii file.
! Data for all GNSS stations are written in the same file.
!
    read(commandline(5:12), '(i8)') date
    write (hour,'(i3.3)') fclen
    read(commandline(5:14), '(a10)') outfile1
    outfile2=trim(outfile1) // '+' // trim(hour) // '_' // pert // '.dat'

    print*,'####### Writing the sky view for station ',stations(i),' to #######'
    print*,outfile2
    !print*,sdelay_nwp(1:360)
    call writeDelays2(resdir, stations(i), lat, lon, hei, sdelay_nwp, outfile2, date, hour_vt, zen, azi, zen_angles, azi_angles)
    do j=1,nobs
      if (sdelay_nwp(j) /= sdelay_nwp(j)) then
        nancount=nancount+1
        !print*,'!!!!!!! NaNs detected !!!!!!!',sdelay_nwp(i)
      endif
    enddo
    print*,'number of NaNs: ',nancount
    nancount=0
    print*,'  ... Done!'
  enddo

CONTAINS


!**************************************************************************
!
! ------------ SUBROUTINES ------------
!
!**************************************************************************

!#########################################################################################
!################# Lauri's subroutines ###################################################
!#########################################################################################

!==============================================================================
! Subroutine GETDOMAIN retrieves the model domain parameters from the given 
! OpenIFS GRIB-file 'gribdir+gfname'.
!
! IN: gribdir, gfname
! OUT: nx, ny, nlev, dlon, dlat, ierr
!
! Variables:
!
! gribdir: path to the directory where the input grib file is
! gfname: name of the input grib file
! filename: locally used name for the grib file, gribdir+gfname
! nx: number of longitudes
! ny: number of latitudes
! nlev: number of vertical levels
! west: longitude of the westernmost point
! south: latitude of the southernmost point
! east: longitude of the easternmost point
! north: latitude of the northernmost point
! dlon: spacing between longitudes in degrees
! dlat: spacing between latitudes in degrees
! ierr: variable for possible error code
! rfile: memory unit where the grib file is opened in
! iret: variable for possible error code
! n: iteration counter
! messages: number of grib messages in a file
! level: model level number read from grib message
! igrib: contains identifiers for all grib messages
! levels: contains all model level information from all messages
!
! Logic:
!
! 1. Open the OpenIFS grib file, and count the number of grib messages in the
! file. Then allocate storage vectors for message identifiers and model level
! information. Set levels to be zero in order to clean possible random data.
!
! 2. Retrieve identifiers for grib messages.
!
! 3. Get horizontal field dimensions from some grib message. The message number
! is arbitrary.
!
! 4. This block retrieves level information from all grib messages, and takes
! the largest value from the model levels. The problem here is how to
! distinguish model levels from pressure levels. Shortcut solution is to
! exclude all large numbers like 1000 that means 1000hPa level. The limit of
! 150 is arbitrary. This should be done in a better way.
!
! 5. Latitude and longitude edges of the model domain are retrieved.
! The message number is again arbitrary.
!
! 6. Longitude and latitude increments are calculated.
!
! 7. Memory is freed, and the file closed.

  subroutine getDomain2(gribdir,gfname,nx,ny,nlev,dlon,dlat,ierr)
    CHARACTER(LEN=90) :: gribdir, gfname
    character(len=110) :: filename
    INTEGER :: nx, ny, nlev
    REAL :: west, south, east, north, dlon, dlat
    INTEGER :: ierr

    integer :: rfile,iret,n,messages,level
    integer,dimension(:),allocatable :: igrib,levels
! 1.
    filename = trim(gribdir)//trim(gfname)
    call codes_open_file(rfile,filename,'r')
    call codes_count_in_file(rfile,messages,ierr)
    allocate(igrib(messages+1))
    allocate(levels(messages+1))
    levels=0
! 2.
    n=1
    LOOP1: do while (iret/= CODES_END_OF_FILE)
      call codes_grib_new_from_file(rfile,igrib(n),iret)
      n=n+1
    end do LOOP1
! 3.
    call codes_get(igrib(1),'Ni',nx)
    call codes_get(igrib(1),'Nj',ny)
! 4.
    n=1
    LOOP2: do while (n<=messages)
      call codes_get(igrib(n),'level',level)
      if (level<150) then
        levels(n)=level
      else
        levels(n)=1
      end if
      n=n+1
    end do LOOP2
    nlev=maxval(levels)
! 5.
    call codes_get(igrib(1),'latitudeOfFirstGridPointInDegrees',north)
    call codes_get(igrib(1),'latitudeOfLastGridPointInDegrees',south)
    call codes_get(igrib(1),'longitudeOfFirstGridPointInDegrees',west)
    call codes_get(igrib(1),'longitudeOfLastGridPointInDegrees',east)
! 6.
    dlon  = (east - west)/real(nx)!-1)
    dlat  = (north-south)/real(ny)!-1)
! 7.
    do n=1,messages
      call codes_release(igrib(n))
    end do
    call codes_close_file(rfile)

  end subroutine getDomain2

!==========================================================================
! Subroutine GETFIELDS retrieves the needed fields of orography, logarithm
! of surface pressure, temperature and specific humidity from the given
! OpenIFS GRIB-file 'gfname'. NOTE: the file must be post-processed into
! regular gaussian grid.
!
! IN: gribdir, gfname, nx, ny, nlev
! OUT: afull, bfull, fis_bg, lnps_bg, t_bg, q_bg, ierr
!
! Variables:
!
! gribdir: path to the directory where the input grib file is
! gfname: name of the input grib file
! filename: locally used name for the grib file, gribdir+gfname
! nx: number of longitudes
! ny: number of latitudes
! nlev: number of vertical levels
! paramId: grib parameter of OpenIFS output field
! afull: hybrid A coefficient, vertical sigma coordinate
! bfull: hybrid B coefficient, vertical sigma coordinate
! fis_bg: 2D geopotential of the Earth's surface
! lnps_bg: 2D ln(pressure at the Earth's surface)
! t_bg: 3D atmospheric temperature
! q_bg: 3D specific humidity
! ierr: variable for possible error code
! field: temporary storage for field data
! jx: longitude counter
! jy: latitude counter
! jz: level counter
! rfile: memory unit where the grib file is opened in
! iret: variable for possible error code
! n: iteration counter
! i: iteration counter
! nb_values: number of data elements in one grib message
! messages: number of grib messages in a file
! coordsfile: contains vertical coordinates of OpenIFS
! igrib: contains identifiers for all grib messages
! avgz: average geopotential
! switch: used for resetting level counter jz
!
! Logic:
!
! 1. Initially set everything to zero in order to clean possible random data.
!
! 2. Open the OpenIFS grib file, and count the number of grib messages in the
! file. Allocate vector for message identifiers.
!
! 3. Retrieve identifiers for grib messages.
!
! 4. Retrieve number of elements in the data in the grib message, and allocate
! temporary storage vector accordingly. The message number is arbitrary.
!
! 5. Read vertical coordinates of OpenIFS from a separate file, and put them
! into afull and bfull. Initial idea was to read also these directly from
! the grib file but it did not start working.
!
! 6. This is the main do loop that reads the data in the grib messages, and
! sorts the data into apprpriate variables. First retrieve data and grib code
! from the grib messages.
! 6.1. Switch system, which controls the vertical level counter, is turned on
! when the first grib message of temperature (grib code 130) or specific
! humidity (grib code 133) is encountered.
! 6.2. Grib codes are used to sort the data into correct matrices. The data,
! which initially is in vector form is also reshaped into 2D or 3D matrices.
! There is one problem in separating surface geopotential from geopotential
! of the first model level; they have identical metadata. Now there is a
! shortcut solution so that average value of geopotential is used to separate
! these two. This value might be different with different grid resolutions.
! 6.3. This is the latter part of the switch system. It will be turned off
! when the level counter hits the max number of levels.
!
! 7. Memory is freed, and the file closed.

  subroutine getFields2(gribdir,gfname,resolution,nx,ny,nlev,afull,bfull,fis_bg,lnps_bg,t_bg,q_bg, ierr)
    CHARACTER(LEN=90) :: gribdir, gfname
    character(len=5) :: resolution
    character(len=110) filename
    INTEGER :: nx, ny, nlev, paramId
    REAL :: afull(nlev+1), bfull(nlev+1)
    REAL :: fis_bg(nx,ny), lnps_bg(nx,ny)
    REAL :: t_bg(nx,ny,nlev), q_bg(nx,ny,nlev)
    ! REAL :: clwc_bg(nx,ny,nlev), ciwc_bg(nx,ny,nlev), cswc_bg(nx,ny,nlev), crwc_bg(nx,ny,nlev)
    INTEGER :: ierr
    REAL, ALLOCATABLE :: field(:)
    INTEGER :: jx, jy, jz

    integer :: rfile,iret,n,i,nb_values,messages,level
    character(len=90) :: coordsfile
    integer,dimension(:),allocatable :: igrib
    real :: avgz
    logical :: switch=.false.
    real :: coords(3,nlev+1)

    if (resolution=='t639') then
      coordsfile='LTT_conf_files/vert_coord_l91.dat'
    elseif (resolution=='t1279') then
      coordsfile='LTT_conf_files/vert_coord_l137.dat'
    else
      print*,'Damn, define the resolution correctly!!!'
    endif
! 1.
    afull=0.0
    bfull=0.0
    fis_bg=0.0
    lnps_bg=0.0
    t_bg=0.0
    q_bg=0.0
    ! clwc_bg=0.0
    ! ciwc_bg=0.0
    ! cswc_bg=0.0
    ! crwc_bg=0.0
! 2.
    filename = trim(gribdir) // trim(gfname)
    call codes_open_file(rfile,filename,'r')
    call codes_count_in_file(rfile,messages,ierr)   
    allocate(igrib(messages+1))
! 3.
    n=1
    LOOP1: do while (iret/= CODES_END_OF_FILE)
      call codes_grib_new_from_file(rfile,igrib(n),iret)
      n=n+1
    end do LOOP1
! 4.    
    call codes_get(igrib(10),"numberOfPoints",nb_values);allocate(field(nb_values))
! 5.
    call read_vert_coords(coordsfile,nlev,coords)
    afull=coords(2,:)
    bfull=coords(3,:)
! 6.
    do i=1,messages
      call codes_get(igrib(i),"values",field)
      call codes_get(igrib(i),"paramId",paramId)
! 6.1.
      if (.not. switch .and. paramId==130) then ! 3-d temp
        switch=.true.
        jz=1
      elseif (.not. switch .and. paramId==133) then ! 3-d spec. humid
        switch=.true.
        jz=1
      elseif (.not. switch .and. paramId==246) then ! 3-d clwc
        switch=.true.
        jz=1
      elseif (.not. switch .and. paramId==247) then ! 3-d ciwc
        switch=.true.
        jz=1
      elseif (.not. switch .and. paramId==76) then ! 3-d cswc
        switch=.true.
        jz=1
      elseif (.not. switch .and. paramId==75) then ! 3-d crwc
        switch=.true.
        jz=1
      end if
! 6.2.
      if (paramId==152) then ! logarithm of surface pressure
        DO jx=1,nx
          DO jy=1, ny
            lnps_bg(jx,jy)=field((jy-1)*nx+jx)
          END DO
        END DO
      elseif (paramId==130) then ! temperature
        DO jx=1,nx
          DO jy=1, ny
            t_bg(jx,jy,jz)=field((jy-1)*nx+jx)
          END DO
       END DO
      elseif (paramId==133) then ! specific humidity
        DO jx=1,nx
          DO jy=1, ny
            q_bg(jx,jy,jz)=field((jy-1)*nx+jx)
          END DO
        END DO
      ! elseif (paramId==246) then ! clwc
      !   DO jx=1,nx
      !     DO jy=1, ny
      !       clwc_bg(jx,jy,jz)=field((jy-1)*nx+jx)
      !     END DO
      !   END DO
      ! elseif (paramId==247) then ! ciwc
      !   DO jx=1,nx
      !     DO jy=1, ny
      !       ciwc_bg(jx,jy,jz)=field((jy-1)*nx+jx)
      !     END DO
      !   END DO
      ! elseif (paramId==76) then ! cswc
      !   DO jx=1,nx
      !     DO jy=1, ny
      !       cswc_bg(jx,jy,jz)=field((jy-1)*nx+jx)
      !     END DO
      !   END DO
      ! elseif (paramId==75) then ! crwc
      !   DO jx=1,nx
      !     DO jy=1, ny
      !       crwc_bg(jx,jy,jz)=field((jy-1)*nx+jx)
      !     END DO
      !   END DO
      elseif (paramId==129) then ! 2D surface geopotential
        !call codes_get(igrib(i),"average",avgz) ! use this trick if the file contains geopotential at surface AND model levels
        !call codes_get(igrib(i),"level",level) !
        !if (level==1 .and. avgz<3720.0) then !
        DO jx=1,nx
          DO jy=1, ny
            fis_bg(jx,jy)=field((jy-1)*nx+jx)
          END DO
        END DO
        !end if !
      end if
! 6.3.
      if (switch) then
        jz=jz+1
      end if
      if (jz>nlev) then
        switch=.false.
      end if
    end do
! 7.
    do n=1,messages
      call codes_release(igrib(n))
    end do
    call codes_close_file(rfile)

  end subroutine getFields2

!==========================================================================
! Subroutine read_vert_coords reads the vertical coordinates of OpenIFS
! from an ASCII file
!
! Variables:
!
! coordsfile: contains vertical coordinates of OpenIFS
! nlev: number of vertical levels
! coords: final coordinates
! coords2: helper variable used to read the coordinate file
!
! Logic:
!
! Read the hybrid sigma vertical coordinates from a separate file. The file
! contains level numbering, and nlev+1 A and B coordinates. Level 0 is the
! space so it is omitted.

  subroutine read_vert_coords(coordsfile,nlev,coords)
    character(len=90) :: coordsfile
    integer :: nlev
    real,dimension(3,nlev+1) :: coords
    ! real,dimension(3,nlev+1) :: coords2
    open(unit=10, file=coordsfile)
    read(10,*) coords
    ! coords=coords2(:,2:nlev+1)
  end subroutine read_vert_coords

!==========================================================================
! Subroutine read_y_coords reads the latitudes of OpenIFS
! from an ASCII file
!
! Variables:
!
! ycoordsfile: file containing the list of latitudes
! nlats: number of latitudes
! lats: latitudes in degrees
!
! Logic:
! Count the number of latitudes and
! then read the latitudes line by line.

  subroutine read_y_coords(ycoordsfile,nlats,lats)
    character(len=90) :: ycoordsfile
    integer :: nlines=0, ierr, io, nlats
    real,dimension(:),allocatable :: lats
    open (1, file = ycoordsfile)
    do 
      read(1,*,iostat=io)
      if (io/=0) exit
      nlines = nlines + 1
    enddo
    close (1)
    nlats=nlines
    allocate(lats(nlines))
    open (2, file = ycoordsfile)
    do i=1,nlines
      read(2,*) lats(i)
    enddo
    close (2)
  end subroutine read_y_coords

!==========================================================================
! Subroutine createSlantGPSfield initializes the slant delay sky view
!
! Variables:
!
! station: four-letter GNSS station identifier
! longitude: longitude of the station
! latitude: latitude of the station
! height: height above mean sea level of the station
! nrays: number of rays in one sky view
! stid: four-letter GNSS identifier assigned to each ray
! lats: GNSS station latitude assigned to each ray
! lons: GNSS station longitude assigned to each ray
! heis: GNSS station height assigned to each ray
! izen: number of zenith angles
! iazi: number of azimuth angles
! aa_geom: azimuth angle of ray
! za_geom: zenith angle of ray
!
! d_az: azimuth angle increment
! d_zen: zenith angle increment
! az: azimuth angle
! zen: zenith angle
! i: loop index
! j: loop index
! counter: helper variable to initialize the angles corrctly
!
! Logic:
! Prepare vectors of station identifier, station latitude, station longitude,
! station height, azimuth angles and zenith angles.
! These vectors are later used element by element in the slant delay calculations.

  subroutine createSlantGPSfield(station, longitude, latitude, height, nrays, stid, lats, lons, heis, izen, iazi, &
             aa_geom, za_geom, zen_angles, azi_angles)
    implicit none
    character(len=4) :: station
    real :: longitude, latitude, height
    integer :: nrays
    CHARACTER(LEN=4), DIMENSION(nrays) :: stid
    REAL, DIMENSION(nrays) :: lats, lons, heis
    integer :: izen, iazi
    REAL, DIMENSION(nrays) :: aa_geom, za_geom
    real :: d_az, d_zen1, d_zen2, d_zen3, az, zen
    integer :: i, j, counter
    real, dimension(izen) :: zen_angles
    real, dimension(iazi) :: azi_angles
    zen=1.0
    d_az=1.0
    d_zen1=1.0!3.0
    d_zen2=1.0
    d_zen3=1.0!0.2
    do i=1,nrays
      stid(i)=station
      lats(i)=latitude
      lons(i)=longitude
      heis(i)=height
    enddo
    aa_geom(1)=1.0
    za_geom(1)=0.0
    counter=2
    do i=1,izen
      az=1.0
      do j=1,iazi
        aa_geom(counter)=az
        za_geom(counter)=zen
        counter=counter+1
        az=az+d_az
      enddo
      !zen=zen+d_zen
      !print*,i,zen
      zen_angles(i)=zen
      if (0==1) then
        if (i<=19) then
          zen=zen+d_zen1
        elseif (i>19 .and. i<=29) then
          zen=zen+d_zen2
        elseif (i>29) then
          zen=zen+d_zen3
        endif
      endif
      if (0==1) then
        zen=zen+0.2
      endif
      if (0==1) then
        zen=zen+0.1
      endif
      if (1==1) then
        zen=zen+1.0
      endif
    enddo
    az=1.0
    if (1==1) then
      do i=1,iazi
        azi_angles(i)=az
        az=az+1.0
      enddo
    endif
    if (0==1) then
      do i=1,iazi
        azi_angles(i)=az
        az=az+0.1
      enddo
    endif
  end subroutine createSlantGPSfield

!==========================================================================
! Subroutine readStations reads the station coordinates from a given text file.
!
! Variables:
!
! stationfile: text file containing the station details
! ierr: 
! stations: a vector containing four letter station identifiers
! longitudes: a vector containing longitudes of the stations
! latitudes: a vector containing latitudes of the stations
! heights:a vector containing altitudes of the stations respect to mean sea level
!
! fname: path and name of the text file of the station coordinates
! nlines: number of lines in the text file
! io: a helper variable to test whether the line exists
! i: loop index
!
! Logic:
! Count the number of lines in the station coordinate file.
! Allocate the vectors for the station details.
! Read the the station coordinate file

  subroutine readStations(stationfile, ierr, stations, longitudes, latitudes, heights)
    implicit none
    character(len=90) :: stationfile
    character(len=120) :: fname
    integer :: nlines, ierr, io, i
    character(len=4),dimension(:),allocatable :: stations
    real,dimension(:),allocatable :: longitudes, latitudes, heights

    fname=trim(stationfile)
    nlines = 0 
    open (1, file = fname)
    do
      read(1,*,iostat=io)
      if (nlines==0) then
        ierr=io
      endif
      if (io/=0) exit
      nlines = nlines + 1
    enddo
    close (1)

    allocate(stations(nlines))
    allocate(longitudes(nlines))
    allocate(latitudes(nlines))
    allocate(heights(nlines))

    open (2, file = fname)
    do i=1,nlines
      read(2,*) stations(i),longitudes(i),latitudes(i),heights(i)
    enddo
    close (2)
  end subroutine readStations

!==========================================================================
! Subroutine writeDelays writes the slant delay sky view to a text file.
!
! Variables:
!
! resdir: output directory
! statid: four-letter station identifier
! zlat: latitude of the station
! zlon: longitude of the station
! zhei: altitude of the station respect to mean sea level
! zstd: the slant delay vector
! outfile: text file for the output
! idate: analysis time
! itime: forecast length
! izen: number of zenith angles
! iazi: number of azimuth angles
!
! ia: loop index over azimuth angles
! iz: loop index over zenith angles
! eofr: marker for end of record
! fname_out: output directory + output file
! counter: helper index to decompose the slant delay vector
! i: loop index
! exist: a logical variable to check whether the output file exists
! istart: start index of the garbage check
! istop: end index of the garbage check
!
! Logic:
! First check if the slant delay vector contains any unmeaningful
! values and then set those to zero.
! Then write some metadata at the beginning of each record and then
! decompose the slant delay vector so that the zenith delay is
! written first and then each azimuth circle beginning from 1 degree.
! The azimuth circles begin from 1 deg and run to 360 deg.
  subroutine writeDelays2(resdir, statid, zlat, zlon, zhei, zstd, outfile, idate, itime, izen, iazi, zen_angles, azi_angles)
    implicit none
    character*4 statid,eofr
    integer :: ia, iz, idate, itime
    integer :: izen
    integer :: iazi
    real,dimension(nobs) :: zlat, zlon, zhei, zstd
    character(len=140) :: fname_out
    character(len=90) :: resdir
    character(len=25) :: outfile
    real,dimension(iazi) :: zstd_a
    integer :: counter,i
    logical :: exist
    integer :: istart, istop
    real, dimension(izen) :: zen_angles
    real, dimension(iazi) :: azi_angles
    counter=2
    eofr='EOR'

    istop=size(zstd)
    istart=istop-2*iazi
    do i=istart,istop
      if (zstd(i)<0.0 .or. zstd(i)>100.0 .or. zstd(i)/=zstd(i)) then
        zstd(i)=0.0
        print*,'WARNING: inappropriate value at',i,', Setting to zero.'
      endif
    enddo
    fname_out=trim(resdir) // trim(outfile)
    inquire(file=fname_out, exist=exist)

    if (exist) then
      open(12, file=fname_out, status="old", position="append", action="write")
    else
      open(12, file=fname_out, status="new", action="write")
    end if
    write(12,'(A4,3(3X,F23.18),3X,I8,3X,I2)') statid,zlat(1),zlon(1),zhei(1),idate,itime
    write(12,'(F11.6)') zstd(1)
    write(12,'(361(F11.6,3X))') -999.9999,azi_angles !###
    !write(12,'(3601(F11.6,3X))') -999.9999,azi_angles !###
    do iz = 1, izen
      do ia = 1, iazi   
        zstd_a(ia)=zstd(counter)
        counter=counter+1      
      enddo
      write(12,'(361(F11.6,3X))') zen_angles(iz),zstd_a
    enddo
    write(12,'(1X,A4)') eofr

  end subroutine writeDelays2

!#########################################################################################
!#########################################################################################
!#########################################################################################

!==========================================================================
! Subroutine CYCLEID returns a string of 10 characters of the
! form yyyymmddhh, which identifies the current analysis
! cycle.
!
  SUBROUTINE cycleid(year, month, day, hour, cycid)
    IMPLICIT NONE
    INTEGER :: year, month, day, hour
    CHARACTER(LEN=10) :: cycid
    IF (year<0 .OR. year>9999 .OR. &
         month<1 .OR. month>12 .OR. &
         day<1 .OR. day>31 .OR. &
         hour<0 .OR. hour>24) THEN
       cycid='-undefined'
    ELSE
       WRITE (cycid,'(I4,3A2)') year, ch2(month), ch2(day), ch2(hour)
    END IF
  END SUBROUTINE cycleid

!==========================================================================
! The subroutine slantGPS_LTT is the main subroutine of the ray tracing
! based slant delay forward model.
!
  SUBROUTINE slantGPS_ltt( &
       maxobs, nobs_tot, &
       stid, lat, lon, hei, aa_geom, za_geom, &
       nx, ny, nlev, dlon, &
       afull, bfull, lnps_bg, fis_bg, t_bg, q_bg, &
       nhor, dvecrad, lood, sdelay_nwp, nm_lats, dvec_lon, dvec_lat, m_lats )
    IMPLICIT NONE
    INTEGER :: maxobs, nobs_tot
    CHARACTER(LEN=4), DIMENSION(maxobs) :: stid
    REAL, DIMENSION(maxobs) :: lat, lon, hei
    REAL, DIMENSION(maxobs) :: aa_geom, za_geom
    INTEGER :: nx, ny, nlev
    REAL :: dlon
    REAL, DIMENSION(nlev+1) :: afull, bfull
    REAL, DIMENSION(nx,ny) :: lnps_bg, fis_bg
    REAL, DIMENSION(nx,ny,nlev) :: t_bg, q_bg
    INTEGER :: nhor
    REAL :: dvec_lon, dvecrad, dzmax
    real,dimension(nm_lats-1) :: dvec_lat
    INTEGER, DIMENSION(maxobs) :: lood
    REAL, DIMENSION(maxobs) :: sdelay_nwp
    INTEGER, DIMENSION(maxobs) :: nhor_used
    REAL, DIMENSION(maxobs,nhor) :: preglat, preglon
    REAL, DIMENSION(maxobs,nlev+1,nhor) :: z_pr, p_pr, t_pr, q_pr, n_pr, rad_pr
    ! REAL, DIMENSION(maxobs,nlev+1,nhor) :: clwc_pr, ciwc_pr, cswc_pr, crwc_pr
    REAL, DIMENSION(maxobs,nhor) :: locrad
    REAL, DIMENSION(nlev+1,nhor) :: refrac, radius
    INTEGER :: jobs, jlev, jhor
    REAL :: theta, rrec, rsat
    LOGICAL, PARAMETER :: timing=.FALSE.
    INTEGER, DIMENSION(8) :: t
    REAL :: t0, t1, t2, t3, t4, t5
    integer :: i, nm_lats
    real,dimension(nm_lats) :: m_lats
!
    IF (nobs_tot<=0) THEN
       RETURN
    END IF

    IF (timing) THEN
       OPEN (101,FILE='timing_slantGPS_ltt',STATUS='replace')

       CALL DATE_AND_TIME(values=t)
       t0 = t(5)*3600.0+t(6)*60.0+t(7)*1.0+t(8)*0.001
    END IF

! Define the number of actually needed vertical profiles for each observation: nhor_used
    IF (afull(1)>999.0) THEN
       dzmax=35000.0
    ELSE
       dzmax=100000.0
    END IF
    sdelay_nwp=0.0
    DO jobs=1,nobs_tot
       !sdelay_nwp(jobs)=0.0
       !nhor_used(jobs) = NINT(dzmax*tandeg(za_geom(jobs)) / dvec)+2
       nhor_used(jobs) = NINT(dzmax*tandeg(za_geom(jobs)) / dvec_lon)+2
       IF (nhor_used(jobs)>nhor) nhor_used(jobs)=nhor
    END DO

    IF (timing) THEN
       CALL DATE_AND_TIME(values=t)
       t1 = t(5)*3600.0+t(6)*60.0+t(7)*1.0+t(8)*0.001
    END IF

! Project the satellite - receiver vector along the surface of the ellipsoid
    CALL pathprojection2( &
         maxobs, nobs_tot, nhor, nhor_used, lat, lon, aa_geom, nm_lats, dvec_lon, dvec_lat, m_lats, &
         preglat, preglon )

    IF (timing) THEN
       CALL DATE_AND_TIME(values=t)
       t2 = t(5)*3600.0+t(6)*60.0+t(7)*1.0+t(8)*0.001
    END IF

! Interpolate the NWP fields to the projected two dimensional grid
    CALL twodimtableinterp( &
         nx, ny, nlev, dlon, &
         afull, bfull, lnps_bg, fis_bg, t_bg, q_bg, &
        !  clwc_bg, ciwc_bg, cswc_bg, crwc_bg, &
         maxobs, nobs_tot, nhor, nhor_used, lood, preglat, preglon, aa_geom, &
         z_pr, p_pr, t_pr, q_pr, locrad, rad_pr, nm_lats, m_lats)

    IF (timing) THEN
       CALL DATE_AND_TIME(values=t)
       t3 = t(5)*3600.0+t(6)*60.0+t(7)*1.0+t(8)*0.001
    END IF

! Calculate the refractivity in the two dimensional grid
    CALL twodimrefractivity( &
         maxobs, nobs_tot, nlev, nhor, nhor_used, p_pr, t_pr, q_pr, n_pr)

    IF (timing) THEN
       CALL DATE_AND_TIME(values=t)
       t4 = t(5)*3600.0+t(6)*60.0+t(7)*1.0+t(8)*0.001
    END IF

! Call to the ray tracer (by Sean Healy); only one observation at a time
    rsat = 26608.293
    DO jobs=1,nobs_tot
       DO jlev=1,nlev+1
          DO jhor=1,nhor_used(jobs)
             refrac(jlev,jhor) = n_pr(jobs,nlev-jlev+2,jhor)
             radius(jlev,jhor) = rad_pr(jobs,nlev-jlev+2,jhor)
          END DO
       END DO
       rrec = 0.001*(locrad(jobs,1)+hei(jobs))
       theta = degtor(za_geom(jobs)) - &
            asin((rrec/rsat)*sin(degtor(za_geom(jobs))))

       CALL calc_slantRK( &
            nlev+1, nhor_used(jobs), dvecrad, &
            degtor(za_geom(jobs)), hei(jobs)-z_pr(jobs,nlev+1,1), &
            z_pr(jobs,1,1), &
            preglat(jobs,1)*(4.0*ATAN(1.0)/180.0), &
            1.0, &
            rsat*1000.0, theta, &
            refrac, radius, sdelay_nwp(jobs))

    END DO
    !print*,'>>>>>>>>>>>>>>>',nlev+1,nhor_used,dvecrad,degtor(za_geom(1)),hei,z_pr(1,nlev+1,1),z_pr(1,1,1),preglat
    !print*,'>>>>>>>>>>>>>>>',rsat,theta,refrac,radius,sdelay_nwp
    !print*,'>>>>>>>>>>>>>>>',rad_pr(:,1,:)

    IF (timing) THEN
       CALL DATE_AND_TIME(values=t)
       t5 = t(5)*3600.0+t(6)*60.0+t(7)*1.0+t(8)*0.001

       WRITE (101,'(A,3(1x,F12.6))') 'Step 1: ', t1-t0, 100*(t1-t0)/(t5-t0)
       WRITE (101,'(A,3(1x,F12.6))') 'Step 2: ', t2-t1, 100*(t2-t1)/(t5-t0)
       WRITE (101,'(A,3(1x,F12.6))') 'Step 3: ', t3-t2, 100*(t3-t2)/(t5-t0)
       WRITE (101,'(A,3(1x,F12.6))') 'Step 4: ', t4-t3, 100*(t4-t3)/(t5-t0)
       WRITE (101,'(A,3(1x,F12.6))') 'Step 5: ', t5-t4, 100*(t5-t4)/(t5-t0)
       CLOSE (101)
    END IF

  END SUBROUTINE slantGPS_ltt

!==========================================================================
! Subroutine PATHPROJECTION calculates the coordinates along the projection
! of the vector connecting a satellite with the receiver. The vector is
! projected on the surface of the ellipsoid representing the Earth. Input
! consists of the receiver coordinates and the satellite azimuth angle; the
! output is a set of coordinates separated by the distance corresponding to
! the grid spacing of the NWP model.
!

  SUBROUTINE pathprojection2( &
       maxobs, nobs, nhor, nhor_used, relat, relon, recazan, nm_lats, dvec_lon, dvec_lat, m_lats, preglat, preglon )
    IMPLICIT NONE
!
! Input (mainly observational) data
!
    INTEGER :: maxobs, nobs, nhor, nm_lats
    INTEGER, DIMENSION(maxobs) :: nhor_used
    REAL, DIMENSION(maxobs) :: relat, relon, recazan
    real, dimension(nm_lats) :: m_lats
    real, dimension(nm_lats-1) :: dvec_lat
    REAL :: dvec_lon
!
! Output: projection of the satellite - receiver geometry at the ellipsoid
!
    REAL, DIMENSION(maxobs,nhor) :: preglat, preglon
!
! Local variables
!
    INTEGER :: jhor, jobs, i
    REAL :: rcurv, dn, de, dlat, dlon
    REAL, DIMENSION(nhor) :: pazloc
!
!
!
! Loop over the observations
!
    DO jobs=1,nobs

       DO jhor=1,nhor ! Initialize
          preglat(jobs,jhor) = relat(jobs)
          preglon(jobs,jhor) = relon(jobs)
          pazloc(jhor) = recazan(jobs)
       END DO
!
! Loop over the projected grid point coordinates
!
       CALL radicurv( preglat(jobs,1), rcurv )
       i=1
       do
         if (m_lats(i)<preglat(jobs,1)) exit
         i=i+1
       enddo

       DO jhor=2,nhor_used(jobs)

          de   = dvec_lon * sindeg(pazloc(jhor-1))
          dn   = dvec_lat(i) * cosdeg(pazloc(jhor-1))
          dlat = rtodeg(dn/rcurv) !local variable
          dlon = rtodeg(de/(rcurv*cosdeg(preglat(jobs,jhor-1))))
          preglat(jobs,jhor) = preglat(jobs,jhor-1)+dlat
          preglon(jobs,jhor) = preglon(jobs,jhor-1)+dlon
          pazloc(jhor) = pazloc(jhor-1) !!!

       END DO

    END DO

  END SUBROUTINE pathprojection2
!==========================================================================
! Subroutine TWODIMTABLEINTERP projects the model fields from the model
! grid to the two-dimensional plane defined by two points, e.g. GPS
! satellite and the receiver coordinates, and the center of the Earth.
!
! Variables:
! *** Subroutine call
! nlon: number of longitudes
! nlat: number of latitudes
! nlev: number of vertical levels
! dlon: longitude increment
! afull: hybrid A coefficient, vertical sigma coordinate
! bfull: hybrid B coefficient, vertical sigma coordinate
! lnps_bg: 2D ln(pressure at the Earth's surface)
! fis_bg: 2D geopotential of the Earth's surface
! t_bg: 3D atmospheric temperature
! q_bg: 3D specific humidity
! c*wc_bg: 3D water in clouds
! c*wc_pr: 2D projection of water in clouds
! maxobs: limit for number of rays
! nobs: actual number of rays
! nhor: maximum number of model columns to be used to interpolate 2D slice
! nhor_used: actual number of model columns used to interpolate 2D slice
! lood
! preglat: latitudes of the model columns
! preglon: longitudes of the model columns
! aa_geom: geometric azimuth angle of the ray
! z_pr: geometric heights of the model levels
! p_pr: pressure at model levels
! nm_lats: number of latitudes in OpenIFS
! m_lats: list of latitudes of OpenIFS
  SUBROUTINE twodimtableinterp( &
       nlon, nlat, nlev, dlon, &
       afull, bfull, lnps_bg, fis_bg, t_bg, q_bg, &
       maxobs, nobs, nhor, nhor_used, lood, preglat, preglon, aa_geom, &
       z_pr, p_pr, t_pr, q_pr, locrad, rad_pr, nm_lats, m_lats)
    IMPLICIT NONE
!
! Model grid specifics
!
    INTEGER :: nlon, nlat, nlev
    REAL :: dlon
    REAL, DIMENSION(nlev+1) :: afull, bfull
!
! Model grid fields
!
    REAL, DIMENSION(nlon,nlat) :: lnps_bg, fis_bg
    REAL, DIMENSION(nlon,nlat,nlev) :: t_bg, q_bg
    ! REAL, DIMENSION(nlon,nlat,nlev) :: clwc_bg, ciwc_bg, cswc_bg, crwc_bg
!
! The observational input specifics
!
    INTEGER :: maxobs, nobs, nhor
    INTEGER, DIMENSION(maxobs) :: nhor_used, lood
    REAL, DIMENSION(maxobs,nhor) :: preglat, preglon
    REAL, DIMENSION(maxobs) :: aa_geom  ! azimuthal angle
!
! The output projected arrays
!
    REAL, DIMENSION(maxobs,nlev+1,nhor) :: z_pr, p_pr, t_pr, q_pr, rad_pr
    ! REAL, DIMENSION(maxobs,nlev+1,nhor) :: clwc_pr, ciwc_pr, cswc_pr, crwc_pr
    REAL, DIMENSION(maxobs,nhor) :: locrad

    integer :: nm_lats, i=1
    real,dimension(nm_lats) :: m_lats
!
! Local variables
!
    INTEGER :: jobs, jhor, jlev, jx, jy, jiter, iiter
    REAL :: gridx, gridy, da, db, dc, dd, gridxpre, gridypre
    REAL :: lnps_ip, fis_ip, zs_ip, tlr
    REAL :: phl, pf, phu, dlnp, zhl, zf, zhu, dz, tm, qm, tvm
    REAL, PARAMETER :: Rd=287.04 ! R_dry = R/M_dry [J/kg/K]
    REAL, PARAMETER :: eps=0.622 ! M_h2o/M_dry
    REAL :: rcurv

    real :: gravity1, gravity2, dz1, dz2, p_upper, p_lower
    ! real :: mult

    locrad=0.0
    z_pr=0.0
!
!
!
! Loop over the observations
!
    obloop: DO jobs=1,nobs
!
! Loop over the 2d array in horizontal
!

! calculate the radius of curvature given lat. and azimuthal angle

       rcurv=radcurv(preglat(jobs,1),aa_geom(jobs))
      !  locrad(:,:)=rcurv
       !print*,'>>>>>>>rcurv',rcurv

       lood(jobs)=0
       !if (nhor_used(jobs) > 30) print*, 'Ncolumns=', nhor_used(jobs), 'jobs=', jobs
       locloop: DO jhor=1,nhor_used(jobs)
!
! Find the grid point indeces corresponding to the geographical coordinates
!
          if (preglon(jobs,jhor)<0.0) then
            gridx = (preglon(jobs,jhor)+360.0)/dlon
          else
            gridx = preglon(jobs,jhor)/dlon!+1.0
          endif
          do
            if (m_lats(i)<preglat(jobs,jhor)) exit
            i=i+1
          enddo
          !if (preglat(jobs,jhor) < -90) print*, 'PregLat=', preglat(jobs,jhor), 'jobs=', jobs, 'jhor=', jhor, 'i=', i
          gridy=float(i-1)+(m_lats(i-1)-preglat(jobs,jhor))/(m_lats(i-1)-m_lats(i))
          i=1
          IF (gridy<1.0 .OR. gridy>nlat*1.0 .OR. gridx<1.0 .OR. gridx>nlon*1.0) THEN
!             WRITE (*,'(A,2(A,I6))') 'The projection goes beyond the grid domain!', &
!                  '   jobs=', jobs, '   jhor=', jhor
             IF (lood(jobs)==0) lood(jobs)=jhor
             IF (jhor==1) THEN
                CYCLE obloop
             ELSE
                gridx=gridxpre
                gridy=gridypre
             END IF
          ELSE
             gridxpre=gridx
             gridypre=gridy
          END IF
          jx=FLOOR(gridx); jy=FLOOR(gridy)
          da=gridx-jx*1.0; db=gridy-jy*1.0; dc=(jx+1)*1.0-gridx; dd=(jy+1)*1.0-gridy
!
! Interpolate the surface variables:
!
! * Surface pressure and its logarithm
!
          lnps_ip = &
               dd*dc*lnps_bg(jx,jy) + &
               dd*da*lnps_bg(jx+1,jy) + &
               db*dc*lnps_bg(jx,jy+1) + &
               db*da*lnps_bg(jx+1,jy+1)
          p_pr(jobs,nlev+1,jhor) = EXP(lnps_ip)
!
! * Surface geopotential
!
          fis_ip = &
               dd*dc*fis_bg(jx,jy) + &
               dd*da*fis_bg(jx+1,jy) + &
               db*dc*fis_bg(jx,jy+1) + &
               db*da*fis_bg(jx+1,jy+1)
!
! * Surface geometric height (above the mean sea level)
!
          zs_ip=0.0

          DO jiter=1,2
            zs_ip = fis_ip/gravaccel(preglat(jobs,jhor),zs_ip)
          END DO
!
! * Local radius
!
          locrad(jobs,jhor)=rcurv
          !locrad(jobs,jhor)=radcurv(preglat(jobs,jhor),aa_geom(jobs))
          !print*,'>>>>>>>>locrad before',locrad
!
! Initialize geometrical heights of the lowest model level and surface
!
          z_pr(jobs,nlev+1,jhor)=zs_ip
          z_pr(jobs,nlev,jhor)  =z_pr(jobs,nlev+1,jhor)+10.0
!
! * Pressure temperature, water content and specific humidity at the model levels
!
          DO jlev=1,nlev
             p_upper=afull(jlev) + bfull(jlev)*EXP(lnps_ip)
             p_lower=afull(jlev+1) + bfull(jlev+1)*EXP(lnps_ip)
             p_pr(jobs,jlev,jhor) = (p_upper+p_lower)/2.0
             t_pr(jobs,jlev,jhor) = &
                  dd*dc*t_bg(jx,jy,jlev) + &
                  dd*da*t_bg(jx+1,jy,jlev) + &
                  db*dc*t_bg(jx,jy+1,jlev) + &
                  db*da*t_bg(jx+1,jy+1,jlev)
            !  clwc_pr(jobs,jlev,jhor) = &
            !       dd*dc*clwc_bg(jx,jy,jlev) + &
            !       dd*da*clwc_bg(jx+1,jy,jlev) + &
            !       db*dc*clwc_bg(jx,jy+1,jlev) + &
            !       db*da*clwc_bg(jx+1,jy+1,jlev)
            !  ciwc_pr(jobs,jlev,jhor) = &
            !       dd*dc*ciwc_bg(jx,jy,jlev) + &
            !       dd*da*ciwc_bg(jx+1,jy,jlev) + &
            !       db*dc*ciwc_bg(jx,jy+1,jlev) + &
            !       db*da*ciwc_bg(jx+1,jy+1,jlev)
            !  cswc_pr(jobs,jlev,jhor) = &
            !       dd*dc*cswc_bg(jx,jy,jlev) + &
            !       dd*da*cswc_bg(jx+1,jy,jlev) + &
            !       db*dc*cswc_bg(jx,jy+1,jlev) + &
            !       db*da*cswc_bg(jx+1,jy+1,jlev)
            !  crwc_pr(jobs,jlev,jhor) = &
            !       dd*dc*crwc_bg(jx,jy,jlev) + &
            !       dd*da*crwc_bg(jx+1,jy,jlev) + &
            !       db*dc*crwc_bg(jx,jy+1,jlev) + &
            !       db*da*crwc_bg(jx+1,jy+1,jlev)
             q_pr(jobs,jlev,jhor) = &
                  dd*dc*q_bg(jx,jy,jlev) + &
                  dd*da*q_bg(jx+1,jy,jlev) + &
                  db*dc*q_bg(jx,jy+1,jlev) + &
                  db*da*q_bg(jx+1,jy+1,jlev)
          END DO
! convert water content units from [kg/kg] to [g/m3]    See: https://nwp-saf.eumetsat.int/site/download/documentation/rtm/docs_rttov12/rttov_gas_cloud_aerosol_units.pdf
          ! DO jlev=1,nlev
          !   mult = 1000*p_pr(jobs,jlev,jhor)/t_pr(jobs,jlev,jhor)/(Rd * (1+(1-eps)/eps*q_pr(jobs,jlev,jhor)))

          !   clwc_pr(jobs,jlev,jhor) = clwc_pr(jobs,jlev,jhor) * mult
          !   ciwc_pr(jobs,jlev,jhor) = ciwc_pr(jobs,jlev,jhor) * mult
          !   cswc_pr(jobs,jlev,jhor) = cswc_pr(jobs,jlev,jhor) * mult
          !   crwc_pr(jobs,jlev,jhor) = crwc_pr(jobs,jlev,jhor) * mult
          ! END DO
!
! Surface specific humidity
!
          q_pr(jobs,nlev+1,jhor) = q_pr(jobs,nlev,jhor)
!
! Iterate surface temperature and geometric height of the lowest model level
!
          tlr=0.0065
          if (1==1) then
            DO iiter=1,2
!
! * Surface temperature
!
              t_pr(jobs,nlev+1,jhor) = &
                  dd*dc*t_bg(jx,jy,nlev) + &
                  dd*da*t_bg(jx+1,jy,nlev) + &
                  db*dc*t_bg(jx,jy+1,nlev) + &
                  db*da*t_bg(jx+1,jy+1,nlev) + &
                  tlr * (z_pr(jobs,nlev,jhor)-z_pr(jobs,nlev+1,jhor))

              zhl = zs_ip
              phl = EXP(lnps_ip)

              pf  = p_pr(jobs,nlev,jhor)
              phu = 2*pf-phl
              tm  = t_pr(jobs,nlev,jhor)
              qm  = q_pr(jobs,nlev,jhor)
              tvm = (1.0-qm+qm/0.622)*tm

              zf  = zhl
              zhu = zhl
              DO jiter=1,2
                dlnp = LOG(pf/phl)
                dz   = -(Rd*tvm/gravaccel(preglat(jobs,jhor),0.5*zhl+0.5*zf))*dlnp
                zf   = zhl+dz
                dlnp = LOG(phu/phl)
                dz   = -(Rd*tvm/gravaccel(preglat(jobs,jhor),0.5*zhl+0.5*zhu))*dlnp
                zhu  = zhl+dz
              END DO

            END DO
          endif

! Calculation of the geometric model level heights (taken from the slfwd_fmi_v2.f90)
          if (1==1) then
            DO jiter=1,3
              DO jlev=nlev-1,1,-1
                ! phl=afull(jlev+1)*1.0 + bfull(jlev+1)*p_pr(jobs,nlev+1,jhor)
                ! phu=afull(jlev  )*1.0 + bfull(jlev  )*p_pr(jobs,nlev+1,jhor)
                phl=p_pr(jobs,jlev+1,jhor)
                phu=p_pr(jobs,jlev,jhor)
                dlnp=LOG(phu/phl)
                tm=0.5*t_pr(jobs,jlev+1,jhor)+0.5*t_pr(jobs,jlev,jhor)
                qm=0.5*q_pr(jobs,jlev+1,jhor)+0.5*q_pr(jobs,jlev,jhor)
                tvm=(1.0-qm+qm/0.622)*tm

                !print*,jobs
                if (jiter==1) then
                  dz=-(Rd*tvm/gravaccel(preglat(jobs,jhor),0.5*z_pr(jobs,jlev+2,jhor)+0.5*z_pr(jobs,jlev+1,jhor)))*dlnp
                else
                  dz=-(Rd*tvm/gravaccel(preglat(jobs,jhor),0.5*z_pr(jobs,jlev+1,jhor)+0.5*z_pr(jobs,jlev,jhor)))*dlnp
                endif

                ! dz=-(Rd*tvm/gravaccel(preglat(jobs,jhor),0.5*z_pr(jobs,jlev+1,jhor)+0.5*z_pr(jobs,jlev,jhor)))*dlnp
                
                z_pr(jobs,jlev,jhor) = z_pr(jobs,jlev+1,jhor)+dz
                !print*,'>>>>',jhor,jiter,jlev,z_pr(jobs,jlev+1,jhor),z_pr(jobs,jlev,jhor),dz
              END DO
            END DO
          endif
!
! * Radius at the model levels
!
          !print*,'>>>>>>>>>>>locrad after',locrad
          DO jlev=1,nlev+1
             rad_pr(jobs,jlev,jhor) = &
                  locrad(jobs,jhor) + z_pr(jobs,jlev,jhor)
             !print*,jlev,rad_pr(jobs,jlev,jhor),locrad(jobs,jhor),z_pr(jobs,jlev,jhor)
          END DO

       END DO locloop

    END DO obloop

  END SUBROUTINE twodimtableinterp


  REAL FUNCTION radcurv(lat_in_deg,azim_in_deg)

    IMPLICIT NONE

    REAL :: lat_in_deg,azim_in_deg

! local 

    REAL, PARAMETER :: R_q=6378138.0
    REAL, PARAMETER :: R_p=6356752.0
    REAL :: lat,azim
    REAL :: R_ns,R_ew

! conver to radians

    lat = lat_in_deg*ASIN(1.0)/90.0
    azim = azim_in_deg*ASIN(1.0)/90.0

! calculate the radius of curvature

    R_ns=(R_q*R_p)**2/((R_q*COS(lat))**2+(R_p*SIN(lat))**2)**1.5
    R_ew=R_q**2/SQRT((R_q*COS(lat))**2+(R_p*SIN(lat))**2)

    radcurv=1.0/(COS(azim)**2/R_ns + SIN(azim)**2/R_ew)


    RETURN

  END FUNCTION radcurv


!==========================================================================
! Subroutine LOCRADIUS calculates the local radius of the Earth as a
! function of latitude.
!
  SUBROUTINE locradius(lat,radius)
    IMPLICIT NONE
    REAL :: lat, radius, pi, sinlat, coslat
    REAL, PARAMETER :: a=6378137.0, f=1.0/298.257223563, b=a*(1.0-f)
    pi=4.0*ATAN(1.0)
    sinlat=SIN(lat*pi/180.0)
    coslat=COS(lat*pi/180.0)
    radius = a*b/SQRT(b**2*coslat**2 + a**2*sinlat**2)
  END SUBROUTINE locradius


!==========================================================================
! Subroutine RADICURV calculates the local radius of Earth curvature as a
! function of latitude.
!
  SUBROUTINE radicurv(phi,radeff)
    IMPLICIT NONE
    REAL :: phi, radeff
    REAL, PARAMETER :: a=6378137.0, f=1.0/298.257223563, pi=3.14159265359
    REAL :: b, e2
    b=a*(1.0-f)
    e2=(a**2-b**2)/a**2
    radeff=a/SQRT(1.0-e2*(SIN(pi*phi/180.0))**2)
  END SUBROUTINE radicurv


!**************************************************************************
! calculate slant delays

  SUBROUTINE calc_slantRK(nlev,        & ! no. of vertical levels
       nhoriz,       & ! no. of horizontal locations
       dsep,         & ! the angular spacing of horizontal grid (radians)
       zenith_angle, & ! angle to zenith in radians 
       z_height,     & ! height of receiver above the surface
       z_top,        & ! height of the model top above the surface
       lat_top,      & ! (approximate) latitude at the model top
       p_top,        & ! pressure at the model top
       r_gps,        & ! radius of gps
       theta_gps,    & ! theta of GPS. (theta=0.0 for receiver)	
       refrac,       & ! 2d refractivity
       radius,       & ! 2d radius values
       slant)          ! slant delay in m

! Description:
!
! History:
!
! Version Date     Comment
! ------- -------- -------
! 1.0     ??/??/?? Original version. Sean Healy
!
! Code Description:
!   Language:		Fortran 90 
!   Software Standards: GTDP 8
!
!
! End of header -------------------------------------------------------		   

    IMPLICIT NONE

!
! subroutine args. 
!

    INTEGER, INTENT(IN)  :: nlev           ! no. of refractivity levels
    INTEGER, INTENT(IN)  :: nhoriz         ! no. of horizontal locations
    REAL,    INTENT(IN)  :: dsep           ! angular spacing of grid
    REAL,    INTENT(IN)  :: zenith_angle   ! Geometrical zenith angle **in radians**
    REAL,    INTENT(IN)  :: z_height       ! height of receiver above surface   
    REAL,    INTENT(IN)  :: z_top          ! height of the model top above the surface
    REAL,    INTENT(IN)  :: lat_top        ! latitude at the model top
    REAL,    INTENT(IN)  :: p_top          ! pressure at the model top
    REAl,    INTENT(IN)  :: r_gps,theta_gps ! m and radians
    REAL,    INTENT(IN)  :: refrac(nlev,nhoriz)   ! refractivity values on levels
    REAL,    INTENT(IN)  :: radius(nlev,nhoriz)   ! radius values
    REAL,    INTENT(OUT) :: slant   ! slant delay in m
            
!
! local variables
!

    INTEGER :: i,j,n,ibot,iray,k,kp1,idum,ilev,ihor
    INTEGER, PARAMETER :: msplit = 10
    REAL, parameter :: pi=3.14159265359
    REAL :: rad
    REAL :: hwt1,hwt2
    REAL :: amult
    REAL :: h,h2,huse
    REAL :: y(4),yt(4)
    REAL :: dydh(4),dydht(4,4)
    REAL :: theta,theta_tan,theta_min,theta_max
    REAL :: dr_max,dr_dtheta,rtan,dr
    REAL :: kval
    REAL :: slant_tmp(10),theta_ray(10)
    REAL :: map_fac,zenith_delay,sdelay_top,alpha,wt,ref
    REAL :: hdc, zf

    REAL :: za_start
    INTEGER, PARAMETER :: r128 = SELECTED_REAL_KIND(r=128)
    REAL(KIND=r128) :: delta_d,       &    ! difference between true and geometrical paths
                        d,tTOA,ro,beta,d_za ! aux. values
    
    REAL :: dndr, dndt, dndr_sum, dndt_sum


!
! the angles
!
    theta_min = 0.0
    theta_tan = 0.0   ! theta_tan is used in the radio occulation problem - set to zero here
    theta_max = REAL(nhoriz-1)*dsep

!
! set the radius value of the starting point. Radius of surface + height above surface of receiver
!

    rtan = radius(1,1) + z_height

    alpha = 0.0   ! bending angle
    za_start = zenith_angle
    
    DO iray = 1,10 ! number of iterations
   
       amult = 1.0
      
!
! initialise vector
!
       za_start = za_start - sqrt(0.5)*alpha

       y(1) = 0.0                   ! height above receiver           
       y(2) = 0.0                   ! theta by definition theta=0.0 at receiver
       y(3) = za_start              ! thi - bending
       y(4) = 0.0                   ! slant_path delay 

       n=0
!   if (iray == 2) write (6,*) n,y(1)+rtan,y(2),y(3) 

       ibot = 1

! which model levels is the receiver between?

       DO 

          IF (rtan < radius(ibot+1,1)) EXIT 

          ibot=ibot+1

       ENDDO

! now integrate ray-path to top of model atmosphere
        
       DO i = ibot,nlev-1		

! estimate the step length
 
          k = INT((y(2) + theta_tan)/dsep)+1
          k = MIN(MAX(1,k),nhoriz)
          dr_max = (radius(i+1,k)- MAX(radius(i,k),rtan))/REAL(msplit)
          h = dr_max/MAX(COS(y(3)),1.0E-10)
          ! h = dr_max*SIN(y(3))
      
            
! limit to horizontal distance between grid points
    
          h = MIN(h,6.371E6*dsep) 
      
          h2 = 0.5*h
!
! now calculate the path-length with a RUNGE-KUTTA
!
    
          DO  j = 1,msplit 
                    
              yt(:) = y(:)  
      
! first calculation of derivs

              CALL pderivs(nlev,nhoriz,i,dsep,theta_min,theta_max,theta_tan,&
                  rtan,amult,refrac,radius,yt,dydht(:,1),dndr,dndt)
    
              yt(:) = y(:) + dydht(:,1)*h2
            !  dndr_sum = dndr_sum + dndr
            !  dndt_sum = dndt_sum + dndt

! second call at new yt
      
              CALL pderivs(nlev,nhoriz,i,dsep,theta_min,theta_max,theta_tan,&
                  rtan,amult,refrac,radius,yt,dydht(:,2),dndr,dndt)
  
              yt(:) = y(:) + dydht(:,2)*h2
            !  dndr_sum = dndr_sum + dndr
            !  dndt_sum = dndt_sum + dndt
      
! third call at new yt

             CALL pderivs(nlev,nhoriz,i,dsep,theta_min,theta_max,theta_tan,&
                  rtan,amult,refrac,radius,yt,dydht(:,3),dndr,dndt)

              yt(:) = y(:) + dydht(:,3)*h
            !  dndr_sum = dndr_sum + dndr
            !  dndt_sum = dndt_sum + dndt

! fourth last call

             CALL pderivs(nlev,nhoriz,i,dsep,theta_min,theta_max,theta_tan,&
                  rtan,amult,refrac,radius,yt,dydht(:,4),dndr,dndt) 
	    	     	 
             dydh(:) = (dydht(:,1)+dydht(:,4)+2.0*(dydht(:,2)+dydht(:,3)))/6.0
	    
             yt(:) = y(:) + dydh(:)*h
            !  dndr_sum = dndr_sum + dndr
            !  dndt_sum = dndt_sum + dndt
        !print*,dydh
!
! check the radius - have we exited the level
!	    
	    
             k = INT((yt(2) + theta_tan)/dsep)+1
             k = MIN(MAX(1,k),nhoriz-1)
             kp1 = k+1
	 
! horizontal weighting factor	       
	    
             IF ( yt(2) < theta_max .AND. yt(2) >= theta_min) THEN
	    	 
                hwt1 = (REAL(k)*dsep - (yt(2)+theta_tan))/dsep    
                hwt2 = 1.0 - hwt1
	    
             ELSE IF (yt(2) < theta_min) THEN
	    	    
                hwt1 = 1.0
                hwt2 = 0.0
	       
             ELSE IF (yt(2) > theta_max) THEN
	    
                hwt1 = 0.0
                hwt2 = 1.0
	       
             ENDIF
	          	    
	    	
! radius of pressure level	     

             rad = hwt1*radius(i+1,k)+hwt2*radius(i+1,kp1)   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            if (j==msplit) write (6,*) i,j,rad-(yt(1)+rtan)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
	    
! if gone over the boundary scale h
	    
             IF ( j == msplit ) THEN
	    
                dr_dtheta = 0.0
	    
                IF (yt(2) < theta_max .AND. yt(2) > theta_min) &	    
                     dr_dtheta = (radius(i+1,kp1)-radius(i+1,k))/dsep
	
                huse = h - (yt(1)+rtan-rad)/(dydh(1) - dr_dtheta*dydh(2))

             ELSE 

                huse = h  

             ENDIF
	 
!
! update the position vector
!	    	 
             y(:) = y(:) + dydh(:)*huse
	    
!
! output the refrac etc. at the model level.
!	    
             
             IF ( j== msplit .and. iray>1) THEN
    
                yt(:) = y(:)
	     	    
                k = INT((yt(2) + theta_tan)/dsep)+1
                k = MIN(MAX(1,k),nhoriz-1)
                kp1 = k+1
	 
! horizontal weighting factor	       
	    
                IF ( yt(2) < theta_max .AND. yt(2) >= theta_min) THEN
	    	 
                   hwt1 = (REAL(k)*dsep - (yt(2)+theta_tan))/dsep    
                   hwt2 = 1.0 - hwt1
	     
                ELSE IF (yt(2) < theta_min) THEN
	    	    
                   hwt1 = 1.0
                   hwt2 = 0.0
	       
                ELSE IF (yt(2) > theta_max) THEN
	    
                   hwt1 = 0.0
                   hwt2 = 1.0
	       
                ENDIF
   
                rad = hwt1*radius(i+1,k)+hwt2*radius(i+1,kp1)
                ref = hwt1*refrac(i+1,k)+hwt2*refrac(i+1,kp1)
	    
             ENDIF
	    
	    	    
! try to maintain roughly the same radial increment by adjusting h

             IF (j < msplit) THEN
	 
                dr = (rad-y(1)-rtan)/REAL(msplit-j)
	 
                h = MIN(h,dr/MAX(COS(y(3)),1.0E-10))
	    
                h2 = 0.5*h
	    
             ENDIF

          ENDDO  ! complete path thru ith layer

!         if (iray == 2) write (6,*) i+1,rtan+y(1),y(2),y(3)
	 
	 
       ENDDO  ! i the layers

!
! estimate the slant path above model top
!      

       kval = LOG(refrac(nlev-1,k)/refrac(nlev,k))/ &
            (radius(nlev,k)-radius(nlev-1,k))
     
! zenith delay above the model top

       hdc=2.2779E-5
       zf = 1.0 - 0.00266*COS(2.0*lat_top) - 0.00000028*z_top
       zenith_delay = (hdc/zf)*p_top

! calculate mapping factor     
               
       map_fac = mapping_factor((y(1)+rtan),y(3),kval)

! slant delay above model top

       sdelay_top = zenith_delay * map_fac

! total slant delay 
     
       slant_tmp(iray) = y(4) + sdelay_top
       !print*,y(1),y(2),y(3),y(4),sdelay_top
               
!
! estimate the theta value at r_gps, assuming a straightline above the model top
!          
        theta_ray(iray) = y(2) + y(3) - ASIN(SIN(y(3))*(y(1)+rtan)/r_gps)
        alpha = theta_ray(iray) - theta_gps

        !alpha = y(2)+(y(3)-za_start) 

!!     write (6,*) 'alpha=',theta_ray(iray) - theta_gps, &
!!     (y(2)+(y(3)-zenith_angle))*180.0/3.14159 

!
! estimate the length difference between geometrical and bended paths
!

        if (abs(zenith_angle) > 1.0E-4) then
        
        d_za = zenith_angle - za_start
        beta = pi/2 - zenith_angle
        ro = rtan * cos(beta)
        tTOA = acos(ro/(y(1)+rtan)) - beta
        d = (y(1)+rtan) * (tTOA-y(2))

        delta_d = 0.5 * abs(d) * abs(d_za) * cos(y(3))

        !  if (zenith_angle*180.0/3.14159 > 85.0) then
        !  print*, delta_d, (tTOA-y(2)), d_za
        !  endif
        
        slant_tmp(iray) = slant_tmp(iray) + delta_d
        
        endif

    !  if (zenith_angle*180.0/3.14159 > 75.0) then
    !   print*, dndr_sum, dndt_sum
    !  endif

    ENDDO ! iray


    slant = 0.0
    do  j = 0,3
      slant = slant + slant_tmp(10-j)
    enddo
    slant = slant / 4.0

    !print*,slant_tmp(1),slant_tmp(2)

    RETURN

          
  END SUBROUTINE calc_slantRK


  SUBROUTINE pderivs(nlev,   & ! no.of observations
       nhoriz, & ! no. of horizontal layers  ODD
       i, &
       dsep, &
       theta_min, & 
       theta_max, & 
       theta_tan, & 
       rtan, &
       amult, &
       refrac, &
       radius, &
       y, &
       dydh, dndr, dndt)
    
    IMPLICIT NONE

		   
    INTEGER, INTENT(IN)  :: nlev           ! no. of refractivity levels
    INTEGER, INTENT(IN)  :: nhoriz         ! no. of horizontal locations
    INTEGER, INTENT(IN)  :: i 
    REAL,    INTENT(IN)  :: dsep           ! angular spacing of grid
    REAL,    INTENT(IN)  :: theta_min
    REAL,    INTENT(IN)  :: theta_max
    REAL,    INTENT(IN)  :: theta_tan
    REAL,    INTENT(IN)  :: rtan
    REAL,    INTENT(IN)  :: amult          !
    REAL,    INTENT(IN)  :: refrac(nlev,nhoriz)   ! refractivity values on levels
    REAL,    INTENT(IN)  :: radius(nlev,nhoriz)   ! radius values
    REAL,    INTENT(IN)  :: y(4)   !  current location
    REAL,    INTENT(OUT) :: dydh(4)
    REAL,    INTENT(OUT) :: dndr, dndt

! local

    INTEGER :: k,kp1
    REAL :: hwt1,hwt2
    REAL :: ref_up,ref_low
    REAL :: rad_up,rad_low
    REAL :: kval
    
! the easy bits

    dydh(1) = COS(y(3))
    dydh(2) = amult*SIN(y(3))/(y(1)+rtan)


    IF ( y(2) >= theta_min .AND. y(2) <= theta_max) THEN

       k = INT((y(2) + theta_tan)/dsep)+1
       k = MIN(MAX(1,k),nhoriz)
       kp1 = MIN(nhoriz,k+1)

! horizontal weighting factor	       
	       
       hwt1 = (REAL(k)*dsep - (y(2)+theta_tan))/dsep   
       hwt2 = 1.0 - hwt1

    ELSE IF (y(2) < theta_min ) THEN

       k = 1
       kp1 = 2
       hwt1 = 1.0
       hwt2 = 0.0
   
    ELSE IF  (y(2) > theta_max) THEN

       k = nhoriz -1 
       kp1 = nhoriz
       hwt1 = 0.0
       hwt2 = 1.0
   
    ENDIF

! now calculate the local radial gradient of n

    ref_up  = hwt1*refrac(i+1,k)+hwt2*refrac(i+1,kp1)
    ref_low = hwt1*refrac(i,k) + hwt2*refrac(i,kp1)

    rad_up = hwt1*radius(i+1,k)+hwt2*radius(i+1,kp1)
    rad_low = hwt1*radius(i,k) + hwt2*radius(i,kp1)

    kval = LOG(ref_low/ref_up)/(rad_up-rad_low)

    dndr = - 1.0E-6*kval*ref_low*EXP(-kval*(y(1)+rtan-rad_low))

    !dndt = 0.0
    !if (dydh(2) > 0.0) then
    if (k == 1) then
      dndt = 1.0E-6*hwt1*(refrac(i,kp1)-refrac(i,k))/dsep
    else
      dndt = 1.0E-6*(hwt1*(refrac(i,kp1)-refrac(i,k)) + hwt2*(refrac(i,k)-refrac(i,k-1)))/dsep
    endif
    !endif

    dydh(3) = -SIN(y(3))*(1.0/(y(1)+rtan) + dndr) + dndt*COS(y(3))/(y(1)+rtan)

    dydh(4) = 1.0E-6*ref_low*EXP(-kval*(y(1)+rtan-rad_low))


    !dndr = -SIN(y(3))*dndr
    !dndt = dndt*COS(y(3))/(y(1)+rtan)

    RETURN

  END SUBROUTINE pderivs


  REAL FUNCTION mapping_factor(rad,phi,kval)

    REAL, INTENT(IN) :: rad       ! radius value
    REAL, INTENT(IN) :: phi       ! angle with local radius
    REAL, INTENT(IN) :: kval      ! inverse scale-height


    INTEGER :: i
    REAL :: alpha,amult,asym_series,aterm

!
! estimate the mapping factor using an asymptotic expansion in "alpha"
! should not be used if angle > 85 degrees
!

    alpha = TAN(phi)**2/(rad*kval)

    asym_series = 1.0
    aterm = 1.0
    amult = -1.0
    
    i = 1

    DO

       aterm = REAL(2*i-1)*alpha*aterm

       asym_series = asym_series + amult*aterm

       amult = -1.0*amult

       i=i+1

       IF (ABS(aterm/asym_series) < 0.01 .OR. i > 5) EXIT

    ENDDO

    mapping_factor = asym_series/COS(phi)
    !print*,mapping_factor
    RETURN

  END FUNCTION mapping_factor


!==========================================================================
! Subroutine TWODIMREFRACTIVITY calculates the values of total recractivity
! in the model level - signal path -intersection points as a function of
! pressure, temperature and specific humidity values.
!
  SUBROUTINE twodimrefractivity( &
       maxobs, nobs_tot, nlev, nhor, nhor_used, p_pr, t_pr, q_pr, n_pr)
    IMPLICIT NONE
    INTEGER :: maxobs, nobs_tot, nlev, nhor
    INTEGER, DIMENSION(maxobs) :: nhor_used
    REAL, DIMENSION(maxobs,nlev+1,nhor) :: p_pr, t_pr, q_pr, n_pr
    REAL, PARAMETER :: k1=77.607E-2, k2=70.4E-2, k3=3.739E3, k4=1.45, eps=0.622
    INTEGER :: jobs, jlev, jhor
    DO jobs=1,nobs_tot
       DO jlev=1,nlev
          DO jhor=1,nhor_used(jobs)
             n_pr(jobs,jlev,jhor) = &
                  k1*p_pr(jobs,jlev,jhor)/t_pr(jobs,jlev,jhor) + &
                  (k2-k1)*p_pr(jobs,jlev,jhor)*q_pr(jobs,jlev,jhor) &
                  / ((eps+0.378*q_pr(jobs,jlev,jhor))*t_pr(jobs,jlev,jhor)) + &
                  k3*p_pr(jobs,jlev,jhor)*q_pr(jobs,jlev,jhor) &
                  / ((eps+0.378*q_pr(jobs,jlev,jhor))*t_pr(jobs,jlev,jhor)**2) !+ &
                ! liquid water content impact: included ONLY liquid
                ! for others (ice, snow and rain) - need to apply scattering models instead
                  ! k4*(clwc_pr(jobs,jlev,jhor))
            !  if (clwc_pr(jobs,jlev,jhor) > 5) print *, 'CLWC EXTREMELY BIG !!', clwc_pr(jobs,jlev,jhor), &
                  ! ' (g/m3) at location jobs=', jobs, 'jlev=', jlev, 'jhor=', jhor
          END DO
       END DO
       do jhor=1,nhor_used(jobs)
        n_pr(jobs, nlev+1, jhor) = n_pr(jobs, nlev, jhor)
       enddo
    END DO
  END SUBROUTINE twodimrefractivity


!**************************************************************************



!**************************************************************************
!
! ------------ FUNCTIONS ------------
!
!**************************************************************************

  FUNCTION ch2(ipar)
    IMPLICIT NONE
    INTEGER :: ipar
    CHARACTER(LEN=2) :: ch2
    IF (ipar< 0  .OR.  ipar> 99) ch2='00'
    IF (ipar>=0  .AND. ipar< 10) ch2='0' // CHAR(48+ipar)
    IF (ipar>=10 .AND. ipar<=99) ch2=CHAR(48+(ipar/10)) // &
         CHAR(48+MOD(ipar,10))
    RETURN
  END FUNCTION ch2

  REAL FUNCTION gravaccel(glat,gz)
    IMPLICIT NONE
    REAL :: glat, gz
    REAL, PARAMETER :: c1=5.2885e-3,c2=-5.9e-6,c3=-3.086e-6, &
         pifac=1.745329e-2,g00=9.780356
    gravaccel = &
         g00*(1.0e0+c1*(sin(glat*pifac))**2 &
         + c2*(sin(2.0e0*glat*pifac))**2)+c3*gz
  END FUNCTION gravaccel

  REAL FUNCTION degtor(x)
    IMPLICIT NONE
    REAL :: x
    REAL, parameter :: pi=3.14159265359
    degtor=x*pi/180.
  END FUNCTION degtor
  
  REAL FUNCTION rtodeg(x)
    IMPLICIT NONE
    REAL :: x
    REAL, parameter :: pi=3.14159265359
    rtodeg=x*180./pi
  END FUNCTION rtodeg
  
  REAL FUNCTION sindeg(x)
    IMPLICIT NONE
    REAL :: x
    sindeg=sin(degtor(x))
  END FUNCTION sindeg
  
  REAL FUNCTION cosdeg(x)
    IMPLICIT NONE
    REAL :: x
    cosdeg=cos(degtor(x))
  END FUNCTION cosdeg
  
  REAL FUNCTION tandeg(x)
    IMPLICIT NONE
    REAL :: x
    tandeg=tan(degtor(x))
  END FUNCTION tandeg
  
  REAL FUNCTION asindeg(x)
    IMPLICIT NONE
    REAL :: x
    asindeg=rtodeg(asin(x))
  END FUNCTION asindeg
  
  REAL FUNCTION acosdeg(x)
    IMPLICIT NONE
    REAL :: x
    acosdeg=rtodeg(acos(x))
  END FUNCTION acosdeg
  
  REAL FUNCTION atandeg(x)
    IMPLICIT NONE
    REAL :: x
    atandeg=rtodeg(atan(x))
  END FUNCTION atandeg

END PROGRAM slfwd_ltt_v8
