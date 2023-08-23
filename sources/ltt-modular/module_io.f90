module module_io

contains

    subroutine read_arguments(arguments)
        use module_data_types, only: argumentsType

        implicit none
        type(argumentsType), intent(in out) :: arguments

        integer n

        ! assure that it is possible to read all arguments
        n = iargc()
        if (n<4) then
            write (*, '(A, I1, A)') "Not enough arguments given!  (Needed 4, while only ", n, " has given)"
            stop
        end if

        ! read arguments and fill in "arguments" structure
        call get_command_argument_Allocated(1, arguments % setupFile)
        call get_command_argument_Allocated(2, arguments % inputFileName)
        call get_command_argument_Allocated(3, arguments % inputFileDir)
        call get_command_argument_Allocated(4, arguments % outputFileDir)

    end subroutine read_arguments


    ! work-around procedure to get command line arguments of unknown length
    subroutine get_command_argument_Allocated(iarg, arg)
        implicit none
        integer iarg
        character(len=:), allocatable :: arg
        integer arglen

        call get_command_argument(iarg,length=arglen)
        allocate(character(arglen) :: arg)
        call get_command_argument(iarg,value=arg)
    end subroutine get_command_argument_Allocated


    subroutine read_parameters(arguments, parameters)
        use module_data_types, only: parametersType, argumentsType

        implicit none
        type(argumentsType), intent(in) :: arguments
        type(parametersType), intent(in out) :: parameters
        integer :: nlines=0, i, io
        character(len=64) :: key, value

        ! include program arguments to parameters 
        parameters%ltt_args = arguments

        open(1, file = trim(arguments%setupFile))
        do
            read(1,*,iostat=io)

            if (nlines==0) then
                if (io/=0) then
                    print *,'Cannot read the setup file!'
                    stop
                endif
            endif

            if (io/=0) exit
            nlines = nlines+1

        enddo
        close(1)

        open(2, file = trim(arguments%setupFile))
        do i=1,nlines
            read(2,*) key, value
            print *, key, 'is ', value
            if (key == 'startStation') then
                read(value,*) parameters%startStation

            elseif (key == 'endStation') then
                read(value,*) parameters%endStation

            elseif (key == 'resolution') then
                parameters%resolution = value

            elseif (key == 'dAzimuth_deg') then
                read(value,*) parameters%dAzimuth_deg

            elseif (key == 'zenAngleLimit_deg') then
                read(value,*) parameters%zenAngleLimit_deg

            elseif (key == 'clwc') then
                if (value=='on') then
                    parameters%include_clwc = .true.
                elseif (value=='off') then
                    parameters%include_clwc = .false.
                endif

            elseif (key == 'use_MSL_heights') then
                if (value=='on') then
                    parameters%use_MSL_heights = .true.
                elseif (value=='off') then
                    parameters%use_MSL_heights = .false.
                endif

            else
                print*,'Error! Unknown parameter ',  key
                stop
            endif
        enddo
        close(2)

    end subroutine read_parameters


    subroutine read_stations(stations)
        use module_data_types, only: stationType, stationsListType

        implicit none
        type(stationsListType), intent(in out) :: stations

        character(len=90) :: file = 'LTT_conf_files/stationCoordinates_extd.txt'
        integer :: nlines=0, i, io

        open(1, file = trim(file))
        do
            read(1,*,iostat=io)

            if (nlines==0) then
                if (io/=0) then
                    print *,'Cannot read the list of stations!'
                    stop
                endif
            endif

            if (io/=0) exit
            nlines = nlines+1
        enddo
        close(1)

        allocate(stations%list(nlines))
        stations%file = trim(file)
        stations%size = nlines

        open(2, file = trim(file))
        do i=1,nlines
            read(2,*) stations%list(i)%name, stations%list(i)%lon, stations%list(i)%lat, stations%list(i)%MSL_height
        enddo
        close(2)

    end subroutine read_stations


    subroutine read_nwpConfiguration(parameters, domain, anTime, foTime)
        use module_data_types, only: nwpDomainType, dateTimeType, parametersType
        use eccodes

        implicit none
        type(nwpDomainType), intent(in out) :: domain
        type(dateTimeType), intent(in out) :: anTime, foTime  ! ANalysis and FOrecast datetimes
        type(parametersType), intent(in) :: parameters

        character(len=150) :: latsFileName, gribFileName, vertCoordFile
        integer :: nlines=0, i, io
        integer :: fclen
        integer, dimension(12), parameter :: &
        dpm=    (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /), &
        dpmleap=(/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
        integer :: rfile,iret,messages,level,ierr,placeholder
        integer, dimension(:), allocatable :: igrib


        ! access used resolution of OpenIFS fields
        if (parameters%resolution == 't639') then
            latsFileName='LTT_conf_files/lats_t639.txt'
            vertCoordFile='LTT_conf_files/vert_coord_l91.dat'
        elseif (parameters%resolution == 't1279') then
            latsFileName='LTT_conf_files/lats_t1279.txt'
            vertCoordFile='LTT_conf_files/vert_coord_l137.dat'
        else
            print*,'Unknown resolution. Choose supported t639 or t1279.'
            stop
        endif


        ! read the latitudes of OpenIFS from an ASCII file
        open(1, file = trim(latsFileName))
        do
            read(1,*,iostat=io)

            if (nlines==0) then
                if (io/=0) then
                    print *,'Cannot read the list of OIFS grid latitudes!'
                    stop
                endif
            endif

            if (io/=0) exit
            nlines = nlines+1

        enddo
        close(1)

        allocate(domain%latsList(nlines))
        open(2, file = trim(latsFileName))
        do i=1,nlines
            read(2,*) domain%latsList(i)
        enddo
        close(2)


        ! Define the time of the analysis and forecast field
        ! analysis from OIFS file name
        read(parameters%ltt_args%inputFileName(5:8), '(i4)') anTime%year
        read(parameters%ltt_args%inputFileName(9:10), '(i2)') anTime%month
        read(parameters%ltt_args%inputFileName(11:12), '(i2)') anTime%day
        read(parameters%ltt_args%inputFileName(13:14), '(i2)') anTime%hour

        ! forecast time = analysis time + ForeCast LENgth
        read(parameters%ltt_args%inputFileName(16:18), '(i3)') fclen
        foTime = anTime
        foTime%hour = foTime%hour + fclen
        do while (foTime%hour >= 24)
            foTime%hour = foTime%hour - 24
            foTime%day = foTime%day + 1
        enddo
        IF ( (MOD(anTime%year,4)==0 .AND. MOD(anTime%year,100)/=0) &
            .OR. MOD(anTime%year,400)==0) THEN
            IF (foTime%day > dpm(anTime%month)) THEN
                foTime%day = foTime%day-dpmleap(foTime%month)
                foTime%month = foTime%month+1
            END IF
        ELSE
            IF (foTime%day > dpm(anTime%month)) THEN
                foTime%day = foTime%day-dpm(foTime%month)
                foTime%month = foTime%month+1
            END IF
        END IF


        ! Retrieve the model domain parameters from the given OpenIFS GRIB-file.
        ! 1. Read OpenIFS file, count messages inside.
        gribFileName = trim(parameters%ltt_args%inputFileDir)//trim(parameters%ltt_args%inputFileName)
        call codes_open_file(rfile, gribFileName, 'r')
        call codes_count_in_file(rfile, messages, ierr)
        allocate(igrib(messages+1))

        ! 2. Retrieve identifiers for grib messages.
        i=1
        do while (iret /= CODES_END_OF_FILE)
            call codes_grib_new_from_file(rfile, igrib(i), iret)
            i=i+1
        enddo

        ! 3. Get horizontal field dimensions from some grib message. The message number is arbitrary.
        call codes_get(igrib(1), 'Ni', domain%gridNLon)
        call codes_get(igrib(1), 'Nj', domain%gridNLat)

        ! 4. Retrieve level count by searching largest level number in grib messages.
        ! Leveling by pressure value (1000hPa, 500hPa, etc.) VS leveling by number..
        ! .. solved by excluding large (>150) numbers..
        ! .. FOR NOW
        i=1
        domain%nLevels = 0
        do while (i<=messages)
            call codes_get(igrib(i), 'level', level)
            if (level<150 .and. domain%nLevels<level) then
                domain%nLevels = level
            endif
            i=i+1
        enddo

        ! 5. Grid spacing in degrees. Longitudal increment is a single number.
        ! Latitudal increment is varying.
        domain%deltaLon = 360./real(domain%gridNLon)

        allocate(domain%deltaLat(domain%gridNLat-1))
        do i=1,domain%gridNLat-1
            domain%deltaLat(i) = domain%latsList(i) - domain%latsList(i+1)
        enddo

        ! 6. Retrieve hybrid A and B coefficients, vertical sigma coordinate
        allocate(domain%afull(domain%nLevels+1))
        allocate(domain%bfull(domain%nLevels+1))
        open(3, file = trim(vertCoordFile))
        do i=1,domain%nLevels+1
            read(3,*) placeholder, domain%afull(i), domain%bfull(i)
        enddo
        close(3)

        ! 7. Memory is freed, and the OIFS file closed.
        do i=1,messages
            call codes_release(igrib(i))
        enddo
        call codes_close_file(rfile)

    end subroutine read_nwpConfiguration


    subroutine read_nwpBackgroundFields(parameters, domain, fields)
        use module_data_types, only: nwpDomainType, parametersType, weatherBackgroundDataType
        use eccodes

        implicit none
        type(parametersType), intent(in) :: parameters
        type(nwpDomainType), intent(in) :: domain
        type(weatherBackgroundDataType), intent(in out) :: fields

        character(len=150) :: gribFileName
        integer :: rfile,iret,messages,level,ierr,nValues,paramId
        integer :: i,jx,jy
        integer, dimension(:), allocatable :: igrib
        double precision, dimension(:), allocatable :: values

        ! 1. Allocate and initiate fields data.
        fields%nLon = domain%gridNLon
        fields%nLat = domain%gridNLat
        fields%nLevels = domain%nLevels

        allocate(fields%fis(fields%nLon, fields%nLat))
        allocate(fields%lnps(fields%nLon, fields%nLat))
        allocate(fields%T(fields%nLon, fields%nLat, fields%nLevels))
        allocate(fields%Q(fields%nLon, fields%nLat, fields%nLevels))
        allocate(fields%clwc(fields%nLon, fields%nLat, fields%nLevels))

        fields%fis=0.0
        fields%lnps=0.0
        fields%T=0.0
        fields%Q=0.0
        fields%clwc=0.0

        ! 2. Read OpenIFS file, count messages inside.
        gribFileName = trim(parameters%ltt_args%inputFileDir)//trim(parameters%ltt_args%inputFileName)
        call codes_open_file(rfile, gribFileName, 'r')
        call codes_count_in_file(rfile, messages, ierr)
        allocate(igrib(messages+1))

        ! 3. Retrieve identifiers for grib messages.
        i=1
        do while (iret /= CODES_END_OF_FILE)
            call codes_grib_new_from_file(rfile, igrib(i), iret)
            i=i+1
        enddo

        ! 4. Retrieve number of elements in the data in the grib message, and allocate
        ! temporary storage vector accordingly. The message number is arbitrary.
        call codes_get(igrib(10),"numberOfPoints",nValues);allocate(values(nValues))

        ! 5. This is the main do loop that reads the data in the grib messages
        do i=1,messages
            call codes_get(igrib(i),"values",values)
            call codes_get(igrib(i),"paramId",paramId)
            call codes_get(igrib(i),"level",level)
            if (level>150) then ! ignore pressure leveling, same way as in "read_nwpConfiguration" subroutine
                cycle
            endif

            select case (paramId)
            case (129) ! 2-d surface geopotential
                do jx=1,fields%nLon
                    do jy=1,fields%nLat
                        fields%fis(jx,jy) = values((jy-1)*fields%nLon+jx)
                    enddo
                enddo
            case (130) ! 3-d temperature
                do jx=1,fields%nLon
                    do jy=1,fields%nLat
                        fields%T(jx,jy,level) = values((jy-1)*fields%nLon+jx)
                    enddo
                enddo
            case (133) ! 3-d specific humidity
                do jx=1,fields%nLon
                    do jy=1,fields%nLat
                        fields%Q(jx,jy,level) = values((jy-1)*fields%nLon+jx)
                    enddo
                enddo
            case (152) ! 2-d logarithm of surface pressure
                do jx=1,fields%nLon
                    do jy=1,fields%nLat
                        fields%lnps(jx,jy) = values((jy-1)*fields%nLon+jx)
                    enddo
                enddo
            case (246) ! 3-d clwc
                do jx=1,fields%nLon
                    do jy=1,fields%nLat
                        fields%clwc(jx,jy,level) = values((jy-1)*fields%nLon+jx)
                    enddo
                enddo
            ! case (someNewId) ! field you wish to add
            !     do jx=1,fields%nLon
            !         do jy=1,fields%nLat
            !             fields%xxxx(jx,jy,level) = values((jy-1)*fields%nLon+jx)
            !         enddo
            !     enddo
            end select
        enddo

        ! 6. Memory is freed, and the OIFS file closed.
        do i=1,messages
            call codes_release(igrib(i))
        enddo
        call codes_close_file(rfile)

    end subroutine read_nwpBackgroundFields


    subroutine writeDelays(parameters, station, delays, forecastTime, skyview)
        use module_data_types

        implicit none
        type(parametersType), intent(in) :: parameters
        type(stationType), intent(in) :: station
        type(SkyDelaysDataType), intent(in out) :: delays
        type(dateTimeType), intent(in) :: forecastTime
        type(SkyDirectionsDataType), intent(in) :: skyview

        integer :: counter,i,fclen,iz,ia,date
        logical :: exist
        character*4 :: eofr='EOR'
        character(len=150) :: outFileName
        character(len=3) :: hour
        character(len=16) :: lineFormat
        double precision, dimension(skyview%nAzimuths) :: dataPerLine

        write(lineFormat,'("(", I0, "(F11.6,3X))")') skyview%nAzimuths+1

        ! check inappropriate values
        do i=1,size(delays%slant)
            if (delays%slant(i)<0.0 .or. delays%slant(i)>100.0 .or. delays%slant(i)/=delays%slant(i)) then
                delays%slant(i) = 0.0
                print*,'WARNING: inappropriate value at',i,', Setting to zero.'
            endif
        enddo
        if (delays%ZTD<0.0 .or. delays%ZTD>100.0 .or. delays%ZTD/=delays%ZTD) then
            print*,'WARNING: inappropriate value at',i,', Setting to zero.'
        endif

        read(parameters%ltt_args%inputFileName(5:12), '(i8)') date

        read(parameters%ltt_args%inputFileName(16:18), '(i3)') fclen
        write(hour,'(i3.3)') fclen
        read(parameters%ltt_args%inputFileName(5:14), '(a10)') outFileName
        outFileName = trim(outFileName) // '+' // trim(hour) // '.dat'
        outFileName = trim(parameters%ltt_args%outputFileDir)//trim(outFileName)

        inquire(file=outFileName, exist=exist)
        if (exist) then
            open(12, file=outFileName, status="old", position="append", action="write")
        else
            open(12, file=outFileName, status="new", action="write")
        end if

        if (parameters%use_MSL_heights) then
            write(12,'(A4,3(3X,F23.18),3X,I8,3X,I2)') &
                    station%name,station%lat,station%lon,station%MSL_height,date,forecastTime%hour
        else
            write(12,'(A4,3(3X,F23.18),3X,I8,3X,I2)') &
                    station%name,station%lat,station%lon,station%elp_height,date,forecastTime%hour
        endif
        write(12,'(F11.6)') delays%ZTD
        write(12,lineFormat) -999.9999,skyview%uniqueAzimuths !###

        counter=1
        do iz = 1, size(skyview%uniqueZenAngles)
            do ia = 1, skyview%nAzimuths   
                dataPerLine(ia) = delays%slant(counter)
                counter=counter+1
            enddo
            write(12,lineFormat) skyview%uniqueZenAngles(iz), dataPerLine
        enddo
        write(12,'(1X,A4)') eofr

    end subroutine writeDelays

end module module_io
