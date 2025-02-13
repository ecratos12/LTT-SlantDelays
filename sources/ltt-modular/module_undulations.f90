module module_undulations

    ! This type represents gridded Geoid undulation data
    ! Geoid undulation is difference between height above mean sea level and ellipsoidal height
    !
    ! N = H - h
    ! N: Geoid undulation
    ! H: Orthometric height (MSL height)
    ! h: Height above the reference ellipsoid
    type, public :: undulationsType

        integer :: nLon, nLat
        double precision, dimension(:,:), allocatable :: values

    end type undulationsType

contains


    subroutine read_undulations(parameters, nwp_domain, undulations)
        use module_data_types, only: parametersType, nwpDomainType

        implicit none
        type(parametersType), intent(in) :: parameters
        type(nwpDomainType), intent(in) :: nwp_domain
        type(undulationsType), intent(out) :: undulations

        character(len=150) :: undulationsFile
        logical :: undFileExists
        double precision, dimension(:), allocatable :: line
        integer :: i,j

        ! access used resolution of NWP fields
        select case (parameters%resolution)
        case ('t639')
            undulationsFile='LTT_conf_files/undul-639-egm08.txt'
        case ('t1279')
            undulationsFile='LTT_conf_files/undul-1279-egm08.txt'
        case default
            print*,'Unknown resolution. Choose supported t639 or t1279.'
            stop
        end select

        undulations%nLon = nwp_domain%gridNLon
        undulations%nLat = nwp_domain%gridNLat
        allocate(undulations%values(undulations%nLon, undulations%nLat))
        undulations%values = 0.0

        inquire(file=undulationsFile, exist=undFileExists)
        if (undFileExists) then
            open (2, file=undulationsFile)
            allocate(line(undulations%nLon))
            do i=1,undulations%nlat
                read (2,*) line
                do j=1,undulations%nLon
                    undulations%values(j,i) = line(j)
                enddo
            enddo
            close(2)
        else
            print *,'Undulations file was not found! Setting undulations to zero.'
        endif

    end subroutine read_undulations


    ! compute station's ellipsoidal height using MSL height and geoid undulation
    subroutine station_undulations(parameters, undulations, nwp_domain, stationsList)
        use module_data_types, only: nwpDomainType, stationsListType, stationType, parametersType

        implicit none
        type(parametersType), intent(in) :: parameters
        type(undulationsType), intent(in) :: undulations
        type(nwpDomainType), intent(in) :: nwp_domain
        type(stationsListType), intent(in out) :: stationsList

        integer :: i, jx,jy, ilat
        double precision :: gridx,gridy, da,db,dc,dd, localUnd

        do i=1,stationsList%size

            ! find the grid point indeces corresponding to the geographical coordinates of the station
            if (stationsList%list(i)%lon < 0.0) then
                gridx = (stationsList%list(i)%lon + 360.0)/nwp_domain%deltaLon + 1.0
            else
                gridx = stationsList%list(i)%lon/nwp_domain%deltaLon + 1.0
            endif

            ilat=1
            do
                if (nwp_domain%latsList(ilat) < stationsList%list(i)%lat) exit
                if (ilat==nwp_domain%gridNLat) exit  ! watch out for grid ending before south pole
                ilat = ilat+1
            enddo
            if (ilat==1) ilat=2     ! watch out for grid ending before north pole
            gridy = float(ilat-1) + (nwp_domain%latsList(ilat-1) - stationsList%list(i)%lat)/ &
                (nwp_domain%latsList(ilat-1) - nwp_domain%latsList(ilat))

            jx = floor(gridx); jy = floor(gridy)
            da = gridx - jx*1.0
            db = gridy - jy*1.0
            dc = (jx+1)*1.0 - gridx
            dd = (jy+1)*1.0 - gridy
            if (jx==0) then
                da = 0.0
                dc = 1.0
                jx = 1
            endif
            if (jy==0) then
                db = 0.0
                dd = 1.0
                jy = 1
            endif

            ! interpolate undulation values to the station location
            localUnd = &
                dd*dc*undulations%values(jx,jy) + &
                dd*da*undulations%values(jx+1,jy) + &
                db*dc*undulations%values(jx,jy+1) + &
                db*da*undulations%values(jx+1,jy+1)
            ! h = H + N
            ! N: Geoid undulation
            ! H: Orthometric height (MSL height)
            ! h: Height above the reference ellipsoid
            if (parameters%use_MSL_heights) then    ! compute elp from MSL
                stationsList%list(i)%elp_height = stationsList%list(i)%MSL_height + localUnd
            else                                    ! compute MSL from elp
                stationsList%list(i)%MSL_height = stationsList%list(i)%elp_height - localUnd
            endif

            ! print *,stationsList%list(i)%name, localUnd, stationsList%list(i)%elp_height
            ! print *,stationsList%list(i)%lat, ',', stationsList%list(i)%lon, ',', stationsList%list(i)%elp_height 

        enddo

    end subroutine station_undulations

end module module_undulations