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
! OpenIFS implementation: Lauri Tuppi and Maksym Vasiuta
! Implementation of dn/dpsi term, geometrical corrections and condensate fields: Maksym Vasiuta
! Modularization and refactoring: Maksym Vasiuta
! date 24.01.2022
!
! Example run call:
! ./slfwd_v9.exe oifs2016120100+000 p000 $inputdir $outputdir 1 66 t1279
!
!--------------------------------------------------------------------------
!
!**************************************************************************


program ltt_v9
    use module_data_types
    use module_io
    use module_skyview
    use module_nwpProjection
    use module_computeDelays
    use module_undulations

    implicit none
    type(parametersType) :: parameters
    type(nwpDomainType) :: oifsDomain
    type(weatherBackgroundDataType) :: oifsFields
    type(SkyDirectionsDataType) :: skyview
    type(stationsListType) :: stations
    type(dateTimeType) :: analysisTime, forecastTime
    type(undulationsType) :: undulations

    ! Other crutial parameters of LTT operator
    !
    ! distance of satellite from centre of the Earth (in m)
    ! finish condition for ray tracer to construct signal path and calculate delay..
    ! ..is r==rsat
    double precision :: rsat=26608293.0
    ! number of unique zenith and azimuth angles to compute and output the delays
    ! if one want to make variable zenith angle resolution, ..
    ! ..provide total number of unique angles here and define them in "constructSkyview" subroutine
    integer :: nZen=85, nAzi=360

    ! Intermediate variables and data
    integer :: ist
    double precision :: sec1, sec2
    integer, dimension(8) :: t


    print*,'###### Welcome to using slfwd slant delay operator. ######'

    
    ! 1. Read command line arguments as LTT parameters
    print*,'Reading the configuration ...'
    call read_parameters(parameters)

    ! 2. Read all observing stations
    print*,'Obtaining list of observing stations ...'
    call read_stations(stations)

    ! 3. Read OpenIFS NWP data domain configuration
    print*,'Reading the GRIB file '//TRIM(parameters%inputFileDir)//TRIM(parameters%inputFileName)//' ...'
    print*,'Obtaining meta information (grid, model levels, etc.) ...'
    call read_nwpConfiguration(parameters, oifsDomain, analysisTime, forecastTime)

    ! 4. Read geoid undulations corresponding the OpenIFS resolution
    print*,'Obtaining geoid undulations ...'
    call read_undulations(parameters, oifsDomain, undulations)
    ! 4.1 Apply geoid undulations to station heights
    call station_undulations(undulations, oifsDomain, stations)

    ! 5. Read OpenIFS NWP data: fields of orography, logarithm of surface pressure, 
    ! temperature, specific humidity and specific cloud liquid water content
    print*,'Obtaining NWP fields ...'
    call read_nwpBackgroundFields(parameters, oifsDomain, oifsFields)

    ! 6. Initialize observations as a set of sky directions
    ! The delays are then computed for signals propagating from a virtual GNSS satellite 
    ! placed at certain direction at certain distance rsat from centre of the Earth (in m).
    print*,'Initialize "observations" as a set of sky directions ...'
    call constructSkyview(nZen, nAzi, rsat, skyview)

    ! 7. Perform slant delays computation
    ! in directions defined by skyview for every station in selected range
    do ist=parameters%startStation,parameters%endStation

        print*,'####### Calculating actual slant delays for station ',stations%list(ist)%name,' #######'
        
        call date_and_time(values=t)
        sec1 = t(5)*3600.0+t(6)*60.0+t(7)*1.0+t(8)*0.001

        call computeAndWriteSkyDelays(&
            parameters, stations%list(ist), skyview, oifsDomain, oifsFields, undulations, forecastTime&
        )

        call date_and_time(values=t)
        sec2 = t(5)*3600.0+t(6)*60.0+t(7)*1.0+t(8)*0.001

        print*,'  ... Done! Time:', sec2-sec1
        
    enddo

    
end program ltt_v9