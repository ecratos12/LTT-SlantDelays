
module module_data_types

    ! "argumentsType" represents LTT program arguments
    type, public :: argumentsType

        ! setupFile: input parameters for LTT operator
        ! inputFileName: grib-file with weather data and domain
        ! outputFileDir: directory with computed LTT slant delays

        character(len=:), allocatable :: setupFile
        character(len=:), allocatable :: inputFileName
        character(len=:), allocatable :: inputFileDir
        character(len=:), allocatable :: outputFileDir

    end type argumentsType


    ! "parametersType" represents essential parameters for LTT operator
    type, public :: parametersType

        ! [start,end]Station: go through station list only within this range of indexes
        ! resolution: pre-defined NWP grid resolution (supports only t639, t1279)
        ! dAzimuth_deg: azimuth angle resolution in skyview
        ! zenAngleLimit_deg: maximum zenith angle in skyview
        ! include_clwc: enables liquid water content in refractivity computation
        ! use_MSL_heights: use MSL heights for stations, if .false. -- use ellipsoid instead
        ! ltt_args: arguments of LTT program

        integer :: startStation
        integer :: endStation
        character(len=:), allocatable :: resolution
        double precision :: dAzimuth_deg
        double precision :: zenAngleLimit_deg
        logical :: include_clwc
        logical :: use_MSL_heights
        type(argumentsType) :: ltt_args

    end type parametersType



    ! "dateTimeType" represents date and time
    type, public :: dateTimeType

        ! ccyc: 10-characters identifier - yyyymmddhh

        integer :: year, month, day, hour
        character(len=10) :: ccyc

    end type dateTimeType



    ! "stationType" represents single observing station
    type, public :: stationType

        ! MSL_height: height of station above mean sea level
        ! elp_height: height of station above reference ellipsoid

        character(len=4) :: name
        double precision :: lat
        double precision :: lon
        double precision :: MSL_height
        double precision :: elp_height
    
    end type stationType



    ! "stationsListType" stores all observing stations in use, and refers to source of this list (stations file)
    type, public :: stationsListType

        type(stationType), dimension(:), allocatable :: list
        integer :: size
        character(len=90) :: file

    end type stationsListType



    ! "nwpDomainType" stores information about a domain of input NWP data - horizontal and vertical resolution and configuration
    type, public :: nwpDomainType

        ! latsList: grid of latitude values (degrees)
        ! lonsList: grid of longitude values (degrees)
        ! afull: hybrid A coefficient, vertical sigma coordinate
        ! bfull: hybrid B coefficient, vertical sigma coordinate
        ! deltaLat: intervals between latitude values (degrees)
        ! deltaLon: interval between longitude values (degrees)
        ! gridSize: number of unique longitudes, number of unique latitudes

        double precision, dimension(:), allocatable :: latsList
        double precision, dimension(:), allocatable :: lonsList

        double precision, dimension(:), allocatable :: deltaLat
        double precision :: deltaLon

        integer :: gridNLon, gridNLat
        integer :: nLevels

        double precision, dimension(:), allocatable :: afull, bfull

    end type nwpDomainType



    ! "weatherBackgroundDataType" stores the NWP model fields from analysis or forecast grib file, given in a certain domain
    type, public :: weatherBackgroundDataType

        ! DIMENSION SHAPE IS => number of longitudes X number of latitudes X number of vertical levels
        ! ALSO => number of longitudes X number of latitudes

        integer :: nLon, nLat, nLevels

        ! fis: 2D geopotential of the Earth's surface
        ! lnps: 2D ln(pressure at the Earth's surface)
        ! T: 3D atmospheric temperature
        ! Q: 3D specific humidity
        ! c*wc: 3D non-gaseous water content (l-liquid, i-ice, s-snow, r-rain)

        double precision, dimension(:,:),   allocatable :: fis
        double precision, dimension(:,:),   allocatable :: lnps
        double precision, dimension(:,:,:), allocatable :: T
        double precision, dimension(:,:,:), allocatable :: Q
        double precision, dimension(:,:,:), allocatable :: clwc, ciwc, cswc, crwc

    end type weatherBackgroundDataType



    ! "weatherInterpolatedDataType" stores projection of the NWP model fields onto set of 2-D planes
    type, public :: weatherInterpolatedDataType

        ! dimension size are => number of azimuths X number of vertical levels+1 X number of path columns
        !                    or number of azimuths X number of path columns
        !
        ! Interpolated fields array are constructed as a projection onto the 2-D path.
        ! Projections are stacked according to number of individual "Path2DType" instances.

        integer :: nPaths, nLevels, nColumns

        ! z: geometric heights of the model levels
        ! P: pressure at model levels
        ! T: interpolated atmospheric temperature
        ! Q: interpolated specific humidity
        ! c*wc: interpolated non-gaseous water content (l-liquid, i-ice, s-snow, r-rain)
        ! localRadius: Earth's curvature radius at path columns        <------------- OSBOLETE, NEED TO REWORK

        double precision, dimension(:,:,:), allocatable :: z
        double precision, dimension(:,:,:), allocatable :: P
        double precision, dimension(:,:,:), allocatable :: T
        double precision, dimension(:,:,:), allocatable :: Q
        double precision, dimension(:,:,:), allocatable :: r
        double precision, dimension(:,:,:), allocatable :: clwc, ciwc, cswc, crwc
        double precision, dimension(:,:),   allocatable :: localRadius

    end type weatherInterpolatedDataType



    type, public :: refractivityDataType

        integer :: nPaths, nLevels, nColumns
        double precision, dimension(:,:,:), allocatable :: values

    end type refractivityDataType

    

    ! According to assumption of this LTT model, ray propagates in a two-dimensional plane
    ! "Path2DType" (a.k.a. 2-D path) represents coordinates along the projection of the vector connecting a satellite with the receiver,
    ! defining the two-dimensional plane
    type, public :: Path2DType

        ! latsList: columns' latitude
        ! lonsList: columns' longitude
        ! station: observing location data, origin of the path

        ! nColumns: number of columns (default = 120)
        ! First column is in opposite direction from ray propagation, second is at station location
        ! Third and so on are placed in direction towards ray propagation

        double precision, dimension(:), allocatable :: latsList
        double precision, dimension(:), allocatable :: lonsList
        integer :: nColumns = 60
        type(stationType) :: station

    end type Path2DType



    ! "SkyPointType" is a point (direction) in the sky
    type, public :: SkyPointType

        double precision :: A
        double precision :: Z

    end type SkyPointType



    ! This LTT model operates with skyviews, i.e. set of sky directions.
    ! The delays are then computed for signals propagating from a virtual GNSS satellite 
    ! placed at certain direction.
    ! "SkyDirectionsDataType" represents a set of sky directions from single observing station.
    type, public :: SkyDirectionsDataType

        ! set: array of sky directions
        ! uniqueAzimuths: array of all possible azimuth angles in the set (used for creating an array of 2-D paths)
        ! rsat: distance of satellite from centre of the Earth (in m)

        double precision :: rsat

        type(SkyPointType), dimension(:), allocatable :: set
        double precision, dimension(:), allocatable :: uniqueAzimuths
        double precision, dimension(:), allocatable :: uniqueZenAngles

        integer :: nPoints
        integer :: nAzimuths

    end type SkyDirectionsDataType



    type, public :: SkyDelaysDataType

        ! slant: signal slant delays, array equal in size to SkyDirectionsDataType.set
        ! ZTD: zenith total delay

        double precision, dimension(:), allocatable :: slant
        double precision :: ZTD

    end type SkyDelaysDataType


contains


end module module_data_types