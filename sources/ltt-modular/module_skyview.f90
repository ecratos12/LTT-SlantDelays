module module_skyview

contains

    subroutine constructSkyview(nZen, nAzi, rsat, skyview)
        use module_data_types, only: SkyPointType, SkyDirectionsDataType

        implicit none
        integer, intent(in) :: nZen, nAzi
        double precision, intent(in) :: rsat
        type(SkyDirectionsDataType), intent(out) :: skyview

        integer :: i,ia,iz
        double precision :: dZen = 1.0, dAzi = 1.0

        ! distance of satellite from centre of the Earth (in m)
        ! approximate, and assuming curcular orbit
        skyview%rsat = rsat

        allocate(skyview%uniqueZenAngles(nZen))
        do i=1,nZen
            skyview%uniqueZenAngles(i) = real(i)*dZen
        enddo

        ! example for variable zenith angle resolution

        ! integer :: nZen1 = 20, nZen2 = 10, nZen3 = 75
        ! double precision :: dZen1 = 3.0, dZen2 = 1.0, dZen3 = 0.2
        ! nZen = nZen1+nZen2+nZen3

        ! allocate(skyview%uniqueZenAngles(nZen))

        ! do i=1,nZen1
        !     skyview%uniqueZenAngles(i) = i*dZen1
        ! enddo
        ! do i=1,nZen2
        !     skyview%uniqueZenAngles(i+nZen1) = skyview%uniqueZenAngles(nZen1) + i*dZen2
        ! enddo
        ! do i=1,nZen3
        !     skyview%uniqueZenAngles(i+nZen1+nZen2) = skyview%uniqueZenAngles(nZen1+nZen2) + i*dZen3
        ! enddo

        skyview%nPoints = nZen*nAzi + 1
        skyview%nAzimuths = nAzi

        allocate(skyview%uniqueAzimuths(nAzi))
        do i=1,nAzi
            skyview%uniqueAzimuths(i) = real(i)*dAzi
        enddo

        allocate(skyview%set(skyview%nPoints))
        do iz=1,nZen
            do ia=1,nAzi
                i = iz*nAzi-nAzi+ia
                skyview%set(i)%A = skyview%uniqueAzimuths(ia)
                skyview%set(i)%Z = skyview%uniqueZenAngles(iz)
            enddo
        enddo

        skyview%set(skyview%nPoints)%Z = 0.0    ! zenith direction is the last element of array

    end subroutine constructSkyview

end module module_skyview