module module_nwpProjection

contains


    !==========================================================================
    ! Subroutine construct_path2D calculates the coordinates along the projection
    ! of the vector connecting a satellite with the receiver. The vector is
    ! projected on the surface of the ellipsoid representing the Earth - as spherical geodesic. 
    ! Input consists of the receiver coordinates and the satellite azimuth angle; the
    ! output is a set of coordinates separated by the distance corresponding to
    ! the grid spacing of the NWP model.
    !
    ! Single path is constructed while executing subroutine
    subroutine construct_path2D(station, domain, azimuth, path)
        use module_data_types, only: Path2DType, stationType, nwpDomainType
        use module_utility

        implicit none
        type(stationType), intent(in) :: station
        type(nwpDomainType), intent(in) :: domain
        double precision, intent(in) :: azimuth
        type(Path2DType), intent(in out) :: path

        integer i
        double precision :: clat0, slat0, sa, ca, ccolat

        allocate(path%latsList(path%nColumns), path%lonsList(path%nColumns))
        path%station = station

        path%latsList(2) = station%lat
        path%lonsList(2) = station%lon

        ! compute ray propagation path, based on spherical geodesic
        ca = cosdeg(azimuth)
        sa = sindeg(azimuth)
        slat0 = sindeg(station%lat)
        clat0 = cosdeg(station%lat)
        do i=1,path%nColumns-2

            ccolat = cosdeg(i*domain%deltaLon)*slat0 + sindeg(i*domain%deltaLon)*clat0*ca
            path%latsList(i+2) = 90.0 - acosdeg(ccolat)

            ! not secure; need to rework sometimes
            if (azimuth==180.0 .or. azimuth==360.0) then 
                path%lonsList(i+2) = station%lon
                cycle
            endif

            if (azimuth < 180.0) then ! towards east
                path%lonsList(i+2) = station%lon + acosdeg((cosdeg(i*domain%deltaLon)- &
                    ccolat*slat0)/clat0/sqrt(1-ccolat**2))
            else ! towards west
                path%lonsList(i+2) = station%lon - acosdeg((cosdeg(i*domain%deltaLon)- &
                    ccolat*slat0)/clat0/sqrt(1-ccolat**2))
            endif

        enddo

        ! compute a backward column
        path%latsList(1) = 90.0 - acosdeg(cosdeg(domain%deltaLon)*slat0 - &
                            sindeg(domain%deltaLon)*clat0*ca)
        path%lonsList(1) = station%lon + asindeg(-sindeg(domain%deltaLon)*sa / &
                            cosdeg(path%latsList(1)))

    end subroutine construct_path2D


    !==========================================================================
    ! Subroutine nwpInterpolation2D projects the model fields from the model
    ! grid to the two-dimensional plane defined by projection path
    !
    ! Model fields projection is done of all possible paths,
    ! that are emerging from set of unique azimuths in skyview
    subroutine nwpInterpolation2D(domain, skyview, paths, bgr_fields, undulations, int_fields)
        use module_data_types
        use module_utility
        use module_undulations

        implicit none
        type(nwpDomainType), intent(in) :: domain
        type(SkyDirectionsDataType), intent(in) :: skyview
        type(Path2DType), dimension(skyview%nAzimuths), intent(in) :: paths
        type(weatherBackgroundDataType), intent(in) :: bgr_fields
        type(undulationsType), intent(in) :: undulations
        type(weatherInterpolatedDataType), intent(in out) :: int_fields

        integer i,c,ilev,ilat, jx,jy, iiter, jiter
        double precision locrad, gridx, gridy, gridxpre, gridypre, da,db,dc,dd
        double precision columnLat, columnLon
        double precision lnps_ip, fis_ip, zs_ip, p_upper, p_lower, undul_ip
        double precision :: phl, pf, phu, dlnp, zhl, zf, zhu, dz, tm, qm, tvm, tlr

        double precision mult
        double precision, PARAMETER :: Rd=287.04 ! R_dry = R/M_dry [J/kg/K]
        double precision, PARAMETER :: eps=0.622 ! M_h2o/M_dry

        double precision, parameter :: gn = 9.80665

        ! this flag tells should model level radius be computed using it's local radius
        ! If true, local radius is computed for each column based on Earth's oblateness
        ! If false, local radius is computed only for station's location..
        ! ..and all columns are equi-distant from Earth's centre (as if Earth is a shpere)
        logical, parameter :: variableRadius = .false.

        int_fields%nPaths = skyview%nAzimuths
        int_fields%nLevels = domain%nLevels+1
        int_fields%nColumns = paths(1)%nColumns

        allocate(int_fields%z(int_fields%nPaths, int_fields%nLevels, int_fields%nColumns))
        allocate(int_fields%P(int_fields%nPaths, int_fields%nLevels, int_fields%nColumns))
        allocate(int_fields%T(int_fields%nPaths, int_fields%nLevels, int_fields%nColumns))
        allocate(int_fields%Q(int_fields%nPaths, int_fields%nLevels, int_fields%nColumns))
        allocate(int_fields%r(int_fields%nPaths, int_fields%nLevels, int_fields%nColumns))
        allocate(int_fields%clwc(int_fields%nPaths, int_fields%nLevels, int_fields%nColumns))
        allocate(int_fields%localRadius(int_fields%nPaths, int_fields%nColumns))

        ! iterate paths
        pathsloop: do i=1,int_fields%nPaths

            ! compute local radius only at the station
            ! call locradius(paths(i)%latsList(2), locrad)
            if (.not. variableRadius) locrad=radcurv(paths(i)%latsList(2), skyview%uniqueAzimuths(i))

            ! iterate path elements (columns)
            columnsloop: do c=1,int_fields%nColumns

                columnLat = paths(i)%latsList(c)
                columnLon = paths(i)%lonsList(c)

                ! print *,i,c,columnLat

                ! find distance from surface to centre of the Earth from surface of ellipsoid
                ! at this column
                ! call locradius(columnLat, locrad)
                if (variableRadius) locrad=radcurv(columnLat, skyview%uniqueAzimuths(i))

                ! find the grid point indeces corresponding to the geographical coordinates
                if (columnLon < 0.0) then
                    gridx = (columnLon + 360.0)/domain%deltaLon
                else
                    gridx = columnLon/domain%deltaLon
                endif

                ilat = 1
                do
                    if (domain%latsList(ilat) < columnLat) exit
                    if (ilat==domain%gridNLat) exit  ! watch out for grid ending before south pole
                    ilat = ilat+1
                enddo
                if (ilat==1) ilat = 2    ! watch out for grid ending before north pole
                gridy = float(ilat-1) + (domain%latsList(ilat-1) - columnLat)/ &
                                    (domain%latsList(ilat-1) - domain%latsList(ilat))
                
                ! ingore columns beyond the grid domain
                if (gridy<0.0 .or. gridy>float(domain%gridNLat) .or. gridx<0.0 .or. gridx>float(domain%gridNLon)) then
                    print *,'The projection goes beyond the grid domain!', &
                        '   azimuth=', skyview%uniqueAzimuths(i), '   columnNumber=', c, 'gridx,gridy', gridx,gridy
                    if (c==1) then
                        cycle pathsloop
                    else
                        gridx = gridxpre
                        gridy = gridypre
                    endif
                else
                    gridxpre = gridx
                    gridypre = gridy
                endif

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

                ! interpolate the surface variables
                ! Surface pressure and its logarithm
                lnps_ip = &
                    dd*dc*bgr_fields%lnps(jx,jy) + &
                    dd*da*bgr_fields%lnps(jx+1,jy) + &
                    db*dc*bgr_fields%lnps(jx,jy+1) + &
                    db*da*bgr_fields%lnps(jx+1,jy+1)
                int_fields%P(i,domain%nLevels+1,c) = exp(lnps_ip)

                ! Surface geopotential
                fis_ip = &
                    dd*dc*bgr_fields%fis(jx,jy) + &
                    dd*da*bgr_fields%fis(jx+1,jy) + &
                    db*dc*bgr_fields%fis(jx,jy+1) + &
                    db*da*bgr_fields%fis(jx+1,jy+1)

                ! Surface geometric height (above the reference ellipsoid)
                zs_ip = 0.0
                do jiter=1,2
                    zs_ip = fis_ip/gravaccel(columnLat, zs_ip)
                enddo

                ! Local radius
                int_fields%localRadius(i,c) = locrad

                ! Initialize geometrical heights of the lowest model level and surface
                int_fields%z(i,domain%nLevels+1,c) = zs_ip
                int_fields%z(i,domain%nLevels,c) = zs_ip + 10.0

                !  interpolate variables at model levels
                do ilev=1,domain%nLevels

                    p_upper = domain%afull(ilev) + domain%bfull(ilev)*exp(lnps_ip)
                    p_lower = domain%afull(ilev+1) + domain%bfull(ilev+1)*exp(lnps_ip)
                    int_fields%P(i,ilev,c) = (p_upper+p_lower)/2.0

                    int_fields%T(i,ilev,c) = &
                        dd*dc*bgr_fields%T(jx,jy,ilev) + &
                        dd*da*bgr_fields%T(jx+1,jy,ilev) + &
                        db*dc*bgr_fields%T(jx,jy+1,ilev) + &
                        db*da*bgr_fields%T(jx+1,jy+1,ilev)

                    int_fields%clwc(i,ilev,c) = &
                        dd*dc*bgr_fields%clwc(jx,jy,ilev) + &
                        dd*da*bgr_fields%clwc(jx+1,jy,ilev) + &
                        db*dc*bgr_fields%clwc(jx,jy+1,ilev) + &
                        db*da*bgr_fields%clwc(jx+1,jy+1,ilev)

                    int_fields%Q(i,ilev,c) = &
                        dd*dc*bgr_fields%Q(jx,jy,ilev) + &
                        dd*da*bgr_fields%Q(jx+1,jy,ilev) + &
                        db*dc*bgr_fields%Q(jx,jy+1,ilev) + &
                        db*da*bgr_fields%Q(jx+1,jy+1,ilev)
                enddo

                ! convert water content units from [kg/kg] to [g/m3]
                ! See: https://nwp-saf.eumetsat.int/site/download/documentation/rtm/docs_rttov12/rttov_gas_cloud_aerosol_units.pdf
                do ilev=1,domain%nLevels
                    mult = 1000*int_fields%P(i,ilev,c)/int_fields%T(i,ilev,c)/ &
                            (Rd * (1+(1-eps)/eps*int_fields%Q(i,ilev,c)))

                    int_fields%clwc(i,ilev,c) = int_fields%clwc(i,ilev,c) * mult
                enddo

                ! Surface water content = 10m water content
                int_fields%clwc(i,domain%nLevels+1,c) = int_fields%clwc(i,domain%nLevels,c)

                ! Surface specific humidity = 10m specific humidity
                int_fields%Q(i,domain%nLevels+1,c) = int_fields%Q(i,domain%nLevels,c)

                ! Iterate surface temperature and geometric height of the lowest model level
                tlr=0.0065
                do iiter=1,2

                    int_fields%T(i,domain%nLevels+1,c) = &
                        dd*dc*bgr_fields%T(jx,jy,domain%nLevels) + &
                        dd*da*bgr_fields%T(jx+1,jy,domain%nLevels) + &
                        db*dc*bgr_fields%T(jx,jy+1,domain%nLevels) + &
                        db*da*bgr_fields%T(jx+1,jy+1,domain%nLevels) + &
                        tlr * (int_fields%z(i,domain%nLevels,c)-int_fields%z(i,domain%nLevels+1,c))

                    zhl = zs_ip
                    phl = EXP(lnps_ip)
    
                    pf  = int_fields%P(i,domain%nLevels,c)
                    phu = 2*pf-phl
                    tm  = int_fields%T(i,domain%nLevels,c)
                    qm  = int_fields%Q(i,domain%nLevels,c)
                    tvm = (1.0-qm+qm/eps)*tm
    
                    zf  = zhl
                    zhu = zhl
                    do jiter=1,2
                        dlnp = LOG(pf/phl)
                        dz   = -(Rd*tvm/gravaccel(columnLat,0.5*zhl+0.5*zf))*dlnp
                        zf   = zhl+dz
                        dlnp = LOG(phu/phl)
                        dz   = -(Rd*tvm/gravaccel(columnLat,0.5*zhl+0.5*zhu))*dlnp
                        zhu  = zhl+dz
                    enddo
                enddo

                ! Calculation of the geometric model level heights (taken from the slfwd_fmi_v2.f90)
                do jiter=1,3
                    do ilev=domain%nLevels-1,1,-1
                      ! phl=domain%afull(ilev+1)*1.0 + domain%bfull(ilev+1)*int_fields%P(i,domain%nLevels+1,c)
                      ! phu=domain%afull(ilev  )*1.0 + domain%bfull(ilev  )*int_fields%P(i,domain%nLevels+1,c)
                      phl=int_fields%P(i,ilev+1,c)
                      phu=int_fields%P(i,ilev,c)
                      dlnp=LOG(phu/phl)
                      tm=0.5*int_fields%T(i,ilev+1,c)+0.5*int_fields%T(i,ilev,c)
                      qm=0.5*int_fields%Q(i,ilev+1,c)+0.5*int_fields%Q(i,ilev,c)
                      tvm=(1.0-qm+qm/0.622)*tm
      
                      !print*,i
                      if (jiter==1) then
                        dz=-(Rd*tvm/gravaccel(columnLat,0.5*int_fields%z(i,ilev+2,c)+0.5*int_fields%z(i,ilev+1,c)))*dlnp
                      else
                        dz=-(Rd*tvm/gravaccel(columnLat,0.5*int_fields%z(i,ilev+1,c)+0.5*int_fields%z(i,ilev,c)))*dlnp
                      endif
      
                      ! dz=-(Rd*tvm/gravaccel(columnLat,0.5*int_fields%z(i,ilev+1,c)+0.5*int_fields%z(i,ilev,c)))*dlnp
                      
                      int_fields%z(i,ilev,c) = int_fields%z(i,ilev+1,c)+dz
                      !print*,'>>>>',c,jiter,ilev,int_fields%z(i,ilev+1,c),int_fields%z(i,ilev,c),dz
                    end do
                end do

                ! Radius at the model levels
                !print*,'>>>>>>>>>>>locrad after', int_fields%localRadius
                do ilev=1,int_fields%nLevels
                    int_fields%r(i,ilev,c) = &
                        int_fields%localRadius(i,c) + int_fields%z(i,ilev,c)
                    !print*,ilev,rad_pr(i,ilev,c),int_fields%localRadius(i,c),int_fields%z(i,ilev,c)
                end do

            enddo columnsloop

        enddo pathsloop

    end subroutine nwpInterpolation2D


    !==========================================================================
    ! Subroutine refractivity2D calculates the values of total recractivity N
    ! in the projected model level - signal path -intersection points.
    ! 
    ! Computation is done of all possible paths,
    ! that are emerging from set of unique azimuths in skyview
    ! int_fields and N have the same dimension
    subroutine refractivity2D(int_fields, N)
        use module_data_types, only: weatherInterpolatedDataType, refractivityDataType
        use module_utility

        implicit none
        type(weatherInterpolatedDataType), intent(in) :: int_fields
        type(refractivityDataType), intent(in out) :: N

        double precision, PARAMETER :: k1=77.607E-2, k2=70.4E-2, k3=3.739E3, k4=1.45, eps=0.622
        integer :: p,l,c
        logical, parameter :: includeLWC = .false.

        ! interpolated fields and N have the same dimension
        N%nPaths = int_fields%nPaths
        N%nLevels = int_fields%nLevels
        N%nColumns = int_fields%nColumns
        allocate(N%values(N%nPaths, N%nLevels, N%nColumns))

        if (includeLWC) then

            do p=1,N%nPaths
                do l=1,N%nLevels
                    do c=1,N%nColumns
                        N%values(p,l,c) = &
                            ! I term (hydrostatic reftactivity)
                            k1*int_fields%P(p,l,c)/int_fields%T(p,l,c) + &
                            ! II term (refractivity of water vapor)
                            (k2-k1)*int_fields%P(p,l,c)*int_fields%Q(p,l,c) &
                            /((eps+0.378*int_fields%Q(p,l,c))*int_fields%T(p,l,c)) + &
                            ! III term (refractivity due to dipole moment of water vapor)
                            k3*int_fields%P(p,l,c)*int_fields%Q(p,l,c) &
                            /((eps+0.378*int_fields%Q(p,l,c))*int_fields%T(p,l,c)**2) + &
                            ! IV term (refractivity due to permittivity of liquid water)
                            ! non-gaseous water content impact: included ONLY liquid
                            ! for others (ice, snow and rain) - need to apply more advanced scattering models instead
                            k4*int_fields%clwc(p,l,c)
                    enddo
                enddo
            enddo

        else

            do p=1,N%nPaths
                do l=1,N%nLevels
                    do c=1,N%nColumns
                        N%values(p,l,c) = &
                            ! I term (hydrostatic reftactivity)
                            k1*int_fields%P(p,l,c)/int_fields%T(p,l,c) + &
                            ! II term (refractivity of water vapor)
                            (k2-k1)*int_fields%P(p,l,c)*int_fields%Q(p,l,c) &
                            /((eps+0.378*int_fields%Q(p,l,c))*int_fields%T(p,l,c)) + &
                            ! III term (refractivity due to dipole moment of water vapor)
                            k3*int_fields%P(p,l,c)*int_fields%Q(p,l,c) &
                            /((eps+0.378*int_fields%Q(p,l,c))*int_fields%T(p,l,c)**2)
                    enddo
                enddo
            enddo

        endif

    end subroutine refractivity2D


end module module_nwpProjection
