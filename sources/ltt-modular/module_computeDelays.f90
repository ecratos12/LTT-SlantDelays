module module_computeDelays


    ! "lttDomainType" represents the domain for performing ray tracing ..
    ! .. - input refractivity and radius values.
    ! The propagation equation is given in 2-d polar coordinate system, 
    ! originated from centre of the Earth (r=0), and observing station location (r=r_st, \psi=0)
    ! Refractivity field N, physical grounding of propagation, is given as 2-d polar grid
    ! located at corresponding r's and \psi's
    ! Upper boundary limit of ray tracing is given as r=rTOA.
    type, public :: lttDomainType

        ! N: refractivity field gridded in 2-d polar coordinates
        ! r: corresponding radius values
        double precision, dimension(:,:), allocatable :: r
        double precision, dimension(:,:), allocatable :: N

        ! dPsi: horizontal (angular) spacing of grid
        ! The refractivity and radius grid has nVert X nHoriz size

        ! Horizontal indexing is according to increase of \psi angle
        ! 1st horizontal index is \psi = -dPsi (backward column)
        ! 2nd horizontal index is at station (\psi=0)
        ! Vertical indexing is from surface to top of atmosphere; 
        ! For example, r(1,2) is surface local radius at station's location
        double precision :: dPsi
        integer :: nHoriz, nVert

        ! zTOA: height of top of atmosphere (where NWP model levels end)
        ! pTOA: pressure at top of atmosphere
        ! latTOA: latitude at top of atmosphere (in radians)
        double precision :: pTOA, zTOA, latTOA

    end type lttDomainType

    
contains


    subroutine computeAndWriteSkyDelays(parameters, station, skyview, domain, fields, undulations, forecastTime)
        use module_data_types
        use module_utility
        use module_io
        use module_nwpProjection
        use module_undulations

        implicit none
        type(parametersType), intent(in) :: parameters
        type(stationType), intent(in) :: station
        type(SkyDirectionsDataType), intent(in) :: skyview
        type(nwpDomainType), intent(in) :: domain
        type(weatherBackgroundDataType), intent(in) :: fields
        type(undulationsType), intent(in) :: undulations
        type(dateTimeType), intent(in) :: forecastTime

        integer :: iaz, i,j, nancount
        type(refractivityDataType) :: N
        type(SkyDelaysDataType) :: delays
        type(weatherInterpolatedDataType) :: interpolatedFields
        type(lttDomainType), dimension(:), allocatable :: lttDomains
        type(Path2DType), dimension(:), allocatable :: pathProjections
        double precision :: stModelHeight, stR, zenAngle, delay, psiSat, azZen
        integer :: iazN, iazE

        double precision :: t0,t1,t2,t3
        integer, dimension(8) :: t
        logical :: print_flag, include_clwc


        call date_and_time(values=t)
        t0 = t(5)*3600.0+t(6)*60.0+t(7)*1.0+t(8)*0.001
        ! Allocate intermediate arrays:
        ! -> two-dimensional planes of signal propagation
        ! -> ray tracer domains given as 2-d polar grid, corresponding to planes
        ! -> slant delays, corresponding to skyview total number of points
        allocate(pathProjections(skyview%nAzimuths))
        allocate(lttDomains(skyview%nAzimuths))
        allocate(delays%slant(skyview%nPoints-1))
        

        call date_and_time(values=t)
        t1 = t(5)*3600.0+t(6)*60.0+t(7)*1.0+t(8)*0.001
        ! Create all possible path projections (a.k.a. 2-D planes)
        ! originated from observing station location in direction of certain unique azimuth
        do iaz=1,skyview%nAzimuths
            call construct_path2D(station, domain, skyview%uniqueAzimuths(iaz), pathProjections(iaz))
        enddo


        ! Project NWP fields value onto 2-D planes
        call nwpInterpolation2D(domain, skyview, pathProjections, fields, undulations, interpolatedFields)


        ! Use projected NWP fields to create refractivity N field
        include_clwc = parameters%include_clwc
        call refractivity2D(interpolatedFields, include_clwc, N)


        ! Define ray tracer domain of each projection plane:
        do iaz=1,skyview%nAzimuths

            ! define grid size: angular separation, number of vertical levels and horizontal columns
            lttDomains(iaz)%dPsi = degtor(domain%deltaLon)
            lttDomains(iaz)%nHoriz = pathProjections(1)%nColumns
            lttDomains(iaz)%nVert = N%nLevels

            ! fill the grid with corresponding radius and refractivity values
            ! counterpart angle \psi values are accessed using angular separation, so no filling of those
            allocate(lttDomains(iaz)%r(lttDomains(iaz)%nVert, lttDomains(iaz)%nHoriz))
            allocate(lttDomains(iaz)%N(lttDomains(iaz)%nVert, lttDomains(iaz)%nHoriz))
            do i=1,lttDomains(iaz)%nVert
                do j=1,lttDomains(iaz)%nHoriz
                    lttDomains(iaz)%r(i,j) = interpolatedFields%r(iaz,lttDomains(iaz)%nVert-i+1,j)
                    lttDomains(iaz)%N(i,j) = N%values(iaz,lttDomains(iaz)%nVert-i+1,j)
                enddo
                !! TESTING: ignore backward column effect by copying station column data into backward column
                !lttDomains(iaz)%r(i,1) = lttDomains(iaz)%r(i,2)
                !lttDomains(iaz)%N(i,1) = lttDomains(iaz)%N(i,2)
            enddo

            ! set upper bourdary condition, where modeled atmosphere ends
            lttDomains(iaz)%pTOA = 1.0 ! Pa
            lttDomains(iaz)%zTOA = interpolatedFields%z(iaz,1,2)
            lttDomains(iaz)%latTOA = degtor(station%lat)

        enddo


        ! Perform LTT to compute slant delay. Loop over every non-zenith point in the skyview.

        ! compute station's model height
        call date_and_time(values=t)
        t2 = t(5)*3600.0+t(6)*60.0+t(7)*1.0+t(8)*0.001
        if (parameters%use_MSL_heights) then
            stModelHeight = station%MSL_height - interpolatedFields%z(1,interpolatedFields%nLevels,2)
        else
            stModelHeight = station%elp_height - interpolatedFields%z(1,interpolatedFields%nLevels,2)
        endif

        print_flag = .false.

        ! Compute slant delays based on skyview
        do i=1,skyview%nPoints-1

            ! get path index corresponding to current direction
            iaz = mod(i, skyview%nAzimuths)
            if (iaz==0) iaz=skyview%nAzimuths

            ! get start and finish (satellite) r-psi coordinates
            ! compute target's psi angle
            ! target's radius is constant - skyview%rsat
            if (parameters%use_MSL_heights) then
                stR = interpolatedFields%localRadius(iaz,2) + station%MSL_height
            else
                stR = interpolatedFields%localRadius(iaz,2) + station%elp_height
            endif
            zenAngle = degtor(skyview%set(i)%Z)
            psiSat = zenAngle - asin((stR/skyview%rsat)*sin(zenAngle))

            ! execute ray-tracer
            call ltt_operator(stModelHeight, skyview%rsat, psiSat, zenAngle, &
                lttDomains(iaz)%nHoriz, lttDomains(iaz)%nVert, lttDomains(iaz)%r, lttDomains(iaz)%N, &
                lttDomains(iaz)%dPsi, lttDomains(iaz)%pTOA, lttDomains(iaz)%zTOA, lttDomains(iaz)%latTOA, &
                print_flag, delay)
            delays%slant(i) = delay

        enddo

        ! North = 360 deg (iaz = nAzimuths)
        ! East = 90 deg (iaz = nAzimuths/4)
        iazN = skyview%nAzimuths
        iazE = skyview%nAzimuths/4

        ! Select plane of propagation for zenith delay case..
        ! .. in direction of horizontal component of refraction gradient:
        ! dN/dx = (Neast - Nwest)/(2*dPsi), dN/dy = (Nnorth - Nsouth)/(2*dPsi)
        ! (approximately)
        azZen = atan2(lttDomains(iazN)%n(1,3) - lttDomains(iazN)%n(1,1), &
                      lttDomains(iazE)%n(1,3) - lttDomains(iazE)%n(1,1))
        azZen = 90 - rtodeg(azZen)
        if (azZen < parameters%dAzimuth_deg) azZen = azZen+360.0
        iaz = floor(azZen / parameters%dAzimuth_deg)

        ! Compute zenith delay
        psiSat = 0.0
        zenAngle = 0.0
        !iaz = 1 ! arbitary azimuth (==dAzimuth deg)

        call ltt_operator(stModelHeight, skyview%rsat, psiSat, zenAngle, &
            lttDomains(iaz)%nHoriz, lttDomains(iaz)%nVert, lttDomains(iaz)%r, lttDomains(iaz)%N, &
            lttDomains(iaz)%dPsi, lttDomains(iaz)%pTOA, lttDomains(iaz)%zTOA, lttDomains(iaz)%latTOA, &
            .false., delay)
        delays%ZTD = delay

        call date_and_time(values=t)
        t3 = t(5)*3600.0+t(6)*60.0+t(7)*1.0+t(8)*0.001

        WRITE (*,'(A,F12.6,A)') 'LTT operator average time: ', 1000*(t3-t2)/real(skyview%nPoints), ' ms'

        ! 7.7. Write skyview in ascii file
        call writeDelays(parameters, station, delays, forecastTime, skyview)

        ! report NaNs
        nancount=0
        do i=1,size(delays%slant)
            if (delays%slant(i) /= delays%slant(i)) nancount = nancount+1
        enddo
        if (nancount/=0) print*,'number of NaNs: ', nancount


        WRITE (*,'(A,F12.6,A)') 'Memory allocation time:    ', 100*(t1-t0)/(t3-t0),' %'
        WRITE (*,'(A,F12.6,A)') 'Fields interpolation time: ', 100*(t2-t1)/(t3-t0),' %'
        WRITE (*,'(A,F12.6,A)') 'LTT operator time:         ', 100*(t3-t2)/(t3-t0),' %'

    end subroutine computeAndWriteSkyDelays


    subroutine ltt_operator(stModelHeight, rSat, psiSat, zenAngleRad, &
                            nHor, nVert, radius, refrac, dPsi, pTOA, zTOA, latTOA, debug_printing, delay)
        use module_utility

        implicit none
        double precision, intent(in) :: stModelHeight  ! in meters above the lowest model level (can be negative!!)
        double precision, intent(in) :: rSat           ! radius at satellite location
        double precision, intent(in) :: psiSat         ! \psi angle at satellite location
        double precision, intent(in) :: zenAngleRad    ! zenith angle in radians
        integer, intent(in) :: nHor, nVert
        double precision, dimension(nVert, nHor), intent(in) :: radius, refrac
        double precision, intent(in) :: dPsi
        double precision, intent(in) :: pTOA, zTOA, latTOA
        logical, intent(in) :: debug_printing
        double precision, intent(out) :: delay

        double precision, parameter :: pi=3.14159265359
        integer, parameter :: msplit=10       ! number of differential steps within one model layer
        integer, parameter :: niterations=8   ! number of iterations towards true starting zenith angle
        double precision, parameter :: zaCorRate = sqrt(0.5) ! rate of correction while iterating

        double precision :: za_start, alpha
        double precision :: psiRay, psiTOA, psiTan !, psiMin, psiMax
        
        double precision :: delta_d, &      ! difference between true and geometrical paths
                            d, ro, beta, d_za

        double precision, dimension(niterations) :: delays_tmp
        double precision :: zDelay_extra, delay_extra, hdc, zf
        double precision :: y(4), yt(4), dydh(4), dydht(4,4)

        double precision :: rad,ref, hwt1,hwt2, amult, h,h2,huse, dr_max,dr_ds,rSt,dr, kval
        integer :: ilev,j,iray,ibot,k,kp1


        ! psiMin = -dPsi
        ! psiMax = dPsi * real(nHor-2)
        psiTan = 0.0   ! psiTan is used in the radio occulation problem - set to zero here
        alpha = 0.0    ! ray bending angle
        za_start = zenAngleRad

        ! set the radius value of the starting point
        ! radius = radius of the lowest model level + station height above the lowest model level
        rSt = radius(1,2) + stModelHeight

        ibot=1
        ! which model levels is the receiver between?
        do
            if (rSt < radius(ibot+1,2)) exit
            ibot = ibot+1
        enddo
        if (debug_printing) print *, ibot, rSt, stModelHeight

        do iray=1,niterations

            amult = 1.0

            if (debug_printing) print *, iray, alpha
            
            za_start = za_start - zaCorRate*alpha

            ! current position vector
            y(1) = 0.0       ! height above receiver
            y(2) = 0.0       ! psi
            y(3) = za_start  ! theta
            y(4) = 0.0       ! delay

            ! integrate ray-path to top of model atmosphere
            do ilev = ibot,nVert-1

                ! estimate the step length
                ! k - column index behind current position y
                ! k+1 - column index in front of current position y
                ! h - horizontal step length (in m)
 
                k = int((y(2) + psiTan)/dPsi)+2
                k = min(k,nHor)
                k = max(k,1)
                dr_max = (radius(ilev+1,k) - MAX(radius(ilev,k),rSt))/REAL(msplit)
                h = dr_max/MAX(COS(y(3)),1.0E-10)

                ! limit to horizontal distance between grid points
                h = MIN(h,6.371E6*dPsi)

                h2 = 0.5*h

                ! now calculate the path-length with a fourth-order RUNGE-KUTTA
                do j=1,msplit

                    yt(:) = y(:)

                    call pderivs(nHor, nVert, radius, refrac, dPsi, ilev, psiTan, rSt, amult, yt, dydht(:,1))
                    yt(:) = y(:) + dydht(:,1)*h2

                    call pderivs(nHor, nVert, radius, refrac, dPsi, ilev, psiTan, rSt, amult, yt, dydht(:,2))
                    yt(:) = y(:) + dydht(:,2)*h2

                    call pderivs(nHor, nVert, radius, refrac, dPsi, ilev, psiTan, rSt, amult, yt, dydht(:,3))
                    yt(:) = y(:) + dydht(:,3)*h

                    call pderivs(nHor, nVert, radius, refrac, dPsi, ilev, psiTan, rSt, amult, yt, dydht(:,4))

                    dydh(:) = (dydht(:,1)+dydht(:,4)+2.0*(dydht(:,2)+dydht(:,3)))/6.0
                    yt(:) = y(:) + dydh(:)*h

                    ! check the radius - have we exited the level
                    k = int((yt(2) + psiTan)/dPsi)+2
                    k = min(k,nHor-1)
                    k = max(k,1)
                    kp1 = k+1

                    ! horizontal weighting factor
                    hwt1 = (REAL(k-1)*dPsi - (yt(2)+psiTan))/dPsi
                    hwt2 = 1.0 - hwt1

                    ! radius of pressure level interpolated between columns
                    rad = hwt1*radius(ilev+1,k)+hwt2*radius(ilev+1,kp1)

                    ! if gone over the boundary scale h
                    if (j==msplit) then
                        ! dr_dpsi = 0.0
                        ! if (yt(2) < psiMax .AND. yt(2) > psiMin) &

                        ! dr/ds = dr/dpsi * dpsi/ds
                        dr_ds = dydh(1) - (radius(ilev+1,kp1)-radius(ilev+1,k))/dPsi*dydh(2)

                        huse = h - (yt(1)+rSt-rad)/dr_ds
                        if (dr_ds < 0) &
                        print *, 'Orography change leads to ray path fall into lower atmospheric levels, dz/ds=', dr_ds

                        if (debug_printing) print *, ilev, huse, yt(1)+rSt-rad, dr_ds
                    else
                        huse = h
                    endif

                    ! update the position vector
                    y(:) = y(:) + dydh(:)*huse
                    ! 
                    if (y(1) < 0) then
                        print *, 'Collision with the ground!!! Set slant delay to 0'
                        delay = 0.0
                        return
                    endif

                    ! output the refrac etc. at the model level.
                    if (j==msplit .and. iray>1) then
                        yt(:) = y(:)
	     	    
                        k = int((yt(2) + psiTan)/dPsi)+2
                        k = min(k,nHor-1)
                        k = max(k,1)
                        kp1 = k+1
                
                        ! horizontal weighting factor
                        hwt1 = (REAL(k-1)*dPsi - (yt(2)+psiTan))/dPsi
                        hwt2 = 1.0 - hwt1

                        if (debug_printing) print *, ilev+1, k, yt
           
                        rad = hwt1*radius(ilev+1,k)+hwt2*radius(ilev+1,kp1)
                        ref = hwt1*refrac(ilev+1,k)+hwt2*refrac(ilev+1,kp1)
                    endif

                    ! try to maintain roughly the same radial increment by adjusting h
                    if (j<msplit) then
                        dr = (rad-y(1)-rSt)/REAL(msplit-j)
                        h = min(h,dr/max(cos(y(3)),1.0E-10))
                        h2 = 0.5*h
                    endif

                enddo

            enddo

            ! estimate the slant path above model top
            kval = LOG(refrac(nVert-1,k)/refrac(nVert,k))/(radius(nVert,k)-radius(nVert-1,k))

            ! delay above the model top
            hdc=2.2779E-5
            zf = 1.0 - 0.00266*COS(2.0*latTOA) - 0.00000028*zTOA
            zDelay_extra = (hdc/zf)*pTOA
            delay_extra = zDelay_extra * mapping_factor((y(1)+rSt),y(3),kval)

            ! estimate the \psi value at r=rSat, assuming a straight line above the model top
            psiRay = y(2)+y(3) - asin(sin(y(3))*(y(1)+rSt)/rSat)
            ! difference between target \psi (psiSat) and calculated is due to bending
            alpha = psiRay - psiSat

            ! estimate the length difference between geometrical and bended paths
            if (abs(zenAngleRad) > 1.0E-4) then
        
                d_za = zenAngleRad - za_start
                beta = pi/2 - zenAngleRad
                ro = rSt * cos(beta)
                psiTOA = acos(ro/(y(1)+rSt)) - beta
                d = (y(1)+rSt) * (psiTOA-y(2))
        
                delta_d = 0.5 * abs(d) * abs(d_za) * cos(y(3))
        
                !  if (zenith_angle*180.0/3.14159 > 85.0) then
                !  print*, delta_d, (psiTOA-y(2)), d_za
                !  endif

            else 
                delta_d = 0.0
            endif

            ! total slant delay on this iteration
            delays_tmp(iray) = y(4) + delay_extra + delta_d

        enddo

        ! calculate total slant delay as average value from latest iterations, 
        ! assuming total slant delay values converge and oscillate near an optimum value.
        !delay = 0.0
        !do j=0,3
        !    delay = delay + delays_tmp(niterations-j)
        !enddo
        !delay = delay / 4.0

        delay = delays_tmp(niterations)

        return

    end subroutine ltt_operator


    subroutine pderivs(nHor, nVert, r, N, dPsi, i, psiTan, rSt, amult, y, dydh)
        ! use module_computeDelays, only: lttDomainType

        implicit none
        integer, intent(in) :: nHor, nVert
        double precision, dimension(nVert, nHor), intent(in) :: r,N
        double precision, intent(in) :: dPsi
        integer, intent(in) :: i                ! below level's number
        double precision, intent(in) :: psiTan
        double precision, intent(in) :: rSt, amult
        double precision, intent(in) :: y(4)
        double precision, intent(out) :: dydh(4)

        integer :: k,kp1
        double precision :: hwt1,hwt2
        double precision :: ref_up,ref_low
        double precision :: rad_up,rad_low
        double precision :: N_frw, N_bkw
        double precision :: kval, localN, localR
        double precision :: dndr, dndpsi

        logical :: preciseGradients
        preciseGradients = .false.

        localR = y(1)+rSt

        ! column index
        k = int((y(2) + psiTan)/dPsi)+2
        k = min(k,nHor)
        k = max(k,1)
        kp1 = min(k+1,nHor)

        ! horizontal weighting factor
        hwt1 = (REAL(k-1)*dPsi - (y(2)+psiTan))/dPsi   
        hwt2 = 1.0 - hwt1

        ! calculate the local radial gradient of n
        ref_up  = hwt1*N(i+1,k)+hwt2*N(i+1,kp1)
        ref_low = hwt1*N(i,k) + hwt2*N(i,kp1)
    
        rad_up = hwt1*r(i+1,k)+hwt2*r(i+1,kp1)
        rad_low = hwt1*r(i,k) + hwt2*r(i,kp1)
    
        kval = LOG(ref_low/ref_up)/(rad_up-rad_low)
        localN = ref_low*EXP(-kval*(localR-rad_low))
    
        dndr = -1.0E-6*kval*localN

        ! calculate the local tangent gradient of n
        if (preciseGradients .eqv. .false.) then
            ! dn/dpsi; dr=const is not guaranted
            if (k/=1) then
                dndpsi = 1.0E-6*(hwt1*(N(i,kp1)-N(i,k)) + hwt2*(N(i,k)-N(i,k-1)))/dPsi
            else
                dndpsi = 1.0E-6*hwt1*(N(i,kp1) - N(i,k))/dPsi
            endif
        else
            N_frw = N(i,k+1)*exp(-log(N(i,k+1)/N(i+1,k+1))*(localR-r(i,k+1))/(r(i+1,k+1)-r(i,k+1)))
            N_bkw = N(i,k)*exp(-log(N(i,k)/N(i+1,k))*(localR-r(i,k))/(r(i+1,k)-r(i,k)))
            dndpsi = 1.0E-6*(N_frw - N_bkw)/dPsi
        endif

        ! compute derivatives for ray-tracing equations
        dydh(1) = cos(y(3))
        dydh(2) = amult*sin(y(3))/localR
        dydh(3) = -SIN(y(3))*(1.0/localR + dndr) + dndpsi*COS(y(3))/localR
        dydh(4) = 1.0E-6*localN

        return

    end subroutine pderivs


    double precision FUNCTION mapping_factor(rad,phi,kval)

        double precision, INTENT(IN) :: rad       ! radius value
        double precision, INTENT(IN) :: phi       ! angle with local radius
        double precision, INTENT(IN) :: kval      ! inverse scale-height

        integer :: i
        double precision :: alpha,amult,asym_series,aterm

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
    
end module module_computeDelays