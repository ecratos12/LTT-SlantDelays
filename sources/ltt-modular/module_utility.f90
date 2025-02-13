module module_utility
    
contains

    double precision function degtor(x)
        implicit none
        double precision :: x
        double precision, parameter :: pi=3.14159265359
        degtor=x*pi/180.
    end function degtor

    double precision function rtodeg(x)
        implicit none
        double precision :: x
        double precision, parameter :: pi=3.14159265359
        rtodeg=x*180./pi
    end function rtodeg

    double precision function sindeg(x)
        implicit none
        double precision :: x
        sindeg=sin(degtor(x))
    end function sindeg

    double precision function cosdeg(x)
        implicit none
        double precision :: x
        cosdeg=cos(degtor(x))
    end function cosdeg

    double precision function tandeg(x)
        implicit none
        double precision :: x
        tandeg=tan(degtor(x))
    end function tandeg

    double precision function asindeg(x)
        implicit none
        double precision :: x
        asindeg=rtodeg(asin(x))
    end function asindeg

    double precision function acosdeg(x)
        implicit none
        double precision :: x
        acosdeg=rtodeg(acos(x))
    end function acosdeg

    double precision function atandeg(x)
        implicit none
        double precision :: x
        atandeg=rtodeg(atan(x))
    end function atandeg
    

    !==========================================================================
    ! Subroutine RADICURV calculates the local radius of Earth curvature as a
    ! function of latitude.
    !
    double precision function radcurv(lat_deg, azi_deg)
        implicit none
        double precision :: lat_deg, azi_deg

        double precision, PARAMETER :: R_q=6378138.0
        double precision, PARAMETER :: R_p=6356752.0
        double precision :: R_ns,R_ew

        R_ns=(R_q*R_p)**2/((R_q*cosdeg(lat_deg))**2+(R_p*sindeg(lat_deg))**2)**1.5
        R_ew=R_q**2/SQRT((R_q*cosdeg(lat_deg))**2+(R_p*sindeg(lat_deg))**2)

        radcurv=1.0/(cosdeg(azi_deg)**2/R_ns + sindeg(azi_deg)**2/R_ew)
    end function radcurv


    !==========================================================================
    ! Subroutine LOCRADIUS calculates the local radius of the Earth as a
    ! function of latitude.
    !
    subroutine locradius(lat,radius)
        implicit none
        double precision :: lat, radius, pi, sinlat, coslat
        double precision, PARAMETER :: a=6378137.0, f=1.0/298.257223563, b=a*(1.0-f)
        pi=4.0*ATAN(1.0)
        sinlat=SIN(lat*pi/180.0)
        coslat=COS(lat*pi/180.0)
        radius = a*b/SQRT(b**2*coslat**2 + a**2*sinlat**2)
    end subroutine locradius


    double precision function gravaccel(glat, gz)
        implicit none
        double precision :: glat, gz
        double precision, PARAMETER :: c1=5.2885e-3,c2=-5.9e-6,c3=-3.086e-6, &
            pifac=1.745329e-2,g00=9.780356
        gravaccel = &
            g00*(1.0e0+c1*(sin(glat*pifac))**2 &
            + c2*(sin(2.0e0*glat*pifac))**2)+c3*gz
    end function gravaccel


    double precision function gravaccelMSL(helmertH)
        implicit none
        double precision :: helmertH
        double precision, parameter :: g = 9.80263

        gravaccelMSL = g

    end function gravaccelMSL


    !==========================================================================
    ! Function microwaveN calculates the refractivity of air
    ! in the microwave band. Uses formula for phase refractive index 
    ! following Eresmaa_Järvinen2006
    ! Optional terms like from water content are not included
    ! 
    ! Computation is done unvectorized (single value is output)
    double precision function microwaveN(T, p, q)

        implicit none
        double precision :: T          ! in K
        double precision :: p          ! in Pa
        double precision :: q          ! unitless (in kg/kg)

        double precision, parameter :: &
        k1=77.607E-2, &
        k2=70.4E-2, &
        k3=3.739E3, &
        eps=0.622

        microwaveN = k1*p/T + (k2-k1)*p*q/((eps+0.378*q)*T) + k3*p*q/((eps+0.378*q)*T**2)

    end function microwaveN
        

    ! Converts geopotential (=dynamic) height "gpH" [m] to orthometric height "gpH2ortH" [m]
    ! using the formula by meteorologist (Kraus, 2004), 
    ! whose formula inherits an approximation approach
    double precision function gpH2ortH(gpH, lat_rad)
        implicit none
        double precision :: gpH, lat_rad

        gpH2ortH = 1/(2*1.57d-7) - &
                sqrt( (1/(2*1.57d-7))**2 - gpH/(1-0.0026373d0*cos(2*lat_rad) + &
                0.0000059d0*(cos(2*lat_rad))**2)/1.57d-7 )
    
    end function gpH2ortH
    
end module module_utility