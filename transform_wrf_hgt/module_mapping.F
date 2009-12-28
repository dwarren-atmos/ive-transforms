      MODULE MODULE_MAPPING
      
      real, parameter ::
     &                  R_earth=6370., ! radius of the earth
     &                  pid=3.1415/180., ! radians to degrees
     &                  missing=-9999.0  ! missing value

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   earth_dist - function that calculates the earth distance between 
!                two points.
!
!    lon1 - longitude of point 1
!    lat1 - latitude of point 1
!    lon2 - longitude of point 2
!    lat2 - latitude of point 2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function earth_dist(lon1, lat1, lon2, lat2)

      real, intent(in) :: lon1, lat1, lon2, lat2
      real :: earth_dist, pro

      pro =    sin(lat1*pid) * sin(lat2*pid) +
     &         cos(lat1*pid) * cos(lat2*pid) * 
     &         cos((lon2-lon1)*pid)

      if ( abs(pro) .gt. 1.001 ) then
        print*,'In earth_dist, pro > 1, returning',pro
        earth_dist = missing
      elseif ( abs(pro) .gt. 1 ) then
        earth_dist = 0.
      else
        earth_dist = R_earth * acos(pro)
      endif

      return
      end function

      END MODULE
