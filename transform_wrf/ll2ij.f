      subroutine ll2ij( x, y, lat, lon, xlat, xlon, nx, ny)
!     this routine takes a lat and lon and returns an x y in 
!     array space.

      integer nx,ny
      real x, y, lat, lon, xlat(nx,ny), xlon(nx,ny)
      integer i,j, besti, ip1, jp1
      real dist, bestdist, lastdist, savdist, sav2dist
      real dx1,dx2,dx3,dx4, dy1,dy2,dy3,dy4

      savdist=sqrt((xlat(nx,ny)-xlat(1,1))**2 +
     &     (xlon(nx,ny)-xlon(1,1))**2)
      sav2dist=savdist
      
      do while(j .le. ny)
         i=1
         lastdist=savdist
         dist = sqrt((lat-xlat(1,j))**2 +
     &     (lon-xlon(1,j))**2)
         if(sav2dist .lt. dist) then
            j=ny+1
         else
            sav2dist = dist
         endif
         do while(i .le. nx)
            dist = sqrt((lat -xlat(i,j))**2 + (lon -xlon(i,j))**2)
            if(bestdist .gt. dist)then 
               bestdist=dist
               x = i
               y = j
               i = i + 1
            else
               if(savdist .gt. dist) then
                  savdist = dist
                  i=nx+1
               else
                  i = nx+1
               endif
            endif
         end do
         j = j + 1
      end do

!     have closest i,j
      if(x .le. lat .and. y .le.lon) then
         i = x; j = y; ip1 = x + 1; jp1 = y + 1;
         if(i .eq. nx)ip1 = nx
         if(j .eq. ny)jp1 = ny
      elseif(x .ge. lat .and. y .ge.lon) then
         i = x - 1; j = y - 1; ip1 = x; jp1 = y;
         if(ip1 .eq. 1)i = 1
         if(jp1 .eq. 1)j = 1
      elseif(x .le. lat .and. y .ge.lon) then
         i = x; j = y - 1; ip1 = x + 1; jp1 = y;
         if(i .eq. nx)ip1 = nx
         if(jp1 .eq. 1)j = 1
      elseif(x .ge. lat .and. y .le.lon) then
         i = x - 1; j = y; ip1 = x; jp1 = y + 1;
         if(ip1 .eq. 1)i = 1
         if(j .eq. ny)jp1 = ny
      else
         write(6,*)'cant get here x.ge.lat.ge.x etc.'
      endif
      dx1 = (lon - xlon(i,j))/(xlon(ip1,j)-xlon(i,j))
      dx2 = 1. - dx1
      dx3 = (lon - xlon(i,jp1))/(xlon(ip1,jp1)-xlon(i,jp1))
      dx4 = 1. - dx3
      dy1 = (lat - xlat(i,j))/(xlat(i,jp1)-xlat(i,j))
      dy2 = 1. - dy1
      dy3 = (lat - xlat(ip1,j))/(xlat(ip1,jp1)-xlat(ip1,j))
      x = i*dx1*dy1 + i*dx3*dy3 + ip1*dx2*dy2 + ip1*dx4*dy4
      y = j*dx1*dy1 + j*dx3*dy3 + jp1*dx2*dy2 + jp1*dx4*dy4
      return
      end subroutine
      
      
