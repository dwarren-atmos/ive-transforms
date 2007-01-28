      subroutine grid
     1  (igrid,reflat,reflon,iref,jref,stdlt1,stdlt2
     2  ,stdlon,delx,delx2,dely,dely2,grdi,grdj,m,n,grdlat
     3  ,grdlon,gcon,f,flat,hx,hy,xpos,ypos
     4  ,distx,disty,dxav1,dxav2,dxm,dxu,dxv
     5  ,dyav1,dyav2,dym,dyu,dyv,grdrot,istr1,istr2,strgrd)

c#include <grid.prol>
c rcs keywords: $RCSfile: grid.f,v $ 
c               $Revision: 1.1.1.1 $ $Date: 2001/04/10 21:59:35 $c
c***********************************************************************
c
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer igrid
      integer iref
      integer istr1
      integer istr2
      integer jref
      integer m
      integer n
c
      real delx
      real delx2
      real dely
      real dely2
      real distx  (m,n)
      real disty  (m,n)
      real dxav1  (m,n)
      real dxav2  (m,n)
      real dxm    (m,n)
      real dxu    (m,n)
      real dxv    (m,n)
      real dyav1  (m,n)
      real dyav2  (m,n)
      real dym    (m,n)
      real dyu    (m,n)
      real dyv    (m,n)
      real f      (m,n)
      real flat
      real gcon
      real grdi   (m,n)
      real grdj   (m,n)
      real grdlat (m,n)
      real grdlon (m,n)
      real hx     (m,n)
      real hy     (m,n)
      real reflat
      real reflon
      real grdrot (m,n)
      real stdlon
      real stdlt1
      real stdlt2
      real strgrd
      real xpos   (m,n)
      real ypos   (m,n)
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
      integer ihem
      integer istr3
      integer istr4
      integer j
      integer jstr1
      integer jstr2
      integer jstr3
      integer jstr4
      integer m1
      integer mn
      integer n1
c
      real angle
      real cn1
      real cn2
      real cn3
      real cn4
      real cn5
      real cn6
      real con1
      real con2
      real d2r
      real deg
      real ogcon
      real omega4
      real onedeg
      real pi
      real pi2
      real      pi4
      real r2d
      real radius
      real rih
      real rr
      real x
      real xih
      real xx
      real y
      real yih
      real yy
 
c
c***********************************************************************
c          subroutine: grid
c
c          purpose: to give latitude, longitude, grid constant,
c                   coriolis force and map factors in x- and y-
c                   directions for 'm x n' points within a fixed
c                   grid.  all latitudes in this routine start
c                   with -90.0 at the south pole and increase
c                   northward to +90.0 at the north pole.  the
c                   longitudes start with 0.0 at the greenwich
c                   meridian and increase to the west, so that
c                   90.0 refers to 90.0w, 180.0 is the inter-
c                   national dateline and 270.0 is 90.0e.
c
c          input variables:
c
c            igrid:  type of grid projection:
c
c                      =1, mercator projection
c                      =2, lambert conformal projection
c                      =3, polar stereographic projection
c                      =4, cartesian coordinates
c                      =5, spherical projection
c
c            reflat: latitude at reference point (iref,jref)
c
c            reflon: longitude at reference point (iref,jref)
c
c            iref:   i-coordinate value of reference point
c
c            jref:   j-coordinate value of reference point
c
c            stdlt1: standard latitude of grid
c
c            stdlt2: second standard latitude of grid (only required
c                    if igrid=2, lambert conformal)
c
c            stdlon: standard longitude of grid (longitude that
c                     points to the north)
c
c            delx:   grid spacing of grid in x-direction
c                    for igrid=1,2,3,4 or 5, delx must be in meters
c
c            dely:   grid spacing (in meters) of grid in y-direction
c                    for igrid=1,2,3,4 or 5, delx must be in meters
c
c            grdi:   i-coordinate(s) that this routine will generate
c                    information for
c
c            grdj:   j-coordinate(s) that this routine will generate
c                    information for
c
c            flat:   latitude for f-plane (cartesian grid only)
c
c            m   :   number of points in x-direction
c
c            n   :   number of points in y-direction
c
c          output variables:
c
c            grdlat: latitude of point (grdi,grdj)
c
c            grdlon: longitude of point (grdi,grdj)
c
c            gcon:   grid constant for grid
c
c            f:      coriolis force at point (grdi,grdj)
c
c            hx:     map factor in x-direction at point (grdi,grdj)
c                    for igrid=1,2,3,4 or 5; grid spacing in meters
c
c            hy:     map factor in y-direction at point (grdi,grdj)
c                    for igrid=1,2,3,4 or 5; grid spacing in meters
c
c            grdrot: Angle difference between grid north and
c                    the true north at grid point.
c
c            xpos:   x-position in meters from the line extending
c                    from the pole through the standard longitue
c                    (reflat)
c
c            ypos:   y-position in meters from the pole
c
c***********************************************************************
c          make sure igrid is an acceptable value
c***********************************************************************
c
      if (igrid.lt.1.or.igrid.gt.5) then
        print 800
        print 805,igrid
        print 800
        return
      endif
c
c***********************************************************************
c          local constants
c***********************************************************************
c
      pi=4.0*atan(1.0)
      pi2=pi/2.0
      pi4=pi/4.0
      d2r=pi/180.0
      r2d=180.0/pi
      radius=6371229.0
      omega4=4.0*pi/86400.0
      onedeg=radius*2.0*pi/360.0
c
      istr3=m-istr2+1
      istr4=m-istr1+1
      mn=m*n
      m1=m-1
      n1=n-1
c
c      print*," istr1=",istr1," istr2=",istr2," iref=",iref
c     1      ," istr3=",istr3," istr4=",istr4
c
c Set grdrot to 0.0
c
      do j=1,n
         do i=1,m
            grdrot(i,j)=0.0
         enddo
      enddo
c
      do i=istr2,istr3-1
        dxu(i,1)=delx2
      enddo
      if (istr2.ge.2) then
        do i=istr2-1,istr1,-1
          dxu(i,1)=dxu(i+1,1)*strgrd
        enddo
        do i=istr3,istr4-1
          dxu(i,1)=dxu(i-1,1)*strgrd
        enddo
        if (istr1.ge.2) then
          do i=istr1-1,1,-1
            dxu(i,1)=dxu(i+1,1)
          enddo
          do i=istr4,m1
            dxu(i,1)=dxu(i-1,1)
          enddo
        endif
      endif
      do j=2,n
        do i=1,m1
          dxu(i,j)=dxu(i,1)
        enddo
      enddo
      do j=1,n
        dxu(m,j)=0.0
      enddo
c
      do i=2,m1
        dxm(i,1)=(dxu(i,1)+dxu(i-1,1))*0.5
      enddo
      dxm(1,1)=dxu( 1,1)
      dxm(m,1)=dxu(m1,1)
      do j=2,n
        do i=1,m
          dxm(i,j)=dxm(i,1)
        enddo
      enddo
c
      do j=1,n1
        do i=1,m
          dxv(i,j)=(dxm(i,j)+dxm(i,j+1))*0.5
        enddo
      enddo
      do i=1,m
        dxv(i,n)=0.0
      enddo
c
c     call qprntn(dxu,'dxu     ','        ',1,1,m,n)
c     call qprntn(dxm,'dxm     ','        ',1,1,m,n)
c     call qprntn(dxv,'dxv     ','        ',1,1,m,n)
c
      jstr1=istr1
      jstr2=istr2
      jstr3=n-jstr2+1
      jstr4=n-jstr1+1
c
      do j=jstr2,jstr3-1
        dyv(1,j)=dely2
      enddo
      if (jstr2.ge.2) then
        do j=jstr2-1,jstr1,-1
          dyv(1,j)=dyv(1,j+1)*strgrd
        enddo
        do j=jstr3,jstr4-1
          dyv(1,j)=dyv(1,j-1)*strgrd
        enddo
        if (jstr1.ge.2) then
          do j=jstr1-1,1,-1
            dyv(1,j)=dyv(1,j+1)
          enddo
          do j=jstr4,n1
            dyv(1,j)=dyv(1,j-1)
          enddo
        endif
      endif
      do j=1,n1
        do i=2,m
          dyv(i,j)=dyv(1,j)
        enddo
      enddo
      do i=1,m
        dyv(i,n)=0.0
      enddo
c
      do j=2,n1
        dym(1,j)=(dyv(1,j)+dyv(1,j-1))*0.5
      enddo
      dym(1,1)=dyv(1, 1)
      dym(1,n)=dyv(1,n1)
      do j=1,n
        do i=2,m
          dym(i,j)=dym(1,j)
        enddo
      enddo
c
      do j=1,n
        do i=1,m1
          dyu(i,j)=(dym(i,j)+dym(i+1,j))*0.5
        enddo
      enddo
      do j=1,n
        dyu(m,j)=0.0
      enddo
c
c     call qprntn(dyv,'dyv     ','        ',1,1,m,n)
c     call qprntn(dym,'dym     ','        ',1,1,m,n)
c     call qprntn(dyu,'dyu     ','        ',1,1,m,n)
c
c***********************************************************************
c          compute weights for averaging winds at mass points
c***********************************************************************
c
      do j=2,n1
        do i=2,m1
          dxav1(i,j)=dxu(i-1,j)/(dxu(i-1,j)+dxu(i,j))
          dxav2(i,j)=dxu(i  ,j)/(dxu(i-1,j)+dxu(i,j))
          dyav1(i,j)=dyv(i,j-1)/(dyv(i,j-1)+dyv(i,j))
          dyav2(i,j)=dyv(i,j  )/(dyv(i,j-1)+dyv(i,j))
        enddo
      enddo
      do i=1,m
        dxav1(i,1)=0.0
        dxav2(i,1)=0.0
        dyav1(i,n)=0.0
        dyav2(i,n)=0.0
      enddo
      do j=2,n1
        dxav1(1,j)=0.0
        dxav2(1,j)=0.0
        dyav1(m,j)=0.0
        dyav2(m,j)=0.0
      enddo
c     call qprntn(dxav1,'dxav1   ','        ',1,1,m,n)
c     call qprntn(dxav2,'dxav2   ','        ',1,1,m,n)
c     call qprntn(dyav1,'dyav1   ','        ',1,1,m,n)
c     call qprntn(dyav2,'dyav2   ','        ',1,1,m,n)
c
c***********************************************************************
c          compute distances on grid
c***********************************************************************
c
      do j=1,n
        distx(1,j)=(grdi(1,j)-iref)*delx
        do i=2,m
          distx(i,j)=distx(i-1,j)+dxu(i-1,j)
        enddo
      enddo
c
      do i=1,m
        disty(i,1)=(grdj(i,1)-jref)*dely
        do j=2,n
          disty(i,j)=disty(i,j-1)+dyv(i,j-1)
        enddo
      enddo
c     call qprntn(distx,'distx   ','        ',1,1,m,n)
c     call qprntn(disty,'disty   ','        ',1,1,m,n)
c
c***********************************************************************
c          mercator projection (igrid=1)
c***********************************************************************
c
      if (igrid.eq.1) then
        gcon=0.0
        deg=abs(stdlt1)*d2r
        con1=cos(deg)
        con2=radius*con1
        deg=reflat*0.5*d2r
        rih=con2*alog(tan(pi4+deg))
        do j=1,n
          do i=1,m
            xpos(i,j)=1.0
            ypos(i,j)=1.0
           rr=rih+(grdj(i,j)-jref)*dely
c            rr=rih+disty(i,j)-disty(i,jref)
            grdlat(i,j)=(2.0*atan(exp(rr/con2))-pi2)*r2d
           grdlon(i,j)=reflon+(grdi(i,j)-iref)*r2d*delx/con2
c            grdlon(i,j)=reflon+(distx(i,j)-distx(iref,j))*r2d/con2
            if (grdlon(i,j).gt.360.0) grdlon(i,j)=grdlon(i,j)-360.0
            if (grdlon(i,j).lt.  0.0) grdlon(i,j)=grdlon(i,j)+360.0
            deg=grdlat(i,j)*d2r
            f(i,j)=omega4*sin(deg)
            hx(i,j)=cos(deg)/con1
            hy(i,j)=hx(i,j)
          enddo
        enddo
        return
c
c***********************************************************************
c          lambert conformal (igrid=2) or
c          polar stereographic (igrid=3)
c***********************************************************************
c
      elseif (igrid.eq.2.or.igrid.eq.3) then
        if (igrid.eq.2) then
          if (stdlt1.eq.stdlt2) then
            gcon=sin(abs(stdlt1)*d2r)
          else
            gcon=(log(sin((90.0-abs(stdlt1))*d2r))
     1           -log(sin((90.0-abs(stdlt2))*d2r)))
     2          /(log(tan((90.0-abs(stdlt1))*0.5*d2r))
     3           -log(tan((90.0-abs(stdlt2))*0.5*d2r)))
          endif
        else
          gcon=1.0
        endif
        ogcon=1.0/gcon
        ihem=nint(abs(stdlt1)/stdlt1)
        deg=(90.0-abs(stdlt1))*d2r
        cn1=sin(deg)
        cn2=radius*cn1*ogcon
        deg=deg*0.5
        cn3=tan(deg)
        deg=(90.0-abs(reflat))*0.5*d2r
        cn4=tan(deg)
        rih=cn2*(cn4/cn3)**gcon
        deg=(reflon-stdlon)*d2r*gcon
        xih= rih*sin(deg)
        yih=-rih*cos(deg)*ihem
        do j=1,n
          do i=1,m
            x=xih+distx(i,j)
            y=yih+disty(i,j)
            xpos(i,j)=x*ihem
            ypos(i,j)=y*ihem
            rr=sqrt(x*x+y*y)
            grdlat(i,j)=r2d*(pi2-2.0*atan(cn3*(rr/cn2)**ogcon))*ihem
            xx= x
            yy=-y*ihem
            if (yy.eq.0.0) then
              if (xx.le.0.0) then
                angle=-90.0
              elseif (xx.gt.0.0) then
                angle=90.0
              endif
            else
              angle=atan2(xx,yy)*r2d
            endif
            grdlon(i,j)=stdlon+angle*ogcon
            deg=grdlat(i,j)*d2r
            f(i,j)=omega4*sin(deg)
            if (grdlon(i,j).gt.360.0) grdlon(i,j)=grdlon(i,j)-360.0
            if (grdlon(i,j).lt.  0.0) grdlon(i,j)=grdlon(i,j)+360.0
            deg=(90.0-grdlat(i,j)*ihem)*d2r
            cn5=sin(deg)
            deg=deg*0.5
            cn6=tan(deg)
            if (igrid.eq.2) then
              hx(i,j)=cn5/cn1*(cn6/cn3)**(-gcon)
            else
              hx(i,j)=(1.0+sin(abs(grdlat(i,j))*d2r))
     1               /(1.0+sin(abs(stdlt1)*d2r))
            endif
            hy(i,j)=hx(i,j)
          enddo
        enddo
c
c  Calculate angle betwen grid north and true north at every grid
c  point.
c 

          do 100 j=1, n
             do 100 i=1, m
                 grdrot(i,j)=(stdlon-grdlon(i,j))*gcon
c
c  account for crossing greenwich meridian
c
c  case 1: standard longitude is east of Meridian, (between 0 and 90E)
c          and the grid longitude is west of Meridian (between 90W and Meridian)
c      
                if (stdlon .ge. 0.0 .and. stdlon .le. 90.0 .and. 
     &             grdlon(i,j) .ge. 270.0 .and. grdlon(i,j) .lt. 360.0) 
     &          then
                   grdrot(i,j)=(stdlon+360.0-grdlon(i,j))*gcon
                endif
c
c  case 2: standard longitude is west of Meridian, (between 90W and Meridian)
c          and the grid longitude is east of Meridian (between 0 and 90E)
                if (stdlon .ge. 270.0 .and. stdlon .lt. 360.0 .and.
     &              grdlon(i,j) .ge. 0.0 .and. grdlon(i,j) .le. 90.0) 
     &          then
                    grdrot(i,j)=(stdlon-grdlon(i,j)-360.0)*gcon
                endif            
100	  continue
        return
c
c***********************************************************************
c          analytic grid (igrid=4)
c***********************************************************************
c
      elseif (igrid.eq.4) then
        gcon=2.0
        cn2=delx/onedeg
        do j=1,n
          do i=1,m
            xpos(i,j)=1.0
            ypos(i,j)=1.0
c            grdlat(i,j)=reflat+(disty(i,j)-disty(i,jref))/onedeg
c            grdlon(i,j)=reflon+(distx(i,j)-distx(iref,j))/onedeg
            grdlat(i,j)=(grdj(i,j)-jref)*dely/onedeg+reflat
            grdlon(i,j)=(grdi(i,j)-iref)*delx/onedeg+reflon
            if (grdlon(i,j).gt.360.0) grdlon(i,j)=grdlon(i,j)-360.0
            if (grdlon(i,j).lt.  0.0) grdlon(i,j)=grdlon(i,j)+360.0
            hx(i,j)=1.0
            hy(i,j)=1.0
            f(i,j)=sin(flat     *d2r)*omega4
          enddo
        enddo
        return
c
c***********************************************************************
c          spherical grid (igrid=5)
c***********************************************************************
c
      elseif (igrid.eq.5) then
        gcon=3.0
        do j=1,n
          do i=1,m
            xpos(i,j)=1.0
            ypos(i,j)=1.0
c            grdlat(i,j)=(disty(i,j)-disty(i,jref))/onedeg+reflat
c            grdlon(i,j)=(distx(i,j)-distx(iref,j))/onedeg+reflon
            grdlat(i,j)=(grdj(i,j)-jref)*dely/onedeg+reflat
            grdlon(i,j)=(grdi(i,j)-iref)*delx/onedeg+reflon
            if (grdlon(i,j).gt.360.0) grdlon(i,j)=grdlon(i,j)-360.0
            if (grdlon(i,j).lt.  0.0) grdlon(i,j)=grdlon(i,j)+360.0
            f(i,j)=sin(grdlat(i,j)*d2r)*omega4
            hx(i,j)=cos(grdlat(i,j)*d2r)
            hy(i,j)=1.0
          enddo
        enddo
        return
      endif
c
c***********************************************************************
c          format statements
c***********************************************************************
c
  800 format(/,' ',72('-'),/)
  805 format(/,' error from subroutine grid:',/
     1        ,'   igrid must be one of the following values:',//
     2        ,'            1: mercator projection',/
     3        ,'            2: lambert conformal projection',/
     4        ,'            3: polar steographic projection',/
     5        ,'            4: cartesian coordinates',/
     6        ,'            5: spherical projection',//
     7        ,' Your entry was:',i6,', Correct and try again',/)
c
c***********************************************************************
c
      return
      end
















