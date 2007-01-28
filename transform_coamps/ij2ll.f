      subroutine ij2ll(igrid,reflat,reflon,iref,jref,stdlt1,stdlt2
     1                ,stdlon,delx,dely,grdi,grdj,npts
     2                ,grdlat,grdlon)

c rcs keywords: $RCSfile: ij2ll.f,v $ 
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
      integer jref
      integer npts
c
      real delx
      real dely
      real grdi   (npts)
      real grdj   (npts)
      real grdlat (npts)
      real grdlon (npts)
      real reflat
      real reflon
      real stdlon
      real stdlt1
      real stdlt2
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
      integer ihem
c
      real angle
      real cn1
      real cn2
      real cn3
      real cn4
      real con1
      real con2
      real d2r
      real deg
      real gcon
      real ogcon
      real omega4
      real onedeg
      real pi
      real pi2
      real pi4
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
c
c          subroutine: ij2ll
c
c          purpose: to compute latitude and longitude of specified
c                   i- and j-points on a grid.  all latitudes in
c                   this routine start with -90.0 at the south pole
c                   and increase northward to +90.0 at the north
c                   pole.  the longitudes start with 0.0 at the
c                   greenwich meridian and increase to the east,
c                   so that 90.0 refers to 90.0e, 180.0 is the
c                   international dateline and 270.0 is 90.0w.
c
c          input variables:
c
c            igrid:  type of grid projection:
c
c                       = 1, mercator projection
c                       = 2, lambert conformal projection
c                       = 3, polar stereographic projection
c                       = 4, cartesian coordinates
c                       = 5, spherical projection
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
c                    if igrid = 2, lambert conformal)
c
c            stdlon: standard longitude of grid (longitude that
c                     points to the north)
c
c            delx:   grid spacing of grid in x-direction
c                    for igrid = 1,2,3 or 4, delx must be in meters
c                    for igrid = 5, delx must be in degrees
c
c            dely:   grid spacing (in meters) of grid in y-direction
c                    for igrid = 1,2,3 or 4, delx must be in meters
c                    for igrid = 5, dely must be in degrees
c
c            grdi:   i-coordinate(s) that this routine will generate
c                    information for
c
c            grdj:   j-coordinate(s) that this routine will generate
c                    information for
c
c            npts:   number of points to find information for
c
c          output variables:
c
c            grdlat: latitude of point (grdi,grdj)
c
c            grdlon: longitude of point (grdi,grdj)
c
c********************************************************************
c
c********************************************************************
c          make sure igrid is an acceptable value
c********************************************************************
c
      if (igrid.lt.1.or.igrid.gt.5) then
        write(6,800)
        write(6,805) igrid
        write(6,800)
        return
      endif
c
c********************************************************************
c          local constants
c********************************************************************
c
      pi=4.0*atan(1.0)
      pi2=pi/2.0
      pi4=pi/4.0
      d2r=pi/180.0
      r2d=180.0/pi
      radius=6371229.0
      omega4=4.0*pi/86400.0
c
c********************************************************************
c          mercator projection (igrid=1)
c********************************************************************
c
      if (igrid.eq.1) then
        gcon=0.0
        deg=abs(stdlt1)*d2r
        con1=cos(deg)
        con2=radius*con1
        deg=reflat*0.5*d2r
        rih=con2*alog(tan(pi4+deg))
        do i=1,npts
          rr=rih +(grdj(i)-jref)*dely
          grdlat(i)=(2.0*atan(exp(rr/con2))-pi2)*r2d
          grdlon(i)=reflon +(grdi(i)-iref)*r2d*delx/con2
          if (grdlon(i).gt.360.0) grdlon(i)=grdlon(i)-360.0
          if (grdlon(i).lt.  0.0) grdlon(i)=grdlon(i)+360.0
        enddo
        return
c
c********************************************************************
c          lambert conformal (igrid=2) or
c          polar stereographic (igrid=3)
c********************************************************************
c
      else if (igrid.eq.2.or.igrid.eq.3) then
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
        do i=1,npts
          x=xih+(grdi(i)-iref)*delx
          y=yih+(grdj(i)-jref)*dely
          rr=sqrt(x*x+y*y)
          grdlat(i)=r2d*(pi2-2.0*atan(cn3*(rr/cn2)**ogcon))*ihem
          xx= x
          yy=-y*ihem
          if (yy.eq.0.0) then
            if (xx.le.0.0) then
              angle=-90.0
            else if (xx.gt.0.0) then
              angle=90.0
            endif
          else
            angle=atan2(xx,yy)*r2d
          endif
          grdlon(i)=stdlon+angle*ogcon
          if (grdlon(i).gt.360.0) grdlon(i)=grdlon(i)-360.0
          if (grdlon(i).lt.  0.0) grdlon(i)=grdlon(i)+360.0
        enddo
        return
c
c********************************************************************
c          analytic grid (igrid=4)
c********************************************************************
c
      else if (igrid.eq.4) then
        onedeg=radius*2.0*pi/360.0
        cn2=delx/onedeg
        do i=1,npts
          grdlat(i)=reflat+(grdj(i)-jref)*cn2
          grdlon(i)=reflon+(grdi(i)-iref)*cn2
          if (grdlon(i).gt.360.0) grdlon(i)=grdlon(i)-360.0
          if (grdlon(i).lt.  0.0) grdlon(i)=grdlon(i)+360.0
        enddo
        return
c
c********************************************************************
c          spherical grid (igrid=5)
c********************************************************************
c
      else if (igrid.eq.5) then
        do i=1,npts
          grdlat(i)= (grdj(i)-float(jref))*dely+reflat
          grdlon(i)= (grdi(i)-float(iref))*delx+reflon
          if (grdlon(i).gt.360.0) grdlon(i)=grdlon(i)-360.0
          if (grdlon(i).lt.  0.0) grdlon(i)=grdlon(i)+360.0
        enddo
        return
      endif
c
c********************************************************************
c          format statements
c********************************************************************
c
  800 format(/,' ',72('-'),/)
  805 format(/,' ERROR from subroutine ij2ll:',/
     1        ,'   igrid must be one of the following values:',//
     2        ,'            1: mercator projection',/
     3        ,'            2: lambert conformal projection',/
     4        ,'            3: polar steographic projection',/
     5        ,'            4: cartesian coordinates',/
     6        ,'            5: spherical projection',//
     7        ,' Your entry was:',i6,', Correct and try again',/)
c
c********************************************************************
c
      return
      end
