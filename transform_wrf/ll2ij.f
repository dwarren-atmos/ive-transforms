      subroutine ll2ij(igrid,reflat,reflon,iref,jref,stdlt1,stdlt2
     1                ,stdlon,delx,dely,grdi,grdj,npts
     2                ,grdlat,grdlon)

c rcs keywords: $RCSfile: ll2ij.f,v $ 
c               $Revision: 1.2 $ $Date: 1998/11/24 23:38:08 $c
c********************************************************************
c********************************************************************
c          SUBROUTINE: ll2ij
c
c          PURPOSE: To compute i- and j-coordinates of a specified
c                   grid given the latitude and longitude points.
c                   All latitudes in this routine start
c                   with -90.0 at the south pole and increase
c                   northward to +90.0 at the north pole.  The
c                   longitudes start with 0.0 at the Greenwich
c                   meridian and increase to the east, so that
c                   90.0 refers to 90.0E, 180.0 is the inter-
c                   national dateline and 270.0 is 90.0W.
c
c          INPUT VARIABLES:
c
c            igrid:  type of grid projection:
c
c                       = 1, mercator projection
c                       = 2, lambert conformal projection
c                       = 3, polar stereographic projection
c                       = 4, cartesian coordinates
c                       = 5, spherical projection
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
c                    if igrid = 2, lambert conformal)
c
c            stdlon: standard longitude of grid (longitude that
c                     points to the north)
c
c            delx:   grid spacing of grid in x-direction
c                    for igrid = 1,2,3 or 4, and 5 delx must be in meters
c
c            dely:   grid spacing (in meters) of grid in y-direction
c                    for igrid = 1,2,3 or 4, and 5 delx must be in meters
c
c            grdlat: latitude of point (grdi,grdj)
c
c            grdlon: longitude of point (grdi,grdj)
c
c            npts:   number of points to find information for
c
c          OUTPUT VARIABLES:
c
c            grdi:   i-coordinate(s) that this routine will generate
c                    information for
c
c            grdj:   j-coordinate(s) that this routine will generate
c                    information for
c
c********************************************************************
c********************************************************************
c
      implicit none
c
c********************************************************************
c
      integer i
      integer igrid
      integer ihem
      integer iref
      integer jref
      integer npts
c
      real alnfix
      real alon
      real check
      real cn1
      real cn2
      real cn3
      real cn4
      real cnx
      real cny
      real con1
      real con2
      real d2r
      real deg
      real delx
      real dely
      real gcon
      real grdi  (npts)
      real grdj  (npts)
      real grdlat(npts)
      real grdlon(npts)
      real ogcon
      real omega4
      real onedeg
      real pi
      real pi2
      real pi4
      real r2d
      real radius
      real reflat
      real reflon
      real rih
      real rrih
      real stdlon
      real stdlt1
      real stdlt2
      real x
      real xih
      real y
      real yih
c
c********************************************************************
c          make sure igrid is an acceptable value
c********************************************************************
c
      if (igrid.lt.1.or.igrid.gt.5) then
        print 800
        print 805,igrid
        print 800
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
      onedeg=radius*2.0*pi/360.0
c
c********************************************************************
c          mercator projection (igrid=1)
c********************************************************************
c
      if (igrid.eq.1) then
        deg=abs(stdlt1)*d2r
        con1=cos(deg)
        con2=radius*con1
        deg=reflat*0.5*d2r
        rih=con2*alog(tan(pi4+deg))
        do i=1,npts
          alon=grdlon(i)+180.0-reflon
          if (alon.lt.  0.0) alon=alon+360.0
          if (alon.gt.360.0) alon=alon-360.0
          grdi(i)=iref+(alon-180.0)*con2/(r2d*delx)
          deg=grdlat(i)*d2r+pi2
          deg=deg*0.5
          grdj(i)=jref+(con2*alog(tan(deg))-rih)/dely
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
          deg=(90.0-grdlat(i)*ihem)*0.5*d2r
          cn4=tan(deg)
          rrih=cn2*(cn4/cn3)**gcon
          check=180.0-stdlon
          alnfix=stdlon+check
          alon=grdlon(i)+check
          if (alon.lt.  0.0) alon=alon+360.0
          if (alon.gt.360.0) alon=alon-360.0
          deg=(alon-alnfix)*gcon*d2r
          x= rrih*sin(deg)
          y=-rrih*cos(deg)*ihem
          grdi(i)=iref+(x-xih)/delx
          grdj(i)=jref+(y-yih)/dely
        enddo
        return
c
c********************************************************************
c          analytic grid (igrid=4)
c********************************************************************
c
      else if (igrid.eq.4) then
        cnx=delx/onedeg
        cny=dely/onedeg
        do i=1,npts
          grdi(i)=iref+(grdlon(i)-reflon)/cnx
          grdj(i)=jref+(grdlat(i)-reflat)/cny
        enddo
        return
c
c********************************************************************
c          spherical grid (igrid=5)
c********************************************************************
c
      else if (igrid.eq.5) then
        cnx=delx/onedeg
        cny=dely/onedeg
        do i=1,npts
          grdi(i)=iref+(grdlon(i)-reflon)/cnx
          grdj(i)=jref+(grdlat(i)-reflat)/cny
        enddo
        return
      endif
c
c********************************************************************
c          format statements
c********************************************************************
c
  800 format(/,' ',72('-'),/)
  805 format(/,' ERROR from subroutine ll2ij:',/
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
