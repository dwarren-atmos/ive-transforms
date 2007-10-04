c -----------------------------------------------------------------------------
c     Transformationen fuer EM, ECMWF Daten auf Druck und Modellflaechen
c -----------------------------------------------------------------------------
c     $Id: trafos.f,v 1.1 1994/11/14 22:37:29 warren Exp $
c     $Log: trafos.f,v $
c Revision 1.1  1994/11/14  22:37:29  warren
c Christoph Schaer's European Model transforms.
c
cRevision 1.39  1994/11/11  09:02:12  schaer
cAenderung fuer IVE-3-2-beta
cVerbessertes Heading
cStartzeit neu auch auf Datenfile
c
cRevision 1.37  1994/05/26  11:56:36  schaer
cAufsplitten von calc_field.f in mehrere Files.
c
cRevision 1.36  1994/05/26  11:51:59  schaer
cNeue Masken-Funktion num1<field<num2.
cErweiterete Funktion DMEAN[field:maske].
c
cRevision 1.35  1994/05/16  06:37:47  schaer
cNeue Funktion XYSHIFT.
c
cRevision 1.34  1994/04/12  09:59:27  schaer
cZwei kleinere Bugs behoben.
c
cRevision 1.33  1994/04/12  08:59:58  schaer
cNeu wird in NEW_FIELD nicht mehr read_var (funktioniert nicht immer)
csondern immer GETVAR aufgerufen.
c
cRevision 1.32  1994/03/22  15:49:40  schaer
cNeue Funktionen: UG,VG,UA,VA, DIV[x_comp:y_comp], SIGMA, M,B,IPV,
c   SETSTACK, INTSTACK, LOWERCASE, UPPERCASE, INTPOL[var], NOMISSDAT[var]
cNeu Theta-Daten auch auf Druck-Flaechen.
c
cRevision 1.31  1994/03/07  08:50:53  schaer
cAenderung bei der Bestimmung des runnames, cleanup.
c
cRevision 1.30  1994/02/18  07:43:32  schaer
cAnpassung an z-Koordinaten im Terrain fuer vertikale Querschnitte.
c
cRevision 1.29  1994/02/17  14:29:58  dani
cChanges for height as vertical coordinate incorporated (Dani Luethi)
c
cRevision 1.27  1993/12/17  08:32:51  schaer
cKleine Aenderung in Subroutine int3dm ohne Einfluss auf Resultat.
c
cRevision 1.26  1993/12/17  07:57:35  schaer
cBut in subroutine int3m corrected.
c
cRevision 1.25  1993/12/14  07:24:10  schaer
cAdded error processing in routine trafos.f
cin routine t_settimes.
c
cRevision 1.24  1993/12/02  13:15:32  schaer
cNeue Version von IVE mit neuen Transformationen.
c
cRevision 1.22  1993/10/07  13:36:18  schaer
cAdded new heading-routines from David Bresch.
c
cRevision 1.21  1993/08/30  08:01:33  schaer
cNeu auch Theta-Koordinaten auf NETcdf-File erlaubt (datatype=3).
cAufnahme von ak und bk in /transcommon/.
c
cRevision 1.20  1993/08/27  11:40:42  schaer
cAb sofort ist derptype=1 der default.
c
cRevision 1.19  1993/08/27  09:05:56  schaer
cAdded 3-point vertical derivation (centered is default) from Andrea Rossa.
c
cRevision 1.18  1993/08/27  07:01:24  schaer
cNeue Funktion 'privat' fuer individuelle Aenderungen.
cNeue flags datatype und disptype zur Kennzeichnung der Daten (ersetzten
cdie alten flags prsrf und levtype, welche aber vorlaeufig noch
cvorhanden sind).
c
cRevision 1.17  1993/08/25  13:17:39  schaer
cKleiner Bug in Version 1.16 behoben.
c
cRevision 1.16  1993/08/23  14:36:34  schaer
cNeue Philosophie mit PS: Falls nicht auf File, wird kein neuer
cButton generiert, sondern Verwaltung des Memories geschieht
cintern in ftr_read.
c
cRevision 1.15  1993/08/20  06:41:54  schaer
cNeue routine freevar, zum expliziten freimachen von temporaeren
cVariablen.
c
cRevision 1.14  1993/08/19  12:41:53  schaer
cGlobale Attribute der vertikalen Koordinate werden neu aus dem
cDatenfile gelesen, anstelle (1050,0) zu sein.
c
cRevision 1.13  1993/06/25  15:11:55  schaer
cErlaube simultane Felder mit nt=1 und nt=ntimes>1.
c
cRevision 1.12  1993/06/03  14:31:37  schaer
cNeue Berechnung von aslay, etc in ftr_read.
c
cRevision 1.11  1993/06/01  07:07:32  schaer
cRevision der Ableitung.
c
cRevision 1.10  1993/05/28  08:48:56  schaer
cErweiterungen fuer Daten auf reinen Druckflaechen
c
cRevision 1.9  1993/05/25  13:40:47  schaer
cSpecial treatment for data without rotated pole.
c
cRevision 1.8  1993/05/25  11:58:12  schaer
cChanges associated with read a var from a file.
c
cRevision 1.7  1993/05/11  11:43:11  schaer
cAttempt to plot surface terrain.
c
cRevision 1.6  1993/05/10  07:53:10  schaer
cAdded isentropic coordinates computations.
c
cRevision 1.5  1993/05/04  09:28:14  dani
cinclude vertical derivatives of AK and BK values
c
cRevision 1.3  1993/04/08  08:20:44  schaer
cBerechnung einer Druckniveau-Tabelle angefuegt.
c
cRevision 1.2  1993/04/07  13:35:53  schaer
cVersion mit Druck-Slicing.
c -----------------------------------------------------------------------------


      subroutine phys_2_lonlat(x, y, lon, lat, npts)
c-----------------------------------------------------------------------------
c     
c     phys_2_lonlat : This routine converts physical coordinates to longitude-
c     latitude coordinates.
c     
c     Arguments:
c     x		real		Physical x coordinate array.
c     y		real		Physical y coordinate array.
c     lon	real		Longitude array (output).
c     lat	real		Latitude array (output).
c     npts	integer		Number of points to convert.
c
      integer npts
      real    x(npts), y(npts), lon(npts), lat(npts)

      integer i

      do i = 1, npts
         call t_ph2ll(y(i),x(i),lat(i),lon(i))
c     identity:
c         lon(i) = x(i)
c         lat(i) = y(i)
      enddo
      return
      end
c
      subroutine lonlat_2_phys(x, y, lon, lat, npts)
c--------------------------------------------------------------------------
c     
c     lonlat_2_phys : This routine converts longitude-latitude coordinates
c     to physical coordinates.
c     
c     Arguments:
c     x		real		Physical x coordinate array (output).
c     y		real		Physical y coordinate array (output).
c     lon	real		Longitude array (output)
c     lat	real		Latitude array (output)
c     npts	integer		Number of points to convert
c
      integer npts
      real    x(npts), y(npts), lon(npts), lat(npts)

      integer i

      do i = 1, npts
         call t_ll2ph(lat(i),lon(i),y(i),x(i))
c     identity:
c         x(i) = lon(i)
c         y(i) = lat(i)
      enddo
      return
      end

      SUBROUTINE T_LL2PH(PHI,LAM,PHIS,LAMS)
c -----------------------------------------------------------------------------
c Purpose:
c     Adaption of the DWD-Functions to convert real geographical coordinates
c     (PHI,LAM) into coordinates in the roteted system (PHIS,LAMS).
c     The Pole is passed trough the common block  /rotpol/. The first two
c     arguments are input, the second two are output. All angles are in 
c     degrees (north>0, east>0)
c History:
c     05/90   D.MAJEWSKI (DWD), G. DE MORSIER (SMA) 
c     03/93   D.BRESCH (ETHZ)


C Common block, containing pole-position:
      include 'rotpol.icl'

c decalration of input/output vars:
      REAL        LAM,PHI,LAMS,PHIS
 
c decalration of internal vars:
      real        zlam,zphi
      real        zarg,zarg1,zarg2

c do case without rotated pole
      if (abs(polphi-90.).lt.1.e-3) then
        phis=phi
        lams=lam
        return
      endif

C first, the conversion of PHI to PHIS:
      ZPHI    = ZPIR18*PHI
      ZLAM    = LAM
      IF(ZLAM.GT.180.0) ZLAM = ZLAM - 360.0
      ZLAM    = ZPIR18*ZLAM
      ZARG    = ZCOSPOL*COS(ZPHI)*COS(ZLAM-ZLAMPOL) + ZSINPOL*SIN(ZPHI)
      PHIS = ZRPI18*ASIN(ZARG)
 
C now, the conversion for LAMS follws: 
      ZPHI    =     ZPIR18*PHI
      ZLAM    = LAM
      IF(ZLAM.GT.180.0) ZLAM = ZLAM - 360.0
      ZLAM    = ZPIR18*ZLAM
      ZARG1   = - SIN(ZLAM-ZLAMPOL)*COS(ZPHI)
      ZARG2   = - ZSINPOL*COS(ZPHI)*COS(ZLAM-ZLAMPOL)+ZCOSPOL*SIN(ZPHI)
      IF (ABS(ZARG2).LT.1.E-30) THEN
        IF (ABS(ZARG1).LT.1.E-30) THEN
          LAMS =   0.0
        ELSEIF (ZARG1.GT.0.) THEN
          LAMS =  90.0
        ELSE
          LAMS = -90.0
        ENDIF
      ELSE
        LAMS = ZRPI18*ATAN2(ZARG1,ZARG2)
      ENDIF
 
      END


      SUBROUTINE T_PH2LL(PHIS,LAMS,PHI,LAM)
c----------------------------------------------------------------------------
c Purpose:
c     Adaption of the DWD-Functions to convert rotated pole coordinates
c     (PHIS,LAMS) into geogrphic coordinates (PHI,LAM). The location of 
c     the rotated pole is passed trough the common block  /rotpol/. The 
c     first two arguments are input, the second two are output. All angles 
c     are in degrees (north>0, east>0)
c History:
c     05/90   D.MAJEWSKI (DWD)
c     03/93   D.BRESCH (ETHZ)

C Common block, containing pole-position:
      include 'rotpol.icl'
 
c declaration of arguments:
      REAL        LAMS,PHIS,PHI,LAM

c declaration of internal vars:
      real        zarg1,zarg2
      real        zphis,zlams,arg

c do case without rotated pole
      if (abs(polphi-90.).lt.1.e-3) then
        phi=phis
        lam=lams
        return
      endif

C first, the conversion of PHIS to PHI:
      ZPHIS  = ZPIR18*PHIS
      ZLAMS  = LAMS
      IF(ZLAMS.GT.180.0) ZLAMS = ZLAMS - 360.0
      ZLAMS  = ZPIR18*ZLAMS
      ARG     = ZCOSPOL*COS(ZPHIS)*COS(ZLAMS) + ZSINPOL*SIN(ZPHIS)
      PHI = ZRPI18*ASIN(ARG)
 
c follows conversion of LAMS to LAM:
      ZPHIS   = ZPIR18*PHIS
      ZLAMS   = LAMS
      IF(ZLAMS.GT.180.0) ZLAMS = ZLAMS - 360.0
      ZLAMS   = ZPIR18*ZLAMS
      ZARG1   = SIN(ZLAMPOL)*(- ZSINPOL*COS(ZLAMS)*COS(ZPHIS)  +
     1                          ZCOSPOL*           SIN(ZPHIS)) -
     2          COS(ZLAMPOL)*           SIN(ZLAMS)*COS(ZPHIS)
      ZARG2   = COS(ZLAMPOL)*(- ZSINPOL*COS(ZLAMS)*COS(ZPHIS)  +
     1                          ZCOSPOL*           SIN(ZPHIS)) +
     2          SIN(ZLAMPOL)*           SIN(ZLAMS)*COS(ZPHIS)
      IF (ABS(ZARG2).LT.1.E-30) THEN
        IF (ABS(ZARG1).LT.1.E-30) THEN
          LAM =   0.0
        ELSEIF (ZARG1.GT.0.) THEN
          LAM =  90.0
        ELSE
          LAM = -90.0
        ENDIF
      ELSE
        LAM = ZRPI18*ATAN2(ZARG1,ZARG2)
      ENDIF
c      if (lam.gt.360.)  lam=360.
c      if (lam.lt.-360.) lam=0.
      END


      subroutine index_2_phys (phypt, compt, iflag, ndims, npts)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutine, given a set of points in the computational 
c        grid, calculates the location of those points in the physical 
c        grid.
c     Arguments:
c        phypt  real  output  an array(ndims,npts) containing points on 
c                             the physical grid axes.
c        compt  real  input   an array(ndims,npts) containing points on 
c                             the computational grid axes.
c        iflag  int   input   a vector of flags that determines whether
c                             a given point in the computational grid 
c                             will be converted to a point in the
c                             physical grid.
c                             iflag(i) = 1  compt(1) will be 
c                                           converted to phypt(1)
c                             otherwise compt(i) will not be converted.
c        ndims  int   input   number of dimensions in both physical and
c                             computational space.
c        npts   int   input   the number of points to be transformed 
c                             in the current call to this routine. 
c                             Normally, npts will be 1, except when 
c                             3D isosurfaces are computed. 
c     History:
c-----------------------------------------------------------------------

      include 'attributes.icl'
      include 'constants.icl'

c     Argument declarations.
      integer  iflag, ndims, npts
      real     compt, phypt
      dimension iflag(ndims)
      dimension compt(ndims,npts), phypt(ndims,npts)

c     Local variable declarations.
      integer       i, jj, jjp1, id
      real          fract, psrf, int4d, getpsi

c     statement-functions for the computation of pressure and potential
c     temperature
      real      pr,theta,zlay,plev
      integer   is
      pr(is)=ak(is)+bk(is)*psrf
      theta(is)=int4d(%val(ptrth),nx,ny,nz,nt,
     &   compt(1,i),compt(2,i),real(is)+trstag(3)+0.5,compt(4,i))
      zlay(is)=int4d(%val(ptrzl),nx,ny,nz,nt,
     &   compt(1,i),compt(2,i),real(is)+trstag(3)+0.5,compt(4,i))
      plev(is)=int4d(%val(ptrpre),nx,ny,nz,nt,
     &   compt(1,i),compt(2,i),real(is)+trstag(3),compt(4,i))

      do 10 i = 1, npts
        do id=1,2
          if (iflag(id).eq.1) then
c           Convert x and y-index to a point in the physical grid.
            if ( trdims(id) .eq. 1 ) then
               phypt(id,i)=trdmin(id)
            else
               phypt(id,i)=(compt(id,i)-1.0+xyshift(id))*trdelt(id)
     &                             +trdmin(id)
c              the xyshift-variable allows to shift the field by some
c              grid-point-increment into the positive x and y-direction
            endif
          endif
        enddo

        if ( iflag(3) .eq. 1 ) then
c         Convert z-index to a point in the physical grid.
          if (trdims(3).eq.1) then
c           --------------------------------------------------------
c           the selected field is one-dimensional in the z-direction
c           --------------------------------------------------------
            phypt(3,i)=trdmin(3)
          else if (disptype.eq.0) then
c           --------------------------------------------------------
c           vertical coordinates is array-index
c           --------------------------------------------------------
            if (datatype.eq.1) then
              phypt(3,i)=real(trdims(3))-compt(3,i)+1
            else
              phypt(3,i)=compt(3,i)
            endif
          else if ((disptype.eq.1).and.(datatype.lt.3)) then
c           --------------------------------------------------------
c           vertical coordinate is pressure
c           --------------------------------------------------------
c           compute surface pressure at current grid-point
            psrf=getpsi(compt(1,i),compt(2,i),compt(4,i))
c           do a linear interpolation
            jj  =min0(int(compt(3,i)),nz-1)
            jjp1=jj+1
            fract=compt(3,i)-float(jj)
            phypt(3,i)=(1.-fract)*pr(jj)+fract*pr(jjp1)
          else if ((disptype.eq.2).and.(datatype.lt.3)) then
c           --------------------------------------------------------
c           pressure coordinates, display in theta
c           --------------------------------------------------------
c           do a linear interpolation
            if (trstag(3).eq.0.) then
              jj  =max0(1,min0(int(compt(3,i)),nz-2))
            else
              jj  =max0(1,min0(int(compt(3,i)),nz-1))
            endif
            jjp1=jj+1
            fract=compt(3,i)-float(jj)
            phypt(3,i)=(1.-fract)*theta(jj)+fract*theta(jjp1)
          else if ((disptype.eq.3).and.(datatype.lt.3)) then
c           --------------------------------------------------------
c           pressure coordinates, display in cartesian coordinates
c           --------------------------------------------------------
c           do a linear interpolation
            if (trstag(3).eq.0.) then
              jj  =max0(1,min0(int(compt(3,i)),nz-2))
            else
              jj  =max0(1,min0(int(compt(3,i)),nz-1))
            endif
            jjp1=jj+1
            fract=compt(3,i)-float(jj)
            phypt(3,i)=(1.-fract)*zlay(jj)+fract*zlay(jjp1)
          else if ((disptype.eq.2).and.(datatype.eq.3)) then
c           --------------------------------------------------------
c           vertical coordinate (data and display) is theta
c           --------------------------------------------------------
            jj=min0(int(compt(3,i)),nz-1)
            fract=compt(3,i)-float(jj)
            phypt(3,i)=(1.-fract)*aklev(jj)+fract*aklev(jj+1)
          else if ((disptype.eq.1).and.(datatype.eq.3)) then
c           --------------------------------------------------------
c           theta coordinates, display in pressure coordinates
c           --------------------------------------------------------
c           do a linear interpolation
            jj  =max0(1,min0(int(compt(3,i)),nz-1))
            jjp1=jj+1
            fract=compt(3,i)-float(jj)
            phypt(3,i)=(1.-fract)*plev(jj)+fract*plev(jjp1)
          else
            print *,'Illegal combination of data and display-type'
            stop
          endif
        endif

        if ( iflag(4) .eq. 1 ) then
c         Convert t-index to a point in the physical grid.
          if (ntimes.eq.0) then
            phypt(4,i) = 0.
          else if (ntimes.eq.1) then
            phypt(4,i) = timeval(1)
          else if (trdims(4).eq.1) then
            phypt(4,i) = timeval(1)
          else
            jj  =max(1,min0(int(compt(4,i)),ntimes-1))
            jjp1=jj+1
            fract=compt(4,i)-float(jj)
            phypt(4,i) = (1.0-fract)*timeval(jj)+fract*timeval(jjp1)
         endif
      endif
 10   continue

      end

      subroutine phys_2_index(phypt, compt, iflag, ndims, npts)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutine, given a set of points in the physical 
c        grid, calculates the location of those points in the 
c        computational grid.
c     Arguments:
c        phypt  real  input   an array(ndims,npts) containing points on 
c                             the physical grid axes.
c        compt  real  output  an array(ndims,npts) containing points on  
c                             the computational grid axes.
c        iflag  int   input   a vector of flags that determines whether
c                             a given point in the physical grid 
c                             will be converted to a point in the
c                             computational grid.
c                             iflag(i) = 1  phypt(1) will be 
c                                           converted to compt(1)
c                             otherwise phypt(i) will not be converted.
c        ndims  int   input   number of dimensions in both physical and
c                             computational space.
c        npts   int   input   the number of points to be transformed 
c                             in the current call to this routine. 
c                             Normally, npts will be 1, except when 
c                             3D isosurfaces are computed. 
c     History:
c-----------------------------------------------------------------------

      include 'attributes.icl'
      include 'constants.icl'

c     Argument declarations.
      integer        iflag, ndims, npts
      real           phypt, compt
      dimension iflag(ndims)
      dimension phypt(ndims,npts), compt(ndims,npts)

c     Local variable declarations.
      integer i, j, jj, jp1, id, kth1,kth2
      real    compt4, fractj, fract, psrf, pbot, ptop, prmore, prless
      real    tbot, ttop, thmore, thless
      real    int4d, getpsi

c     statement-functions for the computation of pressure and potential
c     temperature
      real      pr,theta,zlay,plev
      integer   is
      pr(is)=ak(is)+bk(is)*psrf
      theta(is)=int4d(%val(ptrth),nx,ny,nz,nt,
     &   compt(1,i),compt(2,i),real(is)+trstag(3)+0.5,compt(4,i))
      zlay(is)=int4d(%val(ptrzl),nx,ny,nz,nt,
     &   compt(1,i),compt(2,i),real(is)+trstag(3)+0.5,compt(4,i))
      plev(is)=int4d(%val(ptrpre),nx,ny,nz,nt,
     &   compt(1,i),compt(2,i),real(is)+trstag(3),compt(4,i))


      if ((iflag(3).eq.1).and.
     &         ((disptype.ge.2).and.(datatype.lt.3))) then
c       compute the range of vertical indices which is inside the theta-
c       range.
        if (trstag(3).eq.0.) then
          kth1=1
          kth2=nz-1
        else
          kth1=1
          kth2=nz
        endif
      endif

      do 10 i = 1, npts

c       do conversion of x and y-coordinate
        do id=1,2
          if (iflag(id).eq.1) then
c           Convert phypt(id,i) to a point in the computational grid.
            if ( trdims(id) .eq. 1 ) then
              compt(id,i)=1.0
            else
              compt(id,i)=(phypt(id,i)-trdmin(id))*trdelr(id)+1.0
     &                    -xyshift(id)
c             the xyshift-variable allows to shift the field by some
c             grid-point-increment into the positive x and y-direction
            endif
          endif
        enddo

c       do conversion of t-coordinate
          if (trdims(4).eq.1) then
            compt4=1.
          else if (ntimes.gt.1) then
            do jj = 1,ntimes
              if (timeval(jj).ge.phypt(4,i)) then
                j=jj-1
                go to 8
              endif
            enddo
 8          continue
            j=max0(1,min0(ntimes-1,j))
            jp1 = j+1
            fractj  = (phypt(4,i)-timeval(j))/(timeval(jp1)-timeval(j))
            compt4 = float(j)*(1.0-fractj) + float(jp1)*fractj
          else
            compt4=1.
          endif
          if (iflag(4).eq.1) then
            compt(4,i)=compt4
          endif

c       do conversion of z-coordinate
        if (iflag(3).eq.1) then
          if (trdims(3).eq.1) then
c           --------------------------------------------------------
c           the selected field is one-dimensional in the z-direction
c           --------------------------------------------------------
            compt(3,i)=1.

          else if (disptype.eq.0) then
c           --------------------------------------------------------
c           vertical coordinates is array-index
c           --------------------------------------------------------
            if ((phypt(3,i).lt.1.).or.
     &          (phypt(3,i).gt.real(trdims(3)))) then
              compt(3,i)=-1.
            else
              if (datatype.eq.1) then
                compt(3,i)=real(trdims(3))-phypt(3,i)+1
              else
                compt(3,i)=phypt(3,i)
              endif
            endif

          else if ((disptype.eq.1).and.(datatype.lt.3)) then
c           --------------------------------------------------------
c           vertical coordinate is pressure
c           --------------------------------------------------------
            if (phypt(3,i).eq.presssrf) then
              compt(3,i)=1.
            else if (phypt(3,i).eq.0.) then
              compt(3,i)=real(trdims(3))
            else
c             compute surface pressure at current grid-point
              psrf=getpsi(compt(1,i),compt(2,i),compt4)
c             compute pressure-range at current grid-point
              pbot=pr(1)
              ptop=pr(trdims(3))
c             Check to see if the point in physical space is outside the
c             domain. If it is, return a value < zero.
c             print *,'pr,prtop,prbot: ',phypt(3,i),ptop,pbot
              if ((phypt(3,i).gt.pbot).or.(phypt(3,i).lt.ptop)) then
                compt(3,i)=-1.
              else
c               Do a linear interpolation. First find the level-index which 
c               is just above phypt(3,i).
                jj=0
                do j=2,nz
                  prless=pr(j)
                  if (phypt(3,i).ge.prless) then
                    prmore=pr(j-1)
                    jj=j
                    goto 100
                  endif
                enddo
 100            continue
                if (jj.eq.0) then
                  compt(3,i)=-1.
                else
c                 print *,'pr,prles,prmor: ',phypt(3,i),prless,prmore
                  fract=(phypt(3,i)-prless)/(prmore-prless)
                  compt(3,i)=(1.-fract)*jj+fract*(jj-1)
                endif
              endif
            endif

          else if ((disptype.eq.1).and.(datatype.eq.3)) then
c           --------------------------------------------------------
c           data on theta-levels, vertical coordinate is pressure
c           --------------------------------------------------------
c           compute pressure-range at current grid-point
            tbot=plev(1)
            ttop=plev(nz)
c           Check to see if the point in physical space is outside the
c           domain. If it is, return a value < zero.
            if ((phypt(3,i).gt.tbot).or.(phypt(3,i).lt.ttop)) then
              compt(3,i)=-1.
c             print *,'outside: te,tetop,tebot: ',phypt(3,i),ttop,tbot
            else
c             Do a linear interpolation. First find the level-index which 
c             is just above phypt(3,i). 
              do j=2,nz
                thmore=plev(j)
                if (phypt(3,i).ge.thmore) then
                  thless=plev(j-1)
                  jj=j
                  goto 207
                endif
              enddo
 207          continue
              if (jj.eq.0) then
c               print *,'jj not found'
                compt(3,i)=-1.
              else
c               print *,'th,thles,thmor: ',phypt(3,i),thless,thmore
                fract=(phypt(3,i)-thless)/(thmore-thless)
                compt(3,i)=(1.-fract)*(jj-1)+fract*jj
              endif
            endif

          else if ((disptype.eq.2).and.(datatype.lt.3)) then
c           --------------------------------------------------------
c           vertical coordinate is potential temperature
c           --------------------------------------------------------
c           compute pot.temperature-range at current grid-point
            tbot=theta(kth1)
            ttop=theta(kth2)
c           Check to see if the point in physical space is outside the
c           domain. If it is, return a value < zero.
            if ((phypt(3,i).lt.tbot).or.(phypt(3,i).gt.ttop)) then
              compt(3,i)=-1.
c             print *,'outside: te,tetop,tebot: ',phypt(3,i),ttop,tbot
            else
c             Do a linear interpolation. First find the level-index which 
c             is just above phypt(3,i). For this, a case distinction for 
c             layers and levels is made.
              do j=kth1+1,kth2
                thmore=theta(j)
                if (phypt(3,i).le.thmore) then
                  thless=theta(j-1)
                  jj=j
                  goto 200
                endif
              enddo
 200          continue
              if (jj.eq.0) then
c               print *,'jj not found'
                compt(3,i)=-1.
              else
c                print *,'th,thles,thmor: ',phypt(3,i),thless,thmore
                fract=(phypt(3,i)-thless)/(thmore-thless)
                compt(3,i)=(1.-fract)*(jj-1)+fract*jj
              endif
            endif

          else if ((disptype.eq.3).and.(datatype.lt.3)) then
c           --------------------------------------------------------
c           vertical coordinate is height
c           --------------------------------------------------------
c           compute pot.temperature-range at current grid-point
            tbot=zlay(kth1)
            ttop=zlay(kth2)
c           Check to see if the point in physical space is outside the
c           domain. If it is, return a value < zero.
            if ((phypt(3,i).lt.tbot).or.(phypt(3,i).gt.ttop)) then
              compt(3,i)=-1.
c             print *,'outside: te,tetop,tebot: ',phypt(3,i),ttop,tbot
            else
c             Do a linear interpolation. First find the level-index which 
c             is just above phypt(3,i). For this, a case distinction for 
c             layers and levels is made.
              do j=kth1+1,kth2
                thmore=zlay(j)
                if (phypt(3,i).le.thmore) then
                  thless=zlay(j-1)
                  jj=j
                  goto 201
                endif
              enddo
 201          continue
              if (jj.eq.0) then
c               print *,'jj not found'
                compt(3,i)=-1.
              else
c                print *,'th,thles,thmor: ',phypt(3,i),thless,thmore
                fract=(phypt(3,i)-thless)/(thmore-thless)
                compt(3,i)=(1.-fract)*(jj-1)+fract*jj
              endif
            endif

          else if ((disptype.eq.2).and.(datatype.eq.3)) then
c           --------------------------------------------------------
c           theta-coordinates (data and display)
c           --------------------------------------------------------
c           Check to see if the point in physical space is outside the
c           domain. If it is, return a value < zero.
            if ((phypt(3,i).lt.ak(1)).or.(phypt(3,i).gt.ak(nz))) then
              compt(3,i)=-1.
            else
c             Do a linear interpolation. First find the level-index which 
c             is just below phypt(3,i).
              jj=nz
              do j=2,nz
                if (phypt(3,i).lt.ak(j)) then
                  jj=j-1
                  goto 300
                endif
              enddo
 300          jj=min0(jj,nz-1)
              fract=amin1( 1.,
     &          (phypt(3,i)-ak(jj))/(ak(jj+1)-ak(jj)) )
              compt(3,i)=(1.-fract)*jj+fract*(jj+1)
cc             testoutput
c              print 780,phypt(3,i),compt(3,i),
c     &                level(jj),level(jj+1),fract
c 780          format ('p->c',5f14.2)

            endif

          else
            print *,'Illegal combination of data and display-type'
            stop
          endif
        endif

 10   continue

c          print *,'compt(4,1),phypt(4,1)',compt(4,1),phypt(4,1)
c          kommentar: selbst wenn time=7.1 bestellt wird, wird routine
c          mit phypt(4)=7.0 aufgerufen

      end


      real function int3d(ar,n1,n2,n3,rid,rjd,rkd)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutine interpolates a 3d-array to an arbitrary
c        location within the grid.
c     Arguments:
c        ar        real  input   surface pressure, define as ar(n1,n2,n3)
c        n1,n2,n3  int   input   dimensions of ar
c        ri,rj,rk  real  input   grid location to be interpolated to
c     History:
c-----------------------------------------------------------------------

c     argument declarations
      integer   n1,n2,n3
      real      ar(n1,n2,n3), rid,rjd,rkd

c     local declarations
      integer   i,j,k,ip1,jp1,kp1,ih,jh,kh
      real      frac0i,frac0j,frac0k,frac1i,frac1j,frac1k,ri,rj,rk

c     do linear interpolation
      ri=amax1(1.,amin1(float(n1),rid))
      rj=amax1(1.,amin1(float(n2),rjd))
      rk=amax1(1.,amin1(float(n3),rkd))
      ih=nint(ri)
      jh=nint(rj)
      kh=nint(rk)

c     Check for interpolation in i
      if (abs(float(ih)-ri).lt.1.e-5) then
        i  =ih
        ip1=ih
      else
        i =min0(int(ri),n1-1)
        ip1=i+1
      endif

c     Check for interpolation in j
      if (abs(float(jh)-rj).lt.1.e-5) then
        j  =jh
        jp1=jh
      else
        j =min0(int(rj),n2-1)
        jp1=j+1
      endif

c     Check for interpolation in k
      if (abs(float(kh)-rk).lt.1.e-5) then
        k  =kh
        kp1=kh
      else
        k =min0(int(rk),n3-1)
        kp1=k+1
      endif

      if (k.eq.kp1) then
c       no interpolation in k
        if ((i.eq.ip1).and.(j.eq.jp1)) then
c          no interpolation at all
           int3d=ar(i,j,k)
c          print *,'int3d 00: ',rid,rjd,rkd,int3d
        else
c          horizontal interpolation only
           frac0i=ri-float(i)
           frac0j=rj-float(j)
           frac1i=1.-frac0i
           frac1j=1.-frac0j
           int3d = ar(i  ,j  ,k  ) * frac1i * frac1j
     &           + ar(i  ,jp1,k  ) * frac1i * frac0j
     &           + ar(ip1,j  ,k  ) * frac0i * frac1j
     &           + ar(ip1,jp1,k  ) * frac0i * frac0j
c          print *,'int3d 10: ',rid,rjd,rkd,int3d
        endif
      else 
        frac0k=rk-float(k)
        kp1=k+1
        frac1k=1.-frac0k
        if ((i.eq.ip1).and.(j.eq.jp1)) then
c          vertical interpolation only
           int3d = ar(i  ,j  ,k  ) * frac1k
     &           + ar(i  ,j  ,kp1) * frac0k
c          print *,'int3d 01: ',rid,rjd,rkd,int3d
        else
c          full 3d interpolation
           frac0i=ri-float(i)
           frac0j=rj-float(j)
           ip1=i+1
           jp1=j+1
           frac1i=1.-frac0i
           frac1j=1.-frac0j
           int3d = ar(i  ,j  ,k  ) * frac1i * frac1j * frac1k
     &           + ar(i  ,jp1,k  ) * frac1i * frac0j * frac1k 
     &           + ar(ip1,j  ,k  ) * frac0i * frac1j * frac1k
     &           + ar(ip1,jp1,k  ) * frac0i * frac0j * frac1k
     &           + ar(i  ,j  ,kp1) * frac1i * frac1j * frac0k
     &           + ar(i  ,jp1,kp1) * frac1i * frac0j * frac0k 
     &           + ar(ip1,j  ,kp1) * frac0i * frac1j * frac0k
     &           + ar(ip1,jp1,kp1) * frac0i * frac0j * frac0k
c          print *,'int3d 11: ',rid,rjd,rkd,int3d
        endif
      endif
      end


      real function int3dm(ar,n1,n2,n3,rid,rjd,rkd,misdat)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutine interpolates a 3d-array to an arbitrary
c        location within the grid. The interpolation includes the 
c        testing of the missing data flag 'misdat'. 
c     Arguments:
c        ar        real  input   surface pressure, define as ar(n1,n2,n3)
c        n1,n2,n3  int   input   dimensions of ar
c        ri,rj,rk  real  input   grid location to be interpolated to
c        misdat    real  input   missing data flag (on if misdat<>0)
c     Warning:
c        This routine has not yet been seriously tested
c     History:
c-----------------------------------------------------------------------

c     argument declarations
      integer   n1,n2,n3
      real      ar(n1,n2,n3), rid,rjd,rkd, misdat

c     local declarations
      integer   i,j,k,ip1,jp1,kp1,ih,jh,kh
      real      frac0i,frac0j,frac0k,frac1i,frac1j,frac1k,ri,rj,rk,int3d

c     check if routine without missing data checking can be called instead
      if (misdat.eq.0.) then
        int3dm=int3d(ar,n1,n2,n3,rid,rjd,rkd)
        return
      endif

c     do linear interpolation
      ri=amax1(1.,amin1(float(n1),rid))
      rj=amax1(1.,amin1(float(n2),rjd))
      rk=amax1(1.,amin1(float(n3),rkd))
      ih=nint(ri)
      jh=nint(rj)
      kh=nint(rk)

c     Check for interpolation in i
      if (abs(float(ih)-ri).lt.1.e-5) then
        i  =ih
        ip1=ih
      else
        i =min0(int(ri),n1-1)
        ip1=i+1
      endif

c     Check for interpolation in j
      if (abs(float(jh)-rj).lt.1.e-5) then
        j  =jh
        jp1=jh
      else
        j =min0(int(rj),n2-1)
        jp1=j+1
      endif

c     Check for interpolation in k
      if (abs(float(kh)-rk).lt.1.e-5) then
        k  =kh
        kp1=kh
      else
        k =min0(int(rk),n3-1)
        kp1=k+1
      endif

      if (k.eq.kp1) then
c       no interpolation in k
        if ((i.eq.ip1).and.(j.eq.jp1)) then
c          no interpolation at all
           int3dm=ar(i,j,k)
c          print *,'int3dm 00: ',rid,rjd,rkd,int3dm
        else
c          horizontal interpolation only
           if ((misdat.eq.ar(i  ,j  ,k  )).or.
     &         (misdat.eq.ar(i  ,jp1,k  )).or.
     &         (misdat.eq.ar(ip1,j  ,k  )).or.
     &         (misdat.eq.ar(ip1,jp1,k  ))) then
             frac0i=ri-float(i)
             frac0j=rj-float(j)
             frac1i=1.-frac0i
             frac1j=1.-frac0j
             int3dm = ar(i  ,j  ,k  ) * frac1i * frac1j
     &              + ar(i  ,jp1,k  ) * frac1i * frac0j
     &              + ar(ip1,j  ,k  ) * frac0i * frac1j
     &              + ar(ip1,jp1,k  ) * frac0i * frac0j
c            print *,'int3dm 10: ',rid,rjd,rkd,int3dm
           else
             int3dm=misdat
           endif
        endif
      else 
        frac0k=rk-float(k)
        kp1=k+1
        frac1k=1.-frac0k
        if ((i.eq.ip1).and.(j.eq.jp1)) then
c          vertical interpolation only
           if ((misdat.eq.ar(i  ,j  ,k  )).or.
     &         (misdat.eq.ar(i  ,j  ,kp1))) then
             int3dm=misdat
           else
             int3dm = ar(i  ,j  ,k  ) * frac1k
     &              + ar(i  ,j  ,kp1) * frac0k
c            print *,'int3dm 01: ',rid,rjd,rkd,int3dm
           endif
        else
c          full 3d interpolation
           if ((misdat.eq.ar(i  ,j  ,k  )).or.
     &         (misdat.eq.ar(i  ,jp1,k  )).or.
     &         (misdat.eq.ar(ip1,j  ,k  )).or.
     &         (misdat.eq.ar(ip1,jp1,k  )).or.
     &         (misdat.eq.ar(i  ,j  ,kp1)).or.
     &         (misdat.eq.ar(i  ,jp1,kp1)).or.
     &         (misdat.eq.ar(ip1,j  ,kp1)).or.
     &         (misdat.eq.ar(ip1,jp1,kp1))) then
             int3dm=misdat
           else
             frac0i=ri-float(i)
             frac0j=rj-float(j)
             frac1i=1.-frac0i
             frac1j=1.-frac0j
             int3dm = ar(i  ,j  ,k  ) * frac1i * frac1j * frac1k
     &              + ar(i  ,jp1,k  ) * frac1i * frac0j * frac1k 
     &              + ar(ip1,j  ,k  ) * frac0i * frac1j * frac1k
     &              + ar(ip1,jp1,k  ) * frac0i * frac0j * frac1k
     &              + ar(i  ,j  ,kp1) * frac1i * frac1j * frac0k
     &              + ar(i  ,jp1,kp1) * frac1i * frac0j * frac0k 
     &              + ar(ip1,j  ,kp1) * frac0i * frac1j * frac0k
     &              + ar(ip1,jp1,kp1) * frac0i * frac0j * frac0k
c            print *,'int3dm 11: ',rid,rjd,rkd,int3dm
           endif
        endif
      endif
      end


      real function int4d(ar,n1,n2,n3,n4,rid,rjd,rkd,rld)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutine interpolates a 4d-array to an arbitrary
c        location within the grid.
c     Arguments:
c        ar        real  input   surface pressure, define as ar(n1,n2,n3)
c        n1,..,n4  int   input   dimensions of ar
c        ri,..,rl  real  input   grid location to be interpolated to
c     History:
c-----------------------------------------------------------------------

c     argument declarations
      integer   n1,n2,n3,n4
      real      ar(n1,n2,n3,n4), rid,rjd,rkd,rld

c     local declarations
      integer   l,lp1,lh
      real      frac0l,frac1l,rl,int3d

c     do linear interpolation in l-direction
      rl=amax1(1.,amin1(float(n4),rld))
      lh=nint(rl)

c     Check for interpolation in l
      if (abs(float(lh)-rl).lt.1.e-5) then
        l  =lh
        lp1=lh
      else
        l =min0(int(rl),n4-1)
        lp1=l+1
      endif

      if (l.eq.lp1) then
c       no interpolation in l
        int4d=int3d(ar(1,1,1,l),n1,n2,n3,rid,rjd,rkd)
      else
c       interpolation in l
        frac0l=rl-float(l)
        frac1l=1.-frac0l
        int4d = int3d(ar(1,1,1,l  ),n1,n2,n3,rid,rjd,rkd) * frac1l
     &        + int3d(ar(1,1,1,lp1),n1,n2,n3,rid,rjd,rkd) * frac0l
      endif
      end


      real function int4dm(ar,n1,n2,n3,n4,rid,rjd,rkd,rld,misdat)
c-----------------------------------------------------------------------
c     Purpose:
c        This subroutine interpolates a 4d-array to an arbitrary
c        location within the grid. The interpolation includes the 
c        testing of the missing data flag 'misdat'. 
c     Arguments:
c        ar        real  input   surface pressure, define as ar(n1,n2,n3)
c        n1,..,n4  int   input   dimensions of ar
c        ri,..,rl  real  input   grid location to be interpolated to
c        misdat    real  input   missing data flag (on if misdat<>0)
c     Warning:
c        This routine has not yet been seriously tested.
c     History:
c-----------------------------------------------------------------------

c     argument declarations
      integer   n1,n2,n3,n4
      real      ar(n1,n2,n3,n4), rid,rjd,rkd,rld, misdat

c     local declarations
      integer   l,lp1,lh
      real      frac0l,frac1l,rl,rint0,rint1,int4d,int3dm

c     check whether missing data checking is required
      if (misdat.eq.0.) then
        int4dm=int4d(ar,n1,n2,n3,n4,rid,rjd,rkd,rld)
        return
      endif

c     do linear interpolation in l-direction
      rl=amax1(1.,amin1(float(n4),rld))
      lh=nint(rl)

c     Check for interpolation in l
      if (abs(float(lh)-rl).lt.1.e-5) then
        l  =lh
        lp1=lh
      else
        l =min0(int(rl),n4-1)
        lp1=l+1
      endif

      if (l.eq.lp1) then
c       no interpolation in l
        int4dm = int3dm(ar(1,1,1,l),n1,n2,n3,rid,rjd,rkd,misdat)
      else
c       interpolation in l
        frac0l=rl-float(l)
        frac1l=1.-frac0l
        rint0 = int3dm(ar(1,1,1,l  ),n1,n2,n3,rid,rjd,rkd,misdat)
        rint1 = int3dm(ar(1,1,1,lp1),n1,n2,n3,rid,rjd,rkd,misdat)
        if ((rint0.eq.misdat).or.(rint1.eq.misdat)) then
          int4dm = misdat
        else
          int4dm = rint0*frac1l + rint1*frac0l
        endif
      endif
      end

      subroutine new_field(name, field, cpmax1, cpmax2, cpmax3, cpmax4)
c-----------------------------------------------------------------------
c
c     Purpose:
c        This routine sets up the common block used by the coordinate
c        transform functions. Its purpose is to allow the transform 
c        functions to have access to all of the data attributes of the
c        current variable. 
c        Use of this common block also improves efficiency, because  
c        values computed here do not have to be recomputed each time a 
c        coordinate transform function is called.      
c     
c     Arguments:
c     name	character	The name of the current field.
c     field	real		Values of "field".
c     cpmax1	integer		Size of dimensions in "Fortran"
c     - cpmax4			order (i.e. inmax1 is x, etc.)
c
c     History:
c-----------------------------------------------------------------------

      include 'attributes.icl'
      include 'constants.icl'

c     Argument declarations:
      character*(*) name
      integer       cpmax1, cpmax2, cpmax3, cpmax4
      real          field(cpmax1, cpmax2, cpmax3, cpmax4)

c     Local variable declarations.
      integer       varid,varpt
      real          infinity, missing
      character *(80) vdisplay_units, vdata_units
      integer       i, k, ndims
      logical       local
      character*(80) dimnam(MAXDIMS)

c     External function declarations.
      integer      getvid, getvar, strbeg, strend

c     get the attributes of the variable
      varid=getvid(name)
      if (varid.le.0) then
        print *,'Unknown field ',name(strbeg(name):strend(name)),
     &          ' in NEW_FIELD'
      endif
      varpt=getvar(name,ndims,trdims,trstag,trdmin,trdmax,missing,
     &               vdata_units, vdisplay_units, dimnam, local)

c     fill up dimension array
      if (ndims.lt.4) then
        do i=ndims+1,4
          trdims(i)=1
        enddo
      endif

c     define trdelt
      infinity = 10.0**38
      do i=1,3
        if (trdims(i).ne.1) then 
          trdelt(i)=(trdmax(i)-trdmin(i))/float(trdims(i)-1)
          if(trdelt(i).ne.0.0) then
            trdelr(i) = 1.0 / trdelt(i)
          else
            trdelr(i) = infinity
          endif
        endif
      enddo

c     define ak and bk appropriate for current field 
      if (trstag(3).eq.0.) then
        do k=1,trdims(3)
          ak(k)=aklev(k)
          bk(k)=bklev(k)
        enddo
      else
        do k=1,trdims(3)
          ak(k)=aklay(k)
          bk(k)=bklay(k)
        enddo
      endif

      return
c     test-output
      print *
      print *,'new_field ',name(strbeg(name):strend(name))
      print *,'dims',trdims(1),trdims(2),trdims(3),trdims(4)
      print *,'delt',trdelt(1),trdelt(2),trdelt(3),trdelt(4) 
      print *,'stag',trstag(1),trstag(2),trstag(3),trstag(4)
      print *,'vmin',trdmin(1),trdmin(2),trdmin(3),trdmin(4)
      print *,'vmax',trdmax(1),trdmax(2),trdmax(3),trdmax(4)
      end



      subroutine copyar(a1,a2,n)
c-----------------------------------------------------------------------
c     copy array a1 onto a2
c-----------------------------------------------------------------------
      integer n,i
      real    a1(n), a2(n)
      do i=1,n
        a2(i)=a1(i)
      enddo
      end


      subroutine initar(val,ar,n)
c-----------------------------------------------------------------------
c     initialize array ar with val
c-----------------------------------------------------------------------
      integer n,i
      real    val, ar(n)
      do i=1,n
        ar(i)=val
      enddo
      end


      subroutine new_file(ncid,times,ntim)
c-----------------------------------------------------------------------
c     This routine reads a netCDF file that contains all 
c     information needed to calculate coordinate transforms.
c     (all is the adapted subroutine ftr_read (filnam, errflg), see
c     IVE 3-0-2 transforms about it)
c
c     Arguments:
c     ncid	integer		Id of currently opened netCDF file.
c     times	integer		Pointer to array of time values that
c     				the time window slider should use.
c     ntim	integer		Number of values in times (output).
c     				0 => No values in times => no restrictions
c     				on the window slider (output).
c
c     History:
c-----------------------------------------------------------------------

c     netcdf include
      include '/usr/include/netcdf.inc'

c     common blocks containing level-infos, pole-position, infos on em
      include 'rotpol.icl'
      include 'attributes.icl'
      include 'constants.icl'

c     Argument vars:
      integer        ncid, times, ntim

c     constants for rotpol-common
      data zrpi18,zpir18 / 57.2957795 , 0.0174532925 /      

c     from old arguments:
      logical        errflg

c     new local vars:
      character*(80) datfil
      integer        i
c
c     special var: to save timesteps:
      integer        ptrtimes
      save ptrtimes
      data           ptrtimes / 0 /
      data           xyshift/0.,0./


c     Local variable declarations.
      integer        ndims, dimlst(MAXDIMS), addvar
      real           stag(MAXDIMS), phmin(MAXDIMS), phmax(MAXDIMS)
      real           misdat,psrf
      character*(80) data_units, data_display_units
      logical        newbut
c      integer        corner(3), edglen(3)
      integer        strbeg,strend,ierr,idate(5)
      integer        k,cdfid,rcode,getmem,size,varid
      logical        there,error
      character*(9)  chars
      character*(80) dimnam(MAXDIMS)
      data           derptype / 1 /
      data           headtyp /1/


c     statement-functions for the computation of pressure
      real      prlev,prlay
      integer   is
      prlev(is)=aklev(is)+bklev(is)*psrf
      prlay(is)=aklay(is)+bklay(is)*psrf

c     Set error flag initially to false
      errflg = .false.
      rcode=0

c     initialize pointer to surface pressure field and potential temperature
      ptrps=0
      ptrth=0
      ptrzl=0
      datatype=0
      disptype=0

c     reset the xyshift-variable to zero
      xyshift(1)=0.
      xyshift(2)=0.
 
c     get the filename of the data-file
      call getavar('datfil',datfil,error)
      if (error) then
        print *,'Unable to get datafilename in NEW_FILE'
        return
      endif
c     Cut the filename from datfiln (-> get the path):
      pathnam=' '
      runname=datfil(1:strend(datfil))
      do i=strend(datfil),strbeg(datfil),-1
        if ( (datfil(i:i) .eq. '/').or.(datfil(i:i) .eq. ']') )then
          pathnam=datfil(strbeg(datfil):i)
          runname=datfil(i+1:strend(datfil))
          go to 80
        endif
      enddo
 80   continue

c     extract the constants-filename out of the netcdf-file:
      cstfiln=' '
      call ncagt(ncid,NCGLOBAL,'constants_file_name',cstfiln,error)
      if (error.gt.0) then
         cstfiln=' '                     ! spaeter muss hier ein Fehler checking
         print *,'error reading cstfiln' ! folgen
      endif                         
c     prepend it with the path:
      cstfiln=pathnam(strbeg(pathnam):strend(pathnam))//
     &        cstfiln(strbeg(cstfiln):strend(cstfiln))

c     compute filename of external file
      extfiln=cstfiln(strbeg(cstfiln):(strend(cstfiln)-3)) // 'ext'

c     print out these names
      print *
      print *,'runname is ',runname(strbeg(runname):strend(runname))
      print *,'pathnam is ',pathnam(strbeg(pathnam):strend(pathnam))
      print *,'cstfiln is ',cstfiln(strbeg(cstfiln):strend(cstfiln))
      print *,'extfiln is ',extfiln(strbeg(extfiln):strend(extfiln))

c     initialize presssrf and presstop (taken from global attributes)
      call getrarr('plmin',phmin,4,ierr)
      presssrf=phmin(3)
      call getrarr('plmax',phmax,4,ierr)
      presstop=phmax(3)

c     Open the specified netCDF file for reading.
      inquire ( file = cstfiln(strbeg(cstfiln):strend(cstfiln)),
     &          exist = there )
      if (there) cdfid=ncopn(cstfiln(strbeg(cstfiln):strend(cstfiln)),
     &          NCNOWRIT, rcode)

      if ((there).and.(rcode.eq.NCNOERR)) then
c       the file was successfully accessed
        print *
        print *,'Opened file ',
     &               cstfiln(strbeg(cstfiln):strend(cstfiln)),' '

c       get dimensions
        varid = ncdid (cdfid, 'nx', rcode)
        call ncdinq (cdfid, varid, chars, nx, rcode)
        if (rcode.ne.NCNOERR) goto 950

        varid = ncdid (cdfid, 'ny', rcode)
        call ncdinq (cdfid, varid, chars, ny, rcode)
        if (rcode.ne.NCNOERR) goto 950

        varid = ncdid (cdfid, 'nz', rcode)
        call ncdinq (cdfid, varid, chars, nz, rcode)
        if (rcode.ne.NCNOERR) goto 950

c       get pole postion
        varid = ncvid (cdfid, 'pollon', rcode)
        call ncvgt (cdfid, varid, 1, 1, pollam, rcode)
        if (rcode.ne.NCNOERR) goto 950

        varid = ncvid (cdfid, 'pollat', rcode)
        call ncvgt (cdfid, varid, 1, 1, polphi, rcode)
        if (rcode.ne.NCNOERR) goto 950

c       get coordinate-coefficients for levels (Schichtgrenzen)
        varid = ncvid (cdfid, 'aklev', rcode)
        call ncvgt (cdfid, varid, 1, nz, aklev, rcode)
        if (rcode.ne.NCNOERR) goto 950

        varid = ncvid (cdfid, 'bklev', rcode)
        call ncvgt (cdfid, varid, 1, nz, bklev, rcode)
        if (rcode.ne.NCNOERR) goto 950

c       get coordinate-coefficients for layers (Schichtmitten)
        varid = ncvid (cdfid, 'aklay', rcode)
        call ncvgt (cdfid, varid, 1, nz, aklay, rcode)
        if (rcode.ne.NCNOERR) goto 950

        varid = ncvid (cdfid, 'bklay', rcode)
        call ncvgt (cdfid, varid, 1, nz, bklay, rcode)
        if (rcode.ne.NCNOERR) goto 950

c       get domain
        varid = ncvid (cdfid, 'lonmin', rcode)
        call ncvgt (cdfid, varid, 1, 1, lonmin, rcode)
        if (rcode.ne.NCNOERR) goto 950

        varid = ncvid (cdfid, 'lonmax', rcode)
        call ncvgt (cdfid, varid, 1, 1, lonmax, rcode)
        if (rcode.ne.NCNOERR) goto 950

        varid = ncvid (cdfid, 'latmin', rcode)
        call ncvgt (cdfid, varid, 1, 1, latmin, rcode)
        if (rcode.ne.NCNOERR) goto 950

        varid = ncvid (cdfid, 'latmax', rcode)
        call ncvgt (cdfid, varid, 1, 1, latmax, rcode)
        if (rcode.ne.NCNOERR) goto 950

        varid = ncvid (cdfid, 'dellon', rcode)
        call ncvgt (cdfid, varid, 1, 1, dellon, rcode)
        if (rcode.ne.NCNOERR) goto 950

        varid = ncvid (cdfid, 'dellat', rcode)
        call ncvgt (cdfid, varid, 1, 1, dellat, rcode)
        if (rcode.ne.NCNOERR) goto 950

c         get start from constants-file         
c         varid = ncvid (cdfid, 'starty', rcode)
c         call ncvgt (cdfid, varid, 1, 1, starty, rcode)
c         if (rcode.ne.NCNOERR) goto 950
c         varid = ncvid (cdfid, 'startm', rcode)
c         call ncvgt (cdfid, varid, 1, 1, startm, rcode)
c         if (rcode.ne.NCNOERR) goto 950
c         varid = ncvid (cdfid, 'startd', rcode)
c         call ncvgt (cdfid, varid, 1, 1, startd, rcode)
c         if (rcode.ne.NCNOERR) goto 950
c         varid = ncvid (cdfid, 'starth', rcode)
c         call ncvgt (cdfid, varid, 1, 1, starth, rcode)
c         if (rcode.ne.NCNOERR) goto 950
c         starts=0

c       get start-date of time-axes.
c       try to get the start-date from the data-file.
        call getstart(ncid,idate,rcode)
        if (rcode.eq.0) then
c         blank out button name for starts, etc
          call set_button_name('starty',' ')
          call set_button_name('startm',' ')
          call set_button_name('startd',' ')
          call set_button_name('starth',' ')
          call set_button_name('starts',' ')
        else
c         try to get the start-date from the constants-file.
          call getstart(cdfid,idate,rcode)
          if (rcode.ne.NCNOERR) goto 950
        endif
c       convert to variables in constants-common
        starty=idate(1)
        startm=idate(2)
        startd=idate(3)
        starth=idate(4)
        starts=idate(5)

        varid = ncvid (cdfid, 'dattyp', rcode)
        call ncvgt (cdfid, varid, 1, 1, dattyp, rcode)
        if (rcode.ne.NCNOERR) goto 950

c       define the times-array
        call t_settimes(ncid,timeval,ntimes,ierr)
c        print *,'time stuff------------------'
c        print *,'ntimes: ',ntimes
c        do i=1,ntimes
c            print *,'timeval: ',timeval(i)
c        enddo
c        print *,'time stuff---------------end'
c       do proper defaults
        if (ntimes.lt.1) then
          ntimes=1
          timeval(1)=0.
        endif
c       do the time pointer-stuff:
        if (ptrtimes.ne.0) then
           call freemem(ptrtimes)
        endif
        ptrtimes=getmem(ntimes)
c       copy the content of timeval to ptrtimes
        call copyar(timeval,%val(ptrtimes),ntimes)
        times=ptrtimes
        ntim=ntimes

c       get surface pressure field
        ptrps=addvar ('PS', 'PS', ndims, dimlst, stag, phmin, phmax, 
     &        misdat, data_units, data_display_units, dimnam, newbut)

        if (ptrps.ne.0) then
c         --------------------------
c         surface pressure available
c         --------------------------
          datatype=1
          disptype=1
          if ((nx.ne.dimlst(1)).or.(ny.ne.dimlst(2)).or.
     &        (lonmin.ne.phmin(1)).or.(lonmax.ne.phmax(1)).or.
     &        (latmin.ne.phmin(2)).or.(latmax.ne.phmax(2))) then
            print *
            print *,'WARNING: Surface pressure and constants-file have'
            print *,'         incompatible attributes. Using surface'
            print *,'         pressure attributes.'
            print *
            nx=dimlst(1)
            ny=dimlst(2)
            lonmin=phmin(1)
            lonmax=phmax(1)
            latmin=phmin(2)
            latmax=phmax(2)
          endif
          if (ndims.eq.4) then
            nt=dimlst(4)
          else
            nt=1
          endif
          if (ntim.ne.nt) then
            print *
            print *,'WARNING: Surface pressure and times-array are',
     &              'incompatible.'
            print *
          endif  
        else
        
c         Derive surface pressure. Note that this should only occur if
c         the data is on pressure- or theta-surfaces (otherwise PS should 
c         be on the data file). 

c         check whether the data is on pressure surfaces
          datatype=2
          do k=1,nz
            if ((bklev(k).ne.0.).or.(bklay(k).ne.0.)) datatype=0
          enddo

          if (datatype.ne.0) then
c           data is on pressure or theta levels, check which of the two
            if ((aklev(1).gt.aklev(2)).or.(aklay(1).gt.aklay(2))) then
c             ----------------------------
c             data is on pressure surfaces
c             ----------------------------
              datatype=2
              disptype=1
            else
c             -------------------------
c             data is on theta surfaces
c             -------------------------
              datatype=3
              disptype=2
            endif
          else
c           -------------------
c           unknown coordinates
c           -------------------
            disptype=0
            print *
            print *,'WARNING: Unable to read surface pressure PS!'
            print *,'WARNING: Select vertical levels with array-index'
            print *,'         rather than pressure.'
            print *
          endif
        endif

      else

c       the file could not be accessed
        print *,'WARNING: unable to open file ',
     &                   cstfiln(strbeg(cstfiln):strend(cstfiln))
        print *,'WARNING: Using default pole postion.'
        print *,'WARNING: Select vertical levels with array-index'
        print *,'         rather than pressure.'

c       put filename into constants common block
        cstfiln=' '
        extfiln=' '

c       hardwire default pole postion
        pollam=-170.
        polphi=  32.5

c       set starting date and time
        starty=0
        startm=0
        startd=0
        starth=0
      endif

c     print out some important flags
      print *,'datatype=',datatype,',    disptype=',disptype

c     test-output
      write (6,151) 'nx, ny, nz:           ',nx,ny,nz
      write (6,152) 'lonmin,lonmax,dellon: ',lonmin,lonmax,dellon
      write (6,152) 'latmin,latmax,dellat: ',latmin,latmax,dellat
      write (6,*)
      write (6,155) 'Simulation started on ',
     &                (starty*100+startm)*100+startd,' at ',starth,' h'
      write (6,*)
      write (6,*) '  k(ive)  k(em)      levels      layers',
     &                           '       aklev       bklev'
      psrf=1000.
      write (6,160),0,nz+1,psrf,aklev(k),bklev(k)
      do k=1,nz
        write (6,160) k,nz+1-k,prlev(k),prlay(k),aklev(k),bklev(k)
      enddo

c     compute constants for rotpol-common
      zsinpol = sin(zpir18*polphi)
      zcospol = cos(zpir18*polphi)
      zlampol = zpir18*pollam

c     compute vertical derivatives of ak's and bk's
      do k=2,nz-1
        aslay(k)=(aklay(k-1)-aklay(k+1))/2.
        bslay(k)=(bklay(k-1)-bklay(k+1))/2.
        aslev(k)=(aklev(k-1)-aklev(k+1))/2.
        bslev(k)=(bklev(k-1)-bklev(k+1))/2.
      enddo
      aslay(1 )=aklay(1)-aklay(2)
      bslay(1 )=bklay(1)-bklay(2)
      aslev(1 )=aklev(1)-aklev(2)
      bslev(1 )=bklev(1)-bklev(2)
      aslev(nz)=aklev(nz-1)-aklev(nz)
      bslev(nz)=bklev(nz-1)-bklev(nz)
      aslay(nz)=aklay(nz-1)-aklay(nz)
      bslay(nz)=bklay(nz-1)-bklay(nz)

      return

c     error-exits
 900  continue
      write (6,100)
      return
 950  continue
      write (6,110)
      return

c     formats
 100  format (1x,'Unable to open file ',a,'.',
     &        ' Level information HAS NOT BEEN READ!')
 110  format (1x,'Unable to read from file ',a,'.',
     &        ' Level information HAS NOT BEEN READ!')
 150  format (1x,2i7,4f12.4)
 151  format (1x,a,3i10)
 152  format (1x,a,3f10.3)
 155  format (1x,a,i6,a,i2,a)
 160  format (1x,2i7,4f12.2)
      end


      subroutine getstart(cdfid,idate,ierr)
C------------------------------------------------------------------------
C     Purpose:
C       Get start date for fields on specified NetCDF file. 
C     Arguments:
C       cdfid   int     input   identifier for NetCDF file
C       idate   int     output  array contains date (year,month,day,time,step)
C                               dimensioned as idate(5)
C       ierr    int     output  error flag
C------------------------------------------------------------------------

c     netcdf include
      include '/usr/include/netcdf.inc'

c     variable declarations
      integer	ierr
      integer   idate(5)
      integer	cdfid,ncopts,idvar

c     Get current value of error options, and make sure NetCDF-errors do 
c     not abort execution
      call ncgopt (ncopts)
      call ncpopt (NCVERBOS)

      idvar=ncvid(cdfid,'starty',ierr)
      if (ierr.ne.0) goto 930
      call ncvgt1(cdfid,idvar,1,idate(1),ierr)
      if (ierr.ne.0) goto 920

      idvar=ncvid(cdfid,'startm',ierr)
      if (ierr.ne.0) goto 920
      call ncvgt1(cdfid,idvar,1,idate(2),ierr)
      if (ierr.ne.0) goto 920

      idvar=ncvid(cdfid,'startd',ierr)
      if (ierr.ne.0) goto920
      call ncvgt1(cdfid,idvar,1,idate(3),ierr)
      if (ierr.ne.0) goto 920

      idvar=ncvid(cdfid,'starth',ierr)
      if (ierr.ne.0) goto 920
      call ncvgt1(cdfid,idvar,1,idate(4),ierr)
      if (ierr.ne.0) goto 920

      idvar=ncvid(cdfid,'starts',ierr)
      if (ierr.ne.0) then
        ierr=0
        idate(5)=0
      else
        call ncvgt1(cdfid,idvar,1,idate(5),ierr)
        if (ierr.ne.0) goto 920
      endif

c     normal exit
      call ncpopt (ncopts)
      return

c     error exit
 920  continue
      write (6, *) 'ERROR: An error occurred while attempting to ',
     &             'read the starting-time in subroutine putstart.'
 930  continue
      call ncpopt (ncopts)
      end


      subroutine t_settimes(cdfid,times,ntimes,ierr)
C------------------------------------------------------------------------
C     Purpose:
C        Get all times on the specified NetCDF file
C     Arguments: 
C        cdfid  int  input   identifier for NetCDF file
C        times  real output  array contains all time values on the file,
C                            dimensioned at least times(ntimes)
C        ntimes int  output  number of times on the file
C        error  int  output  errorflag 
C     History:
C        Heini Wernli, ETHZ  
C------------------------------------------------------------------------

      include '/usr/include/netcdf.inc'

c     arguments
      integer   cdfid,ierr
      real      times(*)
      integer   ntimes

c     local variables
      integer   didtim,idtime,i
      integer   ncopts

c     Get current value of error options, and make sure netCDF-errors do 
c     not abort execution
      call ncgopt (ncopts)
      call ncpopt(NCVERBOS)

      didtim=ncdid(cdfid,'time',ierr)   ! inquire id for time dimension
      if (ierr.ne.0) goto 900
      idtime=ncvid(cdfid,'time',ierr)   ! inquire id for time array
      if (ierr.ne.0) goto 900
      call ncdinq(cdfid,didtim,'time',ntimes,ierr)      ! inquire # of times
      if (ierr.ne.0) goto 900
  
      do 10 i=1,ntimes
        call ncvgt1(cdfid,idtime,i,times(i),ierr) ! get times
        if (ierr.ne.0) goto 900
   10 continue
  
c     normal exit
      call ncpopt (ncopts)
      return

c     error exit
 900  ntimes=1
      times(1)=0.
      call ncpopt (ncopts)
      end


      subroutine horiz_ter(topo, nxw, nyw, stagi, stagj, zero, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to fill an array with terrain heights
c        in the windowed domain for horizontal cross sections.
c        NOTE: THIS IS A USER-WRITTEN ROUTINE CALLED BY UWGAP !
c     Arguments:
c        topo   real  output  the array containing terrain heights.
c        nxw    int   input   number of points to be plotted along the 
c                             first dimension.
c        nyw    int   input   number of points to be plotted along the 
c                             second dimension.
c        stagi  real  input   grid staggering along the first dimension.
c        stagj  real  input   grid staggering along the first dimension.
c        zero   real  input   This value is used as the minimum terrain
c                             height.
c        error  int   output  error flag. error = 0  no errors
c                                         error = 1  error
c-----------------------------------------------------------------------

      include 'constants.icl'

c     Argument declarations.
      integer   nxw, nyw
      logical   error
      real      topo, stagi, stagj, zero
      dimension topo (nxw, nyw)

c     variable declarations
      integer   MAXDIM
      parameter (MAXDIM=4)
      integer   getvar,ptrzb,ndims,dimlst(MAXDIM)
      integer   nxwr,nywr,nzwr,iwi,jwi,kwi
      real      stag(MAXDIM), phmin(MAXDIM), phmax(MAXDIM), misdat
      character*80  data_units, data_display_units
      logical   local
      integer   i, j, iflag(MAXDIM)
      real      compt(MAXDIM), phypt(MAXDIM),tmin,tmax,int3d
      character*(80) dimnam(MAXDIM)
      data      iflag / 1,1,0,0 /

c     topography funktioniert nicht mit mapping. Return zero
      do i=1,nxw
      do j=1,nyw
        topo(i,j)=zero
      enddo
      enddo        
      return

c     get pointer to topography-array
      ptrzb=getvar ('ZB', ndims, dimlst, stag, phmin, phmax, 
     &        misdat, data_units, data_display_units, dimnam, local)
      print *,'nxw ,nyw :  ',nxw,nyw
      print *,'nx  ,ny  :  ',dimlst(1),dimlst(2)

c     get the domain size for the topography-domain
      call getdom(iwi,jwi,kwi,nxwr,nywr,nzwr)
      print *,'nxwr,nywr:  ',nxwr,nywr
      print *,'iwi ,jwi :  ',iwi,jwi

c     fill the topography-array (methode 1)
      error = .false.
      tmin=zero+1.e5
      tmax=zero-1.e5
      do i=1,nxw
      do j=1,nyw
        compt(1)=iwi-1+i
        compt(2)=jwi-1+j
        topo(i,j)=amax1( zero, int3d(%val(ptrzb),nx,ny,1,
     &                 compt(1),compt(2),1.) )
        tmin=amin1(tmin,topo(i,j))
        tmax=amax1(tmax,topo(i,j))
      enddo
      enddo
      print *,'zero,tmin,tmax',zero,tmin,tmax
      return

c     fill topography-array (methode 2)
      error = .false.
      tmin=zero+1.e5
      tmax=zero-1.e5
      do i=1,nxw
      do j=1,nyw
c       Find the location of (i, j) in physical space.
        call cpmpxy (0, float(i), float(j), phypt(1), phypt(2))
c       Convert to a point in the data grid.
        call phys_2_index(phypt, compt, iflag, MAXDIM, 1)
c       interpolate topography-array
        topo(i,j)=amax1( zero, int3d(%val(ptrzb),dimlst(1),dimlst(2),1,
     &                 compt(1),compt(2),1.) )
        tmin=amin1(tmin,topo(i,j))
        tmax=amax1(tmax,topo(i,j))
      enddo
      enddo
      print *,'zero,tmin,tmax',zero,tmin,tmax
      end       


      subroutine vert_ter(topoht, topopt, nh, xi, yi, dh, dhsin, dhcos,
     &                   zero)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine finds the terrain height and location under each
c        point in a vertical slab. This routine may have to be modified
c        by the user.
c     Arguments:
c        topoht  real  output  the array containing topography height.
c        topopt  real  output  the array containing the location of each
c                              point in topo.
c        nh      int   input   the number of points in topo and topopt.
c        xi      real  input   the location in physical space of the
c                              beginning of the slice along the first
c                              dimension.
c        yi      real  input   the location in physical space of the
c                              beginning of the slice along the second
c                              dimension.
c        dh      real  input   the interval along the slice.
c        dhsin   real  input   dh * sin(slice angle)
c        dhcos   real  input   dh * cos(slice angle)
c        zero    real  input   the location of the bottom of the plot
c                              window. This value is used as the minimum
c                              terrain height, rather than 0.0.
c     History:
c-----------------------------------------------------------------------

c     Argument declarations.
      integer   nh
      real      topoht, topopt, xi, yi, dh, dhsin, dhcos, zero
      dimension topoht(nh), topopt(nh)

c     include constants file
      include 'constants.icl'

c     Local variable declarations.
      integer   MAXDIMS
      parameter (MAXDIMS=4)
      integer   i, im1, iflag(MAXDIMS)
      integer   getvid,read_var,varid,ptrzb,ndims,dims(MAXDIMS)
      real      int3d
      real      EPSLON, phypt(MAXDIMS), compt(MAXDIMS), getpsi
      real      stag(MAXDIMS),phmin(MAXDIMS),phmax(MAXDIMS),missing
      character*(80) data_units, display_units
      parameter ( EPSLON = 1.0e-10 )


      if (disptype.eq.1) then
        do i = 1, nh
         im1 = float (i - 1)
         phypt(1) = xi + im1 * dhsin
         phypt(2) = yi + im1 * dhcos
         iflag(1)=1  
         iflag(2)=1
         iflag(3)=0
         iflag(4)=0
         call phys_2_index(phypt, compt, iflag, 4, 1)
         topoht(i) = amin1(zero,getpsi(compt(1),compt(2),1.))
         topopt(i) = im1 * dh
c        print *,'VERT_TER output:'
c        print *,i,topoht(i),topopt(i)
        enddo
      else if (disptype.eq.3) then
        varid=getvid('ZB')
        if (varid.gt.0) then
          ptrzb=read_var(varid,ndims,dims,stag,phmin,phmax,missing,
     &               data_units, display_units)
          do i = 1, nh
           im1 = float (i - 1)
           phypt(1) = xi + im1 * dhsin
           phypt(2) = yi + im1 * dhcos
           iflag(1)=1  
           iflag(2)=1
           iflag(3)=0
           iflag(4)=0
           call phys_2_index(phypt, compt, iflag, 4, 1)
           topoht(i) = 0.001 * amax1( zero, 
     &       int3d(%val(ptrzb),dims(1),dims(2),1,compt(1),compt(2),1.) )
           topopt(i) = im1 * dh
c           print *,'VERT_TER output:'
c           print *,i,topoht(i),topopt(i)
          enddo
        else
          print *,'IN ORDER TO GET THE TERRAIN, COMPUTE field=ZB'
          do i=1,nh
           topoht(i) = zero
           topopt(i) = im1 * dh
          enddo
        endif
      else
        do i=1,nh
         topoht(i) = zero
         topopt(i) = im1 * dh
        enddo
      endif
      end


      subroutine intter (ipt, jpt, terhgt)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is used to interpolate the terrain height at the
c        point (ipt, jpt) in the data grid. This routine may 
c        need to be modified by the individual user.
c     Arguments:
c        ipt     real  input  
c        jpt     real  input
c        terhgt  real  output  the interpolated terrain height at point
c                              (ipt, jpt) in the data grid.
c     History:
c-----------------------------------------------------------------------

c     Argument declarations.
      real       ipt, jpt, terhgt
      terhgt = 0.0
      end


      subroutine insert (source, nsource, dest, ndest, idx)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine inserts the value at source(idx), passed by value,
c        into dest(idx).
c     Arguments:
c        source   real  input   the source array, passed by value.
c        nsource  int   input   the length of source.
c        dest     real  in/out  the destination array.
c        ndest    int   input   the length of dest.
c        idx      int   input   the index of the value in source that 
c                               will be copied to the index of dest.
c-----------------------------------------------------------------------

c     Argument declarations.
      real      source, dest
      integer   nsource, ndest, idx
      dimension source(nsource), dest(ndest)

      dest(idx) = source(idx)
      end


      real function getps(i,j,l)
c-----------------------------------------------------------------------------
c     This routine gets the surface pressure at index-location (i,j,1,l). 
c     For data on model surfaces, the value is taken from the the surface 
c     pressure field (stored at pointer-location ptrps). For pressure or
c     other constant-level data, the returned value corresponds to the
c     minimum value of the z-coordinate.
c-----------------------------------------------------------------------------
      integer     i,j,l
      real        getpsh
      include     'constants.icl'
      if (datatype.eq.1) then
        getps=getpsh(%val(ptrps),nx,ny,nt,i,j,l)
      else
        getps=presssrf
      endif
      end

      real function getpsh(ps,nx,ny,nt,i,j,l)
c     Hilfsfunktion fuer getps
      integer nx,ny,nt,i,j,l
      real    ps(nx,ny,nt)
      getpsh=ps(i,j,l)
      end


      real function getpsi(ri,rj,rl)
c-----------------------------------------------------------------------------
c     This routine gets the surface pressure at index-location (ri,rj,1,rl). 
c     The routine is similar to getps, except that an interpolation in
c     index-space is included.
c-----------------------------------------------------------------------------
      include 'attributes.icl'
      real        ri,rj,rl,int3d
      include 'constants.icl'
      if (datatype.eq.1) then
        if ((trdims(1).ne.1).or.(trdims(2).ne.1)) then
          getpsi=int3d(%val(ptrps),nx,ny,nt,ri,rj,rl)
        else
          getpsi=1000.
        endif
      else
        getpsi=presssrf
      endif
c      print *,'getpsi:',getpsi,datatype,trdims(1),trdims(2)
      end
                

      subroutine freevar(var)
c-----------------------------------------------------------------------------
c     This routine is for explicit freeing of a variable by the user.
c     It could later be replaced by an enhanced garbage-collection,
c     which had to be invoked at every level of recursion.
c-----------------------------------------------------------------------------
      character*(*)    var
      integer          fvid, strbeg, strend, getvid
      fvid = getvid(var(strbeg(var):strend(var))//char(0))
      if (fvid.ge.0) call free_var(fvid)
      end








