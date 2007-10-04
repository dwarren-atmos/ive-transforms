c-----------------------------------------------------------------------
c     $Id: constants.icl,v 1.1 1994/11/14 22:36:59 warren Exp $
c     Purpose:
c        This file contains the common block /constants/ used by the 
c        coordinate transform functions and user-written computations
c        in fcalcfld. Its purpose is to allow these units to have
c        access to a range of infromation. Values in this commom block 
c        are set by subroutine ftr_read.
c            This common block is also utilized by the grbcst-routine 
c        runned on the Cray.
c     History:
c        $Log: constants.icl,v $
C Revision 1.1  1994/11/14  22:36:59  warren
C Christoph Schaer's European Model transforms.
C
cRevision 1.16  1994/11/11  09:15:03  schaer
cAenderung fuer IVE-3-2-beta
cVerbessertes Heading
cStartzeit neu auch auf Datenfile
c
cRevision 1.14  1994/05/16  06:37:47  schaer
cNeue Funktion XYSHIFT.
c
cRevision 1.13  1994/03/22  15:49:40  schaer
cNeue Funktionen: UG,VG,UA,VA, DIV[x_comp:y_comp], SIGMA, M,B,IPV,
c   SETSTACK, INTSTACK, LOWERCASE, UPPERCASE, INTPOL[var], NOMISSDAT[var]
cNeu Theta-Daten auch auf Druck-Flaechen.
c
cRevision 1.12  1994/03/07  08:50:53  schaer
cAenderung bei der Bestimmung des runnames, cleanup.
c
cRevision 1.11  1994/02/17  14:26:40  dani
cChanges for height as vertical coordinate incorporated (Dani Luethi)
c
cRevision 1.10  1993/09/02  09:25:48  schaer
cAdded runname in common /chrconsts/.
c
cRevision 1.9  1993/08/30  08:01:33  schaer
cNeu auch Theta-Koordinaten auf NETcdf-File erlaubt (datatype=3).
cAufnahme von ak und bk in /transcommon/.
c
cRevision 1.8  1993/08/27  11:40:42  schaer
cAb sofort ist derptype=1 der default.
c
cRevision 1.7  1993/08/27  07:01:24  schaer
cNeue Funktion 'privat' fuer individuelle Aenderungen.
cNeue flags datatype und disptype zur Kennzeichnung der Daten (ersetzten
cdie alten flags prsrf und levtype, welche aber vorlaeufig noch
cvorhanden sind).
c
cRevision 1.6  1993/08/19  12:41:53  schaer
cGlobale Attribute der vertikalen Koordinate werden neu aus dem
cDatenfile gelesen, anstelle (1050,0) zu sein.
c
cRevision 1.5  1993/05/28  08:48:56  schaer
cErweiterungen fuer Daten auf reinen Druckflaechen
c
cRevision 1.4  1993/05/25  11:58:12  schaer
cChanges associated with read a var from a file.
c
cRevision 1.3  1993/05/10  07:53:10  schaer
cAdded isentropic coordinates computations.
c
cRevision 1.2  1993/05/04  09:30:00  dani
cinclude vertical derivatives of AK and BK values
c
cRevision 1.1  1993/04/07  13:35:07  schaer
cInitial revision
c
c---------------------------------------------------------------------- 

c     data-definition
c     ---------------
c     datatype: flag specifying the vertical structure of the data-arrays
c               datatype=0: unspecified
c               datatype=1: hybrid sigma model levels
c               datatype=2: pressure surfaces
c               datatype=3: theta surfaces
c     disptype: flag specifying the vertical coordinate used for 
c               display (selected by the user)
c               disptype=0: level-index 
c               disptype=1: pressure levels 
c               disptype=2: theta levels
c               disptype=3: cartesian coordinates
c
c     vertical derivative in p-coordinates
c     ------------------------------------
c     derptype: mode of vertical derivation in pressure coordinates
c               derptype=1: 3-point weighted derivative (Andrea Rossa)
c               derptype=2: 2-point centered derivative
c
c     Access of surface pressure
c     --------------------------
c     In all routines, the surface pressure field should be accessed with
c     a call to
c       subroutine getps(i,j,l)
c     or
c       subroutine getpsi(ri,rj,rl)
c     These routines will provide the surface pressue value at a grid-point
c     (for getps) or at any location within the grid (getpsi) irrespective
c     whether the data is on hybrid or pressure levels. The pointer
c     to the surface pressure array, i.e. %val(ptrps), should only be
c     accessed in computations reserved to hybrid coordinates.
c
c     Shifting of the whole field into the positive x and y-direction
c     ---------------------------------------------------------------
c     Provided through the variable xyshift(2). It shifts the field by
c     a grid-point-increment
c
c     Heading-Type
c     ------------
c     Stored in variable headtyp.
 
      integer    MAXK
      parameter  (MAXK=60)
      real       aklev(MAXK),bklev(MAXK),aklay(MAXK),bklay(MAXK)
      real       aslev(MAXK),bslev(MAXK),aslay(MAXK),bslay(MAXK)
      real       lonmin,lonmax,latmin,latmax,dellon,dellat
      real       presssrf,presstop
      real       xyshift(2)
      integer    starty,startm,startd,starth,starts
      integer    nx,ny,nz,nt
      integer    ptrps,ptrth,ptrzl,ptrpre
      integer    datatype,disptype,derptype,headtyp
      integer    dattyp,datver,cstver

      common /constants/   aklev,bklev,aklay,bklay,
     &                     aslev,bslev,aslay,bslay,
     &                     lonmin,lonmax,latmin,latmax,dellat,dellon,
     &                     presssrf,presstop,
     &                     xyshift,
     &                     starty,startm,startd,starth,starts,
     &                     nx,ny,nz,nt,
     &                     ptrps,ptrth,ptrzl,ptrpre,
     &                     datatype,disptype,derptype,headtyp,
     &                     dattyp,datver,cstver

      character*(80)       cstfiln,extfiln,runname,pathnam
      common /chrconsts/   cstfiln,extfiln,runname,pathnam







