c-----------------------------------------------------------------------
c     $Id: rotpol.icl,v 1.1 1994/11/14 22:37:25 warren Exp $
c     Purpose:
c        This file contains all common blocks which carry global
c        variables related to rotated pole-transformations.
c        It is included in transforms_rp.f and transforms_ecl.f
c
c     History:
c	$Log: rotpol.icl,v $
C Revision 1.1  1994/11/14  22:37:25  warren
C Christoph Schaer's European Model transforms.
C
cRevision 1.1  1993/03/26  07:39:12  schaer
cInitial revision
c
c-----------------------------------------------------------------------

c Location of the rotated pole: 
C       POLPHI   REAL   WAHRE GEOGR. BREITE DES NORDPOLS
C       POLLAM   REAL   WAHRE GEOGR. LAENGE DES NORDPOLS
c       (identity: polphi= 90., pollam= -180.)
c Other vars:
c       zrpi18=360/(2*pi), zpir18=2*pi/360
c       zsinpol = sin(zpir18*polphi)
c       zcospol = cos(zpir18*polphi)
c       zlampol = zpir18*pollam

      common  /rotpol/ polphi,pollam,zsinpol,zcospol,zlampol,
     &                 zrpi18,zpir18
      real	       polphi,pollam,zsinpol,zcospol,zlampol,
     &                 zrpi18,zpir18

