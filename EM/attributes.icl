c-----------------------------------------------------------------------
c     $Id: attributes.icl,v 1.1 1994/11/14 22:36:52 warren Exp $
c     Purpose:
c        This file contains the common block used by the coordinate
c        transform functions. Its purpose is to allow the transform 
c        functions to have access to all of the data attributes. 
c        Use of this common block also improves efficiency, because  
c        these values do not have to be recomputed each time a 
c        coordinate transform function is called.
c        Values in this commom block are set by subroutine settrcom.
c     History:
c        $Log: attributes.icl,v $
C Revision 1.1  1994/11/14  22:36:52  warren
C Christoph Schaer's European Model transforms.
C
cRevision 1.2  1993/08/30  08:01:33  schaer
cNeu auch Theta-Koordinaten auf NETcdf-File erlaubt (datatype=3).
cAufnahme von ak und bk in /transcommon/.
c
cRevision 1.1  1993/03/26  07:39:12  schaer
cInitial revision
c
c---------------------------------------------------------------------- 

      integer    MAXDIMS,MAXKA
      parameter  (MAXDIMS=4,MAXKA=60)
      integer    trdims(MAXDIMS)
      real                       trdelt(MAXDIMS),trdelr(MAXDIMS)
      real       trstag(MAXDIMS),trdmin(MAXDIMS),trdmax(MAXDIMS)
      real       ak(MAXKA),bk(MAXKA)

      common /transcommon/ trdims,trdelt,trdelr,
     &                     trstag,trdmin,trdmax,
     &                     ak,bk

      integer MAXTIMES, ntimes
      parameter (MAXTIMES = 100)
      real    timeval
      dimension timeval(MAXTIMES)
      common / transtimes / ntimes, timeval


