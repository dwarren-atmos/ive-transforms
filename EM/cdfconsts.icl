c-----------------------------------------------------------------------
c     $Id: cdfconsts.icl,v 1.1 1994/11/14 22:36:56 warren Exp $
c     Purpose:
c        This file contains the common block /fileconst/ used by the 
c        put-command and the creation of 3d variables on other levels
c        with create.
c     History:
c        $Log: cdfconsts.icl,v $
C Revision 1.1  1994/11/14  22:36:56  warren
C Christoph Schaer's European Model transforms.
C
cRevision 1.1  1994/03/22  15:27:17  schaer
cInitial revision
c
c---------------------------------------------------------------------- 

c     the parameter MAXK is defined on constants.icl
      real       cdfaklev(MAXK),cdfbklev(MAXK),cdfaklay(MAXK),
     &           cdfbklay(MAXK)
      real       cdfaslev(MAXK),cdfbslev(MAXK),cdfaslay(MAXK),
     &           cdfbslay(MAXK)
      integer    cdfnz,cdfdatatype,cdflowcase

      common /fileconst/   cdfaklev,cdfbklev,cdfaklay,cdfbklay,
     &                     cdfaslev,cdfbslev,cdfaslay,cdfbslay,
     &                     cdfnz,cdfdatatype,cdflowcase

      character*(80)       cdfcstfiln
      common /filechrcnst/ cdfcstfiln
      

