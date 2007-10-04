c     $Id: levels.icl,v 1.1 1994/11/14 22:37:17 warren Exp $
c     This common-block contains level-information for communication
c     within 'ive'. 
c       nnz       number of active levels (nnz<=MAXLEV)
c       iasc      = 1: ascending level-values (e.g. height or theta)
c                 =-1: descending level-values (e.g. pressure)
c       level     the level-values
c
c     $Log: levels.icl,v $
C Revision 1.1  1994/11/14  22:37:17  warren
C Christoph Schaer's European Model transforms.
C
cRevision 1.1  1993/03/26  07:39:12  schaer
cInitial revision
c
cRevision 1.1  1992/12/15  16:34:49  schaer
cInitial revision
c
c
      integer   MAXLEV,nnx,nny,nnz,iasc
      parameter (MAXLEV=100)
      real      level(MAXLEV)
      real      xxmin,xxmax,yymin,yymax
      common    /levels/ nnx,nny,nnz,xxmin,xxmax,yymin,yymax,iasc,level

