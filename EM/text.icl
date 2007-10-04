c-----------------------------------------------------------------------
c     $Id: text.icl,v 1.1 1994/11/14 22:37:26 warren Exp $
c     Purpose:
c        This file contains all common blocks which carry global
c        variables related to text appearance on IVE plots from
c        routines in calc_field.f.
c        It is included in privat.f to allow the user to extend
c        plotting and drawing directly.
c
c     History:
c	$Log: text.icl,v $
C Revision 1.1  1994/11/14  22:37:26  warren
C Christoph Schaer's European Model transforms.
C
cRevision 1.1  1994/05/21  19:43:43  schaer
cInitial revision
c
cRevision 1.1  1993/03/26  07:39:12  david
cInitial revision
c
c-----------------------------------------------------------------------
c variables in common block text_com (add my_var here):
      integer         rangecheck !a flag used to control range checking
      integer         defcolor !user defined color
c as in the NCAR manual to plchlq:
      real            tsize,center,orient
c common block for parameter interchange (add my_var here):
      common /text_com / defcolor,rangecheck,center,orient,tsize
