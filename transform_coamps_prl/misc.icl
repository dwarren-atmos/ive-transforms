
      integer nxmax,nymax,nzmax
      parameter (nxmax=201,nymax=201,nzmax=201)

      real rdzc(nzmax),rdze(nzmax)
      real xlbnd(nxmax),xrbnd(nxmax),ylbnd(nymax),yrbnd(nymax)
      real zlbnd(nzmax),zrbnd(nzmax)
      common/t_bnd/rdzc,rdze,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd

