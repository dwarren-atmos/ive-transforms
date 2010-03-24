
      integer nxmax,nymax,nzmax
#ifdef SIXTY_FOUR_bptr
      parameter (nxmax=5001,nymax=1001,nzmax=201)
#else
      parameter (nxmax=5001,nymax=201,nzmax=201)
#endif

      real rdzc(nzmax),rdze(nzmax)
      real xlbnd(nxmax),xrbnd(nxmax),ylbnd(nymax),yrbnd(nymax)
      real zlbnd(nzmax),zrbnd(nzmax)
      common/t_bnd/rdzc,rdze,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd

