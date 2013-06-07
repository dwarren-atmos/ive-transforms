
      integer nxmax,nymax,nzmax
#ifdef SIXTY_FOUR_bptr 
      parameter (nxmax=25001,nymax=25001,nzmax=1001)
#else
      parameter (nxmax=5001,nymax=201,nzmax=201)
#endif

      real rdzc(nzmax),rdze(nzmax)
      real xlbnd(nxmax),xrbnd(nxmax),ylbnd(nymax),yrbnd(nymax)
      real zlbnd(nzmax),zrbnd(nzmax)
      common/t_bnd/rdzc,rdze,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd

