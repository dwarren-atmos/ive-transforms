
      subroutine metric_calc(zbot_p,zbot_u,zbot_v,sgz,wgz,dqdx_xy,
     >                         dqdx_zt,dqdx_zw,dqdy_xy,dqdy_zt,
     >                         dqdy_zw,dqdz,dzdq,ztop,dx,dy,
     >                         nx,ny,nz)

      integer nx,ny,nz
      real zbot_p(nx+1,ny+1),zbot_u(nx,ny+1),zbot_v(nx+1,ny)
      real sgz(nz),wgz(nz+1)
      real dqdx_xy(nx+1,ny+1,3),dqdx_zt(nz),dqdx_zw(nz+1)
      real dqdy_xy(nx+1,ny+1,3),dqdy_zt(nz),dqdy_zw(nz+1)
      real dqdz(nx+1,ny+1,3),dzdq(nx+1,ny+1,3)
      real ztop,dx,dy

      integer i,j,k
      real rdx,rdy

      rdx = 1./dx
      rdy = 1./dy
*
* Coordinate metrics.
*    1: therm pts
*    2: u pts
*    3: v pts
* Start with dqdz, dzdq (q is transformed distance in vertical).
*
      do j=1,ny+1
      do i=1,nx+1
        dqdz(i,j,1) = ztop/(ztop-zbot_p(i,j))
        dzdq(i,j,1) = 1./dqdz(i,j,1)
      enddo
      enddo
      
      do j=1,ny+1
      do i=1,nx
        dqdz(i,j,2) = ztop/(ztop-zbot_u(i,j))
        dzdq(i,j,2) = 1./dqdz(i,j,2)
      enddo
      enddo 

      do j=1,ny
      do i=1,nx+1
        dqdz(i,j,3) = ztop/(ztop-zbot_v(i,j))
        dzdq(i,j,3) = 1./dqdz(i,j,3)
      enddo
      enddo
*
* x and y dependent parts of dqdx.  Boundary, then interior.  We use 
*   all the terrain pts here (i.e., u, v, and therm).  Note that 
*   dqdx_xy(i,j,3) is normally not used during model integrations. 
*
      do j=1,ny+1
        dqdx_xy(1,j,1) = ztop/(ztop-zbot_p(1,j))
     >                * 2.*rdx*(zbot_u(1,j)-zbot_p(1,j))
c
        dqdx_xy(1,j,2) = ztop/(ztop-zbot_u(1,j))
     >                * rdx*(zbot_p(2,j)-zbot_p(1,j))
c
        dqdx_xy(nx,j,2) = ztop/(ztop-zbot_u(nx,j))
     >                * rdx*(zbot_p(nx+1,j)-zbot_p(nx,j))
c
        dqdx_xy(nx+1,j,1) = ztop/(ztop-zbot_p(nx+1,j))
     >                * 2.*rdx*(zbot_u(nx,j)-zbot_p(nx+1,j))
c
        dqdx_xy(nx,j,1) = ztop/(ztop-zbot_p(nx,j))
     >      *( 4./3.*rdx*(zbot_u(nx,j)-zbot_u(nx-1,j))
     >          - 1./6.*rdx*(zbot_p(nx+1,j)-zbot_p(nx-1,j)) )
      enddo

      do j=1,ny
        dqdx_xy(1,j,3) = ztop/(ztop-zbot_v(1,j))
     >                * rdx*(zbot_v(2,j)-zbot_v(1,j))
        dqdx_xy(nx+1,j,3) = ztop/(ztop-zbot_v(nx+1,j))
     >                * rdx*(zbot_v(nx,j)-zbot_v(nx+1,j))
      enddo

      do j=1,ny+1
      do i=2,nx-1
c	 
        dqdx_xy(i,j,1) = ztop/(ztop-zbot_p(i,j))
     >      *0.5*rdx*(zbot_p(i+1,j)-zbot_p(i-1,j))
c
        dqdx_xy(i,j,2) = ztop/(ztop-zbot_u(i,j))
     >      *rdx*(zbot_p(i+1,j)-zbot_p(i,j))
c			 
      enddo
      enddo 

      do j=1,ny
      do i=2,nx
        dqdx_xy(i,j,3) = ztop/(ztop-zbot_v(i,j))
     >            * 0.5*rdx*(zbot_v(i+1,j)-zbot_v(i-1,j))
      enddo
      enddo
*
* x and y dependent parts of dqdy.  Boundary, then interior.  We use 
*   all the terrain pts here (i.e., u, v, and therm).  Note that
*   dqdy_xy(i,j,2) is normally not used during model integrations.
*
      do i=1,nx+1
        dqdy_xy(i,1,1) = ztop/(ztop-zbot_p(i,1))
     >                * 2.*rdy*(zbot_v(i,1)-zbot_p(i,1))
        dqdy_xy(i,1,3) = ztop/(ztop-zbot_v(i,1))
     >                * rdy*(zbot_p(i,2)-zbot_p(i,1))
        dqdy_xy(i,ny,3) = ztop/(ztop-zbot_v(i,ny))
     >                * rdy*(zbot_p(i,ny+1)-zbot_p(i,ny))
        dqdy_xy(i,ny+1,1) = ztop/(ztop-zbot_p(i,ny+1)) 
     >                * 2.*rdy*(zbot_v(i,ny)-zbot_p(i,ny+1))
        dqdy_xy(i,ny,1) = ztop/(ztop-zbot_p(i,ny))
     >      *( 4./3.*rdy*(zbot_v(i,ny)-zbot_v(i,ny-1))
     >           - 1./6.*rdy*(zbot_p(i,ny+1)-zbot_p(i,ny-1)) )
      enddo

      do i=1,nx
        dqdy_xy(i,1,2) = ztop/(ztop-zbot_u(i,1))
     >                * rdy*(zbot_u(i,2)-zbot_u(i,1))
        dqdy_xy(i,ny+1,2) = ztop/(ztop-zbot_u(i,ny))
     >                * rdy*(zbot_u(i,ny)-zbot_u(i,ny+1))
      enddo

      do j=2,ny-1
      do i=1,nx+1
c
        dqdy_xy(i,j,1) = ztop/(ztop-zbot_p(i,j))
     >      *0.5*rdy*(zbot_p(i,j+1)-zbot_p(i,j-1))
c
        dqdy_xy(i,j,3) = ztop/(ztop-zbot_v(i,j))
     >      *rdy*(zbot_p(i,j+1)-zbot_p(i,j))
c
      enddo
      enddo

      do j=2,ny
      do i=1,nx
        dqdy_xy(i,j,2) = ztop/(ztop-zbot_u(i,j))
     >            * 0.5*rdy*(zbot_u(i,j+1)-zbot_u(i,j-1))
      enddo
      enddo
*
* z dependent parts of dqdx and dqdy (therm heights first, then w heights)
*
      do k=1,nz
        dqdx_zt(k) = (sgz(k)-ztop)/ztop
        dqdy_zt(k) = (sgz(k)-ztop)/ztop
      enddo

      do k=1,nz+1
        dqdx_zw(k) = (wgz(k)-ztop)/ztop
        dqdy_zw(k) = (wgz(k)-ztop)/ztop
      enddo

      return
      end

*---------------------------------------------------------------------------

      subroutine bnd_arrays(sgz,wgz,nx,ny,nz)

      use misc_module

      integer nx,ny,nz
      real sgz(nz),wgz(nz+1)

      integer i,j,k

      do i=2,nx
        xlbnd(i) = 1.
        xrbnd(i) = 1.
      enddo
      
      xlbnd(1) = 0.
      xlbnd(nx+1) = 2.
      xrbnd(1) = 2.
      xrbnd(nx+1) = 0.

      do j=2,ny
        ylbnd(j) = 1.
        yrbnd(j) = 1.
      enddo

      ylbnd(1) = 0.
      ylbnd(ny+1) = 2.
      yrbnd(1) = 2.
      yrbnd(ny+1) = 0.

      do k=2,nz-1
        zlbnd(k) = 1.
        zrbnd(k) = 1.
      enddo

      zlbnd(1) = 0.
      zlbnd(nz) = 2.
      zrbnd(1) = 2.
      zrbnd(nz) = 0.
*
* vertical grid spacing
*
      do k=1,nz
        rdzc(k) = 1./(wgz(k+1)-wgz(k))
      enddo

      do k=2,nz
        rdze(k) = 1./(sgz(k)-sgz(k-1))
      enddo
      rdze(1) = 1.
      rdze(nz+1) = 1.

      return
      end
