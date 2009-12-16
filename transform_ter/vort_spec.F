
#include "pointer.icl"

      ive_ptr function s_ive(ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer k,nt
      logical err,new
      ive_ptr getmem,u,v

      nt = ntime
      if(nobottom) then
         s_ive = 0
         return
      endif
      s_ive = getmem(nx*ny*nz*nt)
      if (s_ive.eq.0) then
        write(6,*) 'Error:  could not allocate space for s_ive.'
        return
      end if

      do k=1,nx*ny*nz*nt
        call s_val(%val(s_ive),k,0.)
      enddo

      u = getvar('u',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                       dim_names,new)
      if (u.eq.0) then
        write(6,*) 'Error:  unable to access u.'
        call freemem(s_ive)
        return
      end if

      v = getvar('v',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                       dim_names,new)
      if (v.eq.0) then
        write(6,*) 'Error:  unable to access v.'
        call freemem(s_ive)
        return
      end if

      call s_ive_calc(%val(s_ive),%val(u),%val(v),%val(dqdx_xy),
     >               %val(dqdx_zt),%val(dqdx_zw),%val(dqdy_xy),
     >               %val(dqdy_zt),%val(dqdy_zw),%val(dqdz),
     >               %val(dzdq),dx,dy,nx,ny,nz,nt)

      dims(1) = nx
      dims(2) = ny
      dims(3) = nz
      stag(1) = 0.5
      stag(2) = 0.5
      stag(3) = 0.5
      min(1) = plmin(1) + 0.5*dx
      min(2) = plmin(2) + 0.5*dy
      min(3) = sgz_min 
      max(1) = plmax(1) - 0.5*dx
      max(2) = plmax(2) - 0.5*dy 
      max(3) = sgz_max
      dim_names(1) = 'nx'
      dim_names(2) = 'ny'
      dim_names(3) = 'nz'
      dunits = '1/s'
      disunits = '1/s'

      return
      end

*-----------------------------------------------------------------------

      subroutine s_ive_calc(s_ive,u,v,dqdx_xy,dqdx_zt,dqdx_zw,
     >                      dqdy_xy,dqdy_zt,dqdy_zw,dqdz,dzdq,
     >                      dx,dy,nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real s_ive(nx,ny,nz,nt),u(nx+1,ny,nz,nt),v(nx,ny+1,nz,nt)
      real dqdx_xy(nx+1,ny+1,3),dqdx_zt(nz),dqdx_zw(nz+1)
      real dqdy_xy(nx+1,ny+1,3),dqdy_zt(nz),dqdy_zw(nz+1)
      real dqdz(nx+1,ny+1,3),dzdq(nx+1,ny+1,3)
      real dx,dy

      integer i,j,k,n,m,ip1,im1,jp1,jm1,kp1,km1
      real rdx,rdy

      rdx = 1./dx
      rdy = 1./dy

      do n=1,nt
* 
* north and south boundaries
*
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do m=1,4
        if (m.eq.1) j = 1
        if (m.eq.2) j = 2
        if (m.eq.3) j = ny-1
        if (m.eq.4) j = ny
        jp1 = min0(j+1,ny)
        jm1 = max0(j-1,1)
      do i=1,nx
        ip1 = min0(i+1,nx)
        im1 = max0(i-1,1)

        s_ive(i,j,k,n) =  - 0.25*rdy
     >       *( yrbnd(j)*(u(i,jp1,k,n)-u(i,j,k,n)
     >                        + u(i+1,jp1,k,n)-u(i+1,j,k,n))
     >        + ylbnd(j)*(u(i,j,k,n)-u(i,jm1,k,n)
     >                        + u(i+1,j,k,n)-u(i+1,jm1,k,n)) )
     >                - 0.25*dqdy_xy(i,j,1)
     >       *( rdze(k+1)*dqdy_zw(k+1)*zrbnd(k)*(u(i,j,kp1,n)
     >             - u(i,j,k,n)+u(i+1,j,kp1,n)-u(i+1,j,k,n))
     >          + rdze(k)*dqdy_zw(k)*zlbnd(k)*(u(i,j,k,n)
     >             - u(i,j,km1,n)+u(i+1,j,k,n)-u(i+1,j,km1,n)) ) 
     >                     + 0.25*rdx
     >       *( xrbnd(i)*(v(ip1,j,k,n)-v(i,j,k,n)
     >                        + v(ip1,j+1,k,n)-v(i,j+1,k,n))
     >        + xlbnd(i)*(v(i,j,k,n)-v(im1,j,k,n)
     >                        + v(i,j+1,k,n)-v(im1,j+1,k,n)) )
     >                + 0.25*dqdx_xy(i,j,1)
     >       *( rdze(k+1)*dqdx_zw(k+1)*zrbnd(k)*(v(i,j,kp1,n)
     >             - v(i,j,k,n)+v(i,j+1,kp1,n)-v(i,j+1,k,n))
     >          + rdze(k)*dqdx_zw(k)*zlbnd(k)*(v(i,j,k,n)
     >             - v(i,j,km1,n)+v(i,j+1,k,n)-v(i,j+1,km1,n)) )
 
      enddo
      enddo
      enddo
* 
* east and west boundaries  
*
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do m=1,4
        if (m.eq.1) i = 1
        if (m.eq.2) i = 2
        if (m.eq.3) i = nx-1
        if (m.eq.4) i = nx
        ip1 = min0(i+1,nx)
        im1 = max0(i-1,1)
      do j=3,ny-2
        jp1 = min0(j+1,ny)
        jm1 = max0(j-1,1)

        s_ive(i,j,k,n) =  - 0.25*rdy
     >       *( yrbnd(j)*(u(i,jp1,k,n)-u(i,j,k,n)
     >                        + u(i+1,jp1,k,n)-u(i+1,j,k,n))
     >        + ylbnd(j)*(u(i,j,k,n)-u(i,jm1,k,n)
     >                        + u(i+1,j,k,n)-u(i+1,jm1,k,n)) )
     >                - 0.25*dqdy_xy(i,j,1)
     >       *( rdze(k+1)*dqdy_zw(k+1)*zrbnd(k)*(u(i,j,kp1,n)
     >             - u(i,j,k,n)+u(i+1,j,kp1,n)-u(i+1,j,k,n))
     >          + rdze(k)*dqdy_zw(k)*zlbnd(k)*(u(i,j,k,n)
     >             - u(i,j,km1,n)+u(i+1,j,k,n)-u(i+1,j,km1,n)) ) 
     >                     + 0.25*rdx
     >       *( xrbnd(i)*(v(ip1,j,k,n)-v(i,j,k,n)
     >                        + v(ip1,j+1,k,n)-v(i,j+1,k,n))
     >        + xlbnd(i)*(v(i,j,k,n)-v(im1,j,k,n)
     >                        + v(i,j+1,k,n)-v(im1,j+1,k,n)) )
     >                + 0.25*dqdx_xy(i,j,1)
     >       *( rdze(k+1)*dqdx_zw(k+1)*zrbnd(k)*(v(i,j,kp1,n)
     >             - v(i,j,k,n)+v(i,j+1,kp1,n)-v(i,j+1,k,n))
     >          + rdze(k)*dqdx_zw(k)*zlbnd(k)*(v(i,j,k,n)
     >             - v(i,j,km1,n)+v(i,j+1,k,n)-v(i,j+1,km1,n)) )

      enddo
      enddo
      enddo
* 
* interior
*
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do j=3,ny-2
      do i=3,nx-2

        s_ive(i,j,k,n) = 
     >            - 1./3.*rdy*(u(i,j+1,k,n)-u(i,j-1,k,n)
     >                      + u(i+1,j+1,k,n)-u(i+1,j-1,k,n))
     >            + 1./24.*rdy*(u(i,j+2,k,n)-u(i,j-2,k,n)
     >                      + u(i+1,j+2,k,n)-u(i+1,j-2,k,n))
     >                - 0.25*dqdy_xy(i,j,1)
     >       *( rdze(k+1)*dqdy_zw(k+1)*zrbnd(k)*(u(i,j,kp1,n)
     >             - u(i,j,k,n)+u(i+1,j,kp1,n)-u(i+1,j,k,n))
     >          + rdze(k)*dqdy_zw(k)*zlbnd(k)*(u(i,j,k,n)
     >             - u(i,j,km1,n)+u(i+1,j,k,n)-u(i+1,j,km1,n)) ) 
     >            + 1./3.*rdx*(v(i+1,j,k,n)-v(i-1,j,k,n)
     >                      + v(i+1,j+1,k,n)-v(i-1,j+1,k,n))
     >            - 1./24.*rdx*(v(i+2,j,k,n)-v(i-2,j,k,n)
     >                      + v(i+2,j+1,k,n)-v(i-2,j+1,k,n)) 
     >                + 0.25*dqdx_xy(i,j,1)
     >       *( rdze(k+1)*dqdx_zw(k+1)*zrbnd(k)*(v(i,j,kp1,n)
     >             - v(i,j,k,n)+v(i,j+1,kp1,n)-v(i,j+1,k,n))
     >          + rdze(k)*dqdx_zw(k)*zlbnd(k)*(v(i,j,k,n)
     >             - v(i,j,km1,n)+v(i,j+1,k,n)-v(i,j+1,km1,n)) )

      enddo
      enddo
      enddo
*
* end outer time loop
*
      enddo

      return
      end

*----------------------------------------------------------------------

      ive_ptr function r_ive(ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer k,nt
      logical err,new
      ive_ptr getmem,u,w

      nt = ntime
      if(nobottom) then
         r_ive = 0
         return
      endif

      r_ive = getmem(nx*ny*nz*nt)
      if (r_ive.eq.0) then
        write(6,*) 'Error:  could not allocate space for r_ive.'
        return
      end if

      do k=1,nx*ny*nz*nt
        call s_val(%val(r_ive),k,0.)
      enddo

      u = getvar('u',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                       dim_names,new)
      if (u.eq.0) then
        write(6,*) 'Error:  unable to access u.'
        call freemem(r_ive)
        return
      end if

      w = getvar('w',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                       dim_names,new)
      if (w.eq.0) then
        write(6,*) 'Error:  unable to access w.'
        call freemem(r_ive)
        return
      end if

      call r_ive_calc(%val(r_ive),%val(u),%val(w),%val(dqdx_xy),
     >                %val(dqdx_zt),%val(dqdx_zw),%val(dqdy_xy),
     >                %val(dqdy_zt),%val(dqdy_zw),%val(dqdz),
     >                %val(dzdq),dx,dy,nx,ny,nz,nt)

      dims(1) = nx
      dims(2) = ny
      dims(3) = nz
      stag(1) = 0.5
      stag(2) = 0.5
      stag(3) = 0.5
      min(1) = plmin(1) + 0.5*dx
      min(2) = plmin(2) + 0.5*dy
      min(3) = sgz_min 
      max(1) = plmax(1) - 0.5*dx
      max(2) = plmax(2) - 0.5*dy 
      max(3) = sgz_max
      dim_names(1) = 'nx'
      dim_names(2) = 'ny'
      dim_names(3) = 'nz'
      dunits = '1/s'
      disunits = '1/s'

      return
      end

*-----------------------------------------------------------------------

      subroutine r_ive_calc(r_ive,u,w,dqdx_xy,dqdx_zt,dqdx_zw,
     >                      dqdy_xy,dqdy_zt,dqdy_zw,dqdz,dzdq,
     >                      dx,dy,nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real r_ive(nx,ny,nz,nt),u(nx+1,ny,nz,nt),w(nx,ny,nz+1,nt)
      real dqdx_xy(nx+1,ny+1,3),dqdx_zt(nz),dqdx_zw(nz+1)
      real dqdy_xy(nx+1,ny+1,3),dqdy_zt(nz),dqdy_zw(nz+1)
      real dqdz(nx+1,ny+1,3),dzdq(nx+1,ny+1,3)
      real dx,dy

      integer i,j,k,n,m,ip1,im1,jp1,jm1,kp1,km1
      real rdx,rdy

      rdx = 1./dx
      rdy = 1./dy

      do n=1,nt
* 
* north and south boundaries
*
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do m=1,4
        if (m.eq.1) j = 1
        if (m.eq.2) j = 2
        if (m.eq.3) j = ny-1
        if (m.eq.4) j = ny
        jp1 = min0(j+1,ny)
        jm1 = max0(j-1,1)
      do i=1,nx
        ip1 = min0(i+1,nx)
        im1 = max0(i-1,1)

        r_ive(i,j,k,n) =  - 0.25*rdx
     >       *( xrbnd(i)*(w(ip1,j,k,n)-w(i,j,k,n)
     >                        + w(ip1,j,k+1,n)-w(i,j,k+1,n))
     >        + xlbnd(i)*(w(i,j,k,n)-w(im1,j,k,n)
     >                        + w(i,j,k+1,n)-w(im1,j,k+1,n)) )
     >                - dqdx_xy(i,j,1)*dqdx_zt(k)*rdzc(k)
     >                          * (w(i,j,k+1,n)-w(i,j,k,n))
     >                + 0.25*dqdz(i,j,1)
     >       *( rdze(k+1)*zrbnd(k)*(u(i,j,kp1,n)-u(i,j,k,n)
     >              + u(i+1,j,kp1,n)-u(i+1,j,k,n))
     >          + rdze(k)*zlbnd(k)*(u(i,j,k,n)-u(i,j,km1,n)
     >              + u(i+1,j,k,n)-u(i+1,j,km1,n)) )

 
      enddo
      enddo
      enddo
* 
* east and west boundaries  
*
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do m=1,4
        if (m.eq.1) i = 1
        if (m.eq.2) i = 2
        if (m.eq.3) i = nx-1
        if (m.eq.4) i = nx
        ip1 = min0(i+1,nx)
        im1 = max0(i-1,1)
      do j=3,ny-2
        jp1 = min0(j+1,ny)
        jm1 = max0(j-1,1)

        r_ive(i,j,k,n) =  - 0.25*rdx
     >       *( xrbnd(i)*(w(ip1,j,k,n)-w(i,j,k,n)
     >                        + w(ip1,j,k+1,n)-w(i,j,k+1,n))
     >        + xlbnd(i)*(w(i,j,k,n)-w(im1,j,k,n)
     >                        + w(i,j,k+1,n)-w(im1,j,k+1,n)) )
     >                - dqdx_xy(i,j,1)*dqdx_zt(k)*rdzc(k)
     >                          * (w(i,j,k+1,n)-w(i,j,k,n))
     >                + 0.25*dqdz(i,j,1)
     >       *( rdze(k+1)*zrbnd(k)*(u(i,j,kp1,n)-u(i,j,k,n)
     >              + u(i+1,j,kp1,n)-u(i+1,j,k,n))
     >          + rdze(k)*zlbnd(k)*(u(i,j,k,n)-u(i,j,km1,n)
     >              + u(i+1,j,k,n)-u(i+1,j,km1,n)) )

      enddo
      enddo
      enddo
* 
* interior
*
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do j=3,ny-2
      do i=3,nx-2

          r_ive(i,j,k,n) =     
     >            - 1./3.*rdx*(w(i+1,j,k,n)-w(i-1,j,k,n)
     >                      + w(i+1,j,k+1,n)-w(i-1,j,k+1,n))
     >            + 1./24.*rdx*(w(i+2,j,k,n)-w(i-2,j,k,n)
     >                      + w(i+2,j,k+1,n)-w(i-2,j,k+1,n))
     >                - dqdx_xy(i,j,1)*dqdx_zt(k)*rdzc(k)
     >                        * (w(i,j,k+1,n)-w(i,j,k,n))
     >                + 0.25*dqdz(i,j,1)
     >       *( rdze(k+1)*zrbnd(k)*(u(i,j,kp1,n)-u(i,j,k,n)
     >              + u(i+1,j,kp1,n)-u(i+1,j,k,n))
     >          + rdze(k)*zlbnd(k)*(u(i,j,k,n)-u(i,j,km1,n)
     >              + u(i+1,j,k,n)-u(i+1,j,km1,n)) )

      enddo
      enddo
      enddo
*
* end outer time loop
*
      enddo

      return
      end

*---------------------------------------------------------------------------

      ive_ptr function q_ive(ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer k,nt
      logical err,new
      ive_ptr getmem,v,w

      nt = ntime

      if(nobottom) then
         q_ive = 0
         return
      endif
      q_ive = getmem(nx*ny*nz*nt)
      if (q_ive.eq.0) then
        write(6,*) 'Error:  could not allocate space for q_ive.'
        return
      end if

      do k=1,nx*ny*nz*nt
        call s_val(%val(q_ive),k,0.)
      enddo

      v = getvar('v',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                       dim_names,new)
      if (v.eq.0) then
        write(6,*) 'Error:  unable to access v.'
        call freemem(q_ive)
        return
      end if

      w = getvar('w',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                       dim_names,new)
      if (w.eq.0) then
        write(6,*) 'Error:  unable to access w.'
        call freemem(q_ive)
        return
      end if

      call q_ive_calc(%val(q_ive),%val(v),%val(w),%val(dqdx_xy),
     >                %val(dqdx_zt),%val(dqdx_zw),%val(dqdy_xy),
     >                %val(dqdy_zt),%val(dqdy_zw),%val(dqdz),
     >                %val(dzdq),dx,dy,nx,ny,nz,nt)

      dims(1) = nx
      dims(2) = ny
      dims(3) = nz
      stag(1) = 0.5
      stag(2) = 0.5
      stag(3) = 0.5
      min(1) = plmin(1) + 0.5*dx
      min(2) = plmin(2) + 0.5*dy
      min(3) = sgz_min 
      max(1) = plmax(1) - 0.5*dx
      max(2) = plmax(2) - 0.5*dy 
      max(3) = sgz_max
      dim_names(1) = 'nx'
      dim_names(2) = 'ny'
      dim_names(3) = 'nz'
      dunits = '1/s'
      disunits = '1/s'

      return
      end

*-----------------------------------------------------------------------

      subroutine q_ive_calc(q_ive,v,w,dqdx_xy,dqdx_zt,dqdx_zw,
     >                      dqdy_xy,dqdy_zt,dqdy_zw,dqdz,dzdq,
     >                      dx,dy,nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real q_ive(nx,ny,nz,nt),v(nx,ny+1,nz,nt),w(nx,ny,nz+1,nt)
      real dqdx_xy(nx+1,ny+1,3),dqdx_zt(nz),dqdx_zw(nz+1)
      real dqdy_xy(nx+1,ny+1,3),dqdy_zt(nz),dqdy_zw(nz+1)
      real dqdz(nx+1,ny+1,3),dzdq(nx+1,ny+1,3)
      real dx,dy

      integer i,j,k,n,m,ip1,im1,jp1,jm1,kp1,km1
      real rdx,rdy

      rdx = 1./dx
      rdy = 1./dy

      do n=1,nt
* 
* north and south boundaries
*
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do m=1,4
        if (m.eq.1) j = 1
        if (m.eq.2) j = 2
        if (m.eq.3) j = ny-1
        if (m.eq.4) j = ny
        jp1 = min0(j+1,ny)
        jm1 = max0(j-1,1)
      do i=1,nx
        ip1 = min0(i+1,nx)
        im1 = max0(i-1,1)

        q_ive(i,j,k,n) =  0.25*rdy
     >       *( yrbnd(j)*(w(i,jp1,k,n)-w(i,j,k,n)
     >                        + w(i,jp1,k+1,n)-w(i,j,k+1,n))
     >        + ylbnd(j)*(w(i,j,k,n)-w(i,jm1,k,n)
     >                        + w(i,j,k+1,n)-w(i,jm1,k+1,n)) )
     >                + dqdy_xy(i,j,1)*dqdy_zt(k)*rdzc(k)
     >                        * (w(i,j,k+1,n)-w(i,j,k,n))
     >                - 0.25*dqdz(i,j,1)
     >       *( rdze(k+1)*zrbnd(k)*(v(i,j,kp1,n)-v(i,j,k,n)
     >              + v(i,j+1,kp1,n)-v(i,j+1,k,n))
     >          + rdze(k)*zlbnd(k)*(v(i,j,k,n)-v(i,j,km1,n)
     >              + v(i,j+1,k,n)-v(i,j+1,km1,n)) )
 
      enddo
      enddo
      enddo
* 
* east and west boundaries  
*
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do m=1,4
        if (m.eq.1) i = 1
        if (m.eq.2) i = 2
        if (m.eq.3) i = nx-1
        if (m.eq.4) i = nx
        ip1 = min0(i+1,nx)
        im1 = max0(i-1,1)
      do j=3,ny-2
        jp1 = min0(j+1,ny)
        jm1 = max0(j-1,1)

        q_ive(i,j,k,n) =  0.25*rdy
     >       *( yrbnd(j)*(w(i,jp1,k,n)-w(i,j,k,n)
     >                        + w(i,jp1,k+1,n)-w(i,j,k+1,n))
     >        + ylbnd(j)*(w(i,j,k,n)-w(i,jm1,k,n)
     >                        + w(i,j,k+1,n)-w(i,jm1,k+1,n)) )
     >                + dqdy_xy(i,j,1)*dqdy_zt(k)*rdzc(k)
     >                        * (w(i,j,k+1,n)-w(i,j,k,n))
     >                - 0.25*dqdz(i,j,1)
     >       *( rdze(k+1)*zrbnd(k)*(v(i,j,kp1,n)-v(i,j,k,n)
     >              + v(i,j+1,kp1,n)-v(i,j+1,k,n))
     >          + rdze(k)*zlbnd(k)*(v(i,j,k,n)-v(i,j,km1,n)
     >              + v(i,j+1,k,n)-v(i,j+1,km1,n)) )

      enddo
      enddo
      enddo
* 
* interior
*
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do j=3,ny-2
      do i=3,nx-2

        q_ive(i,j,k,n) =      
     >            + 1./3.*rdy*(w(i,j+1,k,n)-w(i,j-1,k,n)
     >                      + w(i,j+1,k+1,n)-w(i,j-1,k+1,n))
     >            - 1./24.*rdy*(w(i,j+2,k,n)-w(i,j-2,k,n)
     >                      + w(i,j+2,k+1,n)-w(i,j-2,k+1,n))
     >                + dqdy_xy(i,j,1)*dqdy_zt(k)*rdzc(k)    
     >                        * (w(i,j,k+1,n)-w(i,j,k,n))
     >                - 0.25*dqdz(i,j,1)
     >       *( rdze(k+1)*zrbnd(k)*(v(i,j,kp1,n)-v(i,j,k,n)
     >              + v(i,j+1,kp1,n)-v(i,j+1,k,n))
     >          + rdze(k)*zlbnd(k)*(v(i,j,k,n)-v(i,j,km1,n)
     >              + v(i,j+1,k,n)-v(i,j+1,km1,n)) )

      enddo
      enddo
      enddo
*
* end outer time loop
*
      enddo

      return
      end

