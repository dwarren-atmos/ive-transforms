
#include "pointer.icl"

      ive_ptr function tangent(basenm,vecx,vecy,vecz,dir,ndims,dims,
     >                  stag,min,max,miss,dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),basenm,vecx,vecy,vecz
      character*(1) dir

      integer strbeg,strend,nt,ijk,idir,iu,iv,iw
      logical new,err,requ
      ive_ptr getmem,base,vx,vy,vz

      tangent = 0

      nt = ntime

      if (dir.eq.'X') then
        idir = 1
      else if (dir.eq.'Y') then
        idir = 2
      else if (dir.eq.'Z') then
        idir = 3
      else
        write(6,*) 'Error:  do not know how to compute ',dir,
     >                     ' component of tangent.'
        return
      end if

      base = getvar(basenm,ndims,dims,stag,min,max,miss,dunits,
     >                disunits,dim_names,new)
      if (base.eq.0) then
        write(6,*) 'Error:  could not access variable ',
     >          basenm(strbeg(basenm):strend(basenm))
        return
      end if

      if (requ(min(1),plmin(1)).or.requ(min(3),plmin(3)).or.
     >        (requ(min(2),plmin(2)).and.(twod.ne.1))) then
        write(6,*) 'Error:  base field must be at thermo pts.'
        return
      end if

      tangent = getmem(nx*ny*nz*nt)
      if (tangent.eq.0) then
        write(6,*) 'Error:  could not allocate space for tangent.'
        return
      end if

      do ijk=1,nx*ny*nz*nt
        call s_val(%val(tangent),ijk,0.)
      enddo

      vx = getvar(vecx,ndims,dims,stag,min,max,miss,dunits,disunits,
     >                               dim_names,new)
      if (vx.eq.0) then
        write(6,*) 'Error:  could not access ', 
     >                  vecx(strbeg(vecx):strend(vecx))
        call freemem(tangent)
        return
      end if

      iu = 0
      iv = 0
      iw = 0
      if (requ(min(1),plmin(1))) iu = 1
      if (requ(min(2),plmin(2)).and.(twod.ne.1)) iv = 1
      if (requ(min(3),plmin(3))) iw = 1

      miss = -777.

      call tan_calc_I(%val(tangent),%val(base),%val(vx),%val(zbot_p),
     >             %val(zbot_u),%val(zbot_v),%val(sgz),ztop,dx,dy,
     >             miss,nx,ny,nz,nt,idir,iu,iv,iw)

      vy = getvar(vecy,ndims,dims,stag,min,max,miss,dunits,disunits,
     >                               dim_names,new)
      if (vy.eq.0) then
        write(6,*) 'Error:  could not access ',
     >                  vecy(strbeg(vecy):strend(vecy))
        call freemem(tangent)
        return
      end if

      iu = 0
      iv = 0
      iw = 0
      if (requ(min(1),plmin(1))) iu = 1
      if (requ(min(2),plmin(2)).and.(twod.ne.1)) iv = 1
      if (requ(min(3),plmin(3))) iw = 1

      miss = -777.

      call tan_calc_II(%val(tangent),%val(base),%val(vy),%val(zbot_p),
     >             %val(zbot_u),%val(zbot_v),%val(sgz),ztop,dx,dy,
     >             miss,nx,ny,nz,nt,idir,iu,iv,iw)

      vz = getvar(vecz,ndims,dims,stag,min,max,miss,dunits,disunits,
     >                               dim_names,new)
      if (vz.eq.0) then
        write(6,*) 'Error:  could not access ',
     >                  vecz(strbeg(vecz):strend(vecz))
        call freemem(tangent)
        return
      end if

      iu = 0
      iv = 0
      iw = 0
      if (requ(min(1),plmin(1))) iu = 1
      if (requ(min(2),plmin(2)).and.(twod.ne.1)) iv = 1
      if (requ(min(3),plmin(3))) iw = 1

      miss = -777.

      call tan_calc_III(%val(tangent),%val(base),%val(vz),%val(zbot_p),
     >               %val(zbot_u),%val(zbot_v),%val(sgz),ztop,dx,dy,
     >               miss,nx,ny,nz,nt,idir,iu,iv,iw)

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

      return
      end

*----------------------------------------------------------------------

      subroutine tan_calc_I(tan,base,vx,zs_p,zs_u,zs_v,sgz,ztop,dx,dy,
     >                            miss,nx,ny,nz,nt,idir,iu,iv,iw)

      integer nx,ny,nz,nt,idir,iu,iv,iw
      real tan(nx,ny,nz,nt),base(nx,ny,nz,nt)
      real vx(nx+iu,ny+iv,nz+iw,nt),zs_p(nx,ny),zs_u(nx+1,ny)
      real zs_v(nx,ny+1),sgz(nz),ztop,dx,dy,miss

      integer i,j,k,n,ip1,im1,jp1,jm1,kp1,km1
      real dbdzet,dzetdx,dzetdy,dzetdz,gdg,vxave

      integer nxmax,nymax,nzmax
      parameter (nxmax=401,nymax=401,nzmax=201)
      real xbnd(nxmax),ybnd(nymax),rdx,rdy
      real zlbnd(nzmax),zrbnd(nzmax),rdz(nzmax)

      if (nx.gt.nxmax) then
        write(6,*) 'Error:  nx too big for tangent calculation.'
        stop
      else if (ny.gt.nymax) then
        write(6,*) 'Error:  ny too big for tangent calculation.'
        stop
      else if (nz.gt.nzmax) then
        write(6,*) 'Error:  nz too big for tangent calculation.'
        stop
      end if

      do i=2,nx-1
        xbnd(i) = 0.5
      enddo
      xbnd(1) = 1.
      xbnd(nx) = 1.  

      do j=2,ny-1
        ybnd(j) = 0.5
      enddo
      ybnd(1) = 1.
      ybnd(ny) = 1.        

      do k=1,nz-1
        rdz(k) = 1./(sgz(k+1)-sgz(k))
        zlbnd(k) = 1.
        zrbnd(k) = 1.
      enddo
      rdz(nz) = 1.
      zlbnd(nz) = 2.
      zlbnd(1) = 0.
      zrbnd(nz) = 0.
      zrbnd(1) = 2.

      rdx = 1./dx
      rdy = 1./dy

      do n=1,nt
      do k=1,nz
        kp1 = min(k+1,nz)
        km1 = max(k-1,1)
      do j=1,ny
        jp1 = min(j+1,ny)
        jm1 = max(j-1,1)
      do i=1,nx
        ip1 = min(i+1,nx)
        im1 = max(i-1,1)
 
        dbdzet = 0.5*zrbnd(k)*rdz(k)
     >                  *(base(i,j,kp1,n)-base(i,j,k,n))
     >         + 0.5*zlbnd(k)*rdz(km1)
     >                  *(base(i,j,k,n)-base(i,j,km1,n))
        dzetdx = (sgz(k)-ztop)/(ztop-zs_p(i,j))
     >                 *rdx*(zs_u(i+1,j)-zs_u(i,j))
        dzetdy = (sgz(k)-ztop)/(ztop-zs_p(i,j))
     >                 *rdy*(zs_v(i,j+1)-zs_v(i,j))
        dzetdz = ztop/(ztop-zs_p(i,j))

        gdg = (xbnd(i)*rdx*(base(ip1,j,k,n)-base(im1,j,k,n)) 
     >                     + dbdzet*dzetdx)**2
     >      + (ybnd(j)*rdy*(base(i,jp1,k,n)-base(i,jm1,k,n))
     >                     + dbdzet*dzetdy)**2
     >      + (dbdzet*dzetdz)**2

        if (iu+iv+iw.le.1) then
          vxave = 0.5*(vx(i,j,k,n)+vx(i+iu,j+iv,k+iw,n))
        else if ((iv.eq.1).and.(iw.eq.1)) then
          vxave = 0.25*(vx(i,j,k,n)+vx(i,j+1,k,n)
     >                   + vx(i,j+1,k+1,n)+vx(i,j,k+1,n))
        else
          write(6,*) 'Error:  case not recognized in tan_calc_I.'
          write(6,*) 'Stopping.'
          stop
        end if

        if (idir.eq.1) then

          tan(i,j,k,n) = tan(i,j,k,n) + vxave
     >                       - vxave
     >        * (xbnd(i)*rdx*(base(ip1,j,k,n)-base(im1,j,k,n))
     >                + dbdzet*dzetdx)**2 / gdg

        else if (idir.eq.2) then

          tan(i,j,k,n) = tan(i,j,k,n)
     >                       - vxave 
     >        * (xbnd(i)*rdx*(base(ip1,j,k,n)-base(im1,j,k,n))
     >                     + dbdzet*dzetdx)
     >        * (ybnd(j)*rdy*(base(i,jp1,k,n)-base(i,jm1,k,n))
     >                     + dbdzet*dzetdy) / gdg

        else 

          tan(i,j,k,n) = tan(i,j,k,n)
     >                       - vxave 
     >        * (xbnd(i)*rdx*(base(ip1,j,k,n)-base(im1,j,k,n))
     >                     + dbdzet*dzetdx)
     >                * dbdzet*dzetdz / gdg
  
        end if

      enddo
      enddo
      enddo
      enddo

      return
      end

*-----------------------------------------------------------------------

      subroutine tan_calc_II(tan,base,vy,zs_p,zs_u,zs_v,sgz,ztop,
     >                   dx,dy,miss,nx,ny,nz,nt,idir,iu,iv,iw)

      integer nx,ny,nz,nt,idir,iu,iv,iw
      real tan(nx,ny,nz,nt),base(nx,ny,nz,nt)
      real vy(nx+iu,ny+iv,nz+iw,nt),zs_p(nx,ny),zs_u(nx+1,ny)
      real zs_v(nx,ny+1),sgz(nz),ztop,dx,dy,miss
      
      integer i,j,k,n,ip1,im1,jp1,jm1,kp1,km1
      real dbdzet,dzetdx,dzetdy,dzetdz,gdg,vyave

      integer nxmax,nymax,nzmax
      parameter (nxmax=401,nymax=401,nzmax=201)
      real xbnd(nxmax),ybnd(nymax),rdx,rdy
      real zlbnd(nzmax),zrbnd(nzmax),rdz(nzmax)

      if (nx.gt.nxmax) then
        write(6,*) 'Error:  nx too big for tangent calculation.'
        stop
      else if (ny.gt.nymax) then
        write(6,*) 'Error:  ny too big for tangent calculation.'
        stop
      else if (nz.gt.nzmax) then
        write(6,*) 'Error:  nz too big for tangent calculation.'
        stop
      end if

      do i=2,nx-1
        xbnd(i) = 0.5
      enddo
      xbnd(1) = 1.
      xbnd(nx) = 1.  

      do j=2,ny-1
        ybnd(j) = 0.5
      enddo
      ybnd(1) = 1.
      ybnd(ny) = 1.        

      do k=1,nz-1
        rdz(k) = 1./(sgz(k+1)-sgz(k))
        zlbnd(k) = 1.
        zrbnd(k) = 1.
      enddo
      rdz(nz) = 1.
      zlbnd(nz) = 2.
      zlbnd(1) = 0.
      zrbnd(nz) = 0.
      zrbnd(1) = 2.

      rdx = 1./dx
      rdy = 1./dy

      do n=1,nt
      do k=1,nz
        kp1 = min(k+1,nz)
        km1 = max(k-1,1)
      do j=1,ny
        jp1 = min(j+1,ny)
        jm1 = max(j-1,1)
      do i=1,nx
        ip1 = min(i+1,nx)
        im1 = max(i-1,1)

        dbdzet = 0.5*zrbnd(k)*rdz(k)
     >                  *(base(i,j,kp1,n)-base(i,j,k,n))
     >         + 0.5*zlbnd(k)*rdz(km1)
     >                  *(base(i,j,k,n)-base(i,j,km1,n))
        dzetdx = (sgz(k)-ztop)/(ztop-zs_p(i,j))
     >                 *rdx*(zs_u(i+1,j)-zs_u(i,j))
        dzetdy = (sgz(k)-ztop)/(ztop-zs_p(i,j))
     >                 *rdy*(zs_v(i,j+1)-zs_v(i,j))
        dzetdz = ztop/(ztop-zs_p(i,j))

        gdg = (xbnd(i)*rdx*(base(ip1,j,k,n)-base(im1,j,k,n)) 
     >                     + dbdzet*dzetdx)**2
     >      + (ybnd(j)*rdy*(base(i,jp1,k,n)-base(i,jm1,k,n))
     >                     + dbdzet*dzetdy)**2
     >      + (dbdzet*dzetdz)**2

        if (iu+iv+iw.le.1) then
          vyave = 0.5*(vy(i,j,k,n)+vy(i+iu,j+iv,k+iw,n))
        else if ((iu.eq.1).and.(iw.eq.1)) then
          vyave = 0.25*(vy(i,j,k,n)+vy(i+1,j,k,n)
     >                    + vy(i+1,j,k+1,n)+vy(i,j,k+1,n))
        else
          write(6,*) 'Error:  case not recognized in tan_calc_II.'
          write(6,*) 'Stopping.'
          stop
        end if

        if (idir.eq.1) then

          tan(i,j,k,n) = tan(i,j,k,n)
     >                       - vyave 
     >        * (ybnd(j)*rdy*(base(i,jp1,k,n)-base(i,jm1,k,n))
     >                       + dbdzet*dzetdy)
     >        * (xbnd(i)*rdx*(base(ip1,j,k,n)-base(im1,j,k,n))
     >                       + dbdzet*dzetdx) / gdg

        else if (idir.eq.2) then

          tan(i,j,k,n) = tan(i,j,k,n) + vyave
     >                       - vyave
     >        * (ybnd(j)*rdy*(base(i,jp1,k,n)-base(i,jm1,k,n))
     >                + dbdzet*dzetdy)**2 / gdg

        else 

          tan(i,j,k,n) = tan(i,j,k,n)
     >                       - vyave 
     >        * (ybnd(j)*rdy*(base(i,jp1,k,n)-base(i,jm1,k,n))
     >                       + dbdzet*dzetdy)
     >                * dbdzet*dzetdz / gdg

        end if

      enddo
      enddo
      enddo
      enddo

      return
      end    

*------------------------------------------------------------------------

      subroutine tan_calc_III(tan,base,vz,zs_p,zs_u,zs_v,sgz,ztop,
     >                    dx,dy,miss,nx,ny,nz,nt,idir,iu,iv,iw)

      integer nx,ny,nz,nt,idir,iu,iv,iw
      real tan(nx,ny,nz,nt),base(nx,ny,nz,nt)
      real vz(nx+iu,ny+iv,nz+iw,nt),zs_p(nx,ny),zs_u(nx+1,ny)
      real zs_v(nx,ny+1),sgz(nz),ztop,dx,dy,miss

      integer i,j,k,n,ip1,im1,jp1,jm1,kp1,km1
      real dbdzet,dzetdx,dzetdy,dzetdz,gdg,vzave

      integer nxmax,nymax,nzmax
      parameter (nxmax=401,nymax=401,nzmax=201)
      real xbnd(nxmax),ybnd(nymax),rdx,rdy
      real zlbnd(nzmax),zrbnd(nzmax),rdz(nzmax)

      if (nx.gt.nxmax) then
        write(6,*) 'Error:  nx too big for tangent calculation.'
        stop
      else if (ny.gt.nymax) then
        write(6,*) 'Error:  ny too big for tangent calculation.'
        stop
      else if (nz.gt.nzmax) then
        write(6,*) 'Error:  nz too big for tangent calculation.'
        stop
      end if

      do i=2,nx-1
        xbnd(i) = 0.5
      enddo
      xbnd(1) = 1.
      xbnd(nx) = 1.  

      do j=2,ny-1
        ybnd(j) = 0.5
      enddo
      ybnd(1) = 1.
      ybnd(ny) = 1.        

      do k=1,nz-1
        rdz(k) = 1./(sgz(k+1)-sgz(k))
        zlbnd(k) = 1.
        zrbnd(k) = 1.
      enddo
      rdz(nz) = 1.
      zlbnd(nz) = 2.
      zlbnd(1) = 0.
      zrbnd(nz) = 0.
      zrbnd(1) = 2.

      rdx = 1./dx
      rdy = 1./dy

      do n=1,nt
      do k=1,nz
        kp1 = min(k+1,nz)
        km1 = max(k-1,1)
      do j=1,ny
        jp1 = min(j+1,ny)
        jm1 = max(j-1,1)
      do i=1,nx
        ip1 = min(i+1,nx)
        im1 = max(i-1,1)

        dbdzet = 0.5*zrbnd(k)*rdz(k)
     >                  *(base(i,j,kp1,n)-base(i,j,k,n))
     >         + 0.5*zlbnd(k)*rdz(km1)
     >                  *(base(i,j,k,n)-base(i,j,km1,n))
        dzetdx = (sgz(k)-ztop)/(ztop-zs_p(i,j))
     >                 *rdx*(zs_u(i+1,j)-zs_u(i,j))
        dzetdy = (sgz(k)-ztop)/(ztop-zs_p(i,j))
     >                 *rdy*(zs_v(i,j+1)-zs_v(i,j))
        dzetdz = ztop/(ztop-zs_p(i,j))

        gdg = (xbnd(i)*rdx*(base(ip1,j,k,n)-base(im1,j,k,n)) 
     >                     + dbdzet*dzetdx)**2
     >      + (ybnd(j)*rdy*(base(i,jp1,k,n)-base(i,jm1,k,n))
     >                     + dbdzet*dzetdy)**2
     >      + (dbdzet*dzetdz)**2

        if (iu+iv+iw.le.1) then
          vzave = 0.5*(vz(i,j,k,n)+vz(i+iu,j+iv,k+iw,n))
        else if ((iu.eq.1).and.(iv.eq.1)) then
          vzave = 0.25*(vz(i,j,k,n)+vz(i+1,j,k,n)
     >                    + vz(i+1,j+1,k,n)+vz(i,j+1,k,n))
        else
          write(6,*) 'Error:  case not recognized in tan_calc_III.'
          write(6,*) 'Stopping.'
          stop
        end if

        if (idir.eq.1) then

          tan(i,j,k,n) = tan(i,j,k,n)
     >                       - vzave 
     >                    * dbdzet*dzetdz
     >        * (xbnd(i)*rdx*(base(ip1,j,k,n)-base(im1,j,k,n))
     >                       + dbdzet*dzetdx) / gdg

        else if (idir.eq.2) then

          tan(i,j,k,n) = tan(i,j,k,n)
     >                       - vzave
     >                    * dbdzet*dzetdz
     >        * (ybnd(j)*rdy*(base(i,jp1,k,n)-base(i,jm1,k,n))
     >                       + dbdzet*dzetdy) / gdg

        else

          tan(i,j,k,n) = tan(i,j,k,n) + vzave
     >                       - vzave
     >                    * (dbdzet*dzetdz)**2 / gdg

        end if

      enddo
      enddo
      enddo
      enddo

      return
      end
