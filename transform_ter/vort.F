

#include "pointer.icl"

      ive_ptr function pvs(ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer ijk,nt
      logical new,err
      ive_ptr getmem,u,v,w,thet

      nt = ntime

      pvs = getmem(nx*ny*nz*nt) 
      if (pvs.eq.0) then
        write(6,*) 'Error:  could not allocate space for pvs.'
        return
      end if

      do ijk=1,nx*ny*nz*nt
        call s_val(%val(pvs),ijk,0.)
      enddo

      thet = getvar('thet',ndims,dims,stag,min,max,miss,dunits, 
     >                        disunits,dim_names,new)
      if (thet.eq.0) then
        thet = getvar('b',ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new)
        if (thet.eq.0) then
          write(6,*) 'Error:  could not access thet.'
          call freemem(pvs)
          return
        end if
      end if

      u = getvar('u',ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
      if (u.eq.0) then
        write(6,*) 'Error:  could not access u.'
        call freemem(pvs)
        return
      end if

      call pvs_calc_I(%val(pvs),%val(thet),%val(u),%val(zbot_p),
     >             %val(zbot_v),%val(sgz),ztop,dx,dy,nx,ny,nz,nt)

      v = getvar('v',ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
      if (v.eq.0) then
        write(6,*) 'Error:  could not access v.'
        call freemem(pvs)
        return
      end if

      call pvs_calc_II(%val(pvs),%val(thet),%val(v),%val(zbot_p),
     >             %val(zbot_u),%val(sgz),ztop,dx,dy,nx,ny,nz,nt)

      w = getvar('w',ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
      if (w.eq.0) then
        write(6,*) 'Error:  could not access w.'
        call freemem(pvs)
        return
      end if

      call pvs_calc_III(%val(pvs),%val(thet),%val(w),%val(zbot_p),
     >             %val(zbot_u),%val(zbot_v),%val(sgz),%val(wgz),
     >             ztop,dx,dy,nx,ny,nz,nt)

      if (fcor.ne.0.) then
  
        call pvs_calc_IV(%val(pvs),%val(thet),%val(zbot_p),
     >             %val(zbot_v),%val(sgz),ztop,dx,dy,fcor,
     >             nx,ny,nz,nt)

      end if

      dims(3) = dims(3) - 1
      stag(3) = 0.5
      min(3) = sgz_min
      max(3) = sgz_max
      dunits = 'K/m/s'
      disunits = 'K/m/s'
      dim_names(3) = 'nz'

      return
      end

*--------------------------------------------------------------------

      subroutine pvs_calc_I(pvs,thet,u,zs_p,zs_v,sgz,ztop,dx,dy,
     >                            nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real pvs(nx,ny,nz,nt),thet(nx,ny,nz,nt),u(nx+1,ny,nz,nt)
      real zs_p(nx,ny),zs_v(nx,ny+1),sgz(nz),ztop,dx,dy

      integer i,j,k,n,kp1,km1
      real dudzet,dtdzet,dzetdy,dzetdz

      integer nzmax
      parameter (nzmax=201)
      real zlbnd(nzmax),zrbnd(nzmax),dzet(nzmax)

      if (nz.gt.nzmax) then
        write(6,*) 'Error:  nz too big for pvs calculation.'
        stop
      end if

      do k=1,nz-1
        dzet(k) = sgz(k+1) - sgz(k)
        zlbnd(k) = 1.
        zrbnd(k) = 1.
      enddo

      dzet(nz) = 1.
      zlbnd(nz) = 2.
      zlbnd(1) = 0.
      zrbnd(nz) = 0.
      zrbnd(1) = 2.

      do n=1,nt
      do k=1,nz
        kp1 = min(k+1,nz)
        km1 = max(k-1,1)
      do j=1,ny
      do i=1,nx

        dudzet = 0.25*zrbnd(k)
     >       * (u(i,j,kp1,n)+u(i+1,j,kp1,n)
     >              - u(i,j,k,n)-u(i+1,j,k,n))/dzet(k)
     >         + 0.25*zlbnd(k)
     >       * (u(i,j,k,n)+u(i+1,j,k,n)
     >              - u(i,j,km1,n)-u(i+1,j,km1,n))/dzet(km1)
        dtdzet = 0.5*zrbnd(k)*(thet(i,j,kp1,n)-thet(i,j,k,n))
     >                                     /dzet(k)
     >         + 0.5*zlbnd(k)*(thet(i,j,k,n)-thet(i,j,km1,n))
     >                                   /dzet(km1)
        dzetdy = (sgz(k)-ztop)/(ztop-zs_p(i,j))
     >                 *(zs_v(i,j+1)-zs_v(i,j))/dy
        dzetdz = ztop/(ztop-zs_p(i,j))

        if ((j.ne.1).and.(j.ne.ny)) then

          pvs(i,j,k,n) = pvs(i,j,k,n) + dudzet*dzetdz 
     >         *(0.5*(thet(i,j+1,k,n)-thet(i,j-1,k,n))/dy
     >                     + dtdzet*dzetdy)
     >                 - dtdzet*dzetdz
     >         *(0.25*(u(i,j+1,k,n)+u(i+1,j+1,k,n)
     >               - u(i,j-1,k,n)-u(i+1,j-1,k,n))/dy
     >                     + dudzet*dzetdy)

        else if (j.eq.1) then

          pvs(i,j,k,n) = pvs(i,j,k,n) + dudzet*dzetdz
     >          *((thet(i,j+1,k,n)-thet(i,j,k,n))/dy
     >                     + dtdzet*dzetdy)
     >                 - dtdzet*dzetdz
     >          *(0.5*(u(i,j+1,k,n)+u(i+1,j+1,k,n)
     >               - u(i,j,k,n)-u(i+1,j,k,n))/dy
     >                 + dudzet*dzetdy)

        else if (j.eq.ny) then

          pvs(i,j,k,n) = pvs(i,j,k,n) + dudzet*dzetdz
     >          *((thet(i,j,k,n)-thet(i,j-1,k,n))/dy
     >                     + dtdzet*dzetdy)
     >                 - dtdzet*dzetdz
     >          *(0.5*(u(i,j,k,n)+u(i+1,j,k,n)
     >               - u(i,j-1,k,n)-u(i+1,j-1,k,n))/dy
     >                 + dudzet*dzetdy)

        end if

      enddo
      enddo
      enddo
      enddo

      return
      end

*--------------------------------------------------------------------

      subroutine pvs_calc_II(pvs,thet,v,zs_p,zs_u,sgz,ztop,dx,dy,
     >                               nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real pvs(nx,ny,nz,nt),thet(nx,ny,nz,nt),v(nx,ny+1,nz,nt)
      real zs_p(nx,ny),zs_u(nx+1,ny),sgz(nz),ztop,dx,dy

      integer i,j,k,n,kp1,km1
      real dvdzet,dtdzet,dzetdx,dzetdz

      integer nzmax
      parameter (nzmax=201)
      real zlbnd(nzmax),zrbnd(nzmax),dzet(nzmax)

      if (nz.gt.nzmax) then
        write(6,*) 'Error:  nz too big for pvs calculation.'
        stop
      end if

      do k=1,nz-1
        dzet(k) = sgz(k+1) - sgz(k)
        zlbnd(k) = 1.
        zrbnd(k) = 1.
      enddo

      dzet(nz) = 1.
      zlbnd(nz) = 2.
      zlbnd(1) = 0.
      zrbnd(nz) = 0.
      zrbnd(1) = 2.

      do n=1,nt
      do k=1,nz
        kp1 = min(k+1,nz)
        km1 = max(k-1,1)
      do j=1,ny
      do i=1,nx

        dvdzet = 0.25*zrbnd(k)
     >       * (v(i,j,kp1,n)+v(i,j+1,kp1,n)
     >              - v(i,j,k,n)-v(i,j+1,k,n))/dzet(k)
     >         + 0.25*zlbnd(k)
     >       * (v(i,j,k,n)+v(i,j+1,k,n)
     >              - v(i,j,km1,n)-v(i,j+1,km1,n))/dzet(km1)
        dtdzet = 0.5*zrbnd(k)*(thet(i,j,kp1,n)-thet(i,j,k,n))
     >                                     /dzet(k)
     >         + 0.5*zlbnd(k)*(thet(i,j,k,n)-thet(i,j,km1,n))
     >                                   /dzet(km1)
        dzetdx = (sgz(k)-ztop)/(ztop-zs_p(i,j))
     >                 *(zs_u(i+1,j)-zs_u(i,j))/dx
        dzetdz = ztop/(ztop-zs_p(i,j))

        if ((i.ne.1).and.(i.ne.nx)) then

          pvs(i,j,k,n) = pvs(i,j,k,n) - dvdzet*dzetdz
     >         *(0.5*(thet(i+1,j,k,n)-thet(i-1,j,k,n))/dx
     >                     + dtdzet*dzetdx)
     >                 + dtdzet*dzetdz
     >         *(0.25*(v(i+1,j,k,n)+v(i+1,j+1,k,n)
     >               - v(i-1,j,k,n)-v(i-1,j+1,k,n))/dx
     >                     + dvdzet*dzetdx)

        else if (i.eq.1) then

          pvs(i,j,k,n) = pvs(i,j,k,n) - dvdzet*dzetdz
     >         *((thet(i+1,j,k,n)-thet(i,j,k,n))/dx
     >                    + dtdzet*dzetdx)
     >                 + dtdzet*dzetdz
     >         *(0.5*(v(i+1,j,k,n)+v(i+1,j+1,k,n)
     >              - v(i,j,k,n)-v(i,j+1,k,n))/dx
     >                    + dvdzet*dzetdx)

        else if (i.eq.nx) then
   
          pvs(i,j,k,n) = pvs(i,j,k,n) - dvdzet*dzetdz
     >         *((thet(i,j,k,n)-thet(i-1,j,k,n))/dx
     >                    + dtdzet*dzetdx)
     >                 + dtdzet*dzetdz
     >         *(0.5*(v(i,j,k,n)+v(i,j+1,k,n)
     >              - v(i-1,j,k,n)-v(i-1,j+1,k,n))/dx
     >                    + dvdzet*dzetdx)

        end if

      enddo
      enddo
      enddo
      enddo

      return
      end

*----------------------------------------------------------------------

      subroutine pvs_calc_III(pvs,thet,w,zs_p,zs_u,zs_v,sgz,wgz,
     >                          ztop,dx,dy,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real pvs(nx,ny,nz,nt),thet(nx,ny,nz,nt),w(nx,ny,nz+1,nt)
      real zs_p(nx,ny),zs_u(nx+1,ny),zs_v(nx,ny+1)
      real sgz(nz),wgz(nz+1),ztop,dx,dy

      integer i,j,k,n,ip1,im1,jp1,jm1,kp1,km1
      real dwdzet,dtdzet,dzetdx,dzetdy,fx,fy

      integer nzmax
      parameter (nzmax=201)
      real zlbnd(nzmax),zrbnd(nzmax),dzet(nzmax)

      if (nz.gt.nzmax) then
        write(6,*) 'Error:  nz too big for pvs calculation.'
        stop
      end if

      do k=1,nz-1
        dzet(k) = sgz(k+1) - sgz(k)
        zlbnd(k) = 1.
        zrbnd(k) = 1.
      enddo

      dzet(nz) = 1.
      zlbnd(nz) = 2.
      zlbnd(1) = 0.
      zrbnd(nz) = 0.
      zrbnd(1) = 2.

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

        dwdzet = (w(i,j,k+1,n)-w(i,j,k,n))/(wgz(k+1)-wgz(k))
        dtdzet = 0.5*zrbnd(k)*(thet(i,j,kp1,n)-thet(i,j,k,n))
     >                                     /dzet(k)
     >         + 0.5*zlbnd(k)*(thet(i,j,k,n)-thet(i,j,km1,n))
     >                                   /dzet(km1)
        dzetdx = (sgz(k)-ztop)/(ztop-zs_p(i,j))
     >                 *(zs_u(i+1,j)-zs_u(i,j))/dx
        dzetdy = (sgz(k)-ztop)/(ztop-zs_p(i,j))
     >                 *(zs_v(i,j+1)-zs_v(i,j))/dy

        fx = 0.5 
        if ((i.eq.1).or.(i.eq.nx)) fx = 1.0
        fy = 0.5
        if ((j.eq.1).or.(j.eq.ny)) fy = 1.0

        pvs(i,j,k,n) = pvs(i,j,k,n) 
     >           + (0.5*fy*(w(i,jp1,k+1,n)+w(i,jp1,k,n)
     >                  - w(i,jm1,k+1,n)-w(i,jm1,k,n))/dy
     >                       + dwdzet*dzetdy)
     >              *(fx*(thet(ip1,j,k,n)-thet(im1,j,k,n))/dx
     >                       + dtdzet*dzetdx)
     >           - (0.5*fx*(w(ip1,j,k+1,n)+w(ip1,j,k,n)
     >                  - w(im1,j,k+1,n)-w(im1,j,k,n))/dx
     >                       + dwdzet*dzetdx)
     >              *(fy*(thet(i,jp1,k,n)-thet(i,jm1,k,n))/dy
     >                       + dtdzet*dzetdy)

      enddo
      enddo
      enddo
      enddo

      return
      end

*----------------------------------------------------------------------

      subroutine pvs_calc_IV(pvs,thet,zs_p,zs_v,sgz,ztop,dx,dy,fcor,
     >                                 nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real pvs(nx,ny,nz,nt),thet(nx,ny,nz,nt),zs_p(nx,ny)
      real zs_v(nx,ny+1),sgz(nz),ztop,dx,dy,fcor

      integer i,j,k,n,kp1,km1
      real dtdzet,dzetdy,dzetdz,fcorc

      integer nzmax
      parameter (nzmax=201)
      real zlbnd(nzmax),zrbnd(nzmax),dzet(nzmax)

      fcorc = 2.*7.292e-05*cos(asin(fcor/(2.*7.292e-05)))

      if (nz.gt.nzmax) then
        write(6,*) 'Error:  nz too big for pvs calculation.'
        stop
      end if

      do k=1,nz-1
        dzet(k) = sgz(k+1) - sgz(k)
        zlbnd(k) = 1.
        zrbnd(k) = 1.
      enddo

      dzet(nz) = 1.
      zlbnd(nz) = 2.
      zlbnd(1) = 0.
      zrbnd(nz) = 0.
      zrbnd(1) = 2.

      do n=1,nt
      do k=1,nz
        kp1 = min(k+1,nz)
        km1 = max(k-1,1)
      do j=1,ny
      do i=1,nx
        
        dtdzet = 0.5*zrbnd(k)*(thet(i,j,kp1,n)-thet(i,j,k,n))
     >                                     /dzet(k)
     >         + 0.5*zlbnd(k)*(thet(i,j,k,n)-thet(i,j,km1,n))
     >                                   /dzet(km1)
        dzetdy = (sgz(k)-ztop)/(ztop-zs_p(i,j))
     >                 *(zs_v(i,j+1)-zs_v(i,j))/dy
        dzetdz = ztop/(ztop-zs_p(i,j))

        if ((j.ne.1).and.(j.ne.ny)) then
 
          pvs(i,j,k,n) = pvs(i,j,k,n) + fcorc
     >         *(0.5*(thet(i,j+1,k,n)-thet(i,j-1,k,n))/dy
     >                      + dtdzet*dzetdy)
     >                 + fcor*dtdzet*dzetdz

        else if (j.eq.1) then

          pvs(i,j,k,n) = pvs(i,j,k,n) + fcorc
     >         *((thet(i,j+1,k,n)-thet(i,j,k,n))/dy
     >                      + dtdzet*dzetdy)
     >                 + fcor*dtdzet*dzetdz

        else

          pvs(i,j,k,n) = pvs(i,j,k,n) + fcorc
     >         *((thet(i,j,k,n)-thet(i,j-1,k,n))/dy
     >                      + dtdzet*dzetdy)
     >                 + fcor*dtdzet*dzetdz

        end if

      enddo
      enddo
      enddo
      enddo

      return
      end

*--------------------------------------------------------------------------

      ive_ptr function vortz(ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer k,nt
      logical err,new
      ive_ptr getmem,u,v

      nt = ntime

      vortz = getmem((nx+1)*(ny+1)*nz*nt)
      if (vortz.eq.0) then
        write(6,*) 'Error:  could not allocate space for vortz.'
        return
      end if

      do k=1,(nx+1)*(ny+1)*nz*nt
        call s_val(%val(vortz),k,0.)
      enddo

      u = getvar('u',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                       dim_names,new)
      if (u.eq.0) then
        write(6,*) 'Error:  unable to access u.'
        call freemem(vortz)
        return
      end if

      miss = -777.

      call u_vortz(%val(vortz),%val(u),%val(zbot_u),%val(sgz),ztop,
     >                           dy,miss,nx,ny,nz,nt)

      v = getvar('v',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                       dim_names,new)
      if (v.eq.0) then
        write(6,*) 'Error:  unable to access v.'
        call freemem(vortz)
        return
      end if

      miss = -777.

      call v_vortz(%val(vortz),%val(v),%val(zbot_v),%val(sgz),ztop,
     >                           dx,miss,nx,ny,nz,nt)

      dims(1) = nx+1
      dims(2) = ny+1
      stag(1) = 0.0
      stag(2) = 0.0
      min(1) = plmin(1)
      min(2) = plmin(2)  
      max(1) = plmax(1)
      max(2) = plmax(2) 
      dim_names(1) = 'nxp1'
      dim_names(2) = 'nyp1'
      dunits = '1/s'
      disunits = '1/s'

      return
      end

*-----------------------------------------------------------------------

      subroutine v_vortz(vort,v,zs_v,sgz,ztop,dx,miss,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real vort(nx+1,ny+1,nz,nt),v(nx,ny+1,nz,nt),zs_v(nx,ny+1)
      real sgz(nz),ztop,dx,miss

      integer maxnz,i,j,k,n,kp1,km1
      parameter (maxnz=201)
      real zlbnd(maxnz),zrbnd(maxnz),dzet(maxnz),dvdzet,dzetdx,dvdx
      real z0,rdx

      if (nz.gt.maxnz) then
        write(6,*) 'Error:  maxnz too small in v_vortz.'
        write(6,*) 'Stopping.'
        stop
      end if

      do k=1,nz-1
        dzet(k) = sgz(k+1)-sgz(k)
        zlbnd(k) = 1.
        zrbnd(k) = 1.
      enddo

      dzet(nz) = 1.
      zlbnd(1) = 0.
      zlbnd(nz) = 2.
      zrbnd(1) = 2.
      zrbnd(nz) = 0.

      rdx = 1./dx

      do n=1,nt
      do k=1,nz
        kp1 = min(k+1,nz)
        km1 = max(k-1,1)
      do j=2,ny
      do i=2,nx

        z0 = 0.5*(zs_v(i,j)+zs_v(i-1,j))

        dvdzet = 0.25*zrbnd(k)*(v(i,j,kp1,n)-v(i,j,k,n)
     >                    +v(i-1,j,kp1,n)-v(i-1,j,k,n))/dzet(k)
     >         + 0.25*zlbnd(k)*(v(i,j,k,n)-v(i,j,km1,n)
     >                    +v(i-1,j,k,n)-v(i-1,j,km1,n))/dzet(km1)
        dzetdx = (sgz(k)-ztop)/(ztop-z0)
     >                 *rdx*(zs_v(i,j)-zs_v(i-1,j))

        dvdx = rdx*(v(i,j,k,n)-v(i-1,j,k,n))

        vort(i,j,k,n) = vort(i,j,k,n) + dvdx + dvdzet*dzetdx

      enddo
      enddo
      enddo
      enddo

      do n=1,nt
      do k=1,nz
        do i=2,nx
          vort(i,1,k,n) = 2.*vort(i,2,k,n) - vort(i,3,k,n)
          vort(i,ny+1,k,n) = 2.*vort(i,ny,k,n) - vort(i,ny-1,k,n) 
        enddo
        do j=1,ny+1
          vort(1,j,k,n) = 2.*vort(2,j,k,n) - vort(3,j,k,n)
          vort(nx+1,j,k,n) = 2.*vort(nx,j,k,n) - vort(nx-1,j,k,n)
        enddo
      enddo
      enddo
   
      return
      end

*-----------------------------------------------------------------------

      subroutine u_vortz(vort,u,zs_u,sgz,ztop,dy,miss,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real vort(nx+1,ny+1,nz,nt),u(nx+1,ny,nz,nt),zs_u(nx+1,ny)
      real sgz(nz),ztop,dy,miss

      integer maxnz,i,j,k,n,kp1,km1
      parameter (maxnz=201)
      real zlbnd(maxnz),zrbnd(maxnz),dzet(maxnz),dudzet,dzetdy,dudy
      real z0,rdy 

      if (nz.gt.maxnz) then
        write(6,*) 'Error:  maxnz too small in u_vortz.'
        write(6,*) 'Stopping.'
        stop
      end if

      do k=1,nz-1
        dzet(k) = sgz(k+1)-sgz(k)
        zlbnd(k) = 1.
        zrbnd(k) = 1.
      enddo

      dzet(nz) = 1.
      zlbnd(1) = 0.
      zlbnd(nz) = 2.
      zrbnd(1) = 2.
      zrbnd(nz) = 0.

      rdy = 1./dy
 
      do n=1,nt
      do k=1,nz
        kp1 = min(k+1,nz)
        km1 = max(k-1,1)
      do j=2,ny
      do i=2,nx

        z0 = 0.5*(zs_u(i,j)+zs_u(i,j-1))

        dudzet = 0.25*zrbnd(k)*(u(i,j,kp1,n)-u(i,j,k,n)
     >                   +u(i,j-1,kp1,n)-u(i,j-1,k,n))/dzet(k)
     >         + 0.25*zlbnd(k)*(u(i,j,k,n)-u(i,j,km1,n)
     >                   +u(i,j-1,k,n)-u(i,j-1,km1,n))/dzet(km1)
        dzetdy = (sgz(k)-ztop)/(ztop-z0)
     >               *rdy*(zs_u(i,j)-zs_u(i,j-1))

        dudy = rdy*(u(i,j,k,n)-u(i,j-1,k,n))

        vort(i,j,k,n) = vort(i,j,k,n) - dudy - dudzet*dzetdy

      enddo
      enddo
      enddo
      enddo

      do n=1,nt
      do k=1,nz
        do i=2,nx
          vort(i,1,k,n) = 2.*vort(i,2,k,n) - vort(i,3,k,n)
          vort(i,ny+1,k,n) = 2.*vort(i,ny,k,n) - vort(i,ny-1,k,n) 
        enddo
        do j=1,ny+1
          vort(1,j,k,n) = 2.*vort(2,j,k,n) - vort(3,j,k,n)
          vort(nx+1,j,k,n) = 2.*vort(nx,j,k,n) - vort(nx-1,j,k,n)
        enddo
      enddo
      enddo

      return
      end

*----------------------------------------------------------------------

      ive_ptr function vorty(ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer k,nt
      logical err,new
      ive_ptr getmem,u,w

      nt = ntime

      vorty = getmem((nx+1)*ny*(nz+1)*nt)
      if (vorty.eq.0) then
        write(6,*) 'Error:  could not allocate space for vorty.'
        return
      end if

      do k=1,(nx+1)*ny*(nz+1)*nt
        call s_val(%val(vorty),k,0.)
      enddo

      u = getvar('u',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                       dim_names,new)
      if (u.eq.0) then
        write(6,*) 'Error:  unable to access u.'
        call freemem(vorty)
        return
      end if

      miss = -777.
    
      call u_vorty(%val(vorty),%val(u),%val(zbot_u),%val(sgz),ztop,
     >                        miss,nx,ny,nz,nt)

      w = getvar('w',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                       dim_names,new)
      if (w.eq.0) then
        write(6,*) 'Error:  unable to access w.'
        call freemem(vorty)
        return
      end if

      miss = -777.

      call w_vorty(%val(vorty),%val(w),%val(zbot_p),%val(zbot_u),
     >                %val(wgz),ztop,dx,miss,nx,ny,nz,nt)

      dims(1) = nx+1
      dims(3) = nz+1
      stag(1) = 0.
      stag(3) = 0.
      min(1) = plmin(1)
      min(3) = wgz_min
      max(1) = plmax(1)
      max(3) = wgz_max
      dim_names(1) = 'nxp1'
      dim_names(3) = 'nzp1'
      dunits = '1/s'
      disunits = '1/s'

      return
      end

*---------------------------------------------------------------------

      subroutine u_vorty(vort,u,zs_u,sgz,ztop,miss,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real vort(nx+1,ny,nz+1,nt),u(nx+1,ny,nz,nt),zs_u(nx+1,ny)
      real sgz(nz),ztop,miss

      integer maxnz,i,j,k,n
      parameter (maxnz=201)
      real rdz(maxnz),dudzet,dzetdz

      if (nz.gt.maxnz) then
        write(6,*) 'Error:  maxnz too small in u_vorty.'
        write(6,*) 'Stopping.'
        stop
      end if

      do k=1,nz-1
        rdz(k) = 1./(sgz(k+1)-sgz(k))
      enddo

      do n=1,nt
      do k=2,nz
      do j=1,ny
      do i=2,nx

        dudzet = rdz(k-1)*(u(i,j,k,n)-u(i,j,k-1,n))
        dzetdz = ztop/(ztop-zs_u(i,j))

        vort(i,j,k,n) = vort(i,j,k,n) + dudzet*dzetdz

      enddo
      enddo
      enddo
      enddo

      do n=1,nt
      do j=1,ny
        do k=2,nz
          vort(1,j,k,n) = 2.*vort(2,j,k,n) - vort(3,j,k,n)
          vort(nx+1,j,k,n) = 2.*vort(nx,j,k,n) - vort(nx-1,j,k,n)
        enddo
        do i=1,nx+1
          vort(i,j,1,n) = 2.*vort(i,j,2,n) - vort(i,j,3,n)
          vort(i,j,nz+1,n) = 2.*vort(i,j,nz,n) - vort(i,j,nz-1,n)
        enddo
      enddo
      enddo

      return
      end

*------------------------------------------------------------------

      subroutine w_vorty(vort,w,zs_p,zs_u,wgz,ztop,dx,miss,
     >                              nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real vort(nx+1,ny,nz+1,nt),w(nx,ny,nz+1,nt),zs_p(nx,ny)
      real zs_u(nx+1,ny),wgz(nz+1),ztop,dx,miss

      integer maxnz,i,j,k,n
      parameter (maxnz=201)
      real rdx,rdz(maxnz),dwdzet,dzetdx,dwdx

      if (nz.gt.maxnz) then
        write(6,*) 'Error:  maxnz too small in w_vorty.'
        write(6,*) 'Stopping.'
        stop
      end if

      do k=1,nz
        rdz(k) = 1./(wgz(k+1)-wgz(k))
      enddo

      rdx = 1./dx

      do n=1,nt
      do k=2,nz
      do j=1,ny
      do i=2,nx

        dwdzet = 0.25*rdz(k)*(w(i,j,k+1,n)-w(i,j,k,n)
     >                          + w(i-1,j,k+1,n)-w(i-1,j,k,n))
     >         + 0.25*rdz(k-1)*(w(i,j,k,n)-w(i,j,k-1,n)
     >                          + w(i-1,j,k,n)-w(i-1,j,k-1,n))
        dzetdx = (wgz(k)-ztop)/(ztop-zs_u(i,j))
     >                    *rdx*(zs_p(i,j)-zs_p(i-1,j))

        dwdx = rdx*(w(i,j,k,n)-w(i-1,j,k,n))

        vort(i,j,k,n) = vort(i,j,k,n) - dwdx - dwdzet*dzetdx

      enddo
      enddo
      enddo
      enddo

      do n=1,nt
      do j=1,ny
        do k=2,nz
          vort(1,j,k,n) = 2.*vort(2,j,k,n) - vort(3,j,k,n)
          vort(nx+1,j,k,n) = 2.*vort(nx,j,k,n) - vort(nx-1,j,k,n)
        enddo
        do i=1,nx+1
          vort(i,j,1,n) = 2.*vort(i,j,2,n) - vort(i,j,3,n)
          vort(i,j,nz+1,n) = 2.*vort(i,j,nz,n) - vort(i,j,nz-1,n)
        enddo
      enddo
      enddo

      return
      end

*----------------------------------------------------------------------

      ive_ptr function vortx(ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer k,nt
      logical err,new
      ive_ptr getmem,v,w

      nt = ntime

      vortx = getmem(nx*(ny+1)*(nz+1)*nt)
      if (vortx.eq.0) then
        write(6,*) 'Error:  could not allocate space for vortx.'
        return
      end if

      do k=1,nx*(ny+1)*(nz+1)*nt
        call s_val(%val(vortx),k,0.)
      enddo

      v = getvar('v',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                       dim_names,new)
      if (v.eq.0) then
        write(6,*) 'Error:  unable to access v.'
        call freemem(vortx)
        return
      end if

      miss = -777.

      call v_vortx(%val(vortx),%val(v),%val(zbot_v),%val(sgz),ztop,
     >                        miss,nx,ny,nz,nt)

      w = getvar('w',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                       dim_names,new)
      if (w.eq.0) then
        write(6,*) 'Error:  unable to access w.'
        call freemem(vortx)
        return
      end if

      miss = -777.

      call w_vortx(%val(vortx),%val(w),%val(zbot_p),%val(zbot_v),
     >                %val(wgz),ztop,dy,miss,nx,ny,nz,nt)

      dims(2) = ny+1
      dims(3) = nz+1
      stag(2) = 0.
      stag(3) = 0.
      min(2) = plmin(2)
      min(3) = wgz_min
      max(2) = plmax(2)
      max(3) = wgz_max
      dim_names(2) = 'nyp1'
      dim_names(3) = 'nzp1'
      dunits = '1/s'
      disunits = '1/s'

      return
      end

*---------------------------------------------------------------------

      subroutine v_vortx(vort,v,zs_v,sgz,ztop,miss,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real vort(nx,ny+1,nz+1,nt),v(nx,ny+1,nz,nt),zs_v(nx,ny+1)
      real sgz(nz),ztop,miss

      integer maxnz,i,j,k,n
      parameter (maxnz=201)
      real rdz(maxnz),dvdzet,dzetdz

      if (nz.gt.maxnz) then
        write(6,*) 'Error:  maxnz too small in v_vortx.'
        write(6,*) 'Stopping.'
        stop
      end if

      do k=1,nz-1
        rdz(k) = 1./(sgz(k+1)-sgz(k))
      enddo

      do n=1,nt
      do k=2,nz
      do j=2,ny
      do i=1,nx

        dvdzet = rdz(k-1)*(v(i,j,k,n)-v(i,j,k-1,n))
        dzetdz = ztop/(ztop-zs_v(i,j))

        vort(i,j,k,n) = vort(i,j,k,n) - dvdzet*dzetdz

      enddo
      enddo
      enddo
      enddo

      do n=1,nt
      do i=1,nx
        do k=2,nz
          vort(i,1,k,n) = 2.*vort(i,2,k,n) - vort(i,3,k,n)
          vort(i,ny+1,k,n) = 2.*vort(i,ny,k,n) - vort(i,ny-1,k,n)
        enddo
        do j=1,ny+1
          vort(i,j,1,n) = 2.*vort(i,j,2,n) - vort(i,j,3,n)
          vort(i,j,nz+1,n) = 2.*vort(i,j,nz,n) - vort(i,j,nz-1,n)
        enddo
      enddo
      enddo

      return
      end

*------------------------------------------------------------------

      subroutine w_vortx(vort,w,zs_p,zs_v,wgz,ztop,dy,miss,
     >                          nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real vort(nx,ny+1,nz+1,nt),w(nx,ny,nz+1,nt),zs_p(nx,ny)
      real zs_v(nx,ny+1),wgz(nz+1),ztop,dy,miss

      integer maxnz,i,j,k,n
      parameter (maxnz=201)
      real rdy,rdz(maxnz),dwdzet,dzetdy,dwdy

      if (nz.gt.maxnz) then
        write(6,*) 'Error:  maxnz too small in w_vortx.'
        write(6,*) 'Stopping.'
        stop
      end if

      do k=1,nz
        rdz(k) = 1./(wgz(k+1)-wgz(k))
      enddo

      rdy = 1./dy

      do n=1,nt
      do k=2,nz
      do j=2,ny
      do i=1,nx

        dwdzet = 0.25*rdz(k)*(w(i,j,k+1,n)-w(i,j,k,n)
     >                          + w(i,j-1,k+1,n)-w(i,j-1,k,n))
     >         + 0.25*rdz(k-1)*(w(i,j,k,n)-w(i,j,k-1,n)
     >                          + w(i,j-1,k,n)-w(i,j-1,k-1,n))
        dzetdy = (wgz(k)-ztop)/(ztop-zs_v(i,j))
     >                    *rdy*(zs_p(i,j)-zs_p(i,j-1))

        dwdy = rdy*(w(i,j,k,n)-w(i,j-1,k,n)) 

        vort(i,j,k,n) = vort(i,j,k,n) + dwdy + dwdzet*dzetdy

      enddo
      enddo
      enddo
      enddo

      do n=1,nt
      do i=1,nx
        do k=2,nz
          vort(i,1,k,n) = 2.*vort(i,2,k,n) - vort(i,3,k,n)
          vort(i,ny+1,k,n) = 2.*vort(i,ny,k,n) - vort(i,ny-1,k,n)
        enddo
        do j=1,ny+1
          vort(i,j,1,n) = 2.*vort(i,j,2,n) - vort(i,j,3,n)
          vort(i,j,nz+1,n) = 2.*vort(i,j,nz,n) - vort(i,j,nz-1,n)
        enddo
      enddo
      enddo
   
      return
      end
