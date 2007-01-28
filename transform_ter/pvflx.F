
#include "pointer.icl"

      ive_ptr function pvflx(whch,dir,rprndl,ipert,ndims,dims,stag,
     >                   min,max,miss,dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer ipert
      real rprndl
      character*(1) whch,dir

      integer nt
      ive_ptr pvflx_a,pvflx_f_x,pvflx_f_y,pvflx_f_z,pvflx_h
      ive_ptr pvflx1,pvflx2,pvflx3,getmem

      pvflx = 0
      nt = ntime
*
* Advective part
*
      if (whch.eq.'A') then

        pvflx = pvflx_a(dir,ndims,dims,stag,min,max,miss,dunits,
     >                         disunits,dim_names)
*
* Viscous part
*
      else if (whch.eq.'F') then

        if (dir.eq.'X') then
          pvflx = pvflx_f_x(ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names)
        else if (dir.eq.'Y') then 
          pvflx = pvflx_f_y(ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names)
        else if (dir.eq.'Z') then 
          pvflx = pvflx_f_z(ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names)
        else
          write(6,*) 'Error:  cannot compute pvflx for direction ',
     >                          dir(1:1)
          return
        end if
*
* Thermal part
*
      else if (whch.eq.'H') then

        pvflx = pvflx_h(dir,rprndl,ipert,ndims,dims,stag,min,max,
     >                    miss,dunits,disunits,dim_names)
*
* Total (advective + viscous + thermal)
*
      else if (whch.eq.'T') then

        pvflx1 = pvflx_a(dir,ndims,dims,stag,min,max,miss,dunits,
     >                         disunits,dim_names)
        if (pvflx1.eq.0) then
          write(6,*) 'Error:  cannot compute advective flux.'
          return
        end if

        if (dir.eq.'X') then
          pvflx2 = pvflx_f_x(ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names)
        else if (dir.eq.'Y') then 
          pvflx2 = pvflx_f_y(ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names)
        else if (dir.eq.'Z') then 
          pvflx2 = pvflx_f_z(ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names)
        else
          write(6,*) 'Error:  cannot compute pvflx for direction ',
     >                        dir(1:1)
          return
        end if
        if (pvflx2.eq.0) then
          call freemem(pvflx1)
          write(6,*) 'Error:  cannot compute viscous flux.'
          return
        end if
  
        pvflx3 = pvflx_h(dir,rprndl,ipert,ndims,dims,stag,min,max,
     >                    miss,dunits,disunits,dim_names)
        if (pvflx3.eq.0) then
          call freemem(pvflx1)
          call freemem(pvflx2)
          write(6,*) 'Error:  cannot compute thermal flux.'
          return
        end if

        pvflx = getmem(nx*ny*nz*nt)
        if (pvflx.eq.0) then
          call freemem(pvflx1)
          call freemem(pvflx2)
          call freemem(pvflx3)
          write(6,*) 'Error:  unable to allocate space for pvflx.'
          return
        end if

        call pvflx_t_calc(%val(pvflx),%val(pvflx1),%val(pvflx2),
     >                       %val(pvflx3),nx,ny,nz,nt)

        call freemem(pvflx1)
        call freemem(pvflx2)
        call freemem(pvflx3)

      else
 
        write(6,*) 'Error:  do not recognize ',whch(1:1),
     >                  ' part of PV flux.'
        return
 
      end if

      return
      end 

*-----------------------------------------------------------------------

      ive_ptr function pvflx_f_x(ndims,dims,stag,min,max,miss,dunits,
     >                                 disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer ijk,nt
      real r_val,tmpval
      logical new1,new2,new
      ive_ptr getmem,fy,fz,b,tmp1,tmp2

      nt = ntime
      pvflx_f_x = 0

      tmp1 = getvar('VMIX',ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new1) 
      if (tmp1.eq.0) then
        write(6,*) 'Error:  cannot access vmix.'
        return
      end if

      tmp2 = getvar('FILT[V]',ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names,new2)
      if (tmp2.eq.0) then
        write(6,*) 'Error:  cannot access filt[v].'
        return
      end if

      fy = getmem(nx*(ny+1)*nz*nt)
      if (fy.eq.0) then
        write(6,*) 'Error:  could not allocate space for fy.'
        return
      end if 

      do ijk=1,nx*(ny+1)*nz*nt
        tmpval = r_val(%val(tmp1),ijk) + r_val(%val(tmp2),ijk)   
        call s_val(%val(fy),ijk,tmpval)
      enddo

      call f_extrp(%val(fy),nx,ny,nz,nt,0,1,0)

c     if (new1) call driver('free = VMIX',0)
c     if (new2) call driver('free = FILT[V]',0)

      tmp1 = getvar('WMIX',ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new1) 
      if (tmp1.eq.0) then
        call freemem(fy)
        write(6,*) 'Error:  cannot access wmix.'
        return
      end if

      tmp2 = getvar('FILT[W]',ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new2)
      if (tmp2.eq.0) then
        call freemem(fy)
        write(6,*) 'Error:  cannot access filt[w].'
        return
      end if

      fz = getmem(nx*ny*(nz+1)*nt)
      if (fz.eq.0) then
        call freemem(fy)
        write(6,*) 'Error:  could not allocate space for fz.'
        return
      end if

      do ijk=1,nx*ny*(nz+1)*nt
        tmpval = r_val(%val(tmp1),ijk) + r_val(%val(tmp2),ijk)
        call s_val(%val(fz),ijk,tmpval)
      enddo

      call f_extrp(%val(fz),nx,ny,nz,nt,0,0,1)

c     if (new1) call driver('free = WMIX',0)
c     if (new2) call driver('free = FILT[W]',0)

      b = getvar('b',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                       dim_names,new)
      if (b.eq.0) then
        call freemem(fy)
        call freemem(fz)
        write(6,*) 'Error:  cannot access b.'
        return
      end if

      pvflx_f_x = getmem(nx*ny*nz*nt)
      if (pvflx_f_x.eq.0) then
        call freemem(fy)
        call freemem(fz)
        write(6,*) 'Error:  cannot allocate space for pvflx_f_x.'
        return
      end if

      call pvflx_f_x_calc(%val(pvflx_f_x),%val(fy),%val(fz),%val(b),
     >                 %val(dqdy_xy),%val(dqdy_zw),%val(dqdz),
     >                 rdy,nx,ny,nz,nt)

      dunits = 'm/s/s/s/s'
      disunits = 'm/s/s/s/s'

      call freemem(fy)
      call freemem(fz)
 
      return
      end

*------------------------------------------------------------------------

      subroutine pvflx_f_x_calc(pvflx,fy,fz,b,dqdy_xy,dqdy_zw,dqdz, 
     >                                rdy,nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real pvflx(nx,ny,nz,nt),fy(nx,ny+1,nz,nt),fz(nx,ny,nz+1,nt)
      real b(nx,ny,nz,nt),rdy

      real dqdy_xy(nx+1,ny+1,3),dqdy_zw(nz+1)
      real dqdz(nx+1,ny+1,3)

      integer i,j,k,n,m,jp1,jm1,kp1,km1 
      real dbdy,dbdz

      do n=1,nt
*
* north and south boundaries
*
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do m=1,4
        if (m.eq.1) j=1
        if (m.eq.2) j=2
        if (m.eq.3) j=ny-1
        if (m.eq.4) j=ny
        jp1 = min0(j+1,ny)
        jm1 = max0(j-1,1)
      do i=1,nx

        dbdy = 0.5*rdy*( yrbnd(j)*(b(i,jp1,k,n)-b(i,j,k,n))
     >                    + ylbnd(j)*(b(i,j,k,n)-b(i,jm1,k,n)) )
     >                + 0.5*dqdy_xy(i,j,1)
     >         *( zrbnd(k)*rdze(k+1)*dqdy_zw(k+1)
     >                *(b(i,j,kp1,n)-b(i,j,k,n))
     >          + zlbnd(k)*rdze(k)*dqdy_zw(k)
     >                *(b(i,j,k,n)-b(i,j,km1,n)) )

        dbdz = 0.5*dqdz(i,j,1)
     >       *( zrbnd(k)*rdze(k+1)*(b(i,j,kp1,n)-b(i,j,k,n))
     >        + zlbnd(k)*rdze(k)*(b(i,j,k,n)-b(i,j,km1,n)) )

        pvflx(i,j,k,n) = 0.5*(fz(i,j,k,n)+fz(i,j,k+1,n))*dbdy
     >             - 0.5*(fy(i,j,k,n)+fy(i,j+1,k,n))*dbdz

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
      do i=1,nx

        dbdy = 4./6.*rdy*(b(i,j+1,k,n)-b(i,j-1,k,n))
     >           - 1./12.*rdy*(b(i,j+2,k,n)-b(i,j-2,k,n))
     >               + 0.5*dqdy_xy(i,j,1)
     >        *( zrbnd(k)*rdze(k+1)*dqdy_zw(k+1)
     >               *(b(i,j,kp1,n)-b(i,j,k,n))
     >         + zlbnd(k)*rdze(k)*dqdy_zw(k)
     >               *(b(i,j,k,n)-b(i,j,km1,n)) )

        dbdz = 0.5*dqdz(i,j,1)
     >       *( zrbnd(k)*rdze(k+1)*(b(i,j,kp1,n)-b(i,j,k,n))
     >        + zlbnd(k)*rdze(k)*(b(i,j,k,n)-b(i,j,km1,n)) )

        pvflx(i,j,k,n) = 0.5*(fz(i,j,k,n)+fz(i,j,k+1,n))*dbdy
     >             - 0.5*(fy(i,j,k,n)+fy(i,j+1,k,n))*dbdz

      enddo
      enddo
      enddo
*
* end time loop
*
      enddo

      return
      end

*----------------------------------------------------------------------

      ive_ptr function pvflx_f_y(ndims,dims,stag,min,max,miss,dunits,
     >                                 disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer ijk,nt
      real r_val,tmpval
      logical new1,new2,new
      ive_ptr getmem,fx,fz,b,tmp1,tmp2

      nt = ntime
      pvflx_f_y = 0

      tmp1 = getvar('UMIX',ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new1)
      if (tmp1.eq.0) then
        write(6,*) 'Error:  cannot access umix.'
        return
      end if

      tmp2 = getvar('FILT[U]',ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names,new2)
      if (tmp2.eq.0) then
        write(6,*) 'Error:  cannot access filt[u].'
        return
      end if

      fx = getmem((nx+1)*ny*nz*nt) 
      if (fx.eq.0) then
        write(6,*) 'Error:  could not allocate space for fx.'
        return
      end if

      do ijk=1,(nx+1)*ny*nz*nt
        tmpval = r_val(%val(tmp1),ijk) + r_val(%val(tmp2),ijk)
        call s_val(%val(fx),ijk,tmpval)
      enddo

      call f_extrp(%val(fx),nx,ny,nz,nt,1,0,0)

c     if (new1) call driver('free = UMIX',0)
c     if (new2) call driver('free = FILT[U]',0)

      tmp1 = getvar('WMIX',ndims,dims,stag,min,max,miss,dunits,
     >                     disunits,dim_names,new1)
      if (tmp1.eq.0) then
        call freemem(fx)
        write(6,*) 'Error:  cannot access wmix.'
        return
      end if
 
      tmp2 = getvar('FILT[W]',ndims,dims,stag,min,max,miss,dunits,   
     >                        disunits,dim_names,new2) 
      if (tmp2.eq.0) then
        call freemem(fx)
        write(6,*) 'Error:  cannot access filt[w].'
        return
      end if

      fz = getmem(nx*ny*(nz+1)*nt)
      if (fz.eq.0) then
        call freemem(fx)
        write(6,*) 'Error:  could not allocate space for fz.'
        return
      end if

      do ijk=1,nx*ny*(nz+1)*nt
        tmpval = r_val(%val(tmp1),ijk) + r_val(%val(tmp2),ijk)
        call s_val(%val(fz),ijk,tmpval)
      enddo

      call f_extrp(%val(fz),nx,ny,nz,nt,0,0,1)

c     if (new1) call driver('free = WMIX',0)
c     if (new2) call driver('free = FILT[W]',0)

      b = getvar('b',ndims,dims,stag,min,max,miss,dunits,disunits,  
     >                           dim_names,new)
      if (b.eq.0) then
        call freemem(fx)
        call freemem(fz)
        write(6,*) 'Error:  cannot access b.'
        return
      end if

      pvflx_f_y = getmem(nx*ny*nz*nt)
      if (pvflx_f_y.eq.0) then
        call freemem(fx)
        call freemem(fz)
        write(6,*) 'Error:  cannot allocate space for pvflx_f_y.'
        return
      end if

      call pvflx_f_y_calc(%val(pvflx_f_y),%val(fx),%val(fz),%val(b), 
     >                 %val(dqdx_xy),%val(dqdx_zw),%val(dqdz),
     >                 rdx,nx,ny,nz,nt)

      dunits = 'm/s/s/s/s'
      disunits = 'm/s/s/s/s'

      call freemem(fx)
      call freemem(fz)

      return
      end

*------------------------------------------------------------------------  

      subroutine pvflx_f_y_calc(pvflx,fx,fz,b,dqdx_xy,dqdx_zw,dqdz,
     >                                rdx,nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real pvflx(nx,ny,nz,nt),fx(nx+1,ny,nz,nt),fz(nx,ny,nz+1,nt)
      real b(nx,ny,nz,nt),rdx

      real dqdx_xy(nx+1,ny+1,3),dqdx_zw(nz+1)
      real dqdz(nx+1,ny+1,3)

      integer i,j,k,n,m,ip1,im1,kp1,km1
      real dbdx,dbdz

      do n=1,nt
*
* east and west boundaries
*
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do j=1,ny
      do m=1,4
        if (m.eq.1) i=1
        if (m.eq.2) i=2
        if (m.eq.3) i=nx-1
        if (m.eq.4) i=nx
        ip1 = min0(i+1,nx)
        im1 = max0(i-1,1)

        dbdx = 0.5*rdx*( xrbnd(i)*(b(ip1,j,k,n)-b(i,j,k,n))
     >                    + xlbnd(i)*(b(i,j,k,n)-b(im1,j,k,n)) )
     >                + 0.5*dqdx_xy(i,j,1)
     >         *( zrbnd(k)*rdze(k+1)*dqdx_zw(k+1)
     >                *(b(i,j,kp1,n)-b(i,j,k,n))
     >          + zlbnd(k)*rdze(k)*dqdx_zw(k)
     >                *(b(i,j,k,n)-b(i,j,km1,n)) )

        dbdz = 0.5*dqdz(i,j,1)
     >       *( zrbnd(k)*rdze(k+1)*(b(i,j,kp1,n)-b(i,j,k,n))
     >        + zlbnd(k)*rdze(k)*(b(i,j,k,n)-b(i,j,km1,n)) )

        pvflx(i,j,k,n) = 0.5*(fx(i,j,k,n)+fx(i+1,j,k,n))*dbdz
     >             - 0.5*(fz(i,j,k,n)+fz(i,j,k+1,n))*dbdx

      enddo
      enddo
      enddo
*
* interior
*
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do j=1,ny
      do i=3,nx-2

        dbdx = 4./6.*rdx*(b(i+1,j,k,n)-b(i-1,j,k,n))
     >           - 1./12.*rdx*(b(i+2,j,k,n)-b(i-2,j,k,n))
     >               + 0.5*dqdx_xy(i,j,1) 
     >        *( zrbnd(k)*rdze(k+1)*dqdx_zw(k+1)
     >               *(b(i,j,kp1,n)-b(i,j,k,n))
     >         + zlbnd(k)*rdze(k)*dqdx_zw(k)
     >               *(b(i,j,k,n)-b(i,j,km1,n)) )

        dbdz = 0.5*dqdz(i,j,1)
     >       *( zrbnd(k)*rdze(k+1)*(b(i,j,kp1,n)-b(i,j,k,n))
     >        + zlbnd(k)*rdze(k)*(b(i,j,k,n)-b(i,j,km1,n)) )

        pvflx(i,j,k,n) = 0.5*(fx(i,j,k,n)+fx(i+1,j,k,n))*dbdz
     >             - 0.5*(fz(i,j,k,n)+fz(i,j,k+1,n))*dbdx

      enddo
      enddo
      enddo
*
* end time loop
*
      enddo

      return
      end

*------------------------------------------------------------------------

      ive_ptr function pvflx_f_z(ndims,dims,stag,min,max,miss,dunits,
     >                                 disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer ijk,nt
      real r_val,tmpval
      logical new1,new2,new
      ive_ptr getmem,fx,fy,b,tmp1,tmp2

      nt = ntime
      pvflx_f_z = 0

      tmp1 = getvar('UMIX',ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new1)
      if (tmp1.eq.0) then
        write(6,*) 'Error:  cannot access umix.'
        return
      end if

      tmp2 = getvar('FILT[U]',ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names,new2)
      if (tmp2.eq.0) then
        write(6,*) 'Error:  cannot access filt[u].'
        return
      end if

      fx = getmem((nx+1)*ny*nz*nt)
      if (fx.eq.0) then
        write(6,*) 'Error:  cannot allocate space for fx.'
        return
      end if

      do ijk = 1,(nx+1)*ny*nz*nt
        tmpval = r_val(%val(tmp1),ijk) + r_val(%val(tmp2),ijk)
        call s_val(%val(fx),ijk,tmpval)
      enddo

      call f_extrp(%val(fx),nx,ny,nz,nt,1,0,0)

c     if (new1) call driver('free = UMIX',0)
c     if (new2) call driver('free = FILT[U]',0)

      tmp1 = getvar('VMIX',ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new1)
      if (tmp1.eq.0) then
        call freemem(fx)
        write(6,*) 'Error:  cannot access vmix.'
        return
      end if

      tmp2 = getvar('FILT[V]',ndims,dims,stag,min,max,miss,dunits,
     >                        disunits,dim_names,new2)
      if (tmp2.eq.0) then
        call freemem(fx)
        write(6,*) 'Error:  cannot access filt[v].'
        return
      end if

      fy = getmem(nx*(ny+1)*nz*nt)
      if (fy.eq.0) then
        call freemem(fx)
        write(6,*) 'Error:  cannot allocate space for fy.'
        return
      end if

      do ijk=1,nx*(ny+1)*nz*nt
        tmpval = r_val(%val(tmp1),ijk) + r_val(%val(tmp2),ijk)
        call s_val(%val(fy),ijk,tmpval)
      enddo

      call f_extrp(%val(fy),nx,ny,nz,nt,0,1,0)
 
c     if (new1) call driver('free = VMIX',0)
c     if (new2) call driver('free = FILT[V]',0)

      b = getvar('b',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                       dim_names,new)
      if (b.eq.0) then
        call freemem(fx)
        call freemem(fy)
        write(6,*) 'Error:  cannot access b.'
        return
      end if

      pvflx_f_z = getmem(nx*ny*nz*nt)
      if (pvflx_f_z.eq.0) then
        call freemem(fx)
        call freemem(fy)
        write(6,*) 'Error:  cannot allocate space for pvflx_f_z.'
        return
      end if

      call pvflx_f_z_calc(%val(pvflx_f_z),%val(fx),%val(fy),%val(b),
     >                 %val(dqdx_xy),%val(dqdx_zw),%val(dqdy_xy),
     >                 %val(dqdy_zw),rdx,rdy,nx,ny,nz,nt)

      dunits = 'm/s/s/s/s'
      disunits = 'm/s/s/s/s'

      call freemem(fx)
      call freemem(fy)
  
      return
      end 

*-----------------------------------------------------------------------

      subroutine pvflx_f_z_calc(pvflx,fx,fy,b,dqdx_xy,dqdx_zw,
     >                 dqdy_xy,dqdy_zw,rdx,rdy,nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real pvflx(nx,ny,nz,nt),fx(nx+1,ny,nz,nt),fy(nx,ny+1,nz,nt)
      real b(nx,ny,nz,nt),rdx,rdy

      real dqdx_xy(nx+1,ny+1,3),dqdx_zw(nz+1)
      real dqdy_xy(nx+1,ny+1,3),dqdy_zw(nz+1)

      integer i,j,k,n,m,ip1,im1,jp1,jm1,kp1,km1
      real dbdx,dbdy

      do n=1,nt
*
* north and south boundaries
*
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do m=1,4
        if (m.eq.1) j=1
        if (m.eq.2) j=2
        if (m.eq.3) j=ny-1
        if (m.eq.4) j=ny
        jp1 = min0(j+1,ny)
        jm1 = max0(j-1,1)
      do i=1,nx
        ip1 = min0(i+1,nx)
        im1 = max0(i-1,1)

        dbdx = 0.5*rdx*( xrbnd(i)*(b(ip1,j,k,n)-b(i,j,k,n))
     >                    + xlbnd(i)*(b(i,j,k,n)-b(im1,j,k,n)) )
     >                + 0.5*dqdx_xy(i,j,1)
     >         *( zrbnd(k)*rdze(k+1)*dqdx_zw(k+1)
     >                *(b(i,j,kp1,n)-b(i,j,k,n))
     >          + zlbnd(k)*rdze(k)*dqdx_zw(k)
     >                *(b(i,j,k,n)-b(i,j,km1,n)) )

        dbdy = 0.5*rdy*( yrbnd(j)*(b(i,jp1,k,n)-b(i,j,k,n))
     >                    + ylbnd(j)*(b(i,j,k,n)-b(i,jm1,k,n)) )
     >                + 0.5*dqdy_xy(i,j,1)
     >         *( zrbnd(k)*rdze(k+1)*dqdy_zw(k+1)
     >                *(b(i,j,kp1,n)-b(i,j,k,n))
     >          + zlbnd(k)*rdze(k)*dqdy_zw(k)
     >                *(b(i,j,k,n)-b(i,j,km1,n)) )

        pvflx(i,j,k,n) = 0.5*(fy(i,j,k,n)+fy(i,j+1,k,n))*dbdx
     >               - 0.5*(fx(i,j,k,n)+fx(i+1,j,k,n))*dbdy

      enddo
      enddo
      enddo
*
* east and west boundaries
*
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do j=3,ny-2
        jp1 = min0(j+1,ny)
        jm1 = max0(j-1,1)
      do m=1,4
        if (m.eq.1) i = 1
        if (m.eq.2) i = 2
        if (m.eq.3) i = nx-1
        if (m.eq.4) i = nx
        ip1 = min0(i+1,nx)
        im1 = max0(i-1,1)

        dbdx = 0.5*rdx*( xrbnd(i)*(b(ip1,j,k,n)-b(i,j,k,n))
     >                    + xlbnd(i)*(b(i,j,k,n)-b(im1,j,k,n)) )
     >                + 0.5*dqdx_xy(i,j,1)
     >         *( zrbnd(k)*rdze(k+1)*dqdx_zw(k+1)
     >                *(b(i,j,kp1,n)-b(i,j,k,n))
     >          + zlbnd(k)*rdze(k)*dqdx_zw(k)
     >                *(b(i,j,k,n)-b(i,j,km1,n)) )

        dbdy = 0.5*rdy*( yrbnd(j)*(b(i,jp1,k,n)-b(i,j,k,n))
     >                    + ylbnd(j)*(b(i,j,k,n)-b(i,jm1,k,n)) )
     >                + 0.5*dqdy_xy(i,j,1)
     >         *( zrbnd(k)*rdze(k+1)*dqdy_zw(k+1)
     >                *(b(i,j,kp1,n)-b(i,j,k,n))
     >          + zlbnd(k)*rdze(k)*dqdy_zw(k)
     >                *(b(i,j,k,n)-b(i,j,km1,n)) )

        pvflx(i,j,k,n) = 0.5*(fy(i,j,k,n)+fy(i,j+1,k,n))*dbdx
     >               - 0.5*(fx(i,j,k,n)+fx(i+1,j,k,n))*dbdy

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

        dbdx = 4./6.*rdx*(b(i+1,j,k,n)-b(i-1,j,k,n))
     >           - 1./12.*rdx*(b(i+2,j,k,n)-b(i-2,j,k,n))
     >               + 0.5*dqdx_xy(i,j,1) 
     >        *( zrbnd(k)*rdze(k+1)*dqdx_zw(k+1)
     >               *(b(i,j,kp1,n)-b(i,j,k,n))
     >         + zlbnd(k)*rdze(k)*dqdx_zw(k)
     >               *(b(i,j,k,n)-b(i,j,km1,n)) )

        dbdy = 4./6.*rdy*(b(i,j+1,k,n)-b(i,j-1,k,n))
     >           - 1./12.*rdy*(b(i,j+2,k,n)-b(i,j-2,k,n))
     >               + 0.5*dqdy_xy(i,j,1)
     >        *( zrbnd(k)*rdze(k+1)*dqdy_zw(k+1)
     >               *(b(i,j,kp1,n)-b(i,j,k,n))
     >         + zlbnd(k)*rdze(k)*dqdy_zw(k)
     >               *(b(i,j,k,n)-b(i,j,km1,n)) )

        pvflx(i,j,k,n) = 0.5*(fy(i,j,k,n)+fy(i,j+1,k,n))*dbdx
     >               - 0.5*(fx(i,j,k,n)+fx(i+1,j,k,n))*dbdy

      enddo
      enddo
      enddo
*
* end time loop
*
      enddo

      return
      end

*------------------------------------------------------------------------

      subroutine f_extrp(f,nx,ny,nz,nt,iu,iv,iw)

      integer nx,ny,nz,nt,iu,iv,iw
      real f(nx+iu,ny+iv,nz+iw,nt)

      integer i,j,k,n

      if (iu.eq.1) then
        do n=1,nt
        do k=1,nz
        do j=1,ny
          f(1,j,k,n) = 2.*f(2,j,k,n) - f(3,j,k,n)
          f(nx+1,j,k,n) = 2.*f(nx,j,k,n) - f(nx-1,j,k,n)
        enddo
        enddo
        enddo
      else if (iv.eq.1) then
        do n=1,nt
        do k=1,nz
        do i=1,nx
          f(i,1,k,n) = 2.*f(i,2,k,n) - f(i,3,k,n)
          f(i,ny+1,k,n) = 2.*f(i,ny,k,n) - f(i,ny-1,k,n)
        enddo
        enddo
        enddo
      else if (iw.eq.1) then
        do n=1,nt
        do j=1,ny
        do i=1,nx
          f(i,j,1,n) = 2.*f(i,j,2,n) - f(i,j,3,n)
          f(i,j,nz+1,n) = 2.*f(i,j,nz,n) - f(i,j,nz-1,n)
        enddo
        enddo
        enddo
      else
        write(6,*) 'Error:  forcing terms must be at velocity pts.'
        write(6,*) 'Stopping.'
        stop
      end if
  
      return
      end

*-----------------------------------------------------------------------

      ive_ptr function pvflx_h(dir,rprndl,ipert,ndims,dims,stag,min,
     >                    max,miss,dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer ipert
      real rprndl
      character*(1) dir

      integer nt,ijk,iarg(1),strbeg,strend
      real r_val,tmpval,rarg(1)
      character*(80) mkfield,comstr,freestr,carg(1)
      logical new1,new2,newz
      ive_ptr getmem,fb,zet,tmp1,tmp2

      nt = ntime
      pvflx_h = 0

      comstr(1:len(comstr)) = ' '
      iarg(1) = ipert
      rarg(1) = rprndl
      comstr = mkfield('BMIX',carg,iarg,rarg,0,1,1,0)

      WRITE(*,*) 'Calling getvar with ',
     >             comstr(strbeg(comstr):strend(comstr))

      tmp1 = getvar(comstr,ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names,new1) 
      if (tmp1.eq.0) then
        write(6,*) 'Error:  cannot access bmix.'
        return
      end if

      tmp2 = getvar('FILT[B]',ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names,new2)
      if (tmp2.eq.0) then
        write(6,*) 'Error:  cannot access filt[b].'
        return
      end if

      fb = getmem(nx*ny*nz*nt)
      if (fb.eq.0) then
        write(6,*) 'Error:  could not allocate space for fb.'
        return
      end if

      do ijk=1,nx*ny*nz*nt
        tmpval = r_val(%val(tmp1),ijk) + r_val(%val(tmp2),ijk)
        call s_val(%val(fb),ijk,tmpval)
      enddo

      freestr = 'free = '
      freestr = freestr(strbeg(freestr):strend(freestr))
     >                  //comstr(strbeg(comstr):strend(comstr))

      WRITE(*,*) 'Calling driver with ',
     >       freestr(strbeg(freestr):strend(freestr))

c     if (new1) call driver(freestr,0)
c     if (new2) call driver('free = FILT[B]',0)

      if (dir.eq.'X') then
        zet = getvar('Q_IVE',ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,newz)
        if (zet.eq.0) then
          call freemem(fb)
          write(6,*) 'Error:  cannot access q_ive.'
          return
        end if
      else if (dir.eq.'Y') then
        zet = getvar('R_IVE',ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,newz)
        if (zet.eq.0) then
          call freemem(fb)
          write(6,*) 'Error:  cannot access r_ive.'
          return
        end if
      else if (dir.eq.'Z') then
        zet = getvar('S_IVE',ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,newz)
        if (zet.eq.0) then
          call freemem(fb)
          write(6,*) 'Error:  cannot access s_ive.'
          return
        end if
      else
        write(6,*) 'Error: cannot compute PV flux for direction ',
     >                        dir(1:1)
      end if
 
      pvflx_h = getmem(nx*ny*nz*nt)
      if (pvflx_h.eq.0) then
        call freemem(fb)
        write(6,*) 'Error:  cannot allocate space for pvflx_h.'
        return
      end if

      call pvflx_h_calc(%val(pvflx_h),%val(fb),%val(zet),nx,ny,nz,nt)

      dunits = 'm/s/s/s/s'
      disunits = 'm/s/s/s/s'

      call freemem(fb)

      return
      end

*----------------------------------------------------------------------

      subroutine pvflx_h_calc(pvflx,fb,zet,nx,ny,nz,nt)
 
      integer nx,ny,nz,nt
      real pvflx(nx,ny,nz,nt),fb(nx,ny,nz,nt),zet(nx,ny,nz,nt)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx
        pvflx(i,j,k,n) = - zet(i,j,k,n)*fb(i,j,k,n)
      enddo
      enddo
      enddo
      enddo

      return
      end

*-----------------------------------------------------------------------
      
      ive_ptr function pvflx_a(dir,ndims,dims,stag,min,max,miss,
     >                            dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      character*(1) dir

      integer nt,iu,iv,iw
      logical new,newpv
      ive_ptr getmem,pv,vel

      pvflx_a = 0
      nt = ntime

      iu = 0
      iv = 0
      iw = 0
      if (dir.eq.'X') then
        iu = 1
        vel = getvar('u',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)
        if (vel.eq.0) then
          write(6,*) 'Error:  cannot access u.'
          return
        end if
      else if (dir.eq.'Y') then
        iv = 1
        vel = getvar('v',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)
        if (vel.eq.0) then
          write(6,*) 'Error:  cannot access v.'
          return
        end if
      else if (dir.eq.'Z') then
        iw = 1
        vel = getvar('w',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)
        if (vel.eq.0) then
          write(6,*) 'Error:  cannot access w.'
          return
        end if
      else
        write(6,*) 'Error:  cannot compute PV flux for direction ',
     >                         dir(1:1)
        return
      end if

      pv = getvar('PVS',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                              dim_names,newpv)
      if (pv.eq.0) then
        write(6,*) 'Error:  cannot access pvs.'
        return
      end if

      pvflx_a = getmem(nx*ny*nz*nt)
      if (pvflx_a.eq.0) then
        write(6,*) 'Error:  unable to allocate space for pvflx_a.'
        return
      end if

      call pvflx_a_calc(%val(pvflx_a),%val(pv),%val(vel),nx,ny,nz,nt,
     >                               iu,iv,iw)

      dunits = 'm/s/s/s/s'
      disunits = 'm/s/s/s/s'
 
c     if (newpv) call driver('free = PVS',0)

      return
      end

*-------------------------------------------------------------------------

      subroutine pvflx_a_calc(pvflx,pvs,vel,nx,ny,nz,nt,iu,iv,iw)

      integer nx,ny,nz,nt,iu,iv,iw
      real pvflx(nx,ny,nz,nt),pvs(nx,ny,nz,nt)
      real vel(nx+iu,ny+iv,nz+iw,nt)

      integer i,j,k,n

      if (iu.eq.1) then
        do n=1,nt
        do k=1,nz
        do j=1,ny
        do i=1,nx
          pvflx(i,j,k,n) = 0.5*(vel(i,j,k,n)+vel(i+1,j,k,n))
     >                           * pvs(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
      else if (iv.eq.1) then
        do n=1,nt
        do k=1,nz
        do j=1,ny
        do i=1,nx
          pvflx(i,j,k,n) = 0.5*(vel(i,j,k,n)+vel(i,j+1,k,n))
     >                           * pvs(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
      else if (iw.eq.1) then
        do n=1,nt
        do k=1,nz
        do j=1,ny
        do i=1,nx
          pvflx(i,j,k,n) = 0.5*(vel(i,j,k,n)+vel(i,j,k+1,n))
     >                           * pvs(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
      else
        write(6,*) 'Error:  do not know which component to compute.'
        write(6,*) 'Stopping.'
        stop
      end if

      return
      end

*-------------------------------------------------------------------------

      subroutine pvflx_t_calc(pvflx,pvflx1,pvflx2,pvflx3,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real pvflx(nx,ny,nz,nt),pvflx1(nx,ny,nz,nt),pvflx2(nx,ny,nz,nt)
      real pvflx3(nx,ny,nz,nt)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx
        pvflx(i,j,k,n) = pvflx1(i,j,k,n)+pvflx2(i,j,k,n)
     >                                     + pvflx3(i,j,k,n)
      enddo
      enddo
      enddo
      enddo

      return
      end
