
#include "pointer.icl"

      ive_ptr function umix(ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer ijk,nt
      logical new
      ive_ptr def,u,v,w,km,getmem

      nt = ntime

      umix = 0

      km = getvar('km',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                                 dim_names,new)
      if (km.eq.0) then
        write(6,*) 'ERROR: cannot access variable km.'
        return
      end if

      v = getvar('v',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                            dim_names,new)
      if (v.eq.0) then
        write(6,*) 'ERROR:  cannot access variable v.'
        return
      end if

      w = getvar('w',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                            dim_names,new)
      if (w.eq.0) then
        write(6,*) 'ERROR:  cannot access variable w.'
        return
      end if      

      u = getvar('u',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                            dim_names,new)
      if (u.eq.0) then
        write(6,*) 'ERROR:  cannot access variable u.'
        return
      end if

      def = getmem((nx+1)*(ny+1)*(nz+1)*nt)
      if (def.eq.0) then
        write(6,*) 'ERROR:  unable to allocate space for def.'
        return
      end if

      umix = getmem((nx+1)*ny*nz*nt)
      if (umix.eq.0) then
        write(6,*) 'ERROR:  unable to allocate space for umix.'
        call freemem(def)
        return
      end if

      do ijk=1,(nx+1)*ny*nz*nt
        call s_val(%val(umix),ijk,0.) 
      enddo
*
* calculate contribution from def11
*
      call defii_calc(%val(def),%val(u),%val(v),%val(w),%val(dqdx_xy),
     >             %val(dqdx_zt),%val(dqdy_xy),%val(dqdy_zt),
     >             %val(dqdz),rdx,rdy,nx,ny,nz,nt,1)

      call add_km(%val(def),%val(km),nx,ny,nz,nt,0,0,0)

      call umix_d11calc(%val(umix),%val(def),%val(dqdx_xy),
     >           %val(dqdx_zw),rdx,nx,ny,nz,nt)
*
* calculate contribution from def12
*
      call def12_calc(%val(def),%val(u),%val(v),%val(dqdx_xy),
     >           %val(dqdx_zt),%val(dqdy_xy),%val(dqdy_zt),
     >           rdx,rdy,nx,ny,nz,nt)

      call add_km(%val(def),%val(km),nx,ny,nz,nt,1,1,0)

      call umix_d12calc(%val(umix),%val(def),%val(dqdy_xy),
     >           %val(dqdy_zw),rdy,nx,ny,nz,nt)
*
* calculate contribution from def13
*
      call def13_calc(%val(def),%val(u),%val(w),%val(dqdx_xy),
     >          %val(dqdx_zw),%val(dqdz),rdx,nx,ny,nz,nt)

      call add_km(%val(def),%val(km),nx,ny,nz,nt,1,0,1)

      call umix_d13calc(%val(umix),%val(def),%val(dqdz),
     >                          nx,ny,nz,nt)

      dunits = 'm/s/s'
      disunits = 'm/s/s'

      call freemem(def)

      return
      end

*------------------------------------------------------------------------

      subroutine umix_d11calc(umix,def11,dqdx_xy,dqdx_zw,rdx,nx,ny,
     >                              nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real umix(nx+1,ny,nz,nt),def11(nx,ny,nz,nt),rdx

      real dqdx_xy(nx+1,ny+1,3),dqdx_zw(nz+1)

      integer i,j,k,n,kp1,km1

      do n=1,nt
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do j=1,ny
      do i=2,nx

        umix(i,j,k,n) = umix(i,j,k,n)
     >              + rdx*(def11(i,j,k,n)-def11(i-1,j,k,n))
     >                  + 0.25*dqdx_xy(i,j,2)
     >   *( zrbnd(k)*rdze(k+1)*dqdx_zw(k+1)*(def11(i,j,kp1,n)
     >      + def11(i-1,j,kp1,n)-def11(i,j,k,n)-def11(i-1,j,k,n))
     >    + zlbnd(k)*rdze(k)*dqdx_zw(k)*(def11(i,j,k,n)
     >      + def11(i-1,j,k,n)-def11(i,j,km1,n)-def11(i-1,j,km1,n)) )

      enddo
      enddo
      enddo
      enddo

      return
      end

*------------------------------------------------------------------------

      subroutine umix_d12calc(umix,def12,dqdy_xy,dqdy_zw,rdy,nx,ny,
     >                                nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real umix(nx+1,ny,nz,nt),def12(nx+1,ny+1,nz,nt),rdy

      real dqdy_xy(nx+1,ny+1,3),dqdy_zw(nz+1)

      integer i,j,k,n,kp1,km1 

      do n=1,nt
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do j=2,ny-1
      do i=2,nx

        umix(i,j,k,n) = umix(i,j,k,n)
     >            + rdy*(def12(i,j+1,k,n)-def12(i,j,k,n))
     >                + 0.25*dqdy_xy(i,j,2)
     >   *( zrbnd(k)*rdze(k+1)*dqdy_zw(k+1)*(def12(i,j+1,kp1,n)
     >       + def12(i,j,kp1,n)-def12(i,j+1,k,n)-def12(i,j,k,n))
     >    + zlbnd(k)*rdze(k)*dqdy_zw(k)*(def12(i,j+1,k,n)
     >       + def12(i,j,k,n)-def12(i,j+1,km1,n)-def12(i,j,km1,n)) )

      enddo
      enddo
      enddo
      enddo

      return
      end

*-----------------------------------------------------------------------

      subroutine umix_d13calc(umix,def13,dqdz,nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real umix(nx+1,ny,nz,nt),def13(nx+1,ny,nz+1,nt)

      real dqdz(nx+1,ny+1,3)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz-1
      do j=1,ny
      do i=2,nx

        umix(i,j,k,n) = umix(i,j,k,n)
     >    + rdzc(k)*(def13(i,j,k+1,n)-def13(i,j,k,n))*dqdz(i,j,2)

      enddo
      enddo
      enddo
      enddo

      return
      end

*-----------------------------------------------------------------------

      ive_ptr function vmix(ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer ijk,nt
      logical new
      ive_ptr def,u,v,w,km,getmem

      nt = ntime

      vmix = 0

      km = getvar('km',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                                 dim_names,new)
      if (km.eq.0) then
        write(6,*) 'ERROR: cannot access variable km.'
        return
      end if

      w = getvar('w',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                            dim_names,new)
      if (w.eq.0) then
        write(6,*) 'ERROR:  cannot access variable w.'
        return
      end if      

      u = getvar('u',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                            dim_names,new)
      if (u.eq.0) then
        write(6,*) 'ERROR:  cannot access variable u.'
        return
      end if

      v = getvar('v',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                            dim_names,new)
      if (v.eq.0) then
        write(6,*) 'ERROR:  cannot access variable v.'
        return
      end if

      def = getmem((nx+1)*(ny+1)*(nz+1)*nt)
      if (def.eq.0) then
        write(6,*) 'ERROR:  unable to allocate space for def.'
        return
      end if

      vmix = getmem(nx*(ny+1)*nz*nt)
      if (vmix.eq.0) then
        write(6,*) 'ERROR:  unable to allocate space for vmix.'
        call freemem(def)
        return
      end if

      do ijk=1,nx*(ny+1)*nz*nt
        call s_val(%val(vmix),ijk,0.) 
      enddo
*
* calculate contribution from def22
*
      call defii_calc(%val(def),%val(u),%val(v),%val(w),%val(dqdx_xy),
     >             %val(dqdx_zt),%val(dqdy_xy),%val(dqdy_zt),
     >             %val(dqdz),rdx,rdy,nx,ny,nz,nt,2)

      call add_km(%val(def),%val(km),nx,ny,nz,nt,0,0,0)

      call vmix_d22calc(%val(vmix),%val(def),%val(dqdy_xy),
     >           %val(dqdy_zw),rdy,nx,ny,nz,nt)
*
* calculate contribution from def12
*
      call def12_calc(%val(def),%val(u),%val(v),%val(dqdx_xy),
     >           %val(dqdx_zt),%val(dqdy_xy),%val(dqdy_zt),
     >           rdx,rdy,nx,ny,nz,nt)

      call add_km(%val(def),%val(km),nx,ny,nz,nt,1,1,0)

      call vmix_d12calc(%val(vmix),%val(def),%val(dqdx_xy),
     >           %val(dqdx_zw),rdx,nx,ny,nz,nt)
*
* calculate contribution from def23
*
      call def23_calc(%val(def),%val(v),%val(w),%val(dqdy_xy),
     >          %val(dqdy_zw),%val(dqdz),rdy,nx,ny,nz,nt)

      call add_km(%val(def),%val(km),nx,ny,nz,nt,0,1,1)

      call vmix_d23calc(%val(vmix),%val(def),%val(dqdz),
     >                          nx,ny,nz,nt)

      dunits = 'm/s/s'
      disunits = 'm/s/s'

      call freemem(def)

      return
      end

*------------------------------------------------------------------------

      subroutine vmix_d22calc(vmix,def22,dqdy_xy,dqdy_zw,rdy,
     >                              nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real vmix(nx,ny+1,nz,nt),def22(nx,ny,nz,nt),rdy

      real dqdy_xy(nx+1,ny+1,3),dqdy_zw(nz+1)

      integer i,j,k,n,kp1,km1

      do n=1,nt
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do j=2,ny
      do i=1,nx

        vmix(i,j,k,n) = vmix(i,j,k,n)
     >              + rdy*(def22(i,j,k,n)-def22(i,j-1,k,n))
     >                  + 0.25*dqdy_xy(i,j,3)
     >   *( zrbnd(k)*rdze(k+1)*dqdy_zw(k+1)*(def22(i,j,kp1,n)
     >      + def22(i,j-1,kp1,n)-def22(i,j,k,n)-def22(i,j-1,k,n))
     >    + zlbnd(k)*rdze(k)*dqdy_zw(k)*(def22(i,j,k,n)
     >      + def22(i,j-1,k,n)-def22(i,j,km1,n)-def22(i,j-1,km1,n)) )

      enddo
      enddo
      enddo
      enddo

      return
      end      

*------------------------------------------------------------------------

      subroutine vmix_d12calc(vmix,def12,dqdx_xy,dqdx_zw,rdx,
     >                          nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real vmix(nx,ny+1,nz,nt),def12(nx+1,ny+1,nz,nt),rdx

      real dqdx_xy(nx+1,ny+1,3),dqdx_zw(nz+1)

      integer i,j,k,n,kp1,km1

      do n=1,nt
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do j=2,ny
      do i=2,nx-1

        vmix(i,j,k,n) = vmix(i,j,k,n)
     >            + rdx*(def12(i+1,j,k,n)-def12(i,j,k,n))
     >                + 0.25*dqdx_xy(i,j,3)
     >   *( zrbnd(k)*rdze(k+1)*dqdx_zw(k+1)*(def12(i+1,j,kp1,n)
     >       + def12(i,j,kp1,n)-def12(i+1,j,k,n)-def12(i,j,k,n))
     >    + zlbnd(k)*rdze(k)*dqdx_zw(k)*(def12(i+1,j,k,n)
     >       + def12(i,j,k,n)-def12(i+1,j,km1,n)-def12(i,j,km1,n)) )

      enddo
      enddo
      enddo
      enddo

      return
      end

*------------------------------------------------------------------------

      subroutine vmix_d23calc(vmix,def23,dqdz,nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real vmix(nx,ny+1,nz,nt),def23(nx,ny+1,nz+1,nt)

      real dqdz(nx+1,ny+1,3)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz-1
      do j=2,ny
      do i=1,nx

        vmix(i,j,k,n) = vmix(i,j,k,n) 
     >    + rdzc(k)*(def23(i,j,k+1,n)-def23(i,j,k,n))*dqdz(i,j,3)

      enddo
      enddo
      enddo
      enddo

      return
      end
   
*------------------------------------------------------------------------

      ive_ptr function wmix(ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer ijk,nt
      logical new
      ive_ptr def,u,v,w,km,getmem

      nt = ntime

      wmix = 0

      km = getvar('km',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                                 dim_names,new)
      if (km.eq.0) then
        write(6,*) 'ERROR: cannot access variable km.'
        return
      end if

      u = getvar('u',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                            dim_names,new)
      if (u.eq.0) then
        write(6,*) 'ERROR:  cannot access variable u.'
        return
      end if

      v = getvar('v',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                            dim_names,new)
      if (v.eq.0) then
        write(6,*) 'ERROR:  cannot access variable v.'
        return
      end if

      w = getvar('w',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                            dim_names,new)
      if (w.eq.0) then
        write(6,*) 'ERROR:  cannot access variable w.'
        return
      end if      

      def = getmem((nx+1)*(ny+1)*(nz+1)*nt)
      if (def.eq.0) then
        write(6,*) 'ERROR:  unable to allocate space for def.'
        return
      end if

      wmix = getmem(nx*ny*(nz+1)*nt)
      if (wmix.eq.0) then
        write(6,*) 'ERROR:  unable to allocate space for wmix.'
        call freemem(def)
        return
      end if

      do ijk=1,nx*ny*(nz+1)*nt
        call s_val(%val(wmix),ijk,0.) 
      enddo
*
* calculate contribution from def33
*
      call defii_calc(%val(def),%val(u),%val(v),%val(w),%val(dqdx_xy),
     >             %val(dqdx_zt),%val(dqdy_xy),%val(dqdy_zt),
     >             %val(dqdz),rdx,rdy,nx,ny,nz,nt,3)

      call add_km(%val(def),%val(km),nx,ny,nz,nt,0,0,0)

      call wmix_d33calc(%val(wmix),%val(def),%val(dqdz),nx,ny,nz,nt)
*
* calculate contribution from def13
*
      call def13_calc(%val(def),%val(u),%val(w),%val(dqdx_xy),
     >          %val(dqdx_zw),%val(dqdz),rdx,nx,ny,nz,nt)

      call add_km(%val(def),%val(km),nx,ny,nz,nt,1,0,1)

      call wmix_d13calc(%val(wmix),%val(def),%val(dqdx_xy),
     >           %val(dqdx_zt),rdx,nx,ny,nz,nt)
*
* calculate contribution from def23
*
      call def23_calc(%val(def),%val(v),%val(w),%val(dqdy_xy),
     >          %val(dqdy_zw),%val(dqdz),rdy,nx,ny,nz,nt)

      call add_km(%val(def),%val(km),nx,ny,nz,nt,0,1,1)

      call wmix_d23calc(%val(wmix),%val(def),%val(dqdy_xy),
     >           %val(dqdy_zt),rdy,nx,ny,nz,nt)  

      dunits = 'm/s/s'
      disunits = 'm/s/s'

      call freemem(def)

      return
      end

*------------------------------------------------------------------------
      
      subroutine wmix_d33calc(wmix,def33,dqdz,nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real wmix(nx,ny,nz+1,nt),def33(nx,ny,nz,nt)

      real dqdz(nx+1,ny+1,3)

      integer i,j,k,n

      do n=1,nt
      do k=2,nz
      do j=1,ny
      do i=1,nx

        wmix(i,j,k,n) = wmix(i,j,k,n) 
     >    + rdze(k)*(def33(i,j,k,n)-def33(i,j,k-1,n))*dqdz(i,j,1)

      enddo
      enddo
      enddo
      enddo

      return
      end

*------------------------------------------------------------------------

      subroutine wmix_d13calc(wmix,def13,dqdx_xy,dqdx_zt,rdx,
     >                               nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real wmix(nx,ny,nz+1,nt),def13(nx+1,ny,nz+1,nt),rdx

      real dqdx_xy(nx+1,ny+1,3),dqdx_zt(nz)

      integer i,j,k,n,kp1,km1

      do n=1,nt
      do k=2,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do j=1,ny
      do i=2,nx-1

        wmix(i,j,k,n) = wmix(i,j,k,n)
     >            + rdx*(def13(i+1,j,k,n)-def13(i,j,k,n))
     >                   + 0.25*dqdx_xy(i,j,1) 
     >   *( zrbnd(k)*rdzc(k)*dqdx_zt(k)*(def13(i,j,kp1,n)
     >      + def13(i+1,j,kp1,n)-def13(i,j,k,n)-def13(i+1,j,k,n))
     >    + zlbnd(k)*rdzc(k-1)*dqdx_zt(k-1)*(def13(i,j,k,n)
     >      + def13(i+1,j,k,n)-def13(i,j,km1,n)-def13(i+1,j,km1,n)) )

      enddo
      enddo
      enddo
      enddo

      return
      end

*------------------------------------------------------------------------

      subroutine wmix_d23calc(wmix,def23,dqdy_xy,dqdy_zt,rdy,
     >                               nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real wmix(nx,ny,nz+1,nt),def23(nx,ny+1,nz+1,nt),rdy

      real dqdy_xy(nx+1,ny+1,3),dqdy_zt(nz)

      integer i,j,k,n,kp1,km1

      do n=1,nt
      do k=2,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do j=2,ny-1
      do i=1,nx   

        wmix(i,j,k,n) = wmix(i,j,k,n)
     >            + rdy*(def23(i,j+1,k,n)-def23(i,j,k,n))
     >                   + 0.25*dqdy_xy(i,j,1)
     >   *( zrbnd(k)*rdzc(k)*dqdy_zt(k)*(def23(i,j,kp1,n)
     >      + def23(i,j+1,kp1,n)-def23(i,j,k,n)-def23(i,j+1,k,n))
     >    + zlbnd(k)*rdzc(k-1)*dqdy_zt(k-1)*(def23(i,j,k,n)
     >      + def23(i,j+1,k,n)-def23(i,j,km1,n)-def23(i,j+1,km1,n)) )

      enddo
      enddo
      enddo
      enddo

      return
      end 

*------------------------------------------------------------------------

      subroutine defii_calc(defii,u,v,w,dqdx_xy,dqdx_zt,dqdy_xy,
     >              dqdy_zt,dqdz,rdx,rdy,nx,ny,nz,nt,iwhch)

#include "misc.icl"

      integer nx,ny,nz,nt,iwhch
      real defii(nx,ny,nz,nt),u(nx+1,ny,nz,nt),v(nx,ny+1,nz,nt)
      real w(nx,ny,nz+1,nt),rdx,rdy

      real dqdx_xy(nx+1,ny+1,3),dqdx_zt(nz)
      real dqdy_xy(nx+1,ny+1,3),dqdy_zt(nz)
      real dqdz(nx+1,ny+1,3)

      real dudx,dvdy,dwdz
      integer i,j,k,n,kp1,km1

      do n=1,nt
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do j=1,ny
      do i=1,nx

        dudx = rdx*(u(i+1,j,k,n)-u(i,j,k,n)) 
     >          + 0.25*dqdx_xy(i,j,1)*dqdx_zt(k)
     >     *( zrbnd(k)*rdze(k+1)*(u(i,j,kp1,n)+u(i+1,j,kp1,n)
     >                    - u(i,j,k,n)-u(i+1,j,k,n))
     >      + zlbnd(k)*rdze(k)*(u(i,j,k,n)+u(i+1,j,k,n)
     >                    - u(i,j,km1,n)-u(i+1,j,km1,n)) ) 

        dvdy = rdy*(v(i,j+1,k,n)-v(i,j,k,n)) 
     >          + 0.25*dqdy_xy(i,j,1)*dqdy_zt(k)
     >     *( zrbnd(k)*rdze(k+1)*(v(i,j,kp1,n)+v(i,j+1,kp1,n)
     >                    - v(i,j,k,n)-v(i,j+1,k,n))
     >      + zlbnd(k)*rdze(k)*(v(i,j,k,n)+v(i,j+1,k,n)
     >                    - v(i,j,km1,n)-v(i,j+1,km1,n)) )

        dwdz = rdzc(k)*(w(i,j,k+1,n)-w(i,j,k,n))*dqdz(i,j,1)

        if (iwhch.eq.1) then
          defii(i,j,k,n) = 4./3.*dudx - 2./3.*(dvdy+dwdz)
        else if (iwhch.eq.2) then
          defii(i,j,k,n) = 4./3.*dvdy - 2./3.*(dudx+dwdz)
        else if (iwhch.eq.3) then
          defii(i,j,k,n) = 4./3.*dwdz - 2./3.*(dudx+dvdy)
        else
          write(6,*) 'Error:  iwhch not equal to 1, 2, or 3 in ',
     >                           'defii_calc.'
          write(6,*) 'Stopping.'
          stop
        end if             

      enddo
      enddo
      enddo 
      enddo 

      return
      end

*--------------------------------------------------------------------------

      subroutine def12_calc(def12,u,v,dqdx_xy,dqdx_zt,dqdy_xy,
     >              dqdy_zt,rdx,rdy,nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real def12(nx+1,ny+1,nz,nt),u(nx+1,ny,nz,nt),v(nx,ny+1,nz,nt)
      real rdx,rdy

      real dqdx_xy(nx+1,ny+1,3),dqdx_zt(nz)
      real dqdy_xy(nx+1,ny+1,3),dqdy_zt(nz)

      integer i,j,k,n,kp1,km1
      real dudy,dvdx

      do n=1,nt
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do j=2,ny
      do i=2,nx

        dudy = rdy*(u(i,j,k,n)-u(i,j-1,k,n))
     >         + 0.5*(dqdy_xy(i,j,2)+dqdy_xy(i,j-1,2))*dqdy_zt(k)
     >      * 0.25*( zrbnd(k)*rdze(k+1)*(u(i,j,kp1,n)+u(i,j-1,kp1,n)
     >                - u(i,j,k,n)-u(i,j-1,k,n))
     >             + zlbnd(k)*rdze(k)*(u(i,j,k,n)+u(i,j-1,k,n)
     >                - u(i,j,km1,n)-u(i,j-1,km1,n)) )

        dvdx = rdx*(v(i,j,k,n)-v(i-1,j,k,n))
     >         + 0.5*(dqdx_xy(i,j,3)+dqdx_xy(i-1,j,3))*dqdy_zt(k)
     >      * 0.25*( zrbnd(k)*rdze(k+1)*(v(i,j,kp1,n)+v(i-1,j,kp1,n)
     >                - v(i,j,k,n)-v(i-1,j,k,n)) 
     >             + zlbnd(k)*rdze(k)*(v(i,j,k,n)+v(i-1,j,k,n)
     >                - v(i,j,km1,n)-v(i-1,j,km1,n)) )

        def12(i,j,k,n) = dudy + dvdx

      enddo
      enddo
      enddo
      enddo
*
* extend to boundaries
*
      do n=1,nt
      do k=1,nz-1
      do j=2,ny
        def12(1,j,k,n) = def12(2,j,k,n)
        def12(nx+1,j,k,n) = def12(nx,j,k,n)
      enddo
      enddo
      enddo

      do n=1,nt
      do k=1,nz-1
      do i=1,nx+1
        def12(i,1,k,n) = def12(i,2,k,n)
        def12(i,ny+1,k,n) = def12(i,ny,k,n)
      enddo
      enddo
      enddo

      return
      end

*-----------------------------------------------------------------------

      subroutine def13_calc(def13,u,w,dqdx_xy,dqdx_zw,dqdz,rdx,
     >                              nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real def13(nx+1,ny,nz+1,nt),u(nx+1,ny,nz,nt),w(nx,ny,nz+1,nt)
      real rdx

      real dqdx_xy(nx+1,ny+1,3),dqdx_zw(nz+1),dqdz(nx+1,ny+1,3)

      integer i,j,k,n
      real dudz,dwdx

      do n=1,nt
      do k=2,nz
      do j=1,ny
      do i=2,nx

        dudz = rdze(k)*(u(i,j,k,n)-u(i,j,k-1,n))*dqdz(i,j,2)

        dwdx = rdx*(w(i,j,k,n)-w(i-1,j,k,n))
     >          + 0.25*dqdx_xy(i,j,2)*dqdx_zw(k)
     >     *( rdzc(k)*(w(i,j,k+1,n)+w(i-1,j,k+1,n)
     >                     - w(i,j,k,n)-w(i-1,j,k,n))
     >      + rdzc(k-1)*(w(i,j,k,n)+w(i-1,j,k,n)
     >                     - w(i,j,k-1,n)-w(i-1,j,k-1,n)) )

        def13(i,j,k,n) = dudz + dwdx

      enddo
      enddo
      enddo
      enddo
*
* extend to boundaries
*
      do n=1,nt
      do j=1,ny
      do i=2,nx
        def13(i,j,1,n) = 0.
        def13(i,j,nz+1,n) = def13(i,j,nz,n)
      enddo
      enddo
      enddo

      do n=1,nt
      do k=1,nz+1
      do j=1,ny
        def13(1,j,k,n) = def13(2,j,k,n)
        def13(nx+1,j,k,n) = def13(nx,j,k,n)
      enddo
      enddo
      enddo

      return
      end

*-----------------------------------------------------------------------

      subroutine def23_calc(def23,v,w,dqdy_xy,dqdy_zw,dqdz,rdy,
     >                              nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt        
      real def23(nx,ny+1,nz+1,nt),v(nx,ny+1,nz,nt),w(nx,ny,nz+1,nt)
      real rdy

      real dqdy_xy(nx+1,ny+1,3),dqdy_zw(nz+1),dqdz(nx+1,ny+1,3)

      integer i,j,k,n
      real dvdz,dwdy

      do n=1,nt
      do k=2,nz
      do j=2,ny
      do i=1,nx

        dvdz = rdze(k)*(v(i,j,k,n)-v(i,j,k-1,n))*dqdz(i,j,3)

        dwdy = rdy*(w(i,j,k,n)-w(i,j-1,k,n))
     >          + 0.25*dqdy_xy(i,j,3)*dqdy_zw(k)
     >     *( rdzc(k)*(w(i,j,k+1,n)+w(i,j-1,k+1,n)
     >                     - w(i,j,k,n)-w(i,j-1,k,n))
     >      + rdzc(k-1)*(w(i,j,k,n)+w(i,j-1,k,n)
     >                     - w(i,j,k-1,n)-w(i,j-1,k-1,n)) )

        def23(i,j,k,n) = dvdz + dwdy

      enddo
      enddo
      enddo
      enddo
*
* extend to boundaries
*
      do n=1,nt
      do j=2,ny
      do i=1,nx
        def23(i,j,1,n) = 0.
        def23(i,j,nz+1,n) = def23(i,j,nz,n)
      enddo
      enddo
      enddo

      do n=1,nt
      do k=1,nz+1
      do i=1,nx
        def23(i,1,k,n) = def23(i,2,k,n)
        def23(i,ny+1,k,n) = def23(i,ny,k,n)
      enddo
      enddo
      enddo

      return
      end 

*----------------------------------------------------------------------

      subroutine add_km(def,km,nx,ny,nz,nt,iu,iv,iw)

      integer nx,ny,nz,nt,iu,iv,iw
      real def(nx+iu,ny+iv,nz+iw,nt),km(nx,ny,nz,nt)

      integer i,j,k,n,im1,jm1,km1
 
      do n=1,nt
      do k=1,nz
        km1 = max0(k-1,1)
      do j=1,ny
        jm1 = max0(j-1,1)       
      do i=1,nx
        im1 = max0(i-1,1)

        if ((iu.eq.0).and.(iv.eq.0).and.(iw.eq.0)) then
*
* defii case
*
          def(i,j,k,n) = km(i,j,k,n)*def(i,j,k,n)

        else if (iw.eq.0) then
*
* def12 case
*
          def(i,j,k,n) = def(i,j,k,n)
     >       * 0.25*(km(i,j,k,n)+km(im1,j,k,n)+km(im1,jm1,k,n)
     >                         + km(i,jm1,k,n)) 

        else if (iv.eq.0) then
*
* def13 case
*
          def(i,j,k,n) = def(i,j,k,n)
     >       * 0.25*(km(i,j,k,n)+km(im1,j,k,n)+km(im1,j,km1,n)
     >                         + km(i,j,km1,n))

        else if (iu.eq.0) then
* 
* def23 case
*
          def(i,j,k,n) = def(i,j,k,n)
     >       * 0.25*(km(i,j,k,n)+km(i,jm1,k,n)+km(i,jm1,km1,n)
     >                         + km(i,j,km1,n))

        else

          write(6,*) 'Error:  do not recognize stress point in ',
     >                           'add_km.'
          return

        end if    

      enddo
      enddo
      enddo
      enddo

      return
      end

*------------------------------------------------------------------------

      ive_ptr function bmix(rprndl,ipert,ndims,dims,stag,min,max,
     >                        miss,dunits,disunits,dim_names)

#include "default.icl"

      integer ipert
      real rprndl 

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer ijk,nt
      real brl,b0rl,r_val
      logical new
      ive_ptr getmem,b,b0,km       

      bmix = 0

      nt = ntime

      km = getvar('km',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                             dim_names,new)
      if (km.eq.0) then
        write(6,*) 'ERROR:  cannot access variable km.'
        return
      end if

      if (ipert.eq.1) then
        b0 = getvar('b0',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                                   dim_names,new)
        if (b0.eq.0) then
          write(6,*) 'ERROR:  cannot access variable b0.'
          return
        end if
      end if

      b = getvar('b',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                                dim_names,new)
      if (b.eq.0) then
        write(6,*) 'ERROR:  cannot access variable b.'
        return
      end if

      if (ipert.eq.1) call sub_b0(%val(b),%val(b0),nx,ny,nz,nt)

      bmix = getmem(nx*ny*nz*nt)
      if (bmix.eq.0) then
        write(6,*) 'ERROR:  cannot allocate space for bmix.'
        return
      end if

      do ijk=1,nx*ny*nz*nt
        call s_val(%val(bmix),ijk,0.)
      enddo

      call bmix_calc(%val(bmix),%val(b),%val(km),%val(dqdx_xy),
     >           %val(dqdx_zt),%val(dqdx_zw),%val(dqdy_xy),
     >           %val(dqdy_zt),%val(dqdy_zw),%val(dqdz),rdx,rdy,
     >           rprndl,nx,ny,nz,nt) 

      dunits = 'm/s/s/s'
      disunits = 'm/s/s/s' 
 
      if (ipert.eq.1) call add_b0(%val(b),%val(b0),nx,ny,nz,nt)

      return
      end

*-----------------------------------------------------------------------

      subroutine bmix_calc(bmix,b,km,dqdx_xy,dqdx_zt,dqdx_zw,dqdy_xy,
     >                       dqdy_zt,dqdy_zw,dqdz,rdx,rdy,rprndl,
     >                       nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real bmix(nx,ny,nz,nt),b(nx,ny,nz,nt),km(nx,ny,nz,nt)
      real rdx,rdy,rprndl

      real dqdx_xy(nx+1,ny+1,3),dqdx_zt(nz),dqdx_zw(nz+1)
      real dqdy_xy(nx+1,ny+1,3),dqdy_zt(nz),dqdy_zw(nz+1)
      real dqdz(nx+1,ny+1,3)

      integer i,j,k,n,kp1,km1
*
* X mixing (sans vertical term--done below)
*
      do n=1,nt
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do j=1,ny
      do i=2,nx-1

        bmix(i,j,k,n) = bmix(i,j,k,n) + 0.5*rdx*rdx*rprndl *
     >      ( (km(i+1,j,k,n)+km(i,j,k,n))*(b(i+1,j,k,n)-b(i,j,k,n)) 
     >      - (km(i,j,k,n)+km(i-1,j,k,n))*(b(i,j,k,n)-b(i-1,j,k,n)) )

        if (k.eq.1) then         

          bmix(i,j,k,n) = bmix(i,j,k,n) + 0.125*rdx*rprndl*(
     >
     >          (km(i+1,j,k,n)+km(i,j,k,n))*dqdx_xy(i+1,j,2)
     >      *( rdze(k+1)*dqdx_zw(k+1)*(b(i,j,kp1,n)
     >            -b(i,j,k,n)+b(i+1,j,kp1,n)-b(i+1,j,k,n)) )
     >
     >      -   (km(i,j,k,n)+km(i-1,j,k,n))*dqdx_xy(i,j,2)
     >      *( rdze(k+1)*dqdx_zw(k+1)*(b(i,j,kp1,n)
     >            -b(i,j,k,n)+b(i-1,j,kp1,n)-b(i-1,j,k,n)) )
     >
     >                         )

          bmix(i,j,k,n) = bmix(i,j,k,n) +
     >           0.25*rprndl*rdx*dqdx_xy(i,j,1)  
     >      *( rdze(k+1)*dqdx_zw(k+1)*zrbnd(k)*
     >          (km(i,j,kp1,n)*(b(i+1,j,kp1,n)-b(i-1,j,kp1,n))
     >             - km(i,j,k,n)*(b(i+1,j,k,n)-b(i-1,j,k,n)))
     >       + rdze(k)*dqdx_zw(k)*zlbnd(k)*
     >          (km(i,j,k,n)*(b(i+1,j,k,n)-b(i-1,j,k,n))
     >             - km(i,j,km1,n)*(b(i+1,j,km1,n)-b(i-1,j,km1,n))) )

        else

          bmix(i,j,k,n) = bmix(i,j,k,n) + 0.125*rdx*rprndl*(
     >
     >          (km(i+1,j,k,n)+km(i,j,k,n))*dqdx_xy(i+1,j,2)
     >      *( rdze(k+1)*dqdx_zw(k+1)*zrbnd(k)*(b(i,j,kp1,n)
     >            -b(i,j,k,n)+b(i+1,j,kp1,n)-b(i+1,j,k,n))
     >         + rdze(k)*dqdx_zw(k)*zlbnd(k)*(b(i,j,k,n)
     >            -b(i,j,km1,n)+b(i+1,j,k,n)-b(i+1,j,km1,n)) )
     >
     >      -   (km(i,j,k,n)+km(i-1,j,k,n))*dqdx_xy(i,j,2)
     >      *( rdze(k+1)*dqdx_zw(k+1)*zrbnd(k)*(b(i,j,kp1,n)
     >            -b(i,j,k,n)+b(i-1,j,kp1,n)-b(i-1,j,k,n))
     >         + rdze(k)*dqdx_zw(k)*zlbnd(k)*(b(i,j,k,n)
     >            -b(i,j,km1,n)+b(i-1,j,k,n)-b(i-1,j,km1,n)) )
     >
     >                         )

          bmix(i,j,k,n) = bmix(i,j,k,n) +
     >           0.25*rprndl*rdx*dqdx_xy(i,j,1)  
     >      *( rdze(k+1)*dqdx_zw(k+1)*zrbnd(k)*
     >          (km(i,j,kp1,n)*(b(i+1,j,kp1,n)-b(i-1,j,kp1,n))
     >             - km(i,j,k,n)*(b(i+1,j,k,n)-b(i-1,j,k,n)))
     >       + rdze(k)*dqdx_zw(k)*zlbnd(k)*
     >          (km(i,j,k,n)*(b(i+1,j,k,n)-b(i-1,j,k,n))
     >             - km(i,j,km1,n)*(b(i+1,j,km1,n)-b(i-1,j,km1,n))) )
  
        end if

      enddo
      enddo
      enddo
      enddo
*
* Y mixing (sans vertical term--done below)
*
      do n=1,nt
      do k=1,nz
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do j=2,ny-1
      do i=1,nx

        bmix(i,j,k,n) = bmix(i,j,k,n) + 0.5*rdy*rdy*rprndl *
     >      ( (km(i,j+1,k,n)+km(i,j,k,n))*(b(i,j+1,k,n)-b(i,j,k,n)) 
     >      - (km(i,j,k,n)+km(i,j-1,k,n))*(b(i,j,k,n)-b(i,j-1,k,n)) ) 

        if (k.eq.1) then

          bmix(i,j,k,n) = bmix(i,j,k,n) + 0.125*rdy*rprndl*(
     >
     >          (km(i,j+1,k,n)+km(i,j,k,n))*dqdy_xy(i,j+1,3)
     >      *( rdze(k+1)*dqdy_zw(k+1)*(b(i,j,kp1,n)
     >            -b(i,j,k,n)+b(i,j+1,kp1,n)-b(i,j+1,k,n)) )
     >
     >      -   (km(i,j,k,n)+km(i,j-1,k,n))*dqdy_xy(i,j,3)
     >      *( rdze(k+1)*dqdy_zw(k+1)*(b(i,j,kp1,n)
     >            -b(i,j,k,n)+b(i,j-1,kp1,n)-b(i,j-1,k,n)) )
     >
     >                           )

          bmix(i,j,k,n) = bmix(i,j,k,n) + 
     >           0.25*rprndl*rdy*dqdy_xy(i,j,1)
     >      *( rdze(k+1)*dqdy_zw(k+1)*zrbnd(k)*
     >          (km(i,j,kp1,n)*(b(i,j+1,kp1,n)-b(i,j-1,kp1,n))
     >             - km(i,j,k,n)*(b(i,j+1,k,n)-b(i,j-1,k,n)))
     >       + rdze(k)*dqdy_zw(k)*zlbnd(k)*
     >          (km(i,j,k,n)*(b(i,j+1,k,n)-b(i,j-1,k,n))
     >             - km(i,j,km1,n)*(b(i,j+1,km1,n)-b(i,j-1,km1,n))) )

        else

          bmix(i,j,k,n) = bmix(i,j,k,n) + 0.125*rdy*rprndl*(
     >
     >          (km(i,j+1,k,n)+km(i,j,k,n))*dqdy_xy(i,j+1,3)
     >      *( rdze(k+1)*dqdy_zw(k+1)*zrbnd(k)*(b(i,j,kp1,n)
     >            -b(i,j,k,n)+b(i,j+1,kp1,n)-b(i,j+1,k,n))
     >         + rdze(k)*dqdy_zw(k)*zlbnd(k)*(b(i,j,k,n)
     >            -b(i,j,km1,n)+b(i,j+1,k,n)-b(i,j+1,km1,n)) ) 
     >
     >      -   (km(i,j,k,n)+km(i,j-1,k,n))*dqdy_xy(i,j,3)
     >      *( rdze(k+1)*dqdy_zw(k+1)*zrbnd(k)*(b(i,j,kp1,n)
     >            -b(i,j,k,n)+b(i,j-1,kp1,n)-b(i,j-1,k,n))
     >         + rdze(k)*dqdy_zw(k)*zlbnd(k)*(b(i,j,k,n)
     >            -b(i,j,km1,n)+b(i,j-1,k,n)-b(i,j-1,km1,n)) )
     >
     >                           )

          bmix(i,j,k,n) = bmix(i,j,k,n) + 
     >           0.25*rprndl*rdy*dqdy_xy(i,j,1)
     >      *( rdze(k+1)*dqdy_zw(k+1)*zrbnd(k)*
     >          (km(i,j,kp1,n)*(b(i,j+1,kp1,n)-b(i,j-1,kp1,n))
     >             - km(i,j,k,n)*(b(i,j+1,k,n)-b(i,j-1,k,n)))
     >       + rdze(k)*dqdy_zw(k)*zlbnd(k)*
     >          (km(i,j,k,n)*(b(i,j+1,k,n)-b(i,j-1,k,n))
     >             - km(i,j,km1,n)*(b(i,j+1,km1,n)-b(i,j-1,km1,n))) )

        end if

      enddo
      enddo
      enddo
      enddo
*
* Z mixing (also vertical terms from horizontal mixing)
*
      do n=1,nt
      do k=1,nz-1
        kp1 = min0(k+1,nz)
        km1 = max0(k-1,1)
      do j=1,ny
      do i=1,nx

        if (k.eq.1) then

          bmix(i,j,k,n) = bmix(i,j,k,n) + 
     >                   0.5*rprndl*rdzc(k)*
     >             dqdx_zt(k)*dqdx_xy(i,j,1)*dqdx_xy(i,j,1)
     >        *( rdze(k+1)*dqdx_zw(k+1)*(km(i,j,kp1,n)+km(i,j,k,n))
     >                     *(b(i,j,kp1,n)-b(i,j,k,n)) )

          bmix(i,j,k,n) = bmix(i,j,k,n) + 
     >                   0.5*rprndl*rdzc(k)*
     >             dqdy_zt(k)*dqdy_xy(i,j,1)*dqdy_xy(i,j,1)
     >        *( rdze(k+1)*dqdy_zw(k+1)*(km(i,j,kp1,n)+km(i,j,k,n))
     >                     *(b(i,j,kp1,n)-b(i,j,k,n)) )

          bmix(i,j,k,n) = bmix(i,j,k,n) + 
     >           0.5*dqdz(i,j,1)*dqdz(i,j,1)*rprndl*rdzc(k)
     >     * rdze(k+1)*(km(i,j,kp1,n)+km(i,j,k,n))
     >                     *(b(i,j,kp1,n)-b(i,j,k,n))       

        else
   
          bmix(i,j,k,n) = bmix(i,j,k,n) + 
     >                   0.5*rprndl*rdzc(k)*
     >             dqdx_zt(k)*dqdx_xy(i,j,1)*dqdx_xy(i,j,1)
     >        *( rdze(k+1)*dqdx_zw(k+1)*(km(i,j,kp1,n)+km(i,j,k,n))
     >                     *(b(i,j,kp1,n)-b(i,j,k,n))
     >           - rdze(k)*dqdx_zw(k)*(km(i,j,k,n)+km(i,j,km1,n))
     >                     *(b(i,j,k,n)-b(i,j,km1,n)) )

          bmix(i,j,k,n) = bmix(i,j,k,n) + 
     >                   0.5*rprndl*rdzc(k)*
     >             dqdy_zt(k)*dqdy_xy(i,j,1)*dqdy_xy(i,j,1)
     >        *( rdze(k+1)*dqdy_zw(k+1)*(km(i,j,kp1,n)+km(i,j,k,n))
     >                     *(b(i,j,kp1,n)-b(i,j,k,n))
     >           - rdze(k)*dqdy_zw(k)*(km(i,j,k,n)+km(i,j,km1,n))
     >                     *(b(i,j,k,n)-b(i,j,km1,n)) )

          bmix(i,j,k,n) = bmix(i,j,k,n) + 
     >           0.5*dqdz(i,j,1)*dqdz(i,j,1)*rprndl*rdzc(k)
     >    *( rdze(k+1)*(km(i,j,kp1,n)+km(i,j,k,n))
     >                     *(b(i,j,kp1,n)-b(i,j,k,n))
     >       - rdze(k)*(km(i,j,k,n)+km(i,j,km1,n))
     >                     *(b(i,j,k,n)-b(i,j,km1,n)) )

        end if

      enddo
      enddo
      enddo
      enddo

      return 
      end

*------------------------------------------------------------------------

      subroutine sub_b0(b,b0,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real b(nx,ny,nz,nt),b0(nx,ny,nz)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx
        b(i,j,k,n) = b(i,j,k,n) - b0(i,j,k)
      enddo
      enddo
      enddo
      enddo

      return
      end

*--------------------------------------------------------------------------

      subroutine add_b0(b,b0,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real b(nx,ny,nz,nt),b0(nx,ny,nz)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx
        b(i,j,k,n) = b(i,j,k,n) + b0(i,j,k)
      enddo
      enddo
      enddo
      enddo

      return
      end

