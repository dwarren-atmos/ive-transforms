
#include "pointer.icl"

      ive_ptr function pvsrc(whch,rprndl,ipert,ndims,dims,stag,min,
     >                     max,miss,dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer ipert
      real rprndl
      character*(1) whch

      integer nt
      ive_ptr pvsrc_f,pvsrc_h,pvsrc1,pvsrc2,getmem

      nt = ntime

      pvsrc = 0
      pvsrc = getmem(nx*ny*nz*nt)
      if(pvsrc.eq.0) then
        write(6,*) 'YOU FOOL:  You need memory'
        pvsrc=0
        return
      endif

      pvsrc1 = 0
      pvsrc1 = getmem(nx*ny*nz*nt)
      if(pvsrc1.eq.0) then
        write(6,*) 'YOU FOOL:  You need memory'
        pvsrc1=0
        return
      endif

      pvsrc2 = 0
      pvsrc2 = getmem(nx*ny*nz*nt)
      if(pvsrc2.eq.0) then
        write(6,*) 'YOU FOOL:  You need memory'
        pvsrc2=0
        return
      endif

*
* Viscous part
*
      if (whch.eq.'F') then
        
        pvsrc = pvsrc_f(ndims,dims,stag,min,max,miss,dunits,disunits,
     >                                 dim_names)
*
* Thermal part
*
      else if (whch.eq.'H') then

        pvsrc = pvsrc_h(rprndl,ipert,ndims,dims,stag,min,max,miss,
     >                        dunits,disunits,dim_names)
*
* Total (viscous + thermal)
*
      else if (whch.eq.'T') then

        pvsrc1 = pvsrc_f(ndims,dims,stag,min,max,miss,dunits,disunits,
     >                                 dim_names)
        if (pvsrc1.eq.0) then
          write(6,*) 'Error:  cannot compute viscous source.'
          return
        end if

        pvsrc2 = pvsrc_h(rprndl,ipert,ndims,dims,stag,min,max,miss,
     >                        dunits,disunits,dim_names)
        if (pvsrc2.eq.0) then
          call freemem(pvsrc1)
          write(6,*) 'Error:  cannot compute thermal source.' 
          return
        end if

        pvsrc = getmem(nx*ny*nz*nt)
        if (pvsrc.eq.0) then
          call freemem(pvsrc1)
          call freemem(pvsrc2)
          write(6,*) 'Error:  unable to allocate space for pvsrc.'
          return
        end if

        call pvsrc_t_calc(%val(pvsrc),%val(pvsrc1),%val(pvsrc2),
     >                             nx,ny,nz,nt)

        call freemem(pvsrc1)
        call freemem(pvsrc2)

      else

        write(6,*) 'Error:  do not recognize ',whch(1:1),
     >                      ' part of PV source.'
        return

      end if

      return
      end

*-----------------------------------------------------------------------

      ive_ptr function pvsrc_f(ndims,dims,stag,min,max,miss,dunits,
     >                             disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer ijk,nt
      real r_val,tmpval
      logical new,new1,new2
      ive_ptr getmem,fx,fy,fz,b,tmp1,tmp2

      nt = ntime
      pvsrc_f = 0

      tmp1 = getvar('UMIX',ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new1)
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

      if (new1) call driver('free = UMIX',0)
      if (new2) call driver('free = FILT[U]',0)

      tmp1 = getvar('VMIX',ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new1)
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

      if (new1) call driver('free = VMIX',0)
      if (new2) call driver('free = FILT[V]',0)

      tmp1 = getvar('WMIX',ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new1)
      if (tmp1.eq.0) then
        write(6,*) 'Error:  cannot access wmix.'
        return
      end if

      tmp2 = getvar('FILT[W]',ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names,new2)
      if (tmp2.eq.0) then
        write(6,*) 'Error:  cannot access filt[w].'
        return
      end if

      fz = getmem(nx*ny*(nz+1)*nt)
      if (fz.eq.0) then
        write(6,*) 'Error:  could not allocate space for fz.'
        return
      end if

      do ijk=1,nx*ny*(nz+1)*nt
        tmpval = r_val(%val(tmp1),ijk) + r_val(%val(tmp2),ijk)
        call s_val(%val(fz),ijk,tmpval)
      enddo

      call f_extrp(%val(fz),nx,ny,nz,nt,0,0,1)

      if (new1) call driver('free = WMIX',0)
      if (new2) call driver('free = FILT[W]',0)

      b = getvar('b',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                          dim_names,new)
      if (b.eq.0) then
        call freemem(fx)
        call freemem(fy)
        call freemem(fz)
        write(6,*) 'Error:  cannot access b.'
        return
      end if

      pvsrc_f = getmem(nx*ny*nz*nt)
      if (pvsrc_f.eq.0) then
        call freemem(fx)
        call freemem(fy)
        call freemem(fz)
        write(6,*) 'Error:  cannot allocate space for pvsrc_f.'
        return
      end if

      call pvsrc_f_calc(%val(pvsrc_f),%val(fx),%val(fy),%val(fz),
     >              %val(b),%val(dqdx_xy),%val(dqdx_zw),
     >              %val(dqdx_zt),%val(dqdy_xy),%val(dqdy_zw),
     >              %val(dqdy_zt),%val(dqdz),rdx,rdy,nx,ny,nz,nt)

      dunits = '1/s/s/s/s'
      disunits = '1/s/s/s/s'
  
      call freemem(fx)
      call freemem(fy)
      call freemem(fz) 
        
      return
      end

*-----------------------------------------------------------------------

      subroutine pvsrc_f_calc(pvsrc,fx,fy,fz,b,dqdx_xy,dqdx_zw,
     >                dqdx_zt,dqdy_xy,dqdy_zw,dqdy_zt,dqdz,rdx,
     >                rdy,nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real pvsrc(nx,ny,nz,nt),fx(nx+1,ny,nz,nt),fy(nx,ny+1,nz,nt)
      real fz(nx,ny,nz+1,nt),b(nx,ny,nz,nt),rdx,rdy

      real dqdx_xy(nx+1,ny+1,3),dqdx_zw(nz+1),dqdx_zt(nz)
      real dqdy_xy(nx+1,ny+1,3),dqdy_zw(nz+1),dqdy_zt(nz)
      real dqdz(nx+1,ny+1,3)

      integer i,j,k,n,m,ip1,im1,jp1,jm1,kp1,km1
      real curlx,curly,curlz,dbdx,dbdy,dbdz

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

        curlx =               0.25*rdy
     >       *( yrbnd(j)*(fz(i,jp1,k,n)-fz(i,j,k,n)
     >                        + fz(i,jp1,k+1,n)-fz(i,j,k+1,n))
     >        + ylbnd(j)*(fz(i,j,k,n)-fz(i,jm1,k,n)
     >                        + fz(i,j,k+1,n)-fz(i,jm1,k+1,n)) )
     >                + dqdy_xy(i,j,1)*dqdy_zt(k)*rdzc(k)
     >                        * (fz(i,j,k+1,n)-fz(i,j,k,n))
     >                - 0.25*dqdz(i,j,1)
     >       *( rdze(k+1)*zrbnd(k)*(fy(i,j,kp1,n)-fy(i,j,k,n)
     >              + fy(i,j+1,kp1,n)-fy(i,j+1,k,n))
     >          + rdze(k)*zlbnd(k)*(fy(i,j,k,n)-fy(i,j,km1,n)
     >              + fy(i,j+1,k,n)-fy(i,j+1,km1,n)) )

        curly =             - 0.25*rdx
     >       *( xrbnd(i)*(fz(ip1,j,k,n)-fz(i,j,k,n)
     >                        + fz(ip1,j,k+1,n)-fz(i,j,k+1,n))
     >        + xlbnd(i)*(fz(i,j,k,n)-fz(im1,j,k,n)
     >                        + fz(i,j,k+1,n)-fz(im1,j,k+1,n)) )
     >                - dqdx_xy(i,j,1)*dqdx_zt(k)*rdzc(k)
     >                          * (fz(i,j,k+1,n)-fz(i,j,k,n))
     >                + 0.25*dqdz(i,j,1)
     >       *( rdze(k+1)*zrbnd(k)*(fx(i,j,kp1,n)-fx(i,j,k,n)
     >              + fx(i+1,j,kp1,n)-fx(i+1,j,k,n))
     >          + rdze(k)*zlbnd(k)*(fx(i,j,k,n)-fx(i,j,km1,n)
     >              + fx(i+1,j,k,n)-fx(i+1,j,km1,n)) )

        curlz =             - 0.25*rdy
     >       *( yrbnd(j)*(fx(i,jp1,k,n)-fx(i,j,k,n)
     >                        + fx(i+1,jp1,k,n)-fx(i+1,j,k,n))
     >        + ylbnd(j)*(fx(i,j,k,n)-fx(i,jm1,k,n)
     >                        + fx(i+1,j,k,n)-fx(i+1,jm1,k,n)) )
     >                - 0.25*dqdy_xy(i,j,1)
     >       *( rdze(k+1)*dqdy_zw(k+1)*zrbnd(k)*(fx(i,j,kp1,n)
     >             - fx(i,j,k,n)+fx(i+1,j,kp1,n)-fx(i+1,j,k,n))
     >          + rdze(k)*dqdy_zw(k)*zlbnd(k)*(fx(i,j,k,n)
     >             - fx(i,j,km1,n)+fx(i+1,j,k,n)-fx(i+1,j,km1,n)) ) 
     >                     + 0.25*rdx
     >       *( xrbnd(i)*(fy(ip1,j,k,n)-fy(i,j,k,n)
     >                        + fy(ip1,j+1,k,n)-fy(i,j+1,k,n))
     >        + xlbnd(i)*(fy(i,j,k,n)-fy(im1,j,k,n)
     >                        + fy(i,j+1,k,n)-fy(im1,j+1,k,n)) )
     >                + 0.25*dqdx_xy(i,j,1)
     >       *( rdze(k+1)*dqdx_zw(k+1)*zrbnd(k)*(fy(i,j,kp1,n)
     >             - fy(i,j,k,n)+fy(i,j+1,kp1,n)-fy(i,j+1,k,n))
     >          + rdze(k)*dqdx_zw(k)*zlbnd(k)*(fy(i,j,k,n)
     >             - fy(i,j,km1,n)+fy(i,j+1,k,n)-fy(i,j+1,km1,n)) )

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

        dbdz = 0.5*dqdz(i,j,1)
     >       *( zrbnd(k)*rdze(k+1)*(b(i,j,kp1,n)-b(i,j,k,n))
     >        + zlbnd(k)*rdze(k)*(b(i,j,k,n)-b(i,j,km1,n)) )

        pvsrc(i,j,k,n) = dbdx*curlx + dbdy*curly + dbdz*curlz

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

        curlx =               0.25*rdy
     >       *( yrbnd(j)*(fz(i,jp1,k,n)-fz(i,j,k,n)
     >                        + fz(i,jp1,k+1,n)-fz(i,j,k+1,n))
     >        + ylbnd(j)*(fz(i,j,k,n)-fz(i,jm1,k,n)
     >                        + fz(i,j,k+1,n)-fz(i,jm1,k+1,n)) )
     >                + dqdy_xy(i,j,1)*dqdy_zt(k)*rdzc(k)
     >                        * (fz(i,j,k+1,n)-fz(i,j,k,n))
     >                - 0.25*dqdz(i,j,1)
     >       *( rdze(k+1)*zrbnd(k)*(fy(i,j,kp1,n)-fy(i,j,k,n)
     >              + fy(i,j+1,kp1,n)-fy(i,j+1,k,n))
     >          + rdze(k)*zlbnd(k)*(fy(i,j,k,n)-fy(i,j,km1,n)
     >              + fy(i,j+1,k,n)-fy(i,j+1,km1,n)) )

        curly =             - 0.25*rdx
     >       *( xrbnd(i)*(fz(ip1,j,k,n)-fz(i,j,k,n)
     >                        + fz(ip1,j,k+1,n)-fz(i,j,k+1,n))
     >        + xlbnd(i)*(fz(i,j,k,n)-fz(im1,j,k,n)
     >                        + fz(i,j,k+1,n)-fz(im1,j,k+1,n)) )
     >                - dqdx_xy(i,j,1)*dqdx_zt(k)*rdzc(k)
     >                          * (fz(i,j,k+1,n)-fz(i,j,k,n))
     >                + 0.25*dqdz(i,j,1)
     >       *( rdze(k+1)*zrbnd(k)*(fx(i,j,kp1,n)-fx(i,j,k,n)
     >              + fx(i+1,j,kp1,n)-fx(i+1,j,k,n))
     >          + rdze(k)*zlbnd(k)*(fx(i,j,k,n)-fx(i,j,km1,n)
     >              + fx(i+1,j,k,n)-fx(i+1,j,km1,n)) )

        curlz =             - 0.25*rdy
     >       *( yrbnd(j)*(fx(i,jp1,k,n)-fx(i,j,k,n)
     >                        + fx(i+1,jp1,k,n)-fx(i+1,j,k,n))
     >        + ylbnd(j)*(fx(i,j,k,n)-fx(i,jm1,k,n)
     >                        + fx(i+1,j,k,n)-fx(i+1,jm1,k,n)) )
     >                - 0.25*dqdy_xy(i,j,1)
     >       *( rdze(k+1)*dqdy_zw(k+1)*zrbnd(k)*(fx(i,j,kp1,n)
     >             - fx(i,j,k,n)+fx(i+1,j,kp1,n)-fx(i+1,j,k,n))
     >          + rdze(k)*dqdy_zw(k)*zlbnd(k)*(fx(i,j,k,n)
     >             - fx(i,j,km1,n)+fx(i+1,j,k,n)-fx(i+1,j,km1,n)) ) 
     >                     + 0.25*rdx
     >       *( xrbnd(i)*(fy(ip1,j,k,n)-fy(i,j,k,n)
     >                        + fy(ip1,j+1,k,n)-fy(i,j+1,k,n))
     >        + xlbnd(i)*(fy(i,j,k,n)-fy(im1,j,k,n)
     >                        + fy(i,j+1,k,n)-fy(im1,j+1,k,n)) )
     >                + 0.25*dqdx_xy(i,j,1)
     >       *( rdze(k+1)*dqdx_zw(k+1)*zrbnd(k)*(fy(i,j,kp1,n)
     >             - fy(i,j,k,n)+fy(i,j+1,kp1,n)-fy(i,j+1,k,n))
     >          + rdze(k)*dqdx_zw(k)*zlbnd(k)*(fy(i,j,k,n)
     >             - fy(i,j,km1,n)+fy(i,j+1,k,n)-fy(i,j+1,km1,n)) )

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

        dbdz = 0.5*dqdz(i,j,1)
     >       *( zrbnd(k)*rdze(k+1)*(b(i,j,kp1,n)-b(i,j,k,n))
     >        + zlbnd(k)*rdze(k)*(b(i,j,k,n)-b(i,j,km1,n)) )

        pvsrc(i,j,k,n) = dbdx*curlx + dbdy*curly + dbdz*curlz

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

        curlx =      
     >            + 1./3.*rdy*(fz(i,j+1,k,n)-fz(i,j-1,k,n)
     >                      + fz(i,j+1,k+1,n)-fz(i,j-1,k+1,n))
     >            - 1./24.*rdy*(fz(i,j+2,k,n)-fz(i,j-2,k,n)
     >                      + fz(i,j+2,k+1,n)-fz(i,j-2,k+1,n))
     >                + dqdy_xy(i,j,1)*dqdy_zt(k)*rdzc(k)    
     >                        * (fz(i,j,k+1,n)-fz(i,j,k,n))
     >                - 0.25*dqdz(i,j,1)
     >       *( rdze(k+1)*zrbnd(k)*(fy(i,j,kp1,n)-fy(i,j,k,n)
     >              + fy(i,j+1,kp1,n)-fy(i,j+1,k,n))
     >          + rdze(k)*zlbnd(k)*(fy(i,j,k,n)-fy(i,j,km1,n)
     >              + fy(i,j+1,k,n)-fy(i,j+1,km1,n)) )

        curly =     
     >            - 1./3.*rdx*(fz(i+1,j,k,n)-fz(i-1,j,k,n)
     >                      + fz(i+1,j,k+1,n)-fz(i-1,j,k+1,n))
     >            + 1./24.*rdx*(fz(i+2,j,k,n)-fz(i-2,j,k,n)
     >                      + fz(i+2,j,k+1,n)-fz(i-2,j,k+1,n))
     >                - dqdx_xy(i,j,1)*dqdx_zt(k)*rdzc(k)
     >                        * (fz(i,j,k+1,n)-fz(i,j,k,n))
     >                + 0.25*dqdz(i,j,1)
     >       *( rdze(k+1)*zrbnd(k)*(fx(i,j,kp1,n)-fx(i,j,k,n)
     >              + fx(i+1,j,kp1,n)-fx(i+1,j,k,n))
     >          + rdze(k)*zlbnd(k)*(fx(i,j,k,n)-fx(i,j,km1,n)
     >              + fx(i+1,j,k,n)-fx(i+1,j,km1,n)) )

        curlz = 
     >            - 1./3.*rdy*(fx(i,j+1,k,n)-fx(i,j-1,k,n)
     >                      + fx(i+1,j+1,k,n)-fx(i+1,j-1,k,n))
     >            + 1./24.*rdy*(fx(i,j+2,k,n)-fx(i,j-2,k,n)
     >                      + fx(i+1,j+2,k,n)-fx(i+1,j-2,k,n))
     >                - 0.25*dqdy_xy(i,j,1)
     >       *( rdze(k+1)*dqdy_zw(k+1)*zrbnd(k)*(fx(i,j,kp1,n)
     >             - fx(i,j,k,n)+fx(i+1,j,kp1,n)-fx(i+1,j,k,n))
     >          + rdze(k)*dqdy_zw(k)*zlbnd(k)*(fx(i,j,k,n)
     >             - fx(i,j,km1,n)+fx(i+1,j,k,n)-fx(i+1,j,km1,n)) ) 
     >            + 1./3.*rdx*(fy(i+1,j,k,n)-fy(i-1,j,k,n)
     >                      + fy(i+1,j+1,k,n)-fy(i-1,j+1,k,n))
     >            - 1./24.*rdx*(fy(i+2,j,k,n)-fy(i-2,j,k,n)
     >                      + fy(i+2,j+1,k,n)-fy(i-2,j+1,k,n)) 
     >                + 0.25*dqdx_xy(i,j,1)
     >       *( rdze(k+1)*dqdx_zw(k+1)*zrbnd(k)*(fy(i,j,kp1,n)
     >             - fy(i,j,k,n)+fy(i,j+1,kp1,n)-fy(i,j+1,k,n))
     >          + rdze(k)*dqdx_zw(k)*zlbnd(k)*(fy(i,j,k,n)
     >             - fy(i,j,km1,n)+fy(i,j+1,k,n)-fy(i,j+1,km1,n)) )

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

        dbdz = 0.5*dqdz(i,j,1)
     >       *( zrbnd(k)*rdze(k+1)*(b(i,j,kp1,n)-b(i,j,k,n))
     >        + zlbnd(k)*rdze(k)*(b(i,j,k,n)-b(i,j,km1,n)) )

        pvsrc(i,j,k,n) = dbdx*curlx + dbdy*curly + dbdz*curlz

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

      ive_ptr function pvsrc_h(rprndl,ipert,ndims,dims,stag,min,max,
     >                       miss,dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer ipert
      real rprndl

      integer ijk,nt,iarg(1),strbeg,strend
      real r_val,tmpval,rarg(1)
      character*(80) mkfield,comstr,freestr,carg(1)
      logical new,new1,new2
      ive_ptr getmem,fb,q,r,s,tmp1,tmp2

      pvsrc_h = 0
      nt = ntime

      comstr(1:len(comstr)) = ' '
      iarg(1) = ipert
      rarg(1) = rprndl
      comstr = mkfield('BMIX',carg,iarg,rarg,0,1,1,0)

      WRITE(*,*) 'Calling getvar with ',
     >                  comstr(strbeg(comstr):strend(comstr))

      tmp1 = getvar(comstr,ndims,dims,stag,min,max,miss,dunits,
     >                    disunits,dim_names,new1)
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
        write(6,*) 'Error:  cannot allocate space for fb.'
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

      if (new1) call driver(freestr,0)
      if (new2) call driver('free = FILT[B]',0)

      q = getvar('Q_IVE',ndims,dims,stag,min,max,miss,dunits,
     >                   disunits,dim_names,new)
      if (q.eq.0) then
        call freemem(fb)
        write(6,*) 'Error:  cannot access q_ive.'
        return
      end if

      r = getvar('R_IVE',ndims,dims,stag,min,max,miss,dunits,
     >                   disunits,dim_names,new)
      if (r.eq.0) then
        call freemem(fb) 
        write(6,*) 'Error:  cannot access r_ive.'
        return
      end if

      s = getvar('S_IVE',ndims,dims,stag,min,max,miss,dunits,
     >                   disunits,dim_names,new)
      if (s.eq.0) then
        call freemem(fb)
        write(6,*) 'Error:  cannot access s_ive.'
        return
      end if

      pvsrc_h = getmem(nx*ny*nz*nt)
      if (pvsrc_h.eq.0) then
        call freemem(fb)
        write(6,*) 'Error:  cannot allocate space for pvsrc_h.'
        return
      end if

      call pvsrc_h_calc(%val(pvsrc_h),%val(fb),%val(q),%val(r),
     >                  %val(s),%val(dqdx_xy),%val(dqdx_zw),
     >                  %val(dqdy_xy),%val(dqdy_zw),%val(dqdz),
     >                  rdx,rdy,nx,ny,nz,nt)

      dunits = '1/s/s/s/s'
      disunits = '1/s/s/s/s'
 
      call freemem(fb)

      return
      end

*-----------------------------------------------------------------------

      subroutine pvsrc_h_calc(pvsrc,fb,q,r,s,dqdx_xy,dqdx_zw,
     >               dqdy_xy,dqdy_zw,dqdz,rdx,rdy,nx,ny,nz,nt)

#include "misc.icl"

      integer nx,ny,nz,nt
      real pvsrc(nx,ny,nz,nt),fb(nx,ny,nz,nt),q(nx,ny,nz,nt)
      real r(nx,ny,nz,nt),s(nx,ny,nz,nt),rdx,rdy

      real dqdx_xy(nx+1,ny+1,3),dqdx_zw(nz+1)
      real dqdy_xy(nx+1,ny+1,3),dqdy_zw(nz+1)
      real dqdz(nx+1,ny+1,3)

      integer i,j,k,n,m,ip1,im1,jp1,jm1,kp1,km1
      real dfbdx,dfbdy,dfbdz
     
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

        dfbdx = 0.5*rdx*( xrbnd(i)*(fb(ip1,j,k,n)-fb(i,j,k,n))
     >                    + xlbnd(i)*(fb(i,j,k,n)-fb(im1,j,k,n)) )
     >                + 0.5*dqdx_xy(i,j,1)
     >         *( zrbnd(k)*rdze(k+1)*dqdx_zw(k+1)
     >                *(fb(i,j,kp1,n)-fb(i,j,k,n))
     >          + zlbnd(k)*rdze(k)*dqdx_zw(k)
     >                *(fb(i,j,k,n)-fb(i,j,km1,n)) )

        dfbdy = 0.5*rdy*( yrbnd(j)*(fb(i,jp1,k,n)-fb(i,j,k,n))
     >                    + ylbnd(j)*(fb(i,j,k,n)-fb(i,jm1,k,n)) )
     >                + 0.5*dqdy_xy(i,j,1)
     >         *( zrbnd(k)*rdze(k+1)*dqdy_zw(k+1)
     >                *(fb(i,j,kp1,n)-fb(i,j,k,n))
     >          + zlbnd(k)*rdze(k)*dqdy_zw(k)
     >                *(fb(i,j,k,n)-fb(i,j,km1,n)) )

        dfbdz = 0.5*dqdz(i,j,1)
     >       *( zrbnd(k)*rdze(k+1)*(fb(i,j,kp1,n)-fb(i,j,k,n))
     >        + zlbnd(k)*rdze(k)*(fb(i,j,k,n)-fb(i,j,km1,n)) )

        pvsrc(i,j,k,n) = q(i,j,k,n)*dfbdx + r(i,j,k,n)*dfbdy
     >                      + s(i,j,k,n)*dfbdz

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

        dfbdx = 0.5*rdx*( xrbnd(i)*(fb(ip1,j,k,n)-fb(i,j,k,n))
     >                    + xlbnd(i)*(fb(i,j,k,n)-fb(im1,j,k,n)) )
     >                + 0.5*dqdx_xy(i,j,1)
     >         *( zrbnd(k)*rdze(k+1)*dqdx_zw(k+1)
     >                *(fb(i,j,kp1,n)-fb(i,j,k,n))
     >          + zlbnd(k)*rdze(k)*dqdx_zw(k)
     >                *(fb(i,j,k,n)-fb(i,j,km1,n)) )

        dfbdy = 0.5*rdy*( yrbnd(j)*(fb(i,jp1,k,n)-fb(i,j,k,n))
     >                    + ylbnd(j)*(fb(i,j,k,n)-fb(i,jm1,k,n)) )
     >                + 0.5*dqdy_xy(i,j,1)
     >         *( zrbnd(k)*rdze(k+1)*dqdy_zw(k+1)
     >                *(fb(i,j,kp1,n)-fb(i,j,k,n))
     >          + zlbnd(k)*rdze(k)*dqdy_zw(k)
     >                *(fb(i,j,k,n)-fb(i,j,km1,n)) )

        dfbdz = 0.5*dqdz(i,j,1)
     >       *( zrbnd(k)*rdze(k+1)*(fb(i,j,kp1,n)-fb(i,j,k,n))
     >        + zlbnd(k)*rdze(k)*(fb(i,j,k,n)-fb(i,j,km1,n)) )

        pvsrc(i,j,k,n) = q(i,j,k,n)*dfbdx + r(i,j,k,n)*dfbdy
     >                      + s(i,j,k,n)*dfbdz

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

        dfbdx = 4./6.*rdx*(fb(i+1,j,k,n)-fb(i-1,j,k,n))
     >           - 1./12.*rdx*(fb(i+2,j,k,n)-fb(i-2,j,k,n))
     >               + 0.5*dqdx_xy(i,j,1) 
     >        *( zrbnd(k)*rdze(k+1)*dqdx_zw(k+1)
     >               *(fb(i,j,kp1,n)-fb(i,j,k,n))
     >         + zlbnd(k)*rdze(k)*dqdx_zw(k)
     >               *(fb(i,j,k,n)-fb(i,j,km1,n)) )

        dfbdy = 4./6.*rdy*(fb(i,j+1,k,n)-fb(i,j-1,k,n))
     >           - 1./12.*rdy*(fb(i,j+2,k,n)-fb(i,j-2,k,n))
     >               + 0.5*dqdy_xy(i,j,1)
     >        *( zrbnd(k)*rdze(k+1)*dqdy_zw(k+1)
     >               *(fb(i,j,kp1,n)-fb(i,j,k,n))
     >         + zlbnd(k)*rdze(k)*dqdy_zw(k)
     >               *(fb(i,j,k,n)-fb(i,j,km1,n)) )

        dfbdz = 0.5*dqdz(i,j,1)
     >       *( zrbnd(k)*rdze(k+1)*(fb(i,j,kp1,n)-fb(i,j,k,n))
     >        + zlbnd(k)*rdze(k)*(fb(i,j,k,n)-fb(i,j,km1,n)) )

        pvsrc(i,j,k,n) = q(i,j,k,n)*dfbdx + r(i,j,k,n)*dfbdy
     >                       + s(i,j,k,n)*dfbdz

      enddo
      enddo
      enddo
*
* end time loop
*        
      enddo

      return
      end

*-----------------------------------------------------------------------

      subroutine pvsrc_t_calc(pvsrc,pvsrc1,pvsrc2,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real pvsrc(nx,ny,nz,nt),pvsrc1(nx,ny,nz,nt),pvsrc2(nx,ny,nz,nt)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx
        pvsrc(i,j,k,n) = pvsrc1(i,j,k,n) + pvsrc2(i,j,k,n)
      enddo
      enddo
      enddo
      enddo

      return
      end
        
