
#include "pointer.icl"

      ive_ptr function deriv(varnm,dir,ifac,ndims,dims,stag,min,max,
     >                        miss,dunits,disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4),ifac
      real stag(4),min(4),max(4),miss
      character*(*) varnm,dir,dunits,disunits,dim_names(4)
      
      integer idir,iu,iv,iw,it,strbeg,strend,k,nt
      real r_val
      character*(80) tmp
      character*(2) dimu
      ive_ptr var,getmem
      logical new,err,requ

      deriv = 0

      nt = ntime

      if (dir(1:1).eq.'X') then
        idir = 1
      else if (dir(1:1).eq.'Y') then
        idir = 2
      else if (dir(1:1).eq.'Z') then
        idir = 3
      else if (dir(1:1).eq.'T') then
        idir = 4
      else
        write(6,*) 'Error:  do not know how to compute derivative ',
     >                'for direction ',dir(strbeg(dir):strend(dir))
        return
      end if

      if ((idir.eq.2).and.(twod.eq.1)) then
        write(6,*) 'Error:  cannot take y derivative for 2d run.'
        return
      end if

      var = getvar(varnm,ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names,new)
      if (var.eq.0) then
        write(6,*) 'Error:  cannot access variable ',
     >               varnm(strbeg(varnm):strend(varnm))
        return
      end if

      if (dims(idir).eq.1) then
        write(6,*) 'Error:  cannot take derivative of trivial ',
     >                 'dimension ',dir(strbeg(dir):strend(dir))
        return
      end if

      deriv = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (deriv.eq.0) then
        write(6,*) 'Error:  cannot allocate space for deriv.'
        return
      end if

      it = 0
      iu = 0
      iv = 0
      iw = 0
      if ((dims(1).eq.1).and.(dims(2).eq.1).and.(dims(3).eq.1)) then 
        it = 1
      else if (requ(plmin(1),min(1))) then
        iu = 1
      else if (requ(plmin(2),min(2)).and.(twod.ne.1)) then
        iv = 1
      else if (requ(plmin(3),min(3))) then
        iw = 1
      end if

      if (ifac.eq.1) then

        if (it.eq.1) then
          call diff1_t(%val(var),%val(deriv),%val(time),nt)
        else if (iu.eq.1) then
          call diff1_upts(%val(var),%val(deriv),%val(dqdx_xy),
     >              %val(dqdx_zw),%val(dqdy_xy),%val(dqdy_zw),
     >              %val(dqdz),dx,dy,%val(time),nx,ny,nz,nt,idir)
        else if (iv.eq.1) then
          call diff1_vpts(%val(var),%val(deriv),%val(dqdx_xy),
     >              %val(dqdx_zw),%val(dqdy_xy),%val(dqdy_zw),
     >              %val(dqdz),dx,dy,%val(time),nx,ny,nz,nt,idir)
        else if (iw.eq.1) then
          call diff1_wpts(%val(var),%val(deriv),%val(dqdx_xy),
     >              %val(dqdx_zt),%val(dqdy_xy),%val(dqdy_zt),
     >              %val(dqdz),dx,dy,%val(time),nx,ny,nz,nt,idir)
        else 
          call diff1_spts(%val(var),%val(deriv),%val(dqdx_xy),
     >              %val(dqdx_zw),%val(dqdy_xy),%val(dqdy_zw),
     >              %val(dqdz),dx,dy,%val(time),nx,ny,nz,nt,idir)
        end if

        if (idir.ne.4) then
          dimu = '/m'
        else
          dimu = '/s'
        end if

      else

        write(6,*) 'Error:  second derivative not yet implemented.'  
        call freemem(deriv)
        return

      end if

      tmp = dunits(strbeg(dunits):strend(dunits))//dimu
      call units_smpfy(tmp,dunits,err)

      tmp = disunits(strbeg(disunits):strend(disunits))//dimu
      call units_smpfy(tmp,disunits,err)

      return
      end

*--------------------------------------------------------------------

      subroutine diff1_upts(var,deriv,dqdx_xy,dqdx_zw,dqdy_xy,
     >                 dqdy_zw,dqdz,dx,dy,time,nx,ny,nz,nt,idir)

      include 'misc.icl'

      integer nx,ny,nz,nt,idir
      real var(nx+1,ny,nz,nt),deriv(nx+1,ny,nz,nt)
      real dqdx_xy(nx+1,ny+1,3),dqdx_zw(nz+1)
      real dqdy_xy(nx+1,ny+1,3),dqdy_zw(nz+1)
      real dqdz(nx+1,ny+1,3),dx,dy,time(nt)

      integer i,j,k,n,m,ip1,im1,jp1,jm1,kp1,km1
      real rdx,rdy,da,db,ta,tb,dta,dtb      

      rdx = 1./dx
      rdy = 1./dy
*
* x derivative
*
      if (idir.eq.1) then

        do n=1,nt
*
* east and west boundaries
*
        do k=1,nz
          kp1 = min0(k+1,nz)
          km1 = max0(k-1,1)
        do j=1,ny

          i = 1
          deriv(i,j,k,n) = rdx*(var(i+1,j,k,n)-var(i,j,k,n))
     >                   + 0.5*dqdx_xy(i,j,2)
     >          *( rdze(k+1)*dqdx_zw(k+1)*zrbnd(k)
     >                       * (var(i,j,kp1,n)-var(i,j,k,n))
     >             + rdze(k)*dqdx_zw(k)*zlbnd(k)
     >                       * (var(i,j,k,n)-var(i,j,km1,n)) )
          i = nx+1
          deriv(i,j,k,n) = rdx*(var(i,j,k,n)-var(i-1,j,k,n))
     >                   + 0.5*dqdx_xy(i,j,2)
     >          *( rdze(k+1)*dqdx_zw(k+1)*zrbnd(k)
     >                       * (var(i,j,kp1,n)-var(i,j,k,n))
     >             + rdze(k)*dqdx_zw(k)*zlbnd(k)
     >                       * (var(i,j,k,n)-var(i,j,km1,n)) )

        enddo
        enddo

        do k=1,nz
          kp1 = min0(k+1,nz)
          km1 = max0(k-1,1)
        do m=1,2
          if (m.eq.1) i = 2
          if (m.eq.2) i = nx
        do j=1,ny

          deriv(i,j,k,n) = 0.5*rdx*(var(i+1,j,k,n)-var(i-1,j,k,n))
     >                   + 0.5*dqdx_xy(i,j,2)
     >          *( rdze(k+1)*dqdx_zw(k+1)*zrbnd(k)
     >                       * (var(i,j,kp1,n)-var(i,j,k,n))
     >             + rdze(k)*dqdx_zw(k)*zlbnd(k)
     >                       * (var(i,j,k,n)-var(i,j,km1,n)) )

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
        do i=3,nx-1

          deriv(i,j,k,n) = 4./6.*rdx*(var(i+1,j,k,n)-var(i-1,j,k,n))
     >                 - 1./12.*rdx*(var(i+2,j,k,n)-var(i-2,j,k,n))
     >                   + 0.5*dqdx_xy(i,j,2)
     >          *( rdze(k+1)*dqdx_zw(k+1)*zrbnd(k)
     >                       * (var(i,j,kp1,n)-var(i,j,k,n))
     >             + rdze(k)*dqdx_zw(k)*zlbnd(k)
     >                       * (var(i,j,k,n)-var(i,j,km1,n)) )

        enddo
        enddo
        enddo
*
* end time loop
*
        enddo
*
* y derivative
*
      else if (idir.eq.2) then

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
        do i=1,nx+1
          
          deriv(i,j,k,n) = 0.5*rdy
     >       *( yrbnd(j)*(var(i,jp1,k,n)-var(i,j,k,n))
     >        + ylbnd(j)*(var(i,j,k,n)-var(i,jm1,k,n)) )
     >                   + 0.5*dqdy_xy(i,j,2)
     >          *( rdze(k+1)*dqdy_zw(k+1)*zrbnd(k)
     >                       * (var(i,j,kp1,n)-var(i,j,k,n))
     >             + rdze(k)*dqdy_zw(k)*zlbnd(k)
     >                       * (var(i,j,k,n)-var(i,j,km1,n)) )
        
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
        do i=1,nx+1

          deriv(i,j,k,n) = 4./6.*rdy*(var(i,j+1,k,n)-var(i,j-1,k,n))
     >                 - 1./12.*rdy*(var(i,j+2,k,n)-var(i,j-2,k,n))
     >                   + 0.5*dqdy_xy(i,j,2)
     >          *( rdze(k+1)*dqdy_zw(k+1)*zrbnd(k)
     >                       * (var(i,j,k+1,n)-var(i,j,k,n))
     >             + rdze(k)*dqdy_zw(k)*zlbnd(k)
     >                       * (var(i,j,k,n)-var(i,j,k-1,n)) )

        enddo
        enddo
        enddo
*
* end time loop
*
        enddo
*
* z derivative
*
      else if (idir.eq.3) then
        
        do n=1,nt
        do k=1,nz
          kp1 = min0(k+1,nz)
          km1 = max0(k-1,1)
        do j=1,ny
        do i=1,nx+1

          deriv(i,j,k,n) = 0.5*dqdz(i,j,2)
     >      *( rdze(k+1)*zrbnd(k)*(var(i,j,k+1,n)-var(i,j,k,n))
     >         + rdze(k)*zlbnd(k)*(var(i,j,k,n)-var(i,j,k-1,n)) )

        enddo
        enddo
        enddo
        enddo
*
* t derivative
*
      else if (idir.eq.4) then

        do k=1,nz
        do j=1,ny
        do i=1,nx+1 

          deriv(i,j,k,1) = (var(i,j,k,2)-var(i,j,k,1))
     >                        / (time(2)-time(1))
          deriv(i,j,k,nt) = (var(i,j,k,nt)-var(i,j,k,nt-1))
     >                        / (time(nt)-time(nt-1))

          do n=2,nt-1
            dta = time(n+1)-time(n)
            dtb = time(n)-time(n-1)
            deriv(i,j,k,n) = dtb/(dta+dtb)
     >               * (var(i,j,k,n+1)-var(i,j,k,n))/dta
     >                     + dta/(dta+dtb)
     >               * (var(i,j,k,n)-var(i,j,k,n-1))/dtb
          enddo

        enddo
        enddo
        enddo

      end if

      return
      end

*------------------------------------------------------------------------
        
      subroutine diff1_vpts(var,deriv,dqdx_xy,dqdx_zw,dqdy_xy,
     >                 dqdy_zw,dqdz,dx,dy,time,nx,ny,nz,nt,idir)

      include 'misc.icl'

      integer nx,ny,nz,nt,idir
      real var(nx,ny+1,nz,nt),deriv(nx,ny+1,nz,nt)
      real dqdx_xy(nx+1,ny+1,3),dqdx_zw(nz+1)
      real dqdy_xy(nx+1,ny+1,3),dqdy_zw(nz+1)
      real dqdz(nx+1,ny+1,3),dx,dy,time(nt)

      integer i,j,k,n,m,ip1,im1,jp1,jm1,kp1,km1
      real rdx,rdy,da,db,ta,tb,dta,dtb      

      rdx = 1./dx
      rdy = 1./dy
*
* x derivative
*
      if (idir.eq.1) then

        do n=1,nt
*
* west and east boundaries
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
        do j=1,ny+1
          
          deriv(i,j,k,n) = 0.5*rdx
     >       *( xrbnd(i)*(var(ip1,j,k,n)-var(i,j,k,n))
     >        + xlbnd(i)*(var(i,j,k,n)-var(im1,j,k,n)) )
     >                   + 0.5*dqdx_xy(i,j,3)
     >          *( rdze(k+1)*dqdx_zw(k+1)*zrbnd(k)
     >                       * (var(i,j,kp1,n)-var(i,j,k,n))
     >             + rdze(k)*dqdx_zw(k)*zlbnd(k)
     >                       * (var(i,j,k,n)-var(i,j,km1,n)) )
        
        enddo
        enddo
        enddo
*
* interior
*
        do k=1,nz
          kp1 = min0(k+1,nz)
          km1 = max0(k-1,1)
        do j=1,ny+1
        do i=3,nx-2

          deriv(i,j,k,n) = 4./6.*rdx*(var(i+1,j,k,n)-var(i-1,j,k,n))
     >                 - 1./12.*rdx*(var(i+2,j,k,n)-var(i-2,j,k,n))
     >                   + 0.5*dqdx_xy(i,j,3)
     >          *( rdze(k+1)*dqdx_zw(k+1)*zrbnd(k)
     >                       * (var(i,j,kp1,n)-var(i,j,k,n))
     >             + rdze(k)*dqdx_zw(k)*zlbnd(k)
     >                       * (var(i,j,k,n)-var(i,j,km1,n)) )

        enddo
        enddo
        enddo
*
* end time loop
*
        enddo
*
* y derivative
*
      else if (idir.eq.2) then

        do n=1,nt
*
* north and south boundaries
*
        do k=1,nz
          kp1 = min0(k+1,nz)
          km1 = max0(k-1,1)
        do i=1,nx

          j = 1
          deriv(i,j,k,n) = rdy*(var(i,j+1,k,n)-var(i,j,k,n))
     >                   + 0.5*dqdy_xy(i,j,3)
     >          *( rdze(k+1)*dqdy_zw(k+1)*zrbnd(k)
     >                       * (var(i,j,kp1,n)-var(i,j,k,n))
     >             + rdze(k)*dqdy_zw(k)*zlbnd(k)
     >                       * (var(i,j,k,n)-var(i,j,km1,n)) )
          j = ny+1
          deriv(i,j,k,n) = rdy*(var(i,j,k,n)-var(i,j-1,k,n))
     >                   + 0.5*dqdy_xy(i,j,3)
     >          *( rdze(k+1)*dqdy_zw(k+1)*zrbnd(k)
     >                       * (var(i,j,kp1,n)-var(i,j,k,n))
     >             + rdze(k)*dqdy_zw(k)*zlbnd(k)
     >                       * (var(i,j,k,n)-var(i,j,km1,n)) )

        enddo
        enddo

        do k=1,nz
          kp1 = min0(k+1,nz)
          km1 = max0(k-1,1)
        do m=1,2
          if (m.eq.1) j = 2
          if (m.eq.2) j = ny
        do i=1,nx

          deriv(i,j,k,n) = 0.5*rdy*(var(i,j+1,k,n)-var(i,j-1,k,n))
     >                   + 0.5*dqdy_xy(i,j,3)
     >          *( rdze(k+1)*dqdy_zw(k+1)*zrbnd(k)
     >                       * (var(i,j,kp1,n)-var(i,j,k,n))
     >             + rdze(k)*dqdy_zw(k)*zlbnd(k)
     >                       * (var(i,j,k,n)-var(i,j,km1,n)) )

        enddo
        enddo
        enddo
*
* interior
*
        do k=1,nz      
          kp1 = min0(k+1,nz)
          km1 = max0(k-1,1)
        do j=3,ny-1
        do i=1,nx

          deriv(i,j,k,n) = 4./6.*rdy*(var(i,j+1,k,n)-var(i,j-1,k,n))
     >                 - 1./12.*rdy*(var(i,j+2,k,n)-var(i,j-2,k,n))
     >                   + 0.5*dqdy_xy(i,j,3)
     >          *( rdze(k+1)*dqdy_zw(k+1)*zrbnd(k)
     >                       * (var(i,j,kp1,n)-var(i,j,k,n))
     >             + rdze(k)*dqdy_zw(k)*zlbnd(k)
     >                       * (var(i,j,k,n)-var(i,j,km1,n)) )

        enddo
        enddo
        enddo
*
* end time loop
*
        enddo
*
* z derivative
*
      else if (idir.eq.3) then
        
        do n=1,nt
        do k=1,nz
          kp1 = min0(k+1,nz)
          km1 = max0(k-1,1)
        do j=1,ny+1
        do i=1,nx

          deriv(i,j,k,n) = 0.5*dqdz(i,j,3)
     >      *( rdze(k+1)*zrbnd(k)*(var(i,j,k+1,n)-var(i,j,k,n))
     >         + rdze(k)*zlbnd(k)*(var(i,j,k,n)-var(i,j,k-1,n)) )

        enddo
        enddo
        enddo
        enddo
*
* t derivative
*
      else if (idir.eq.4) then

        do k=1,nz
        do j=1,ny+1
        do i=1,nx 

          deriv(i,j,k,1) = (var(i,j,k,2)-var(i,j,k,1))
     >                        / (time(2)-time(1))
          deriv(i,j,k,nt) = (var(i,j,k,nt)-var(i,j,k,nt-1))
     >                        / (time(nt)-time(nt-1))

          do n=2,nt-1
            dta = time(n+1)-time(n)
            dtb = time(n)-time(n-1)
            deriv(i,j,k,n) = dtb/(dta+dtb)
     >               * (var(i,j,k,n+1)-var(i,j,k,n))/dta
     >                     + dta/(dta+dtb)
     >               * (var(i,j,k,n)-var(i,j,k,n-1))/dtb
          enddo

        enddo
        enddo
        enddo

      end if

      return
      end
    
*------------------------------------------------------------------------
        
      subroutine diff1_wpts(var,deriv,dqdx_xy,dqdx_zt,dqdy_xy,
     >                 dqdy_zt,dqdz,dx,dy,time,nx,ny,nz,nt,idir)

      include 'misc.icl'

      integer nx,ny,nz,nt,idir
      real var(nx,ny,nz+1,nt),deriv(nx,ny,nz+1,nt)
      real dqdx_xy(nx+1,ny+1,3),dqdx_zt(nz+1)
      real dqdy_xy(nx+1,ny+1,3),dqdy_zt(nz+1)
      real dqdz(nx+1,ny+1,3),dx,dy,time(nt)

      integer i,j,k,n,m,ip1,im1,jp1,jm1,kp1,km1
      real rdx,rdy,da,db,ta,tb,dta,dtb      

      rdx = 1./dx
      rdy = 1./dy
*
* x derivative
*
      if (idir.eq.1) then

        do n=1,nt
*
* west and east boundaries
*
        do m=1,4
          if (m.eq.1) i = 1
          if (m.eq.2) i = 2
          if (m.eq.3) i = nx-1
          if (m.eq.4) i = nx
          ip1 = min0(i+1,nx)
          im1 = max0(i-1,1)
        do j=1,ny

          k = 1 
          deriv(i,j,k,n) = 0.5*rdx
     >       *( xrbnd(i)*(var(ip1,j,k,n)-var(i,j,k,n))
     >        + xlbnd(i)*(var(i,j,k,n)-var(im1,j,k,n)) )
     >                   + dqdx_xy(i,j,1)
     >      * rdzc(k)*dqdx_zt(k)*(var(i,j,k+1,n)-var(i,j,k,n))

          k = nz+1 
          deriv(i,j,k,n) = 0.5*rdx
     >       *( xrbnd(i)*(var(ip1,j,k,n)-var(i,j,k,n))
     >        + xlbnd(i)*(var(i,j,k,n)-var(im1,j,k,n)) )
     >                   + dqdx_xy(i,j,1)
     >      * rdzc(k-1)*dqdx_zt(k-1)*(var(i,j,k,n)-var(i,j,k-1,n))

          do k=2,nz
            deriv(i,j,k,n) = 0.5*rdx
     >         *( xrbnd(i)*(var(ip1,j,k,n)-var(i,j,k,n))
     >          + xlbnd(i)*(var(i,j,k,n)-var(im1,j,k,n)) )
     >                     + 0.5*dqdx_xy(i,j,1)
     >        *( rdzc(k)*dqdx_zt(k)*(var(i,j,k+1,n)-var(i,j,k,n))
     >       + rdzc(k-1)*dqdx_zt(k-1)*(var(i,j,k,n)-var(i,j,k-1,n)) )
          enddo

        enddo
        enddo
*
* interior
*
        do j=1,ny
        do i=3,nx-2

          k = 1 
          deriv(i,j,k,n) = 4./6.*rdx*(var(i+1,j,k,n)-var(i-1,j,k,n))
     >                 - 1./12.*rdx*(var(i+2,j,k,n)-var(i-2,j,k,n))
     >                   + dqdx_xy(i,j,1)
     >      * rdzc(k)*dqdx_zt(k)*(var(i,j,k+1,n)-var(i,j,k,n))

          k = nz+1 
          deriv(i,j,k,n) = 4./6.*rdx*(var(i+1,j,k,n)-var(i-1,j,k,n))
     >                 - 1./12.*rdx*(var(i+2,j,k,n)-var(i-2,j,k,n))
     >                   + dqdx_xy(i,j,1)
     >      * rdzc(k-1)*dqdx_zt(k-1)*(var(i,j,k,n)-var(i,j,k-1,n))

          do k=2,nz
            deriv(i,j,k,n) = 4./6.*rdx*(var(i+1,j,k,n)-var(i-1,j,k,n))
     >                   - 1./12.*rdx*(var(i+2,j,k,n)-var(i-2,j,k,n))
     >                     + 0.5*dqdx_xy(i,j,1)
     >        *( rdzc(k)*dqdx_zt(k)*(var(i,j,k+1,n)-var(i,j,k,n))
     >       + rdzc(k-1)*dqdx_zt(k-1)*(var(i,j,k,n)-var(i,j,k-1,n)) )
          enddo

        enddo
        enddo
*
* end time loop
*
        enddo
*
* y derivative
*
      else if (idir.eq.2) then

        do n=1,nt
*
* north and south boundaries
*
        do m=1,4
          if (m.eq.1) j = 1
          if (m.eq.2) j = 2
          if (m.eq.3) j = ny-1
          if (m.eq.4) j = ny
          jp1 = min0(j+1,ny)
          jm1 = max0(j-1,1)
        do i=1,nx
         
          k = 1 
          deriv(i,j,k,n) = 0.5*rdy
     >       *( yrbnd(j)*(var(i,jp1,k,n)-var(i,j,k,n))
     >        + ylbnd(j)*(var(i,j,k,n)-var(i,jm1,k,n)) )
     >                   + dqdy_xy(i,j,1)
     >      * rdzc(k)*dqdy_zt(k)*(var(i,j,k+1,n)-var(i,j,k,n))

          k = nz+1 
          deriv(i,j,k,n) = 0.5*rdy
     >       *( yrbnd(j)*(var(i,jp1,k,n)-var(i,j,k,n))
     >        + ylbnd(j)*(var(i,j,k,n)-var(i,jm1,k,n)) )
     >                   + dqdy_xy(i,j,1)
     >      * rdzc(k-1)*dqdy_zt(k-1)*(var(i,j,k,n)-var(i,j,k-1,n))

          do k=2,nz
            deriv(i,j,k,n) = 0.5*rdy
     >       *( yrbnd(j)*(var(i,jp1,k,n)-var(i,j,k,n))
     >        + ylbnd(j)*(var(i,j,k,n)-var(i,jm1,k,n)) )
     >                     + 0.5*dqdy_xy(i,j,1)
     >        *( rdzc(k)*dqdy_zt(k)*(var(i,j,k+1,n)-var(i,j,k,n))
     >       + rdzc(k-1)*dqdy_zt(k-1)*(var(i,j,k,n)-var(i,j,k-1,n)) )
          enddo
     
        enddo
        enddo
*
* interior
*
        do j=3,ny-2
        do i=1,nx

          k = 1 
          deriv(i,j,k,n) = 4./6.*rdy*(var(i,j+1,k,n)-var(i,j-1,k,n))
     >                 - 1./12.*rdy*(var(i,j+2,k,n)-var(i,j-2,k,n))
     >                   + dqdy_xy(i,j,1)
     >      * rdzc(k)*dqdy_zt(k)*(var(i,j,k+1,n)-var(i,j,k,n))

          k = nz+1 
          deriv(i,j,k,n) = 4./6.*rdy*(var(i,j+1,k,n)-var(i,j-1,k,n))
     >                 - 1./12.*rdy*(var(i,j+2,k,n)-var(i,j-2,k,n))
     >                   + dqdy_xy(i,j,1)
     >      * rdzc(k-1)*dqdy_zt(k-1)*(var(i,j,k,n)-var(i,j,k-1,n))

          do k=2,nz
            deriv(i,j,k,n) = 4./6.*rdy*(var(i,j+1,k,n)-var(i,j-1,k,n))
     >                   - 1./12.*rdy*(var(i,j+2,k,n)-var(i,j-2,k,n))
     >                     + 0.5*dqdy_xy(i,j,1)
     >        *( rdzc(k)*dqdy_zt(k)*(var(i,j,k+1,n)-var(i,j,k,n))
     >       + rdzc(k-1)*dqdy_zt(k-1)*(var(i,j,k,n)-var(i,j,k-1,n)) )
          enddo

        enddo
        enddo
*
* end time loop
*
        enddo
*
* z derivative
*
      else if (idir.eq.3) then
        
        do n=1,nt
        do j=1,ny
        do i=1,nx

          k = 1
          deriv(i,j,k,n) = dqdz(i,j,1)
     >               * rdzc(k)*(var(i,j,k+1,n)-var(i,j,k,n))

          k = nz+1
          deriv(i,j,k,n) = dqdz(i,j,1)
     >             * rdzc(k-1)*(var(i,j,k,n)-var(i,j,k-1,n))

          do k=2,nz
            deriv(i,j,k,n) = 0.5*dqdz(i,j,1)
     >         *( rdzc(k)*(var(i,j,k+1,n)-var(i,j,k,n))
     >        + rdzc(k-1)*(var(i,j,k,n)-var(i,j,k-1,n)) )      
          enddo

        enddo
        enddo
        enddo
*
* t derivative
*
      else if (idir.eq.4) then

        do k=1,nz
        do j=1,ny
        do i=1,nx 

          deriv(i,j,k,1) = (var(i,j,k,2)-var(i,j,k,1))
     >                        / (time(2)-time(1))
          deriv(i,j,k,nt) = (var(i,j,k,nt)-var(i,j,k,nt-1))
     >                        / (time(nt)-time(nt-1))

          do n=2,nt-1
            dta = time(n+1)-time(n)
            dtb = time(n)-time(n-1)
            deriv(i,j,k,n) = dtb/(dta+dtb)
     >               * (var(i,j,k,n+1)-var(i,j,k,n))/dta
     >                     + dta/(dta+dtb)
     >               * (var(i,j,k,n)-var(i,j,k,n-1))/dtb
          enddo

        enddo
        enddo
        enddo

      end if

      return
      end

*------------------------------------------------------------------------
        
      subroutine diff1_spts(var,deriv,dqdx_xy,dqdx_zw,dqdy_xy,
     >                 dqdy_zw,dqdz,dx,dy,time,nx,ny,nz,nt,idir)

      include 'misc.icl'

      integer nx,ny,nz,nt,idir
      real var(nx,ny,nz,nt),deriv(nx,ny,nz,nt)
      real dqdx_xy(nx+1,ny+1,3),dqdx_zw(nz+1)
      real dqdy_xy(nx+1,ny+1,3),dqdy_zw(nz+1)
      real dqdz(nx+1,ny+1,3),dx,dy,time(nt)

      integer i,j,k,n,m,ip1,im1,jp1,jm1,kp1,km1
      real rdx,rdy,da,db,ta,tb,dta,dtb      

      rdx = 1./dx
      rdy = 1./dy
*
* x derivative
*
      if (idir.eq.1) then

        do n=1,nt
*
* west and east boundaries
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
        do j=1,ny
          
          deriv(i,j,k,n) = 0.5*rdx
     >       *( xrbnd(i)*(var(ip1,j,k,n)-var(i,j,k,n))
     >        + xlbnd(i)*(var(i,j,k,n)-var(im1,j,k,n)) )
     >                   + 0.5*dqdx_xy(i,j,1)
     >          *( rdze(k+1)*dqdx_zw(k+1)*zrbnd(k)
     >                       * (var(i,j,kp1,n)-var(i,j,k,n))
     >             + rdze(k)*dqdx_zw(k)*zlbnd(k)
     >                       * (var(i,j,k,n)-var(i,j,km1,n)) )
        
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

          deriv(i,j,k,n) = 4./6.*rdx*(var(i+1,j,k,n)-var(i-1,j,k,n))
     >                 - 1./12.*rdx*(var(i+2,j,k,n)-var(i-2,j,k,n))
     >                   + 0.5*dqdx_xy(i,j,1)
     >          *( rdze(k+1)*dqdx_zw(k+1)*zrbnd(k)
     >                       * (var(i,j,kp1,n)-var(i,j,k,n))
     >             + rdze(k)*dqdx_zw(k)*zlbnd(k)
     >                       * (var(i,j,k,n)-var(i,j,km1,n)) )

        enddo
        enddo
        enddo
*
* end time loop
*
        enddo
*
* y derivative
*
      else if (idir.eq.2) then

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
          
          deriv(i,j,k,n) = 0.5*rdy
     >       *( yrbnd(j)*(var(i,jp1,k,n)-var(i,j,k,n))
     >        + ylbnd(j)*(var(i,j,k,n)-var(i,jm1,k,n)) )
     >                   + 0.5*dqdy_xy(i,j,1)
     >          *( rdze(k+1)*dqdy_zw(k+1)*zrbnd(k)
     >                       * (var(i,j,kp1,n)-var(i,j,k,n))
     >             + rdze(k)*dqdy_zw(k)*zlbnd(k)
     >                       * (var(i,j,k,n)-var(i,j,km1,n)) )
        
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

          deriv(i,j,k,n) = 4./6.*rdy*(var(i,j+1,k,n)-var(i,j-1,k,n))
     >                 - 1./12.*rdy*(var(i,j+2,k,n)-var(i,j-2,k,n))
     >                   + 0.5*dqdy_xy(i,j,1)
     >          *( rdze(k+1)*dqdy_zw(k+1)*zrbnd(k)
     >                       * (var(i,j,kp1,n)-var(i,j,k,n))
     >             + rdze(k)*dqdy_zw(k)*zlbnd(k)
     >                       * (var(i,j,k,n)-var(i,j,km1,n)) )

        enddo
        enddo
        enddo
*
* end time loop
*
        enddo
*
* z derivative
*
      else if (idir.eq.3) then
        
        do n=1,nt
        do k=1,nz
          kp1 = min0(k+1,nz)
          km1 = max0(k-1,1)
        do j=1,ny
        do i=1,nx

          deriv(i,j,k,n) = 0.5*dqdz(i,j,1)
     >      *( rdze(k+1)*zrbnd(k)*(var(i,j,kp1,n)-var(i,j,k,n))
     >         + rdze(k)*zlbnd(k)*(var(i,j,k,n)-var(i,j,km1,n)) )

        enddo
        enddo
        enddo
        enddo
*
* t derivative
*
      else if (idir.eq.4) then

        do k=1,nz
        do j=1,ny
        do i=1,nx 

          deriv(i,j,k,1) = (var(i,j,k,2)-var(i,j,k,1))
     >                        / (time(2)-time(1))
          deriv(i,j,k,nt) = (var(i,j,k,nt)-var(i,j,k,nt-1))
     >                        / (time(nt)-time(nt-1))

          do n=2,nt-1
            dta = time(n+1)-time(n)
            dtb = time(n)-time(n-1)
            deriv(i,j,k,n) = dtb/(dta+dtb)
     >               * (var(i,j,k,n+1)-var(i,j,k,n))/dta
     >                     + dta/(dta+dtb)
     >               * (var(i,j,k,n)-var(i,j,k,n-1))/dtb
          enddo

        enddo
        enddo
        enddo

      end if

      return
      end

*-------------------------------------------------------------------------
   
      subroutine diff1_t(var,deriv,time,nt)

      integer nt
      real var(nt),deriv(nt),time(nt)
  
      integer n 
      real dta,dtb     

      deriv(1) = (var(2)-var(1))/(time(2)-time(1))
      deriv(nt) = (var(nt)-var(nt-1))/(time(nt)-time(nt-1))

      do n=2,nt-1
        dta = time(n+1)-time(n)
        dtb = time(n)-time(n-1)
        deriv(n) = dtb/(dta+dtb)*(var(n+1)-var(n))/dta
     >                   + dta/(dta+dtb)*(var(n)-var(n-1))/dtb
      enddo

      return
      end
 


        
