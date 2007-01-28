
#include "pointer.icl"

      ive_ptr function def2(ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      logical new
      ive_ptr getmem,u,v,w,def

      def2 = 0

      u = getvar('u',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                                 dim_names,new)
      if (u.eq.0) then
        write(6,*) 'Error:  unable to get variable u.'
        return
      end if

      v = getvar('v',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                                 dim_names,new)
      if (v.eq.0) then
        write(6,*) 'Error:  unable to get variable v.'
        return
      end if

      w = getvar('w',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                                 dim_names,new)
      if (w.eq.0) then
        write(6,*) 'Error:  unable to get variable w.'
        return
      end if

      dims(3) = dims(3) - 1
      stag(3) = 0.5
      min(3) = sgz_min
      max(3) = sgz_max
      dim_names(3) = 'nz'
      dunits = '1/s/s'
      disunits = '1/s/s'

      def = getmem((dims(1)+1)*(dims(2)+1)*(dims(3)+1))
      if (def.eq.0) then
        write(6,*) 'Error:  unable to allocate space for def.'
        return
      end if

      def2 = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (def2.eq.0) then
        write(6,*) 'Error:  unable to allocate space for def2.'
        call freemem(def)
        return
      end if

      call def2_calc(%val(def2),%val(def),%val(u),%val(v),%val(w),
     >          %val(zbot_p),%val(zbot_u),%val(zbot_v),%val(sgz),
     >          %val(wgz),dx,dy,ztop,dims(1),dims(2),dims(3),dims(4))

      call freemem(def)

      return 
      end

*------------------------------------------------------------------------

      subroutine def2_calc(def2,def,u,v,w,zs_p,zs_u,zs_v,sgz,wgz,
     >                           dx,dy,ztop,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real def2(nx,ny,nz,nt),def(nx+1,ny+1,nz+1),u(nx+1,ny,nz,nt)
      real v(nx,ny+1,nz,nt),w(nx,ny,nz+1,nt),zs_p(nx,ny)
      real zs_u(nx+1,ny),zs_v(nx,ny+1),sgz(nz),wgz(nz+1),dx,dy,ztop

      integer i,j,k,n,kp1,km1,nzmax
      parameter(nzmax=201)
      real dzets(nzmax),dzetw(nzmax),zlbnd(nzmax),zrbnd(nzmax)
      real dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
      real dudzet,dvdzet,dwdzet,dzetdx,dzetdy,dzetdz,rdx,rdy

      if (nz+1.gt.nzmax) then
        write(6,*) 'Error:  nz too big!  Increase nzmax in mix.F'
        stop
      end if

      rdx = 1./dx
      rdy = 1./dy

      do k=1,nz-1
        dzets(k) = sgz(k+1) - sgz(k)
        dzetw(k) = wgz(k+1) - wgz(k)
        zlbnd(k) = 1.
        zrbnd(k) = 1.
      enddo

      dzets(nz) = 1.
      dzetw(nz) = wgz(nz+1) - wgz(nz)
      zlbnd(nz) = 2.
      zlbnd(1) = 0.
      zrbnd(nz) = 0.
      zrbnd(1) = 2.

      do k=1,nz
        kp1 = min(k+1,nz)
        km1 = max(k-1,1)
      do j=1,ny
      do i=1,nx

        dzetdx = (sgz(k)-ztop)/(ztop-zs_p(i,j))
     >                   *rdx*(zs_u(i+1,j)-zs_u(i,j))
        dzetdy = (sgz(k)-ztop)/(ztop-zs_p(i,j))
     >                   *rdy*(zs_v(i,j+1)-zs_v(i,j))
        dzetdz = ztop/(ztop-zs_p(i,j))

        do n=1,nt

          dudzet = 0.25*zrbnd(k)*(u(i,j,kp1,n)-u(i,j,k,n)
     >                   + u(i+1,j,kp1,n)-u(i+1,j,k,n))/dzets(k)
     >           + 0.25*zlbnd(k)*(u(i,j,k,n)-u(i,j,km1,n)
     >                   + u(i+1,j,k,n)-u(i+1,j,km1,n))/dzets(km1)
          dvdzet = 0.25*zrbnd(k)*(v(i,j,kp1,n)-v(i,j,k,n)
     >                   + v(i,j+1,kp1,n)-v(i,j+1,k,n))/dzets(k)
     >           + 0.25*zlbnd(k)*(v(i,j,k,n)-v(i,j,km1,n)
     >                   + v(i,j+1,k,n)-v(i,j+1,km1,n))/dzets(km1)
          dwdzet = (w(i,j,k+1,n)-w(i,j,k,n))/dzetw(k)

          dudx = rdx*(u(i+1,j,k,n)-u(i,j,k,n)) + dudzet*dzetdx
          dvdy = rdy*(v(i,j+1,k,n)-v(i,j,k,n)) + dvdzet*dzetdy
          dwdz = dwdzet*dzetdz

          def2(i,j,k,n) = 0.5*(4.*dudx/3. - 2.*(dvdy+dwdz)/3.)**2
     >                  + 0.5*(4.*dvdy/3. - 2.*(dudx+dwdz)/3.)**2
     >                  + 0.5*(4.*dwdz/3. - 2.*(dudx+dvdy)/3.)**2

        enddo

      enddo
      enddo
      enddo

      do n=1,nt

        do k=1,nz
          kp1 = min(k+1,nz)
          km1 = max(k-1,1)
        do j=2,ny
        do i=2,nx

          dzetdx = (sgz(k)-ztop)/(ztop-0.5*(zs_v(i,j)+zs_v(i-1,j)))
     >                     *rdx*(zs_v(i,j)-zs_v(i-1,j))
          dzetdy = (sgz(k)-ztop)/(ztop-0.5*(zs_u(i,j)+zs_u(i,j-1)))
     >                     *rdy*(zs_u(i,j)-zs_u(i,j-1))

          dudzet = 0.25*zrbnd(k)*(u(i,j,kp1,n)-u(i,j,k,n)
     >                  + u(i,j-1,kp1,n)-u(i,j-1,k,n))/dzets(k)
     >           + 0.25*zlbnd(k)*(u(i,j,k,n)-u(i,j,km1,n)
     >                  + u(i,j-1,k,n)-u(i,j-1,km1,n))/dzets(km1)
          dvdzet = 0.25*zrbnd(k)*(v(i,j,kp1,n)-v(i,j,k,n)
     >                  + v(i-1,j,kp1,n)-v(i-1,j,k,n))/dzets(k)
     >           + 0.25*zlbnd(k)*(v(i,j,k,n)-v(i,j,km1,n)
     >                  + v(i-1,j,k,n)-v(i-1,j,km1,n))/dzets(km1)

          dudy = rdy*(u(i,j,k,n)-u(i,j-1,k,n)) + dudzet*dzetdy
          dvdx = rdx*(v(i,j,k,n)-v(i-1,j,k,n)) + dvdzet*dzetdx

          def(i,j,k) = dudy + dvdx

        enddo
        enddo
        enddo

        do k=1,nz
        do i=2,nx
          def(i,1,k) = def(i,2,k)
          def(i,ny+1,k) = def(i,ny,k)
        enddo
        enddo
 
        do k=1,nz
        do j=1,ny+1
          def(1,j,k) = def(2,j,k) 
          def(nx+1,j,k) = def(nx,j,k)
        enddo
        enddo

        do k=1,nz
        do j=1,ny
        do i=1,nx
          def2(i,j,k,n) = def2(i,j,k,n)  
     >             + (0.25*(def(i,j,k)+def(i+1,j,k)
     >                  + def(i,j+1,k)+def(i+1,j+1,k)))**2
        enddo
        enddo
        enddo

        do k=2,nz
        do j=1,ny
        do i=2,nx

          dzetdx = (wgz(k)-ztop)/(ztop-zs_u(i,j))
     >                     *rdx*(zs_p(i,j)-zs_p(i-1,j))
          dzetdz = ztop/(ztop-zs_u(i,j))

          dwdzet = 0.25*(w(i,j,k+1,n)-w(i,j,k,n)
     >                      + w(i-1,j,k+1,n)-w(i-1,j,k,n))/dzetw(k)
     >           + 0.25*(w(i,j,k,n)-w(i,j,k-1,n)
     >                      + w(i-1,j,k,n)-w(i-1,j,k-1,n))/dzetw(k-1)
          dudzet = (u(i,j,k,n)-u(i,j,k-1,n))/dzets(k-1)

          dwdx = rdx*(w(i,j,k,n)-w(i-1,j,k,n)) + dwdzet*dzetdx
          dudz = dudzet*dzetdz

          def(i,j,k) = dwdx + dudz

        enddo
        enddo
        enddo

        do k=2,nz
        do j=1,ny
          def(1,j,k) = def(2,j,k)
          def(nx+1,j,k) = def(nx,j,k)
        enddo
        enddo

        do j=1,ny
        do i=1,nx+1
          def(i,j,1) = 0.
          def(i,j,nz+1) = def(i,j,nz)
        enddo
        enddo

        do k=1,nz
        do j=1,ny
        do i=1,nx
          def2(i,j,k,n) = def2(i,j,k,n)
     >             + (0.25*(def(i,j,k)+def(i+1,j,k)
     >                  + def(i,j,k+1)+def(i+1,j,k+1)))**2
        enddo
        enddo
        enddo

        do k=2,nz
        do j=2,ny
        do i=1,nx

          dzetdy = (sgz(k)-ztop)/(ztop-zs_v(i,j))
     >                     *rdy*(zs_p(i,j)-zs_p(i,j-1))
          dzetdz = ztop/(ztop-zs_v(i,j))

          dwdzet = 0.25*(w(i,j,k+1,n)-w(i,j,k,n)
     >                      + w(i,j-1,k+1,n)-w(i,j-1,k,n))/dzetw(k)
     >           + 0.25*(w(i,j,k,n)-w(i,j,k-1,n)
     >                      + w(i,j-1,k,n)-w(i,j-1,k-1,n))/dzetw(k-1)
          dvdzet = (v(i,j,k,n)-v(i,j,k-1,n))/dzets(k-1)

          dwdy = rdy*(w(i,j,k,n)-w(i,j-1,k,n)) + dwdzet*dzetdy
          dvdz = dvdzet*dzetdz

          def(i,j,k) = dwdy + dvdz

        enddo
        enddo
        enddo

        do k=2,nz
        do i=1,nx
          def(i,1,k) = def(i,2,k)
          def(i,ny+1,k) = def(i,ny,k)
        enddo
        enddo
 
        do j=1,ny+1
        do i=1,nx
          def(i,j,1) = 0.
          def(i,j,nz+1) = def(i,j,nz)
        enddo
        enddo

        do k=1,nz
        do j=1,ny
        do i=1,nx

          def2(i,j,k,n) = def2(i,j,k,n)
     >             + (0.25*(def(i,j,k)+def(i,j+1,k)
     >                  + def(i,j,k+1)+def(i,j+1,k+1)))**2

         enddo
         enddo
         enddo

       enddo

       return
       end

        
