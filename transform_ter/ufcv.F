
#include "pointer.icl"

      ive_ptr function ufcv(ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names)

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer ndims2,dims2(4)
      real stag2(4),min2(4),max2(4),miss2
      character*(80) dunits2,disunits2,dim_names2(4)

      integer k,twod
      real dom_min(4),dom_max(4),dx,dy,ztop,r_val
      logical err,new
      ive_ptr getmem,getvar,u,v,w,rho,zs_p,zs_u,zs_v,sgz,wgz

      ufcv = 0

      call getrarr('plmin',dom_min,4,err)
      if (err) then
        write(*,*) 'Error:  unable to access physical space ',
     >                         'domain minimums.'
        return
      end if
      call getrarr('plmax',dom_max,4,err)
      if (err) then
        write(*,*) 'Error:  unable to access physical space ',
     >                        'domain maximums.'
        return
      end if

      if (dom_min(2).eq.dom_max(2)) then
        twod = 1
      else 
        twod = 0
      end if

      zs_p = getvar('zbot_p',ndims,dims,stag,min,max,miss,dunits,
     >                    disunits,dim_names,new)
      if (zs_p.eq.0) then
        write(*,*) 'Error:  unable to access field zbot_p.'
        return
      end if

      zs_u = getmem((dims(1)+1)*dims(2))
      if (zs_u.eq.0) then
        write(*,*) 'Error:  unable to allocate memory for zs_u.'
        return
      end if
      call zs_2_upts(%val(zs_u),%val(zs_p),dims(1),dims(2))

      if (twod.ne.1) then
        zs_v = getmem(dims(1)*(dims(2)+1))
        if (zs_v.eq.0) then
          write(*,*) 'Error:  unable to allocate memory for zs_v.'
          call freemem(zs_u)
          return
        end if
        call zs_2_vpts(%val(zs_v),%val(zs_p),dims(1),dims(2))
      end if

      dx = (max(1)-min(1))/(dims(1)-1)
      if (twod.eq.1) then
        dy = 1.
      else
        dy = (max(2)-min(2))/(dims(2)-1)
      end if

      wgz = getvar('wgz',ndims,dims,stag,min,max,miss,dunits,
     >                   disunits,dim_names,new)
      if (wgz.eq.0) then
        write(*,*) 'Error:  unable to access field wgz.'
        call freemem(zs_u)
        if (twod.ne.1) call freemem(zs_v)
        return
      end if
      ztop = r_val(%val(wgz),dims(1))

      sgz = getvar('sgz',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)
      if (sgz.eq.0) then
        write(*,*) 'Error:  unable to access field sgz.'
        call freemem(zs_u)
        if (twod.ne.1) call freemem(zs_v)
        return
      end if

      rho = getvar('RHO',ndims,dims,stag,min,max,miss,dunits,
     >                disunits,dim_names,new)
      if (rho.eq.0) then
        write(*,*) 'Error:  unable to access mass density.'
        call freemem(zs_u)
        if (twod.ne.1) call freemem(zs_v)
        return
      end if

      ufcv = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (ufcv.eq.0) then
        write(*,*) 'Error:  unable to allocate space for ufcv.'
        call freemem(zs_u)
        if (twod.ne.1) call freemem(zs_v)
        return
      end if

      do k=1,dims(1)*dims(2)*dims(3)*dims(4)
        call s_val(%val(ufcv),k,0.)
      enddo

      u = getvar('u',ndims2,dims2,stag2,min2,max2,miss2,dunits2,
     >               disunits2,dim_names2,new)
      if (u.eq.0) then
        write(*,*) 'Error:  unable to access u.'
        call freemem(zs_u)
        if (twod.ne.1) call freemem(zs_v)
        return
      end if
    
      call u_ufcv(%val(ufcv),%val(u),%val(rho),%val(zs_u),%val(zs_p),
     >             %val(wgz),%val(sgz),ztop,dx,dims(1),dims(2),
     >             dims(3),dims(4))

      w = getvar('w',ndims2,dims2,stag2,min2,max2,miss2,dunits2,
     >               disunits2,dim_names2,new)
      if (w.eq.0) then
        write(*,*) 'Error:  unable to access w.'
        call freemem(zs_u)
        if (twod.ne.1) call freemem(zs_v)
      end if

      call w_ufcv(%val(ufcv),%val(u),%val(w),%val(rho),%val(zs_p),
     >             %val(wgz),%val(sgz),ztop,dims(1),dims(2),dims(3),
     >             dims(4))

      call driver('free = w',0)

      if (twod.ne.1) then

        v = getvar('v',ndims2,dims2,stag2,min2,max2,miss2,dunits2,
     >                disunits2,dim_names2,new)
        if (v.eq.0) then
          write(*,*) 'Error:  unable to access v.'
          call freemem(zs_u)
          call freemem(zs_v)
          return
        end if

        call v_ufcv(%val(ufcv),%val(u),%val(v),%val(rho),%val(zs_p),
     >               %val(zs_v),%val(wgz),%val(sgz),ztop,dy,dims(1),
     >               dims(2),dims(3),dims(4))

        call driver('free = v',0)

      end if

      call flip_sign(%val(ufcv),dims(1),dims(2),dims(3),dims(4))

      dunits = 'kg/m/m/s/s'
      disunits = 'kg/m/m/s/s'

      call freemem(zs_u)
      if (twod.ne.1) call freemem(zs_v)

      return
      end

*--------------------------------------------------------------------------

      subroutine u_ufcv(ufcv,u,r,zs_u,zs_p,wgz,sgz,ztop,dx,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real ufcv(nx,ny,nz,nt),u(nx+1,ny,nz,nt),r(nx,ny,nz,nt)
      real zs_u(nx+1,ny),zs_p(nx,ny),wgz(nz+1),sgz(nz),ztop,dx

      integer i,j,k,n,ip1,im1,NZMAX
      parameter (NZMAX=201)
      real H,H1,H2,dzsdx,Gx1,Gx2,term1,term2,rdzet(NZMAX)

      do k=1,nz-1
        rdzet(k) = 1/(wgz(k+1)-wgz(k))
      enddo
      rdzet(nz) = 1/(sgz(nz)-sgz(nz-1))

      do n=1,nt
      do j=1,ny
      do i=1,nx
        ip1 = min(i+1,nx)
        im1 = max(i-1,1)

        H = ztop/(ztop-zs_p(i,j))
        H1 = ztop/(ztop-zs_u(i,j))
        H2 = ztop/(ztop-zs_u(i+1,j))
        dzsdx = (zs_u(i+1,j)-zs_u(i,j))/dx

        Gx2 = (wgz(2)-ztop)/(ztop-zs_p(i,j))*dzsdx

        term1 = 0.5*((r(ip1,j,1,n)+r(i,j,1,n))*u(i+1,j,1,n) 
     >                   *u(i+1,j,1,n)/H2 -
     >               (r(i,j,1,n)+r(im1,j,1,n))*u(i,j,1,n) 
     >                   *u(i,j,1,n)/H1)/dx

        term2 = 0.5*(r(i,j,1,n)+r(i,j,2,n))
     >      *(0.25*(u(i,j,1,n)+u(i+1,j,1,n)
     >                         +u(i,j,2,n)+u(i+1,j,2,n)))**2
     >                  *Gx2*rdzet(1)/H

        ufcv(i,j,1,n) = ufcv(i,j,1,n) + H*(term1+term2)

        do k=2,nz-1

          Gx2 = (wgz(k+1)-ztop)/(ztop-zs_p(i,j))*dzsdx
          Gx1 = (wgz(k)-ztop)/(ztop-zs_p(i,j))*dzsdx

          term1 = 0.5*((r(ip1,j,k,n)+r(i,j,k,n))*u(i+1,j,k,n)
     >                      *u(i+1,j,k,n)/H2 -
     >                 (r(i,j,k,n)+r(im1,j,k,n))*u(i,j,k,n)
     >                      *u(i,j,k,n)/H1)/dx

          term2 = 0.5*((r(i,j,k,n)+r(i,j,k+1,n))*Gx2
     >      *(0.25*(u(i,j,k,n)+u(i+1,j,k,n)
     >                    +u(i,j,k+1,n)+u(i+1,j,k+1,n)))**2
     >               - (r(i,j,k-1,n)+r(i,j,k,n))*Gx1
     >      *(0.25*(u(i,j,k-1,n)+u(i+1,j,k-1,n)
     >                    +u(i,j,k,n)+u(i+1,j,k,n)))**2)
     >                         *rdzet(k)/H

          ufcv(i,j,k,n) = ufcv(i,j,k,n) + H*(term1+term2)

        enddo

        Gx2 = (sgz(nz)-ztop)/(ztop-zs_p(i,j))*dzsdx
        Gx1 = (sgz(nz-1)-ztop)/(ztop-zs_p(i,j))*dzsdx
        
        term1 = 0.5*((r(ip1,j,nz,n)+r(i,j,nz,n))*u(i+1,j,nz,n)
     >                       *u(i+1,j,nz,n)/H2 -
     >               (r(i,j,nz,n)+r(im1,j,nz,n))*u(i,j,nz,n)
     >                       *u(i,j,nz,n)/H1)/dx

        term2 = (r(i,j,nz,n)*Gx2
     >       *(0.5*(u(i,j,nz,n)+u(i+1,j,nz,n)))**2
     >         - r(i,j,nz-1,n)*Gx1
     >       *(0.5*(u(i,j,nz-1,n)+u(i+1,j,nz-1,n)))**2)
     >                     *rdzet(nz)/H

        ufcv(i,j,nz,n) = ufcv(i,j,nz,n) + H*(term1+term2)

      enddo
      enddo
      enddo

      return
      end

*-------------------------------------------------------------------------

      subroutine v_ufcv(ufcv,u,v,r,zs_p,zs_v,wgz,sgz,ztop,dy,
     >                        nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real ufcv(nx,ny,nz,nt),u(nx+1,ny,nz,nt),v(nx,ny+1,nz,nt)
      real r(nx,ny,nz,nt),zs_p(nx,ny),zs_v(nx,ny+1),wgz(nz+1),sgz(nz)
      real ztop,dy

      integer i,j,k,n,jp1,jm1,MAXNZ
      parameter (MAXNZ=201)
      real H,H1,H2,dzsdy,Gy1,Gy2,term1,term2,rdzet(MAXNZ)

      do k=1,nz-1
        rdzet(k) = 1/(wgz(k+1)-wgz(k))
      enddo
      rdzet(nz) = 1/(sgz(nz)-sgz(nz-1))

      do n=1,nt
      do j=1,ny
        jp1 = min(j+1,ny)
        jm1 = max(j-1,1)
      do i=1,nx

        H = ztop/(ztop-zs_p(i,j))
        H1 = ztop/(ztop-zs_v(i,j))
        H2 = ztop/(ztop-zs_v(i,j+1))
        dzsdy = (zs_v(i,j+1)-zs_v(i,j))/dy

        Gy2 = (wgz(2)-ztop)/(ztop-zs_p(i,j))*dzsdy

        term1 = 0.125*((r(i,jp1,1,n)+r(i,j,1,n))*v(i,j+1,1,n)
     >    *(u(i,j,1,n)+u(i+1,j,1,n)+u(i,jp1,1,n)+u(i+1,jp1,1,n))/H2 
     >               - (r(i,j,1,n)+r(i,jm1,1,n))*v(i,j,1,n)
     >    *(u(i,jm1,1,n)+u(i+1,jm1,1,n)+u(i,j,1,n)+u(i+1,j,1,n))/H1)
     >                               /dy

        term2 = 0.5*(r(i,j,2,n)+r(i,j,1,n))*Gy2
     >    *0.25*(u(i,j,1,n)+u(i+1,j,1,n)+u(i,j,2,n)+u(i+1,j,2,n))
     >    *0.25*(v(i,j,1,n)+v(i,j+1,1,n)+v(i,j,2,n)+v(i,j+1,2,n))
     >                      *rdzet(1)/H

        ufcv(i,j,1,n) = ufcv(i,j,1,n) + H*(term1+term2)

        do k=2,nz-1

          Gy2 = (wgz(k+1)-ztop)/(ztop-zs_p(i,j))*dzsdy
          Gy1 = (wgz(k)-ztop)/(ztop-zs_p(i,j))*dzsdy

          term1 = 0.125*((r(i,j,k,n)+r(i,jp1,k,n))*v(i,j+1,k,n)
     >     *(u(i,j,k,n)+u(i+1,j,k,n)+u(i,jp1,k,n)+u(i+1,jp1,k,n))/H2 
     >                -(r(i,jm1,k,n)+r(i,j,k,n))*v(i,j,k,n)
     >     *(u(i,jm1,k,n)+u(i+1,jm1,k,n)+u(i,j,k,n)+u(i+1,j,k,n))/H1)
     >                            /dy

          term2 = 0.5*((r(i,j,k,n)+r(i,j,k+1,n))*Gy2
     >     *0.25*(u(i,j,k,n)+u(i+1,j,k,n)+u(i,j,k+1,n)+u(i+1,j,k+1,n))
     >     *0.25*(v(i,j,k,n)+v(i,j+1,k,n)+v(i,j,k+1,n)+v(i,j+1,k+1,n))
     >               - (r(i,j,k-1,n)+r(i,j,k,n))*Gy1
     >     *0.25*(u(i,j,k-1,n)+u(i+1,j,k-1,n)+u(i,j,k,n)+u(i+1,j,k,n))
     >     *0.25*(v(i,j,k-1,n)+v(i,j+1,k-1,n)+v(i,j,k,n)+v(i,j+1,k,n)))
     >                         *rdzet(k)/H

          ufcv(i,j,k,n) = ufcv(i,j,k,n) + H*(term1+term2)

        enddo

        Gy2 = (sgz(nz)-ztop)/(ztop-zs_p(i,j))*dzsdy
        Gy1 = (sgz(nz-1)-ztop)/(ztop-zs_p(i,j))*dzsdy

        term1 = 0.125*((r(i,jp1,nz,n)+r(i,j,nz,n))*v(i,j+1,nz,n)
     >    *(u(i,j,nz,n)+u(i+1,j,nz,n)
     >                     +u(i,jp1,nz,n)+u(i+1,jp1,nz,n))/H2 
     >               - (r(i,j,nz,n)+r(i,jm1,nz,n))*v(i,j,nz,n)
     >    *(u(i,jm1,nz,n)+u(i+1,jm1,nz,n)
     >                     +u(i,j,nz,n)+u(i+1,j,nz,n))/H1)
     >                         /dy

        term2 = 0.25*(r(i,j,nz,n)*Gy2
     >                 *(u(i,j,nz,n)+u(i+1,j,nz,n))
     >                 *(v(i,j,nz,n)+v(i,j+1,nz,n))
     >              - r(i,j,nz-1,n)*Gy1
     >                 *(u(i,j,nz-1,n)+u(i+1,j,nz-1,n))
     >                 *(v(i,j,nz-1,n)+v(i,j+1,nz-1,n)))
     >                       *rdzet(nz)/H

        ufcv(i,j,nz,n) = ufcv(i,j,nz,n) + H*(term1+term2)

      enddo
      enddo
      enddo

      return
      end

*------------------------------------------------------------------------

      subroutine w_ufcv(ufcv,u,w,r,zs_p,wgz,sgz,ztop,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real ufcv(nx,ny,nz,nt),u(nx+1,ny,nz,nt),w(nx,ny,nz+1,nt)
      real r(nx,ny,nz,nt),zs_p(nx,ny),wgz(nz+1),sgz(nz),ztop

      
      integer i,j,k,n,MAXNZ
      parameter (MAXNZ=201)
      real H,rdzet(MAXNZ)

      do k=1,nz-1
        rdzet(k) = 1/(wgz(k+1)-wgz(k))
      enddo
      rdzet(nz) = 1/(sgz(nz)-sgz(nz-1))

      do n=1,nt
      do j=1,ny
      do i=1,nx

        H = ztop/(ztop-zs_p(i,j))

        ufcv(i,j,1,n) = ufcv(i,j,1,n) +
     >       H*0.5*(r(i,j,1,n)+r(i,j,2,n))*w(i,j,2,n)
     >   *0.25*(u(i,j,1,n)+u(i+1,j,1,n)+u(i,j,2,n)+u(i+1,j,2,n))
     >                *rdzet(1)

        do k=2,nz-1

          ufcv(i,j,k,n) = ufcv(i,j,k,n) +
     >       H*0.5*((r(i,j,k,n)+r(i,j,k+1,n))*w(i,j,k+1,n)
     >   *0.25*(u(i,j,k,n)+u(i+1,j,k,n)+u(i,j,k+1,n)+u(i+1,j,k+1,n))
     >            - (r(i,j,k-1,n)+r(i,j,k,n))*w(i,j,k,n)
     >   *0.25*(u(i,j,k-1,n)+u(i+1,j,k-1,n)+u(i,j,k,n)+u(i+1,j,k,n)))
     >                  *rdzet(k)

        enddo

        ufcv(i,j,nz,n) = ufcv(i,j,nz,n) +
     >       H*0.25*(r(i,j,nz,n)*(w(i,j,nz+1,n)+w(i,j,nz,n))
     >            *(u(i,j,nz,n)+u(i+1,j,nz,n)) 
     >             - r(i,j,nz-1,n)*(w(i,j,nz,n)+w(i,j,nz-1,n))
     >            *(u(i,j,nz-1,n)+u(i+1,j,nz-1,n)))
     >                 *rdzet(nz)

      enddo
      enddo
      enddo

      return
      end

*------------------------------------------------------------------------

      subroutine flip_sign(var,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real var(nx,ny,nz,nt)

      integer i,j,k,n

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx

        var(i,j,k,n) = -var(i,j,k,n)

      enddo
      enddo
      enddo
      enddo

      return
      end
