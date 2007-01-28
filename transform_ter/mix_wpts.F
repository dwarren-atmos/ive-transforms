
#include "pointer.icl"

      ive_ptr function mix_wpts(which,ndims,dims,stag,min,max,
     >                    miss,dunits,disunits,dim_names)

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)
      character*(10) which

      ive_ptr getmem,getvar,uptr,vptr,wptr,tptr,zbot_u,zbot_v,zbot_p,
     >            sgz,wgz,d11,d22,d33,d12,d13,d23
      integer iwhch,varsp
      real delx,dely,r_val
      logical new


      if (which(1:4).eq.'def2') then
        iwhch = 1
      else if (which(1:7).eq.'richinv') then
        iwhch = 2
      else if (which(1:2).eq.'km') then
        iwhch = 3
      else
        mix_wpts = 0
        write(*,*) 'Error: do not know which field to calculate'
        return
      end if

      sgz = getvar('sgz',ndims,dims,stag,min,max,miss,dunits,
     >                   disunits,dim_names,new)
      wgz = getvar('wgz',ndims,dims,stag,min,max,miss,dunits,
     >                   disunits,dim_names,new)
      zbot_p = getvar('zbot_p',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new) 
      zbot_u = getvar('zbot_u',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)
      zbot_v = getvar('zbot_v',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)
      uptr = getvar('u',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)
      vptr = getvar('v',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)
      wptr = getvar('w',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)
      tptr = getvar('thet',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)

      varsp = dims(1)*dims(2)*dims(3)*dims(4)
      mix_wpts = getmem(varsp)
      if (mix_wpts.eq.0) then
        write(*,*) 'Error:  could not allocate space for mix_wpts.'
        return
      end if
      d11 = getmem(6*varsp)
      if (d11.eq.0) then
        write(*,*) 'Error:  could not allocate space for deformations.'
        call freemem(mix_wpts)
        mix_wpts = 0
        return
      end if
      d22 = d11 + varsp
      d33 = d22 + varsp
      d12 = d33 + varsp
      d13 = d12 + varsp
      d23 = d13 + varsp

      miss = -9999.
      if (iwhch.eq.1) then
        dunits = '1/(s*s)'
        disunits = '1/(s*s)'
      else if (iwhch.eq.2) then
        dunits = ' '
        disunits = ' '
      else if (iwhch.eq.3) then
        dunits = '(m*m)/s'
        disunits = '(m*m)/s'
      end if
      min(3) = r_val(%val(wgz),1)
      max(3) = r_val(%val(wgz),dims(3))

      delx = (max(1)-min(1))/(dims(1)-1)
      dely = (max(2)-min(2))/(dims(2)-1) 

      call mix_wpts_calc(%val(mix_wpts),%val(d11),%val(d22),
     >           %val(d33),%val(d12),%val(d13),%val(d23),%val(uptr),
     >           %val(vptr),%val(wptr),%val(tptr),%val(zbot_u),
     >           %val(zbot_v),%val(zbot_p),%val(sgz),%val(wgz),
     >           miss,delx,dely,dims(1),dims(2),dims(3),
     >           dims(4),iwhch)

      call freemem(d11)

      return
      end

*************************************************************************

      subroutine mix_wpts_calc(mix,d11,d22,d33,d12,d13,d23,u,v,w,
     >                 thet,teru,terv,terp,sgz,wgz,miss,delx,
     >                 dely,nx,ny,nz,nt,iwhch)

      integer nx,ny,nz,nt,iwhch
      real mix(nx,ny,nz,nt)
      real d11(nx,ny,nz,nt),d22(nx,ny,nz,nt),d33(nx,ny,nz,nt),
     >     d12(nx,ny,nz,nt),d13(nx,ny,nz,nt),d23(nx,ny,nz,nt),
     >     u(nx+1,ny,nz,nt),v(nx,ny+1,nz,nt),w(nx,ny,nz+1,nt),
     >     thet(nx,ny,nz,nt),teru(nx+1,ny),terv(nx,ny+1),
     >     terp(nx,ny),sgz(nz),wgz(nz+1),miss,delx,dely
      
      integer i,j,k,n,ip1,im1,jp1,jm1,kp1,km1,NZMAX
      PARAMETER (NZMAX = 200)
      real dzet(NZMAX),dzetw(NZMAX),delta(NZMAX),g,rprandl,Cm,
     >       ter,ztop
      PARAMETER (g = 9.806, Cm = 0.21, rprandl = 3.)
      real dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,dtdz,
     >       dudzet,dvdzet,dwdzet,dtdzet,dzetdx,dzetdy,dzetdz,
     >       aved12,aved13,aved23,aved11,aved22,aved33,N2


      do k=1,nz-1
        dzet(k) = sgz(k+1) - sgz(k)
        dzetw(k) = wgz(k+1) - wgz(k)
        delta(k) = (wgz(k+1) - wgz(k))*delx*(Cm**2.)
      enddo

      dzet(nz) = 1.
      dzetw(nz) = wgz(nz+1) - wgz(nz)
      ztop = wgz(nz+1)


      do n=1,nt

      do k=1,nz
      do j=1,ny
      do i=1,nx

        kp1 = min(k+1,nz)
        km1 = max(k-1,1)
        jp1 = min(j+1,ny)
        jm1 = max(j-1,1)
        ip1 = min(i+1,nz)
        im1 = max(i-1,1)
*
*  d11,d22,d33 at thermodynamic pts
*
        dzetdx = (sgz(k)-ztop)/(ztop-terp(i,j))
     >               * (teru(i+1,j)-teru(i,j))/delx

        dzetdy = (sgz(k)-ztop)/(ztop-terp(i,j))
     >               * (terv(i,j+1)-terv(i,j))/dely

        dzetdz = ztop/(ztop-terp(i,j))

        dudzet = (u(i,j,kp1,n)+u(i+1,j,kp1,n)-u(i,j,k,n)-u(i+1,j,k,n))
     >                      / (4.*dzet(k))
     >         + (u(i,j,k,n)+u(i+1,j,k,n)-u(i,j,km1,n)-u(i+1,j,km1,n))
     >                      / (4.*dzet(km1))

        dvdzet = (v(i,j,kp1,n)+v(i,j+1,kp1,n)-v(i,j,k,n)-v(i,j+1,k,n))
     >                      / (4.*dzet(k))
     >         + (v(i,j,k,n)+v(i,j+1,k,n)-v(i,j,km1,n)-v(i,j+1,km1,n))
     >                      / (4.*dzet(km1)) 

        dwdzet = (w(i,j,k+1,n)-w(i,j,k,n))/(wgz(k+1)-wgz(k))

        dudx = (u(i+1,j,k,n)-u(i,j,k,n))/delx + dudzet*dzetdx

        dvdy = (v(i,j+1,k,n)-v(i,j,k,n))/dely + dvdzet*dzetdy

        dwdz = dwdzet*dzetdz

        d11(i,j,k,n) = 4.*dudx/3. - 2.*(dvdy+dwdz)/3.
        d22(i,j,k,n) = 4.*dvdy/3. - 2.*(dudx+dwdz)/3.
        d33(i,j,k,n) = 4.*dwdz/3. - 2.*(dudx+dvdy)/3.
*
*  d12 at z-vorticity pts
*
        ter = (teru(i,j)+teru(i,jm1)+terv(i,j)+terv(im1,j))/4.

        dzetdx = (sgz(k)-ztop)/(ztop-ter)
     >                 *(terv(i,j)-terv(im1,j))/delx

        dzetdy = (sgz(k)-ztop)/(ztop-ter)
     >                 *(teru(i,j)-teru(i,jm1))/dely

        dudzet = (u(i,j,kp1,n)+u(i,jm1,kp1,n)-u(i,j,k,n)-u(i,jm1,k,n))
     >                      /(4.*dzet(k))
     >         + (u(i,j,k,n)+u(i,jm1,k,n)-u(i,j,km1,n)-u(i,jm1,km1,n))
     >                      /(4.*dzet(km1))

        dvdzet = (v(i,j,kp1,n)+v(im1,j,kp1,n)-v(i,j,k,n)-v(im1,j,k,n))
     >                      /(4.*dzet(k))
     >         + (v(i,j,k,n)+v(im1,j,k,n)-v(i,j,km1,n)-v(im1,j,km1,n))
     >                      /(4.*dzet(km1))

        dudy = (u(i,j,k,n)-u(i,jm1,k,n))/dely + dudzet*dzetdy

        dvdx = (v(i,j,k,n)-v(im1,j,k,n))/delx + dvdzet*dzetdx

        d12(i,j,k,n) = dudy + dvdx
*
*  d13 at y-vorticity pts
* 
        dzetdx = (wgz(k)-ztop)/(ztop-teru(i,j))
     >                 *(terp(i,j)-terp(im1,j))/delx

        dzetdz = ztop/(ztop-teru(i,j))

        dudzet = (u(i,j,k,n)-u(i,j,km1,n))/dzet(km1)

        dwdzet = (w(i,j,k+1,n)+w(im1,j,k+1,n)-w(i,j,k,n)-w(im1,j,k,n))
     >                    /(4.*dzetw(k))
     >         + (w(i,j,k,n)+w(im1,j,k,n)-w(i,j,km1,n)-w(im1,j,km1,n))
     >                    /(4.*dzetw(km1))

        dwdx = (w(i,j,k,n)-w(im1,j,k,n))/delx + dwdzet*dzetdx

        dudz = dudzet*dzetdz

        d13(i,j,k,n) = dudz + dwdx
*
*  d23 at x-vorticity pts
*
        dzetdy = (wgz(k)-ztop)/(ztop-terv(i,j))
     >                 *(terp(i,j)-terp(i,jm1))/dely

        dzetdz = ztop/(ztop-terv(i,j))

        dvdzet = (v(i,j,k,n)-v(i,j,km1,n))/dzet(km1)

        dwdzet = (w(i,j,k+1,n)+w(i,jm1,k+1,n)-w(i,j,k,n)-w(i,jm1,k,n))
     >                        /(4.*dzetw(k))
     >         + (w(i,j,k,n)+w(i,jm1,k,n)-w(i,j,km1,n)-w(i,jm1,km1,n))
     >                        /(4.*dzetw(km1))

        dvdz = dvdzet*dzetdz

        dwdy = (w(i,j,k,n)-w(i,jm1,k,n))/dely + dwdzet*dzetdy

        d23(i,j,k,n) = dvdz + dwdy

      enddo
      enddo
      enddo


      do k=2,nz  
      do j=2,ny-1
      do i=2,nx-1
*
*  deformation squared at w pts
*
        aved12 = (d12(i,j,k,n)+d12(i+1,j,k,n)
     >                           +d12(i,j+1,k,n)+d12(i+1,j+1,k,n))/8.
     >         + (d12(i,j,k-1,n)+d12(i+1,j,k-1,n)
     >                       +d12(i,j+1,k-1,n)+d12(i+1,j+1,k-1,n))/8.

        aved13 = (d13(i,j,k,n)+d13(i+1,j,k,n))/2.

        aved23 = (d23(i,j,k,n)+d23(i,j+1,k,n))/2.

        aved11 = (d11(i,j,k,n)+d11(i,j,k-1,n))/2.

        aved22 = (d22(i,j,k,n)+d22(i,j,k-1,n))/2.

        aved33 = (d33(i,j,k,n)+d33(i,j,k-1,n))/2.

        mix(i,j,k,n) = (aved11**2. + aved22**2. + aved33**2.)/2.
     >          + aved12**2. + aved13**2. + aved23**2.

        if (iwhch.ne.1) then

          dzetdz = ztop/(ztop-terp(i,j))

          dtdzet = (thet(i,j,k,n)-thet(i,j,k-1,n))/dzet(k)

          dtdz = dtdzet*dzetdz

          N2 = 2.*g*dtdz/(thet(i,j,k,n)+thet(i,j,k-1,n))

          if (iwhch.eq.2) then
 
            mix(i,j,k,n) = mix(i,j,k,n)/N2

          else if ((rprandl*N2/mix(i,j,k,n)).lt.1.) then

            mix(i,j,k,n) = delta(k)*sqrt(mix(i,j,k,n) - rprandl*N2)

          else

            mix(i,j,k,n) = 0.
 
          end if

        end if 

      enddo
      enddo
      enddo
*
*  Fill in missing values around the edge
*
      do i=1,nx
      do j=1,ny
        mix(i,j,1,n) = miss
      enddo
      enddo

      do i=1,nx
      do k=1,nz
        mix(i,1,k,n) = miss
        mix(i,ny,k,n) = miss
      enddo
      enddo

      do j=1,ny
      do k=1,nz
        mix(1,j,k,n) = miss
        mix(nx,j,k,n) = miss
      enddo
      enddo
*
*  Done with this record.  Get next time.
*
      enddo

      return
      end


              


