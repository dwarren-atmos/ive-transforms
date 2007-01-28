
#include "pointer.icl"

      ive_ptr function div(axis,ndims,dims,stag,min,max,miss,
     >               dunits,disunits,dim_names)
*
*  div calculates the divergence of the velocity field.  Options are:
*    1)  dw/dz + dv/dy
*    2)  dw/dz + du/dx
*    3)  du/dx + dv/dy 
*    0)  du/dx + dv/dy + dw/dz
*
      integer axis,ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer ndims2,dims2(4)
      real stag2(4),min2(4),max2(4),miss2
      character*(80) dunits2,disunits2,dim_names2(4)

      integer twod,k
      real dom_min(4),dom_max(4),dx,dy,ztop,r_val,min_p(2),max_p(2)
      logical err,new
      ive_ptr getmem,getvar,uptr,vptr,wptr,zbot_p,zbot_u,zbot_v,
     >                       sgz,wgz

      div = 0

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

      zbot_p = getvar('zbot_p',ndims,dims,stag,min,max,miss,
     >                  dunits,disunits,dim_names,new)
      if (zbot_p.eq.0) then
        write(*,*) 'Error:  unable to access field zbot_p.'
        return
      end if

      do k=1,2
        min_p(k) = min(k)
        max_p(k) = max(k)
      enddo
      dx = (max(1) - min(1))/(dims(1) - 1)
      if (twod.ne.1) then
        dy = (max(2) - min(2))/(dims(2) - 1)
      else
        dy = 1.
      end if

      zbot_u = getmem((dims(1)+1)*dims(2))
      if (zbot_u.eq.0) then
        write(*,*) 'Error:  unable to allocate space for zbot_u.'
        return
      end if
      call zs_2_upts(%val(zbot_u),%val(zbot_p),dims(1),dims(2))

      if (twod.ne.1) then
        zbot_v = getmem(dims(1)*(dims(2)+1))
        if (zbot_v.eq.0) then
          write(*,*) 'Error:  unable to allocate space for zbot_v.'
          call freemem(zbot_u)
          return
        end if
        call zs_2_vpts(%val(zbot_v),%val(zbot_p),dims(1),dims(2))
      else
        zbot_v = 0
      end if

      wgz = getvar('wgz',ndims,dims,stag,min,max,miss,dunits,
     >                   disunits,dim_names,new)
      if (wgz.eq.0) then
        write(*,*) 'Error:  unable to access field wgz.'
        call freemem(zbot_u)
        call freemem(zbot_v)
        return
      end if
      ztop = r_val(%val(wgz),dims(1))

      sgz = getvar('sgz',ndims,dims,stag,min,max,miss,dunits,
     >                      disunits,dim_names,new)
      if (sgz.eq.0) then
        write(*,*) 'Error:  unable to access field sgz.'
        call freemem(zbot_u)
        call freemem(zbot_v)
        return
      end if

      if (axis.eq.1) then

        wptr = getvar('w',ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
        if (wptr.eq.0) then
          write(*,*) 'Error:  unable to access vertical velocity.'
          call freemem(zbot_u)
          call freemem(zbot_v)
          return
        end if
        dims(3) = dims(3)-1
        min(3) = r_val(%val(sgz),1)
        max(3) = r_val(%val(sgz),dims(3))
        stag(3) = 0.5
        dim_names(3) = 'nz'

        div = getmem(dims(1)*dims(2)*dims(3)*dims(4))
        if (div.eq.0) then
          write(*,*) 'Error:  unable to allocate space for div.'
          call freemem(zbot_u)
          call freemem(zbot_v)
        end if
        do k=1,dims(1)*dims(2)*dims(3)*dims(4)
          call s_val(%val(div),k,0.)
        enddo

        call w_div(%val(div),%val(wptr),%val(zbot_p),%val(wgz),ztop,
     >               dims(1),dims(2),dims(3),dims(4))

        if (twod.ne.1) then

          vptr = getvar('v',ndims2,dims2,stag2,min2,max2,miss2,
     >               dunits2,disunits2,dim_names2,new)
          if (vptr.eq.0) then
            write(*,*) 'Error:  unable to find v field.'
            call freemem(zbot_u)
            call freemem(zbot_v)
            return
          end if
          
          call v_div(%val(div),%val(vptr),%val(zbot_p),%val(zbot_v),
     >             %val(sgz),ztop,dy,dims(1),dims(2),dims(3),dims(4))

        end if

      else if (axis.eq.2) then

        wptr = getvar('w',ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
        if (wptr.eq.0) then
          write(*,*) 'Error:  unable to access vertical velocity.'
          call freemem(zbot_u)
          call freemem(zbot_v)
          return
        end if
        dims(3) = dims(3)-1
        min(3) = r_val(%val(sgz),1)
        max(3) = r_val(%val(sgz),dims(3))
        stag(3) = 0.5
        dim_names(3) = 'nz'

        div = getmem(dims(1)*dims(2)*dims(3)*dims(4))
        if (div.eq.0) then
          write(*,*) 'Error:  unable to allocate space for div.'
          call freemem(zbot_u)
          call freemem(zbot_v)
        end if
        
        do k=1,dims(1)*dims(2)*dims(3)*dims(4)
          call s_val(%val(div),k,0.)
        enddo

        call w_div(%val(div),%val(wptr),%val(zbot_p),%val(wgz),ztop,
     >               dims(1),dims(2),dims(3),dims(4))

        uptr = getvar('u',ndims2,dims2,stag2,min2,max2,miss2,dunits2,
     >                       disunits2,dim_names2,new)
        if (uptr.eq.0) then
          write(*,*) 'Error:  unable to access u field.'
          call freemem(zbot_u)
          call freemem(zbot_v)
          return
        end if

        call u_div(%val(div),%val(uptr),%val(zbot_p),%val(zbot_u),
     >               %val(sgz),ztop,dx,dims(1),dims(2),dims(3),
     >                           dims(4))

      else if (axis.eq.3) then

        uptr = getvar('u',ndims,dims,stag,min,max,miss,dunits,
     >                   disunits,dim_names,new)
        if (uptr.eq.0) then
          write(*,*) 'Error:  unable to access u field.'
          call freemem(zbot_u)
          call freemem(zbot_v)
          return
        end if
        dims(1) = dims(1)-1
        min(1) = min_p(1)
        max(1) = max_p(1)
        stag(1) = 0.5
        dim_names(1) = 'nx'

        div = getmem(dims(1)*dims(2)*dims(3)*dims(4))
        if (div.eq.0) then
          write(*,*) 'Error:  unable to allocate space for div.'
          call freemem(zbot_u)
          call freemem(zbot_v)
        end if
        do k=1,dims(1)*dims(2)*dims(3)*dims(4)
          call s_val(%val(div),k,0.)
        enddo

        call u_div(%val(div),%val(uptr),%val(zbot_p),%val(zbot_u),
     >               %val(sgz),ztop,dx,dims(1),dims(2),dims(3),
     >                               dims(4))

        if (twod.ne.1) then

          vptr = getvar('v',ndims2,dims2,stag2,min2,max2,miss2,
     >               dunits2,disunits2,dim_names2,new)
          if (vptr.eq.0) then
            write(*,*) 'Error:  unable to find v field.'
            call freemem(zbot_u)
            call freemem(zbot_v)
            return
          end if
          
          call v_div(%val(div),%val(vptr),%val(zbot_p),%val(zbot_v),
     >             %val(sgz),ztop,dy,dims(1),dims(2),dims(3),dims(4))

        end if

      else

        uptr = getvar('u',ndims,dims,stag,min,max,miss,dunits,
     >                   disunits,dim_names,new)
        if (uptr.eq.0) then
          write(*,*) 'Error:  unable to access u field.'
          call freemem(zbot_u)
          call freemem(zbot_v)
          return
        end if
        dims(1) = dims(1)-1
        min(1) = min_p(1)
        max(1) = max_p(1)
        stag(1) = 0.5
        dim_names(1) = 'nx'

        div = getmem(dims(1)*dims(2)*dims(3)*dims(4))
        if (div.eq.0) then
          write(*,*) 'Error:  unable to allocate space for div.'
          call freemem(zbot_u)
          call freemem(zbot_v)
        end if
        do k=1,dims(1)*dims(2)*dims(3)*dims(4)
          call s_val(%val(div),k,0.)
        enddo

        call u_div(%val(div),%val(uptr),%val(zbot_p),%val(zbot_u),
     >              %val(sgz),ztop,dx,dims(1),dims(2),dims(3),
     >                               dims(4)) 

        wptr = getvar('w',ndims2,dims2,stag2,min2,max2,miss2,
     >            dunits2,disunits2,dim_names2,new)
        if (wptr.eq.0) then
          write(*,*) 'Error:  unable to access vertical velocity.'
          call freemem(zbot_u)
          call freemem(zbot_v)
          return
        end if

        call w_div(%val(div),%val(wptr),%val(zbot_p),%val(wgz),ztop,
     >                dims(1),dims(2),dims(3),dims(4))
 
        if (twod.ne.1) then

          vptr = getvar('v',ndims2,dims2,stag2,min2,max2,miss2,
     >               dunits2,disunits2,dim_names2,new)
          if (vptr.eq.0) then
            write(*,*) 'Error:  unable to find v field.'
            call freemem(zbot_u)
            call freemem(zbot_v)
            return
          end if
          
          call v_div(%val(div),%val(vptr),%val(zbot_p),%val(zbot_v),
     >             %val(sgz),ztop,dy,dims(1),dims(2),dims(3),dims(4))

        end if

      end if
  
      miss = 8675309.
      dunits = '1/s'
      disunits = '1/s'

      call freemem(zbot_u)
      call freemem(zbot_v)

      return
      end

*------------------------------------------------------------------------

      subroutine w_div(div,w,zs,wgz,ztop,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real div(nx,ny,nz,nt),w(nx,ny,nz+1,nt),zs(nx,ny),wgz(nz+1)
      real ztop

      integer MAXNZ,i,j,k,n

      parameter (MAXNZ = 251)

      real dzet(MAXNZ),dzetdz

      if (nz.gt.MAXNZ) then
        write(*,*) 'Error:  nz greater than MAXNZ in w_div.'
        write(*,*) 'Stopping ive.'
        stop
      end if

      do k=1,nz
        dzet(k) = wgz(k+1)-wgz(k)
      enddo

      do n=1,nt
      do j=1,ny
      do i=1,nx

        dzetdz = ztop/(ztop-zs(i,j))

        do k=1,nz

          div(i,j,k,n) = div(i,j,k,n) 
     >           + (w(i,j,k+1,n) - w(i,j,k,n))*dzetdz/dzet(k)

        enddo

      enddo
      enddo
      enddo

      return
      end      

*----------------------------------------------------------------------

      subroutine v_div(div,v,zbot_p,zbot_v,sgz,ztop,dy,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real div(nx,ny,nz,nt),v(nx,ny+1,nz,nt),zbot_p(nx,ny)
      real zbot_v(nx,ny+1),sgz(nz),ztop,dy

      integer MAXNZ,i,j,k,n,kp1,km1

      parameter (MAXNZ=251)

      real zlbnd(MAXNZ),zrbnd(MAXNZ),dvdzet,dzetdy,dzet(MAXNZ)

      if (nz.gt.MAXNZ) then
        write(*,*) 'Error:  nz greater than MAXNZ in v_div.'
        write(*,*) 'Stopping ive.'
        stop
      end if

      do k=1,nz-1
        dzet(k) = sgz(k+1)-sgz(k)
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

        dvdzet = 0.25*zrbnd(k)*(v(i,j,kp1,n)-v(i,j,k,n)
     >            + v(i,j+1,kp1,n)-v(i,j+1,k,n))/dzet(k)  +
     >           0.25*zlbnd(k)*(v(i,j,k,n)-v(i,j,km1,n)
     >            + v(i,j+1,k,n)-v(i,j+1,km1,n))/dzet(km1)

        dzetdy = (sgz(k)-ztop)/(ztop-zbot_p(i,j))
     >                 * (zbot_v(i,j+1)-zbot_v(i,j))/dy

        div(i,j,k,n) = div(i,j,k,n) + 
     >       (v(i,j+1,k,n)-v(i,j,k,n))/dy + dvdzet*dzetdy
 
      enddo
      enddo
      enddo
      enddo

      return
      end

*------------------------------------------------------------------------

      subroutine u_div(div,u,zbot_p,zbot_u,sgz,ztop,dx,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real div(nx,ny,nz,nt),u(nx+1,ny,nz,nt),zbot_p(nx,ny)
      real zbot_u(nx+1,ny),sgz(nz),ztop,dx

      integer MAXNZ,i,j,k,n,kp1,km1

      parameter (MAXNZ=251)

      real zlbnd(MAXNZ),zrbnd(MAXNZ),dudzet,dzetdx,dzet(MAXNZ)

      if (nz.gt.MAXNZ) then
        write(*,*) 'Error:  nz greater than MAXNZ in u_div.'
        write(*,*) 'Stopping ive.'
        stop
      end if

      do k=1,nz-1
        dzet(k) = sgz(k+1)-sgz(k)
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

        dudzet = 0.25*zrbnd(k)*(u(i,j,kp1,n)-u(i,j,k,n) 
     >          + u(i+1,j,kp1,n)-u(i+1,j,k,n))/dzet(k)  +
     >           0.25*zlbnd(k)*(u(i,j,k,n)-u(i,j,km1,n) 
     >          + u(i+1,j,k,n)-u(i+1,j,km1,n))/dzet(km1)

        dzetdx = (sgz(k)-ztop)/(ztop-zbot_p(i,j)) *
     >               (zbot_u(i+1,j)-zbot_u(i,j))/dx

        div(i,j,k,n) = div(i,j,k,n) 
     >          + (u(i+1,j,k,n)-u(i,j,k,n))/dx + dudzet*dzetdx

      enddo
      enddo
      enddo
      enddo

      return
      end
