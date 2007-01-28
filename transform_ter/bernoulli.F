

#include "pointer.icl"

      ive_ptr function bernoulli(ndims,dims,stag,min,max,miss,dunits,
     >                                 disunits,dim_names)

#include "default.icl"

      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)

      integer ijk,nt
      logical new,err
      real N,N2
      parameter (N=0.012)
      ive_ptr getmem,u,v,w,p,b,b0

      nt = ntime
      N2 = N*N

      bernoulli = getmem(nx*ny*nz*nt)
      if (bernoulli.eq.0) then
        write(6,*) 'Error:  could not allocate space for bernoulli.'
        return
      end if

      do ijk=1,nx*ny*nz*nt
        call s_val(%val(bernoulli),ijk,0.)
      enddo

      u = getvar('u',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                          dim_names,new)
      if (u.eq.0) then
        write(6,*) 'Error:  could not access u.'
        call freemem(bernoulli)
        return
      end if

      call bernoulli_calc_I(%val(bernoulli),%val(u),nx,ny,nz,nt,
     >                                 1,0,0)
      call driver('free = u',0)

      v = getvar('v',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                          dim_names,new)
      if (v.eq.0) then
        write(6,*) 'Error:  could not access v.'
        call freemem(bernoulli)
        return
      end if

      call bernoulli_calc_I(%val(bernoulli),%val(v),nx,ny,nz,nt,
     >                                 0,1,0)
      call driver('free = v',0)

      w = getvar('w',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                          dim_names,new)
      if (w.eq.0) then
        write(6,*) 'Error:  could not access w.'
        call freemem(bernoulli)
        return
      end if

      call bernoulli_calc_I(%val(bernoulli),%val(w),nx,ny,nz,nt,
     >                                 0,0,1)
      call driver('free = w',0)

      p = getvar('p',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                          dim_names,new)
      if (p.eq.0) then
        write(6,*) 'Error:  could not access p.'
        call freemem(bernoulli)
        return
      end if

      call bernoulli_calc_I(%val(bernoulli),%val(p),nx,ny,nz,nt,
     >                                 0,0,0)
      call driver('free = p',0)

      b0 = getvar('b0',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                            dim_names,new)
      if (b0.eq.0) then
        write(6,*) 'Error:  could not access b0.'
        call freemem(bernoulli)
        return
      end if

      b = getvar('b',ndims,dims,stag,min,max,miss,dunits,disunits,
     >                          dim_names,new)
      if (b.eq.0) then
        write(6,*) 'Error:  could not access b.'
        call freemem(bernoulli)
        return
      end if

      call bernoulli_calc_II(%val(bernoulli),%val(b),%val(b0),
     >                         %val(zbot_p),%val(sgz),ztop,N2,
     >                         nx,ny,nz,nt)

      call driver('free = b',0)
      call driver('free = b0',0)

      dunits = 'm*m/s/s'
      disunits = 'm*m/s/s'
      
      return
      end

*-------------------------------------------------------------------------

      subroutine bernoulli_calc_I(bernoulli,var,nx,ny,nz,nt,iu,iv,iw)

      integer nx,ny,nz,nt,iu,iv,iw
      real bernoulli(nx,ny,nz,nt),var(nx+iu,ny+iv,nz+iw,nt)

      integer i,j,k,n

      if (iu+iv+iw.ne.0) then

        do n=1,nt
        do k=1,nz
        do j=1,ny
        do i=1,nx

          bernoulli(i,j,k,n) = bernoulli(i,j,k,n) 
     >           + 0.125*(var(i,j,k,n)+var(i+iu,j+iv,k+iw,n))**2

        enddo
        enddo
        enddo
        enddo

      else

        do n=1,nt
        do k=1,nz
        do j=1,ny
        do i=1,nx

          bernoulli(i,j,k,n) = bernoulli(i,j,k,n) + var(i,j,k,n)

        enddo
        enddo
        enddo
        enddo

      end if

      return
      end

*--------------------------------------------------------------------------

      subroutine bernoulli_calc_II(bernoulli,b,b0,zs_p,sgz,ztop,N2,
     >                                  nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real bernoulli(nx,ny,nz,nt),b(nx,ny,nz,nt),b0(nx,ny,nz)
      real zs_p(nx,ny),sgz(nz),ztop,N2

      integer i,j,k,n
      real z0,z,tfac

      do n=1,nt
      do k=1,nz
      do j=1,ny
      do i=1,nx

        z0 = zs_p(i,j)
        tfac = (ztop-z0)/ztop
        z = z0 + tfac*sgz(k)

        bernoulli(i,j,k,n) = bernoulli(i,j,k,n) 
     >           - (b(i,j,k,n)-b0(i,j,k))*z - 0.5*N2*z**2

      enddo
      enddo
      enddo
      enddo

      return
      end

