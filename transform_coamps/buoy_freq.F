
#include "pointer.icl"

      ive_ptr function buoy_freq(fld,ndims,dims,stag,min,max,
     >                 miss,dunits,disunits,dim_names,isqr)

      use default_module

      integer ndims,dims(4),isqr
      real stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4),fld      

      logical new
      ive_ptr thet

      buoy_freq = 0

      thet = getvar(fld,ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
      if (thet.eq.0) then
        write(6,*) 'Error:  unable to access variable ',fld,'.'
        return
      end if

      buoy_freq = getmem(dims(1)*dims(2)*dims(3)*dims(4))
      if (buoy_freq.eq.0) then
        write(6,*) 'Error:  unable to allocate memory for ',
     >                             'buoy_freq.'
        return
      end if

      miss = 0.060458

      call buoy_freq_calc(%val(buoy_freq),%val(thet),%val(zbot_p),
     >                     %val(sgz),ztop,miss,dims(1),dims(2),
     >                             dims(3),dims(4),isqr) 

      dunits = '1/s'
      disunits = '1/s'
 
      return
      end
*------------------------------------------------------------------------

      subroutine buoy_freq_calc(bvf,thet,zs_p,sgz,ztop,miss,
     >                               nx,ny,nz,nt,isqr)

      integer nx,ny,nz,nt,isqr
      real bvf(nx,ny,nz,nt),thet(nx,ny,nz,nt),zs_p(nx,ny),sgz(nz)
      real ztop,miss

      integer i,j,k,n,kp1,km1
      real dtdzet,dzetdz,N2

      integer nzmax
      parameter (nzmax=201)
      real zlbnd(nzmax),zrbnd(nzmax),dzet(nzmax)

      real g
      parameter (g=9.806)

      if (nz.gt.nzmax) then
        write(6,*) 'Error:  nz too big for buoy_freq calculation.'
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

      do n=1,nt ; do k=1,nz
        kp1 = min(k+1,nz)
        km1 = max(k-1,1)
        do j=1,ny ; do i=1,nx


          dtdzet = 0.5*zrbnd(k)*(thet(i,j,kp1,n)-thet(i,j,k,n))
     >                                      /dzet(k)
     >           + 0.5*zlbnd(k)*(thet(i,j,k,n)-thet(i,j,km1,n))
     >                                      /dzet(km1)
          dzetdz = ztop/(ztop-zs_p(i,j))


          N2 = g*dtdzet*dzetdz/( thet(i,j,k,n) )

          if (N2.le.0) then
            if (isqr.ne.0) then
              bvf(i,j,k,n) = N2
            else
c            bvf(i,j,k,n) = miss
              bvf(i,j,k,n) = 0.0
            end if
          else 
            if (isqr.ne.0) then
              bvf(i,j,k,n) = N2
            else
              bvf(i,j,k,n) = sqrt(N2)
            end if
          end if

        enddo ; enddo
      enddo ; enddo

      return
      end
