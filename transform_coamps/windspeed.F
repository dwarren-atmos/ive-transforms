
#include "pointer.icl"
c
c  Calculates the windspeed at mass points
c
      ive_ptr function windspeed(xwind,ywind,
     >                  ndims,dims,stag,min,max,miss,dunits,
     >                  disunits,dim_names)
      implicit none
c
#include "default.icl"
c
      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) xwind,ywind,dunits,disunits,dim_names(4)
c      
      integer strbeg,strend,nxp1,nyp1,nt
      logical new
      ive_ptr uvar,vvar
c
      nxp1 = nx+1;
      nyp1 = ny+1;
c
c Get wind fields
c
      uvar = 0; vvar=0;
      uvar = getvar(xwind,ndims,dims,stag,min,max,miss,dunits,
     >                           disunits,dim_names,new)
      vvar = getvar(ywind,ndims,dims,stag,min,max,miss,dunits,
     >                          disunits,dim_names,new)
      if (uvar.eq.0.or.vvar.eq.0) then
        write(6,*) 'Error:  cannot access variable xwind or ywind'
        return
      end if
c
      nt = dims(4)
c
      windspeed = 0
      windspeed = getmem(nxp1*nyp1*nz*nt)
      if (windspeed.eq.0) then
        write(6,*) 'Error:  cannot allocate space for windspeed.'
        return
      end if
c
         call calc_wspd(%val(windspeed),%val(uvar),%val(vvar),
     >    nx,ny,nz,nt)
c
      dims(1) = nxp1
      dims(2) = nyp1
      dims(3) = nz
      dims(4) = nt
      dim_names(1) = 'nx'
      dim_names(2) = 'ny'
      dim_names(3) = 'nz'
      dim_names(4) = 'time'
c      min(3) = 0
c      max(3) = 0
      dunits = 'm/s'
      disunits = 'm/s'
c
      return
      end
c
      subroutine calc_wspd(wspd,uwind,vwind,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real uwind(nx,ny+1,nz,nt),vwind(nx,ny+1,nz,nt)
      real wspd(nx+1,ny+1,nz,nt)
c
      real u_t(nx+1,ny+1,nz,nt), v_t(nx+1,ny+1,nz,nt)
      integer i,j,k,n
c

c      do n=1,nt; do k=1,nz; do i=1,nx+1; do j=1,ny+1;
c      print *,i,j,k,n,uwind(i,j,k,n)
c      end do; end do; end do; end do;
c
      u_t = 1; v_t = 1;
      call var_at_t(u_t,uwind,1,0,0,nx,ny,nz,nt)
      call var_at_t(v_t,vwind,0,1,0,nx,ny,nz,nt)

      do n=1,nt; do k=1,nz; do i=1,nx+1; do j=1,ny+1;
         wspd(i,j,k,n) = ( u_t(i,j,k,n)*u_t(i,j,k,n) 
     >                 +   v_t(i,j,k,n)*v_t(i,j,k,n) )**(0.5)
      end do; end do; end do; end do;

      return
      end

