
#include "pointer.icl"

      subroutine sound_out

#include "default.icl"
c
c This function output a vertical sounding for a 2d dry run
c The sounding variable output are u and theta and are
c horizontal averages over the domain.
c
      integer ndims,dims(4),isqr
      real stag(4),min(4),max(4),miss
      character*(80) dunits,disunits,dim_names(4)      

	integer k
	real uavg(nz),thetavg(nz),foo,ths,ps,us
	character*12 prefix

	real r_val

      logical new
      ive_ptr u,thet,getmem
c
      thet = getvar('thet',ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
      if (thet.eq.0) then
        write(6,*) 'Error:  unable to access variable thet.'
        return
      end if
c
      u = getvar('u',ndims,dims,stag,min,max,miss,dunits,
     >                       disunits,dim_names,new)
      if (u.eq.0) then
        write(6,*) 'Error:  unable to access variable u.'
        return
      end if
c
      call domain_avg(%val(thet),%val(u),thetavg,uavg,nx,ny,nz,ntime) 
c
40    FORMAT(4(F10.2,2x),F10.2)
c
	prefix = '200404210144' 
	print *,prefix
	open (unit=45,file=prefix//'.sound',status='old')
c
	read(45,*) ps
	read(45,*) foo,ths,foo,us,foo
c
	close(45)
	  
      open (unit=50,file=prefix//'_500hr.sound',status='new')

	write(50,*) ps
	write(50,40) 0.,ths,0.,us,0.
c
	do k=1,nz
	 write(50,40) r_val(%val(sgz),k),thetavg(k),0.,uavg(k),0.
	end do
c
	close(50)
      return
      end
  
*------------------------------------------------------------------------

      subroutine domain_avg(thet,u,thetavg,uavg,nx,ny,nz,nt)

      integer nx,ny,nz,nt
      real thet(nx,ny,nz,nt),u(nx,ny,nz,nt)
	real uavg(nz),thetavg(nz)

      integer i,j,k
c
	do k = 1,nz
	  uavg(k) = 0
	  thetavg(k) = 0
	  do j = 1,ny
	    do i = 1,nx
     	     print *,i,j,k,u(i,j,k,nt),thet(i,j,k,nt)

	      thetavg(k) = thet(i,j,k,nt) + thetavg(k)
	      uavg(k) = u(i,j,k,nt) + uavg(k)
          end do
        end do
c
	  uavg(k) = uavg(k)/(nx*ny)
	  thetavg(k) = thetavg(k)/(nx*ny)
c
      end do

      return
      end

         
      

