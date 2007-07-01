#include "pointer.icl"
c
c  Calculates the Geopotential height of a passed pressure surface
c
      ive_ptr function geopot(p0,
     &  ndims,dims,stag,min,max,miss,dunits,
     &  disunits,dim_names)
      implicit none
c
#include "default.icl"
c
      integer ndims,dims(4)
      real p0,stag(4),min(4),max(4),miss
      character*(*) dunits,disunits,dim_names(4)
c      
      integer strbeg,strend,nxp1,nyp1,nt,len
      logical new
      ive_ptr pvar
c
      nxp1 = nx+1;
      nyp1 = ny+1;
c
c Get pressure field
c
      pvar = 0
      pvar = getvar('ppp',ndims,dims,stag,min,max,miss,dunits,
     &                          disunits,dim_names,new)
      if (pvar.eq.0) then
        write(6,*) 'Error:  cannot access variable PPP'
        return
      end if
c
      nt = dims(4)
c
      geopot = 0
      geopot = getmem(nxp1*nyp1*nt)
      if (geopot.eq.0) then
        write(6,*) 'Error:  cannot allocate space for geopot.'
        return
      end if
c
      len = nxp1*nyp1
      call calc_geopot(%val(geopot),p0,%val(pvar),%val(sgz),
     &                 %val(zbot_p),ztop,miss,len,nz,nt)
c
      dims(1) = nxp1 ; dims(2) = nyp1
      dims(3) = 1 ; dims(4) = nt
      dim_names(1) = 'nx' ; dim_names(2) = 'ny'
      dim_names(3) = 'one' ; dim_names(4) = 'time'
      min(3) = 0 ; max(3) = 0
      dunits = 'm' ; disunits = 'm'
c
      return
      end
c
c
      subroutine calc_geopot(fout,f0,fin,sgz,zsurf,ztop,missing,
     &                       len,nz,nt)

     integer len, nz, nt
	 real f0, fin(len,nz,nt), fout(len,nt), sgz(nz), zsurf(len)
    &     ztop, missing
c
     integer i,k,n,k1,k2
	 real temp, zfrac
c
	do n=1,nt ! outer time loop
c
c 	Assign missing date 
c
	  do i=1,len
	    if(f0.gt.fin(i,1,n).or.f0.lt.fin(i,nz,n)) then
	 	fout(i,n) = missing
	    end if
   	  end do
c
	  do k=2,nz
	    do i=1,len
	      if ( f0.lt.fin(i,k-1,n).and.f0.ge.fin(i,k,n)) then
	    	  zfrac = (ztop - zsurf(i))/ztop
		  k1 = k-1
		  k2 = k
		  temp =  ( alog(f0) - alog(fin(i,k1,n)) )/
     >	          ( alog(fin(i,k2,n)) - alog(fin(i,k1,n)) )
	 	  temp = sgz(k1) + temp*( sgz(k2) - sgz(k1) ) 
	        fout(i,n) = temp*zfrac + zsurf(i)
  		end if
	    end do
 	  end do
	end do

	return
	end

