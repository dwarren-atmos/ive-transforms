
#include "pointer.icl"

      ive_ptr function deriv_nsew(varnm,dir,ndims,dims,stag,min,max,
     >			miss,dunits,disunits,dim_names)
      implicit none
c
#include "default.icl"
c
c  Calculates the South-North and West-East (NSEW) derivatives for 
c  from the model coordinates
c
      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) varnm,dir,dunits,disunits,dim_names(4)
c
	real r_val
      integer iu,iv,iw,it,ntotal,mn
	logical requ
      ive_ptr deriv,xderiv,yderiv,deriv_tmp,foo
c
      xderiv = 0
	yderiv = 0
	deriv_tmp = 0
c
	ntotal = (nx+1)*(ny+1)*nz*ntime
c
c  Get x and interpoloate to mass points
c
      deriv_tmp = deriv(varnm,'X',1,ndims,dims,stag,min,max,
     >		    miss,dunits,disunits,dim_names)
      if (deriv_tmp.eq.0) then
        write(*,*) 'Error:  cannot calculate x derivative'
        return
      end if
c
c  Set up stagering
c
	it = 0
	iu = 0
	iv = 0
	iw = 0
	if ((dims(1).eq.1).and.(dims(2).eq.1).and.(dims(3).eq.1)) then
	  it = 1
	else if (requ(plmin(1),min(1)-dx/2.)) then
	  iu = 1
	else if (requ(plmin(2),min(2)-dy/2.).and.(twod.ne.1)) then
	  iv = 1
	else if (requ(plmin(3),min(3))) then
	  iw = 1
	end if
c
c  Get memory and Interpolate to Mass points
c
	xderiv = getmem(ntotal)
      if (xderiv.eq.0) then
        write(*,*) 'Error:  get memory for xderiv'
        return
      end if
c
	call varatm(%val(xderiv),%val(deriv_tmp),
     >		nx,ny,dims(3),dims(4),iu,iv,iw)
	call freemem(deriv_tmp)
 	deriv_tmp = 0
c
c  Get y direction derivative
c
      deriv_tmp = deriv(varnm,'Y',1,ndims,dims,stag,min,max,
     >     		      miss,dunits,disunits,dim_names)
      if (deriv_tmp.eq.0) then
        write(*,*) 'Error:  cannot calculate y derivative'
        return
      end if
c
c  Get memory and Interpolate to Mass points
c
	yderiv = getmem(ntotal)
      if (yderiv.eq.0) then
        write(*,*) 'Error:  get memory for yderiv'
        return
      end if
c
	call varatm(%val(yderiv),%val(deriv_tmp),
     >		nx,ny,dims(3),dims(4),iu,iv,iw)
	call freemem(deriv_tmp)
c
      foo = getmem(ntotal)
	deriv_nsew = getmem(ntotal)
	mn = (nx+1)*(ny+1)
c
c  Rotate to N-S or E-W directions
c
	if(dir(1:1).eq.'N'.or.dir(1:1).eq.'S'.or.
     >			dir(1:1).eq.'Y') then
	  call uvg2uv(%val(xderiv),%val(yderiv),dims(3),mn, 
     >  		  %val(grdrot),%val(foo),%val(deriv_nsew))
	else if(dir(1:1).eq.'E'.or.dir(1:1).eq.'W'.or.
     >			dir(1:1).eq.'X') then
	  call uvg2uv(%val(xderiv),%val(yderiv),dims(3),mn, 
     >  		  %val(grdrot),%val(deriv_nsew),%val(foo))
      end if
c
	call freemem(foo)
	call freemem(xderiv)
	call freemem(yderiv)

	dims(1) = nx+1;
	dims(2) = ny+1;

      return
      end
