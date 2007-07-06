

      integer(ive_ptr) function geowind(wind,ndims,dims,stag,min,max,
     &                        miss,dunits,disunits,dim_names)
      use default_module
      implicit none

c
c  Calculates the South-North and West-East (NSEW) derivatives for 
c  from the model coordinates
c
      integer ndims,dims(4)
      real stag(4),min(4),max(4),miss
      character*(*) wind,dunits,disunits,dim_names(4)
c
	  real r_val,pm,fm,rho_tmp,deriv_tmp
	  real slope,intercept
      integer ntotal,i,k,n,j
	  logical requ
      character*(1) dir
      integer(ive_ptr) rho,deriv_PPP,deriv,rho_calc
c
      rho = 0
      deriv_PPP = 0
      geowind = 0
c
      ntotal = (nx+1)*(ny+1)*nz*ntime
c
c  Set Direction of Derivative  
c
      if(wind(1:1).eq.'U') then
	  dir = 'Y'
	  pm = -1.0
      else if(wind(1:1).eq.'V') then
	  dir = 'X'
        pm = 1.0
      else
	  write(*,*) 'Error: Not a valid geostrophic component'
	  return
      end if
	  
      deriv_PPP = deriv('PPP',dir,1,ndims,dims,stag,min,max,
     >		    miss,dunits,disunits,dim_names)
      if (deriv_PPP.eq.0) then
        write(*,*) 'Error:  cannot calculate PPP derivative'
        return
      end if
c
c  Get RHO
c
	rho = rho_calc(ndims,dims,stag,min,max,miss,dunits,
     >               disunits,dim_names)
      if (rho.eq.0) then
        write(*,*) 'Error:  cannot calculate rho'
        return
      end if
c
 	geowind = getmem(ntotal)
      if (geowind.eq.0) then
        write(*,*) 'Error:  cannot get memory for geowind'
        return
      end if
c
c	call convert_units('millibars','Pascals',slope,intercept)
	slope = 100.0 
	intercept = 0.0
	i = 0
	do j=1,ntotal
	  i = i + 1
	   
        if(i.gt.(nx+1)*(ny+1)) i = i - (nx+1)*(ny+1)
	  rho_tmp = r_val(%val(rho),j)
	  deriv_tmp = r_val(%val(deriv_PPP),j)
	  fm = r_val(%val(fcor),i)
c	  fm = 0.00001

      call s_val(%val(geowind),j,
     >  real((pm*deriv_tmp)/(rho_tmp*fm)*slope + intercept)) 
      end do
c
	call freemem(rho)
	call freemem(deriv_PPP)

	dims(1) = nx+1
	dims(2) = ny+1
	dims(3) = nz
	dims(4) = ntime
	disunits = 'm/s'
	dunits = 'm/s'

      return
      end
