	
      integer(ive_ptr) function test(ndims,dims,stag,min,max,miss,
     >                      dunits,disunits,dim_names)
      use default_module
c
      implicit none
c
	integer ndims,dims(4)
	real stag(4),min(4),max(4),miss
	character*(*) dunits,disunits,dim_names(4)
c
	integer ntotal,i
	logical new
	integer(ive_ptr) uvar
c
c
c  Get U variable
c
	uvar = 0
	uvar = getvar('UTR',ndims,dims,stag,min,max,miss,dunits,
     >			disunits,dim_names,new)
	if (uvar.eq.0) then
	  write(6,*) 'Error:  Could not get variable U'
	  return
	endif
c
	ntotal = (nx+1)*(ny+1)*(nz+1)*ntime
c
	test = 0
	test = getmem(ntotal)
	if(test.eq.0) then
	  write(6,*) 'Error:  Cannot allocate memory for test'
	  return
	end if
c
	do i=1,ntotal
	  call s_val(%val(test),i,1.)
	enddo
c
	call var_at_w(%val(test),%val(uvar),%val(sgz),%val(wgz),
     >			  1,0,0,nx,ny,nz,ntime)
c
	dims(1) = nx+1
	dims(2) = ny+1
	dims(3) = nz+1
c
	dunits = 'm/s'	
	disunits = 'm/s'	
c
	return
	end
c			    
c#################################################################
